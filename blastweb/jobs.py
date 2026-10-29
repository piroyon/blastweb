"""Background BLAST jobs stored on disk.

Each job is a directory under job_dir named by a random job ID:
    query.fa     query sequences
    params.json  search parameters (for display)
    status.json  {"state": "running"|"done"|"error", "started", "finished", "error"}
    archive.asn  BLAST archive (outfmt 11), reformatted with blast_formatter
    result.json  outfmt 15 JSON used to render the result page

State lives on disk so it is shared between gunicorn workers.
"""
import json
import os
import re
import shutil
import subprocess
import threading
import time
import uuid

JOB_ID_RE = re.compile(r"^[0-9a-f]{32}$")

# name -> (outfmt, file extension, content type)
DOWNLOAD_FORMATS = {
    "pairwise": ("0", "txt", "text/plain"),
    "tabular": ("6", "tsv", "text/tab-separated-values"),
    "tabular-header": ("7", "tsv", "text/tab-separated-values"),
    "csv": ("10", "csv", "text/csv"),
    "xml": ("5", "xml", "application/xml"),
    "json": ("15", "json", "application/json"),
    "asn": ("11", "asn", "text/plain"),
}


def _write_json(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f)
    os.replace(tmp, path)


def _read_json(path):
    with open(path) as f:
        return json.load(f)


def job_path(job_dir, job_id):
    if not JOB_ID_RE.match(job_id):
        return None
    path = os.path.join(job_dir, job_id)
    return path if os.path.isdir(path) else None


def start_job(job_dir, query_fasta, params, search_cmd, formatter_cmd, timeout):
    """Create a job directory and run BLAST in a background thread. Returns the job ID."""
    job_id = uuid.uuid4().hex
    path = os.path.join(job_dir, job_id)
    os.makedirs(path)

    query_path = os.path.join(path, "query.fa")
    with open(query_path, "w") as f:
        f.write(query_fasta)
    _write_json(os.path.join(path, "params.json"), params)
    _write_json(os.path.join(path, "status.json"), {"state": "running", "started": time.time(), "timeout": timeout})

    cmd = search_cmd + ["-query", query_path, "-out", os.path.join(path, "archive.asn"), "-outfmt", "11"]
    thread = threading.Thread(target=_run_job, args=(path, cmd, formatter_cmd, timeout), daemon=True)
    thread.start()
    return job_id


def _run_job(path, cmd, formatter_cmd, timeout):
    status_path = os.path.join(path, "status.json")
    status = _read_json(status_path)
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
        subprocess.run(
            [formatter_cmd, "-archive", os.path.join(path, "archive.asn"),
             "-outfmt", "15", "-out", os.path.join(path, "result.json")],
            check=True, capture_output=True, text=True, timeout=timeout
        )
        status["state"] = "done"
    except subprocess.CalledProcessError as e:
        status["state"] = "error"
        status["error"] = f"BLAST error: {e.stderr.strip() or e}"
    except subprocess.TimeoutExpired:
        status["state"] = "error"
        status["error"] = f"BLAST timed out after {timeout} seconds"
    except Exception as e:
        status["state"] = "error"
        status["error"] = f"Unexpected error: {e}"
    status["finished"] = time.time()
    _write_json(status_path, status)


def get_status(path):
    status = _read_json(os.path.join(path, "status.json"))
    # A job still "running" long after its timeout lost its worker (e.g. server restart)
    if status["state"] == "running" and time.time() - status["started"] > status.get("timeout", 600) * 2 + 60:
        status["state"] = "error"
        status["error"] = "The job was interrupted. Please search again."
    return status


def get_params(path):
    return _read_json(os.path.join(path, "params.json"))


def get_result(path):
    return _read_json(os.path.join(path, "result.json"))


def count_running(job_dir):
    if not os.path.isdir(job_dir):
        return 0
    count = 0
    for name in os.listdir(job_dir):
        path = job_path(job_dir, name)
        if path:
            try:
                if get_status(path)["state"] == "running":
                    count += 1
            except (OSError, ValueError, KeyError):
                pass
    return count


def cleanup_old_jobs(job_dir, retention_days):
    if not os.path.isdir(job_dir):
        return
    cutoff = time.time() - retention_days * 86400
    for name in os.listdir(job_dir):
        path = job_path(job_dir, name)
        if path and os.path.getmtime(path) < cutoff:
            shutil.rmtree(path, ignore_errors=True)


def format_result(path, fmt, formatter_cmd, timeout):
    """Return (bytes, extension, content type) for a download format."""
    outfmt, ext, content_type = DOWNLOAD_FORMATS[fmt]
    archive = os.path.join(path, "archive.asn")
    if fmt == "asn":
        with open(archive, "rb") as f:
            return f.read(), ext, content_type
    if fmt == "json":
        with open(os.path.join(path, "result.json"), "rb") as f:
            return f.read(), ext, content_type
    result = subprocess.run(
        [formatter_cmd, "-archive", archive, "-outfmt", outfmt],
        check=True, capture_output=True, timeout=timeout
    )
    return result.stdout, ext, content_type
