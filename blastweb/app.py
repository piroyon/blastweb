from flask import request, render_template, jsonify, current_app, make_response, redirect, url_for, abort
import os
import tempfile
import time
import yaml

from . import jobs
from .blast import (PROGRAMS, MATRIX_PROGRAMS, MATRICES, TASKS, STRAND_PROGRAMS, DEFAULT_TIMEOUT,
                    BlastError, list_blast_databases, build_search_command, build_query_fasta,
                    run_tabular_search, fetch_subject, get_blast_command)
from .report import build_report, score_legend

PARAM_KEYS = ["program", "database", "evalue", "max_target_seqs", "matrix", "task",
              "word_size", "filter", "strand", "extra_args"]


def register_routes(app):

    def form_context(**kwargs):
        return dict(
            default_extra_args=app.config.get("default_extra_args", ""),
            db_choices=list_blast_databases(),
            programs=PROGRAMS, matrix_programs=MATRIX_PROGRAMS, matrices=MATRICES,
            tasks=TASKS, strand_programs=STRAND_PROGRAMS,
            **kwargs
        )

    def submit_job(sequence, params):
        """Validate and start a background job. Returns the job ID or raises BlastError."""
        search_cmd = build_search_command(params)
        job_dir = app.config["job_dir"]
        os.makedirs(job_dir, exist_ok=True)
        jobs.cleanup_old_jobs(job_dir, app.config["job_retention_days"])
        if jobs.count_running(job_dir) >= app.config["max_running_jobs"]:
            raise BlastError("The server is busy. Please try again later.")
        return jobs.start_job(job_dir, build_query_fasta(sequence), params, search_cmd,
                              get_blast_command("blast_formatter"), app.config["timeout"])

    def find_job(job_id):
        path = jobs.job_path(app.config["job_dir"], job_id)
        if path is None:
            abort(404)
        return path

    @app.route("/", methods=["GET", "POST"])
    def index():
        if request.method == "POST":
            sequence = request.form.get("sequence", "").strip()
            upload = request.files.get("query_file")
            if upload and upload.filename:
                sequence = upload.read().decode("utf-8", errors="replace").strip()
            params = {key: request.form.get(key, "") for key in PARAM_KEYS}

            if not params["database"]:
                return render_template("index.html", **form_context(error="No database"))
            if not sequence:
                return render_template("index.html", **form_context(error="No query sequence"))

            try:
                job_id = submit_job(sequence, params)
            except BlastError as e:
                return render_template("index.html", **form_context(error=str(e)))
            return redirect(url_for("job", job_id=job_id))

        return render_template("index.html", **form_context())


    @app.route("/job/<job_id>")
    def job(job_id):
        path = find_job(job_id)
        status = jobs.get_status(path)
        params = jobs.get_params(path)
        if status["state"] != "done":
            return render_template("job.html", job_id=job_id, status=status, params=params,
                                   elapsed=int((status.get("finished") or time.time()) - status["started"]))

        queries = build_report(jobs.get_result(path))
        return render_template("result.html", job_id=job_id, params=params, queries=queries,
                               legend=score_legend(), download_formats=jobs.DOWNLOAD_FORMATS)


    @app.route("/job/<job_id>/download/<fmt>")
    def download(job_id, fmt):
        path = find_job(job_id)
        if fmt not in jobs.DOWNLOAD_FORMATS or jobs.get_status(path)["state"] != "done":
            abort(404)
        content, ext, content_type = jobs.format_result(path, fmt, get_blast_command("blast_formatter"),
                                                        app.config["timeout"])
        response = make_response(content)
        response.headers["Content-Disposition"] = f"attachment; filename=blast_{job_id[:8]}.{ext}"
        response.headers["Content-Type"] = content_type
        return response


    @app.route("/subject/<db>/<path:subject_id>")
    def get_subject_sequence(db, subject_id):
        text_headers = {"Content-Type": "text/plain; charset=utf-8"}
        try:
            return fetch_subject(db, subject_id), 200, text_headers
        except BlastError as e:
            return f"{e}\n", 404, text_headers


    def api_params(data):
        params = {key: data.get(key, "") for key in PARAM_KEYS}
        params["program"] = params["program"] or "blastn"
        params["evalue"] = params["evalue"] or "1e-3"
        if "extra_args" not in data:
            params["extra_args"] = current_app.config.get("default_extra_args", "")
        return params

    @app.route("/api/blast", methods=["POST"])
    def api_blast():
        """Run a search synchronously and return tabular (outfmt 6) rows."""
        data = request.get_json(silent=True) or {}
        sequence = data.get("sequence", "").strip()
        if not sequence:
            return jsonify({"error": "No sequence provided"}), 400
        if not data.get("database"):
            return jsonify({"error": "No database"}), 400

        try:
            result_lines = run_tabular_search(sequence, api_params(data))
        except BlastError as e:
            return jsonify({"error": str(e)}), 400
        return jsonify({"results": result_lines})


    @app.route("/api/jobs", methods=["POST"])
    def api_submit_job():
        """Start a background search. Poll the returned status_url until state is done."""
        data = request.get_json(silent=True) or {}
        sequence = data.get("sequence", "").strip()
        if not sequence:
            return jsonify({"error": "No sequence provided"}), 400
        if not data.get("database"):
            return jsonify({"error": "No database"}), 400

        try:
            job_id = submit_job(sequence, api_params(data))
        except BlastError as e:
            return jsonify({"error": str(e)}), 400
        return jsonify({
            "job_id": job_id,
            "status_url": url_for("api_job", job_id=job_id, _external=True),
            "result_url": url_for("job", job_id=job_id, _external=True),
        }), 202


    @app.route("/api/jobs/<job_id>")
    def api_job(job_id):
        """Job state; when done, includes the BLAST JSON (outfmt 15) as blast_output."""
        path = find_job(job_id)
        status = jobs.get_status(path)
        body = {"job_id": job_id, "state": status["state"], "params": jobs.get_params(path)}
        if status["state"] == "error":
            body["error"] = status["error"]
        elif status["state"] == "done":
            body["blast_output"] = jobs.get_result(path)
        return jsonify(body)

    return app


def load_config(config_path=None):
    if config_path:
        path = os.path.abspath(config_path)
    else:
        path = os.path.join(os.getcwd(), "blast.yaml")

    if not os.path.exists(path):
        raise FileNotFoundError(f"blast.yaml not found: {path}")

    with open(path) as f:
        config = yaml.safe_load(f) or {}

    if "blast_db" not in config:
        raise ValueError("Need to set 'blast_db' in blast.yaml")
    if not os.path.isdir(config["blast_db"]):
        raise ValueError(f"blast_db directory not found: {config['blast_db']}")
    if not config.get("blast_path"):
        config["blast_path"] = ""
    if not config.get("default_extra_args"):
        config["default_extra_args"] = ""
    config.setdefault("timeout", DEFAULT_TIMEOUT)
    config["job_dir"] = config.get("job_dir") or os.path.join(tempfile.gettempdir(), "blastweb_jobs")
    config.setdefault("job_retention_days", 7)
    config.setdefault("max_running_jobs", 4)
    config["MAX_CONTENT_LENGTH"] = int(config.get("max_upload_mb", 10)) * 1024 * 1024

    return config
