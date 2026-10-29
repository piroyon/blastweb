"""Running BLAST+ commands: databases, parameter validation and command building."""
from flask import current_app
import subprocess
import tempfile
import os
import shlex
import re

# program -> (query type, database type)
PROGRAMS = {
    'blastn': ('nucl', 'nucl'),
    'blastp': ('prot', 'prot'),
    'blastx': ('nucl', 'prot'),
    'tblastn': ('prot', 'nucl'),
    'tblastx': ('nucl', 'nucl'),
}
MATRIX_PROGRAMS = ['blastp', 'blastx', 'tblastn', 'tblastx']
MATRICES = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']
TASKS = {
    'blastn': ['megablast', 'dc-megablast', 'blastn', 'blastn-short'],
    'blastp': ['blastp', 'blastp-fast', 'blastp-short'],
    'blastx': ['blastx', 'blastx-fast'],
    'tblastn': ['tblastn', 'tblastn-fast'],
    'tblastx': [],
}
# Programs searching with a nucleotide query accept -strand
STRAND_PROGRAMS = ['blastn', 'blastx', 'tblastx']
DEFAULT_TIMEOUT = 600

# Options users may pass via "Additional Options".
# Anything touching files, output format, remote search or threads is excluded.
ALLOWED_EXTRA_OPTIONS = {
    '-task', '-word_size', '-gapopen', '-gapextend', '-reward', '-penalty',
    '-strand', '-dust', '-seg', '-soft_masking', '-lcase_masking',
    '-ungapped', '-comp_based_stats', '-threshold', '-perc_identity',
    '-qcov_hsp_perc', '-max_hsps', '-culling_limit', '-window_size',
    '-xdrop_ungap', '-xdrop_gap', '-xdrop_gap_final', '-searchsp',
    '-sum_stats', '-best_hit_overhang', '-best_hit_score_edge',
    '-min_raw_gapped_score', '-template_type', '-template_length',
    '-query_gencode', '-db_gencode', '-query_loc', '-matrix',
    '-max_intron_length', '-use_sw_tback', '-no_greedy', '-subject_besthit',
}


class BlastError(Exception):
    pass


def get_blast_command(program):
    blast_dir = current_app.config["blast_path"]
    if blast_dir:
        full_path = os.path.join(blast_dir, program)
        if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
            return full_path

    raise FileNotFoundError(
        f"BLAST+ command '{program}' not found\n"
        f"Specify blast_path: /your/blast/bin in 'blast.yaml.'"
    )


def list_blast_databases():
    """Return [{"name", "type" ('nucl'|'prot'), "title"}] for databases in blast_db."""
    db_dir = current_app.config["blast_db"]
    try:
        result = subprocess.run(
            [get_blast_command("blastdbcmd"), "-list", db_dir, "-list_outfmt", "%f\t%p\t%t"],
            check=True, capture_output=True, text=True, timeout=60
        )
    except (subprocess.SubprocessError, OSError):
        return _list_databases_from_files(db_dir)

    dbs = []
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        path, dbtype, title = parts[0], parts[1], parts[2]
        dbs.append({
            "name": os.path.basename(path),
            "type": "prot" if dbtype.lower().startswith("prot") else "nucl",
            "title": title,
        })
    return sorted(dbs, key=lambda d: d["name"])


def _list_databases_from_files(db_dir):
    """Fallback when blastdbcmd -list fails: guess databases from file names."""
    pattern = re.compile(r"^(.+?)(\.\d{2,})?\.([np])(in|sq|hr|al)$")
    dbs = {}
    for f in os.listdir(db_dir):
        match = pattern.match(f)
        if match:
            name = match.group(1)
            dbs[name] = {"name": name, "type": "prot" if match.group(3) == "p" else "nucl", "title": name}
    return sorted(dbs.values(), key=lambda d: d["name"])


def find_database(name):
    for db in list_blast_databases():
        if db["name"] == name:
            return db
    return None


def parse_extra_args(extra_args):
    """Split extra_args and reject options not in ALLOWED_EXTRA_OPTIONS.

    default_extra_args from blast.yaml is set by the admin and trusted as-is.
    """
    try:
        tokens = shlex.split(extra_args)
    except ValueError as e:
        raise BlastError(f"Failed to parse additional options: {e}")

    if extra_args.strip() == current_app.config.get("default_extra_args", "").strip():
        return tokens

    for token in tokens:
        # Negative numbers (e.g. "-penalty -3") are values, not options
        if re.match(r"^-[A-Za-z]", token) and token not in ALLOWED_EXTRA_OPTIONS:
            raise BlastError(f"Option not allowed: {token}")
    return tokens


def build_search_command(params):
    """Validate search parameters and return the BLAST command without -query/-out/-outfmt.

    params keys: program, database, evalue, max_target_seqs, matrix, task,
    word_size, filter ('', 'yes', 'no'), strand ('', 'plus', 'minus'), extra_args
    """
    program = params.get("program") or "blastn"
    if program not in PROGRAMS:
        raise BlastError(f"Unknown program: {program}")

    db = find_database(params.get("database") or "")
    if db is None:
        raise BlastError(f"Unknown database: {params.get('database')}")
    if db["type"] != PROGRAMS[program][1]:
        expected = "protein" if PROGRAMS[program][1] == "prot" else "nucleotide"
        raise BlastError(f"{program} needs a {expected} database, but '{db['name']}' is not")

    evalue = str(params.get("evalue") or "10")
    max_target_seqs = str(params.get("max_target_seqs") or "50")
    try:
        float(evalue)
        if int(max_target_seqs) < 1:
            raise ValueError
    except ValueError:
        raise BlastError("E-value and Max target seqs must be numbers")

    cmd = [
        get_blast_command(program),
        "-db", os.path.join(current_app.config["blast_db"], db["name"]),
        "-evalue", evalue,
        "-max_target_seqs", max_target_seqs,
    ]

    task = params.get("task") or ""
    if task:
        if task not in TASKS[program]:
            raise BlastError(f"Task {task} is not available for {program}")
        cmd += ["-task", task]

    matrix = params.get("matrix") or ""
    if matrix and program in MATRIX_PROGRAMS:
        if matrix not in MATRICES:
            raise BlastError(f"Unknown matrix: {matrix}")
        cmd += ["-matrix", matrix]

    word_size = str(params.get("word_size") or "")
    if word_size:
        if not word_size.isdigit():
            raise BlastError("Word size must be a number")
        cmd += ["-word_size", word_size]

    low_complexity = params.get("filter") or ""
    if low_complexity:
        if low_complexity not in ("yes", "no"):
            raise BlastError("Low complexity filter must be 'yes' or 'no'")
        cmd += ["-dust" if program == "blastn" else "-seg", low_complexity]

    strand = params.get("strand") or ""
    if strand and program in STRAND_PROGRAMS:
        if strand not in ("plus", "minus", "both"):
            raise BlastError(f"Unknown strand: {strand}")
        cmd += ["-strand", strand]

    cmd += parse_extra_args(params.get("extra_args") or "")
    return cmd


def build_query_fasta(sequence):
    """Use the input as-is if it is already FASTA (possibly multi-FASTA)."""
    sequence = sequence.strip()
    if sequence.startswith(">"):
        return sequence + "\n"
    return ">query\n" + sequence + "\n"


def run_tabular_search(sequence, params):
    """Run a search synchronously and return outfmt 6 rows (used by /api/blast)."""
    cmd = build_search_command(params)
    timeout = current_app.config.get("timeout", DEFAULT_TIMEOUT)

    with tempfile.TemporaryDirectory() as tmpdir:
        query_path = os.path.join(tmpdir, "query.fa")
        out_path = os.path.join(tmpdir, "result.out")
        with open(query_path, "w") as query_f:
            query_f.write(build_query_fasta(sequence))

        cmd += ["-query", query_path, "-out", out_path, "-outfmt", "6"]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
        except subprocess.CalledProcessError as e:
            raise BlastError(f"BLAST error: {e.stderr.strip() or e}")
        except subprocess.TimeoutExpired:
            raise BlastError(f"BLAST timed out after {timeout} seconds")

        with open(out_path) as out_f:
            return [line.rstrip("\n").split("\t") for line in out_f if line.strip()]


def fetch_subject(db_name, subject_id):
    db = find_database(db_name)
    if db is None:
        raise BlastError(f"Unknown database: {db_name}")
    try:
        result = subprocess.run(
            [get_blast_command("blastdbcmd"),
             "-db", os.path.join(current_app.config["blast_db"], db["name"]),
             "-entry", subject_id],
            check=True, capture_output=True, text=True, timeout=60
        )
    except subprocess.CalledProcessError as e:
        raise BlastError(f"Failed to fetch entry '{subject_id}' from database '{db_name}': {e.stderr}")
    except subprocess.TimeoutExpired:
        raise BlastError(f"Timed out fetching entry '{subject_id}'")
    return result.stdout
