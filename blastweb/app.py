from flask import request, render_template, jsonify, current_app, make_response
import subprocess
import tempfile
import os
import yaml
import shlex
import re

PROGRAMS = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
MATRIX_PROGRAMS = ['blastp', 'blastx', 'tblastn', 'tblastx']
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


def register_routes(app):

    DEFAULT_PROGRAM = 'blastn'

    @app.route("/", methods=["GET", "POST"])
    def index():
        default_extra_args = app.config.get("default_extra_args", "")
        db_dir = app.config.get("blast_db")
        url_prefix = app.config.get("url_prefix", "")
        available_dbs = list_blast_databases(db_dir)

        def form_error(message):
            return render_template("index.html", error=message,
                                   default_extra_args=default_extra_args, db_choices=available_dbs)

        if request.method == "POST":
            sequence = request.form.get("sequence", "").strip()
            program = request.form.get("program", DEFAULT_PROGRAM)
            db = request.form.get("database", '')
            evalue = request.form.get("evalue", "1e-5")
            max_target_seqs = request.form.get("max_target_seqs", "50")
            matrix = request.form.get("matrix", "")
            extra_args = request.form.get("extra_args", "")

            if not db:
                return form_error("No database")

            if not sequence:
                return form_error("No query sequence")

            result_lines, error = run_blast(sequence, program, db, evalue, max_target_seqs, matrix, extra_args)
            if error:
                return form_error(error)
            result_text = "\n".join("\t".join(row) for row in result_lines)

            return render_template("result.html", r_text=result_text, results=result_lines, selected_db=db, url_prefix=url_prefix)

        return render_template("index.html", default_extra_args=default_extra_args, db_choices=available_dbs)


    @app.route("/download", methods=["POST"])
    def download():
        content = request.form.get("result_text", "")
        response = make_response(content)
        response.headers["Content-Disposition"] = "attachment; filename=blast_result.txt"
        response.headers["Content-Type"] = "text/plain"
        return response



    @app.route("/subject/<db>/<subject_id>")
    def get_subject_sequence(db, subject_id):
        text_headers = {"Content-Type": "text/plain; charset=utf-8"}
        if db not in list_blast_databases(app.config.get("blast_db")):
            return f"Unknown database: {db}\n", 404, text_headers

        cmd_path = get_blast_command("blastdbcmd")
        db_path = os.path.join(app.config.get("blast_db"), db)

        try:
            result = subprocess.run(
                [cmd_path, "-db", db_path, "-entry", subject_id],
                check=True, capture_output=True, text=True, timeout=60
            )
        except subprocess.CalledProcessError as e:
            return f"Failed to fetch entry '{subject_id}' from database '{db}': {e.stderr}", 500, text_headers
        except subprocess.TimeoutExpired:
            return f"Timed out fetching entry '{subject_id}'\n", 504, text_headers

        return result.stdout, 200, text_headers


    @app.route("/api/blast", methods=["POST"])
    def api_blast():
        data = request.get_json(silent=True) or {}
        sequence = data.get("sequence", "").strip()
        program = data.get("program", DEFAULT_PROGRAM)
        db = data.get("database", '')
        evalue = data.get("evalue", "1e-3")
        max_target_seqs = data.get("max_target_seqs", "50")
        matrix = data.get("matrix", "")
        extra_args = data.get("extra_args", current_app.config.get("default_extra_args", ""))

        if not sequence:
            return jsonify({"error": "No sequence provided"}), 400

        if not db:
            return jsonify({"error": "No database"}), 400

        result_lines, error = run_blast(sequence, program, db, evalue, max_target_seqs, matrix, extra_args)
        if error:
            return jsonify({"error": error}), 400

        return jsonify({"results": result_lines})

    return app


def build_query_fasta(sequence):
    """Use the input as-is if it is already FASTA (possibly multi-FASTA)."""
    if sequence.startswith(">"):
        return sequence + "\n"
    return ">query\n" + sequence + "\n"


def parse_extra_args(extra_args):
    """Split extra_args and reject options not in ALLOWED_EXTRA_OPTIONS.

    default_extra_args from blast.yaml is set by the admin and trusted as-is.
    """
    try:
        tokens = shlex.split(extra_args)
    except ValueError as e:
        return None, f"Failed to parse additional options: {e}"

    if extra_args.strip() == current_app.config.get("default_extra_args", "").strip():
        return tokens, None

    for token in tokens:
        # Negative numbers (e.g. "-penalty -3") are values, not options
        if re.match(r"^-[A-Za-z]", token) and token not in ALLOWED_EXTRA_OPTIONS:
            return None, f"Option not allowed: {token}"
    return tokens, None


def run_blast(sequence, program, db_name, evalue="1e-5", max_target_seqs="50", matrix="", extra_args=""):
    if program not in PROGRAMS:
        return [], f"Unknown program: {program}"

    db_dir = current_app.config["blast_db"]
    if db_name not in list_blast_databases(db_dir):
        return [], f"Unknown database: {db_name}"

    try:
        float(evalue)
        int(max_target_seqs)
    except (TypeError, ValueError):
        return [], "E-value and Max target seqs must be numbers"

    extra_tokens, error = parse_extra_args(extra_args or "")
    if error:
        return [], error

    blast_path = get_blast_command(program)
    db_path = os.path.join(db_dir, db_name)
    timeout = current_app.config.get("timeout", DEFAULT_TIMEOUT)

    with tempfile.TemporaryDirectory() as tmpdir:
        query_path = os.path.join(tmpdir, "query.fa")
        out_path = os.path.join(tmpdir, "result.out")
        with open(query_path, "w") as query_f:
            query_f.write(build_query_fasta(sequence))

        cmd = [
            blast_path,
            "-query", query_path,
            "-db", db_path,
            "-out", out_path,
            "-outfmt", "6",
            "-evalue", str(evalue),
            "-max_target_seqs", str(max_target_seqs),
        ]

        if program in MATRIX_PROGRAMS and matrix:
            cmd += ["-matrix", matrix]

        cmd += extra_tokens

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
        except subprocess.CalledProcessError as e:
            return [], f"BLAST error: {e.stderr.strip() or e}"
        except subprocess.TimeoutExpired:
            return [], f"BLAST timed out after {timeout} seconds"

        with open(out_path) as out_f:
            return [line.rstrip("\n").split("\t") for line in out_f if line.strip()], None


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

    return config


def get_blast_command(program, config=None):
    blast_dir = current_app.config["blast_path"]
    if blast_dir:
        full_path = os.path.join(blast_dir, program)
        if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
            return full_path

    raise FileNotFoundError(
        f"BLAST+ command '{program}' not found\n"
        f"Specify blast_path: /your/blast/bin in 'blast.yaml.'"
    )

def remove_trailing_digits(s):
    return re.sub(r'\.\d{2}$', '', s)

def list_blast_databases(db_dir):
    files = os.listdir(db_dir)
    pattern = re.compile(r"^(.+)\.(nin|nsq|nhr|pin|psq|phr)$")

    db_names = set()
    for f in files:
        match = pattern.match(f)
        if match:
            db_names.add(remove_trailing_digits(match.group(1)))

    return sorted(db_names)
