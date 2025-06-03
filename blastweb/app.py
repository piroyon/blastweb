from flask import request, render_template, jsonify, current_app, make_response
import subprocess
import tempfile
import os
import logging
import yaml
import shlex
import re

def register_routes(app):

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    DEFAULT_PROGRAM = 'blastn'

    @app.route("/", methods=["GET", "POST"])
    def index():
        #return "<h1>Hello from Flask</h1>"
        default_extra_args = app.config.get("default_extra_args", "")
        db_dir = app.config.get("blast_db")
        url_prefix = app.config.get("url_prefix", "")
        available_dbs = list_blast_databases(db_dir)

        if request.method == "POST":
            sequence = request.form.get("sequence", "").strip()
            program = request.form.get("program", DEFAULT_PROGRAM)
            db = request.form.get("database", '')
            evalue = request.form.get("evalue", "1e-5")
            max_target_seqs = request.form.get("max_target_seqs", "50")
            matrix = request.form.get("matrix", "")
            extra_args = request.form.get("extra_args", "")

            if not db:
                return render_template("index.html", error="No database")

            if not sequence:
                return render_template("index.html", error="No query sequence")

            result_lines, error = run_blast(sequence, program, db, evalue, max_target_seqs, matrix, extra_args)
            result_text = "\n".join("\t".join(row) for row in result_lines)
            if error:
                return render_template("index.html", error=error)

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
        program = "blastdbcmd"
        cmd_path = get_blast_command(program)
        db_dir = app.config.get("blast_db")
        db_path = os.path.join(db_dir, db)

        try:
            result = subprocess.run(
                [cmd_path, "-db", db_path, "-entry", subject_id],
                check=True, capture_output=True, text=True
            )
            fasta = result.stdout
        except subprocess.CalledProcessError as e:
            return f"<p>Failed to fetch entry '{subject_id}' from database '{db}': {e.stderr}</p>", 500

        return f"<pre>{fasta}</pre>"


    @app.route("/api/blast", methods=["POST"])
    def api_blast():
        data = request.get_json()
        sequence = data.get("sequence", "").strip()
        program = data.get("program", DEFAULT_PROGRAM)
        db = data.get("database", '')
        evalue = data.get("evalue", "1e-3")
        max_target_seqs = data.get("max_target_seqs", "50")
        matrix = data.get("matrix", "")
        extra_args = data.get("extra_args", "")
        blast_config = current_app.config["BLAST_CONFIG"]

        if not sequence:
            return jsonify({"error": "No sequence provided"}), 400

        if not db:
            return jsonify({"error": "No database"}), 400

        result_lines, error = run_blast(sequence, program, db, evalue, max_target_seqs, matrix, extra_args)
        if error:
            return jsonify({"error": error}), 500

        return jsonify({"results": result_lines})

    return app


def run_blast(sequence, program, db_name, evalue="1e-5", max_target_seqs="50", matrix="", extra_args=""):
    blast_path = get_blast_command(program)
    db_dir = current_app.config["blast_db"]
    db_path = os.path.join(db_dir, db_name)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as query_f:
        query_f.write(">query\n" + sequence)
        query_f.flush()

        with tempfile.NamedTemporaryFile(mode="r", suffix=".out", delete=False) as out_f:
            cmd = [
                blast_path,
                "-query", query_f.name,
                "-db", db_path,
                "-out", out_f.name,
                "-outfmt", "6",
                "-evalue", evalue,
                "-max_target_seqs", max_target_seqs,
            ]

            if program in ["blastp", "blastx"] and matrix:
                cmd += ["-matrix", matrix]

            if extra_args:
                try:
                    cmd += shlex.split(extra_args)
                except Exception as e:
                    return [], f"Failed to parse additional options.: {e}"
            else:
                extra_args = current_app.config["default_extra_args"]

            try:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                lines = [line.strip().split("\t") for line in out_f]
                return lines, None
            except subprocess.CalledProcessError as e:
                return [], f"BLAST error: {e}"
            finally:
                os.unlink(query_f.name)
                os.unlink(out_f.name)


def load_config(config_path=None):
    if config_path:
        path = os.path.abspath(config_path)
    else:
        path = os.path.join(os.getcwd(), "blast.yaml")

    if not os.path.exists(path):
        raise FileNotFoundError(f"blast.yaml not found: {path}")

    with open(path) as f:
        config = yaml.safe_load(f)

    if "blast_db" not in config:
        raise ValueError("Need to set 'blast_db' in blast.yaml")
    if "blast_path" not in config:
        config["blast_path"] = ""
    if "default_extra_args" not in config:
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
    exts = ('.nin', '.nsq', '.nhr', '.pin', '.psq', '.phr')
    pattern = re.compile(r"^(.+)\.(nin|nsq|nhr|pin|psq|phr)$")

    db_names = set()
    for f in files:
        match = pattern.match(f)
        if match:
            db_names.add(remove_trailing_digits(match.group(1)))

    return sorted(db_names)
