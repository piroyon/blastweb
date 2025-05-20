from flask import Flask, request, render_template, jsonify
import subprocess
import tempfile
import os
import logging
import yaml
import shlex
import re

app = Flask(__name__)

DEFAULT_DB = "/data/blastdb"
DEFAULT_PROGRAM = "blastn"

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

print("現在地:{}".format(os.getcwd()))


@app.route("/", methods=["GET", "POST"])
def index():
    #return "<h1>Hello from Flask</h1>"
    config = load_config()
    default_extra_args = config.get("default_extra_args", "")
    db_dir = config["blast_db"]
    available_dbs = list_blast_databases(db_dir)

    if request.method == "POST":
        sequence = request.form.get("sequence", "").strip()
        program = request.form.get("program", DEFAULT_PROGRAM)
        db = request.form.get("database", DEFAULT_DB)
        evalue = request.form.get("evalue", "1e-5")
        max_target_seqs = request.form.get("max_target_seqs", "50")
        matrix = request.form.get("matrix", "")
        extra_args = request.form.get("extra_args", "")

        if not sequence:
            return render_template("index.html", error="配列が入力されていません")

        result_lines, error = run_blast(sequence, program, db, evalue, max_target_seqs, matrix, extra_args)
        if error:
            return render_template("index.html", error=error)

        return render_template("result.html", results=result_lines)

    return render_template("index.html", default_extra_args=default_extra_args, db_choices=available_dbs)


@app.route("/api/blast", methods=["POST"])
def api_blast():
    data = request.get_json()
    sequence = data.get("sequence", "").strip()
    program = data.get("program", DEFAULT_PROGRAM)
    db = data.get("database", DEFAULT_DB)
    evalue = data.get("evalue", "1e-5")
    max_target_seqs = data.get("max_target_seqs", "50")
    matrix = data.get("matrix", "")
    extra_args = data.get("extra_args", "")
    config = data.get("config", "")

    if not sequence:
        return jsonify({"error": "No sequence provided"}), 400

    result_lines, error = run_blast(sequence, program, db, evalue, max_target_seqs, matrix, extra_args, config)
    if error:
        return jsonify({"error": error}), 500

    return jsonify({"results": result_lines})


import shlex

def run_blast(sequence, program, db_name, evalue="1e-5", max_target_seqs="50", matrix="", extra_args="", config=None):
    config = config or load_config()
    blast_path = get_blast_command(program, config)
    db_dir = config["blast_db"]
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
                    return [], f"追加オプションの解析に失敗しました: {e}"
            else:
                extra_args = config.get("default_extra_args", "")

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
        module_dir = os.path.dirname(__file__)
        path = os.path.join(module_dir, "blast.yaml")

    if not os.path.exists(path):
        raise FileNotFoundError(f"blast.yaml が見つかりません: {path}")

    with open(path) as f:
        config = yaml.safe_load(f)

    if "blast_db" not in config:
        raise ValueError("blast.yaml に 'blast_db' の設定が必要です。")
    if "blast_path" not in config:
        config["blast_path"] = ""
    if "default_extra_args" not in config:
        config["default_extra_args"] = ""

    return config


def get_blast_command(program, config=None):
    config = config or load_config()
    blast_dir = config.get("blast_path", "")
    if blast_dir:
        full_path = os.path.join(blast_dir, program)
        if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
            return full_path

    raise FileNotFoundError(
        f"BLASTコマンド '{program}' が見つかりません。\n"
        f"blast.yaml に blast_path: /your/blast/bin を指定してください。"
    )

def list_blast_databases(db_dir):
    files = os.listdir(db_dir)
    # よくあるBLAST DBの拡張子（nucleotide/protein共通）
    exts = ('.nin', '.nsq', '.nhr', '.pin', '.psq', '.phr')
    pattern = re.compile(r"^(.+)\.(nin|nsq|nhr|pin|psq|phr)$")

    db_names = set()
    for f in files:
        match = pattern.match(f)
        if match:
            db_names.add(match.group(1))

    return sorted(db_names)
