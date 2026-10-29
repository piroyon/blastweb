import os
import random
import shutil
import subprocess
import time

import pytest

from blastweb import create_app

BLASTN = shutil.which("blastn")

random.seed(1)
NUCL_SEQS = {f"seq{i}": "".join(random.choice("ACGT") for _ in range(500)) for i in range(3)}
PROT_SEQS = {f"prot{i}": "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(300)) for i in range(3)}


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _make_db(directory, name, seqs, dbtype):
    fasta = os.path.join(directory, f"{name}.fa")
    with open(fasta, "w") as f:
        for seq_id, seq in seqs.items():
            f.write(f">{seq_id} test {seq_id}\n{seq}\n")
    subprocess.run(["makeblastdb", "-in", fasta, "-dbtype", dbtype, "-parse_seqids",
                    "-out", os.path.join(directory, name), "-title", f"{name} title"],
                   check=True, capture_output=True)


@pytest.fixture(scope="session")
def blast_dbs(tmp_path_factory):
    if BLASTN is None:
        pytest.skip("BLAST+ is not installed")
    db_dir = tmp_path_factory.mktemp("db")
    _make_db(str(db_dir), "nucldb", NUCL_SEQS, "nucl")
    _make_db(str(db_dir), "protdb", PROT_SEQS, "prot")
    return db_dir


@pytest.fixture
def app(blast_dbs, tmp_path):
    config = tmp_path / "blast.yaml"
    config.write_text(
        f"blast_path: {os.path.dirname(BLASTN)}\n"
        f"blast_db: {blast_dbs}\n"
        f"job_dir: {tmp_path / 'jobs'}\n"
        "default_extra_args: ''\n"
    )
    return create_app(str(config), use_dispatcher=False)


@pytest.fixture
def client(app):
    return app.test_client()


def wait_for_job(client, job_id, timeout=60):
    """Poll the API until the job leaves the running state."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = client.get(f"/api/jobs/{job_id}").get_json()
        if body["state"] != "running":
            return body
        time.sleep(0.2)
    raise AssertionError("job did not finish")
