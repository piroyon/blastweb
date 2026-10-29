import os
import time

from blastweb import jobs
from conftest import NUCL_SEQS, PROT_SEQS, revcomp, wait_for_job

QUERY = NUCL_SEQS["seq0"][50:250]


def submit(client, **form):
    data = {"sequence": QUERY, "program": "blastn", "database": "nucldb", "evalue": "1e-5",
            "max_target_seqs": "50", "extra_args": ""}
    data.update(form)
    return client.post("/", data=data)


def test_index_lists_databases_with_types(client):
    r = client.get("/")
    assert r.status_code == 200
    assert b'value="nucldb" data-type="nucl"' in r.data
    assert b'value="protdb" data-type="prot"' in r.data
    assert b"nucldb title" in r.data


def test_search_flow_shows_overview_table_and_alignments(client):
    r = submit(client, sequence=f">plus\n{QUERY}\n>minus\n{revcomp(QUERY)}")
    assert r.status_code == 302
    job_id = r.headers["Location"].rsplit("/", 1)[-1]
    assert wait_for_job(client, job_id)["state"] == "done"

    page = client.get(f"/job/{job_id}").data.decode()
    assert "<svg" in page and "Graphical overview" in page
    assert 'href="#q0-h0"' in page and 'id="q1-h0"' in page
    assert "Query: plus" in page and "Query: minus" in page
    assert "Strand = Plus/Minus" in page
    assert "Query  1" in page


def test_downloads(client):
    job_id = submit(client).headers["Location"].rsplit("/", 1)[-1]
    wait_for_job(client, job_id)
    for fmt in jobs.DOWNLOAD_FORMATS:
        r = client.get(f"/job/{job_id}/download/{fmt}")
        assert r.status_code == 200, fmt
        assert r.data, fmt
    assert client.get(f"/job/{job_id}/download/tabular").data.startswith(b"query\tseq0")
    assert client.get(f"/job/{job_id}/download/nope").status_code == 404


def test_blastp_against_protein_db(client):
    r = submit(client, program="blastp", database="protdb", sequence=PROT_SEQS["prot1"][10:110], matrix="BLOSUM62")
    job_id = r.headers["Location"].rsplit("/", 1)[-1]
    assert wait_for_job(client, job_id)["state"] == "done"
    page = client.get(f"/job/{job_id}").data.decode()
    assert "prot1" in page and "Positives" in page


def test_fasta_file_upload(client):
    from io import BytesIO
    r = submit(client, sequence="", query_file=(BytesIO(f">up\n{QUERY}\n".encode()), "q.fa"))
    job_id = r.headers["Location"].rsplit("/", 1)[-1]
    assert wait_for_job(client, job_id)["state"] == "done"
    assert b"Query: up" in client.get(f"/job/{job_id}").data


def test_form_errors_keep_database_choices(client):
    for form, message in [
        ({"database": "../nucldb"}, b"Unknown database"),
        ({"database": "protdb"}, b"needs a nucleotide database"),
        ({"extra_args": "-out /tmp/x"}, b"Option not allowed: -out"),
        ({"evalue": "abc"}, b"must be numbers"),
        ({"task": "blastp"}, b"not available"),
        ({"sequence": ""}, b"No query sequence"),
    ]:
        r = submit(client, **form)
        assert r.status_code == 200, form
        assert message in r.data, form
        assert b'value="nucldb"' in r.data


def test_blast_error_is_shown_on_job_page(client):
    job_id = submit(client, extra_args="-word_size 1").headers["Location"].rsplit("/", 1)[-1]
    body = wait_for_job(client, job_id)
    assert body["state"] == "error" and "BLAST error" in body["error"]
    assert b"Search failed" in client.get(f"/job/{job_id}").data


def test_unknown_job_is_404(client):
    assert client.get("/job/0123456789abcdef0123456789abcdef").status_code == 404
    assert client.get("/job/..%2F..%2Fetc").status_code == 404


def test_api_blast_sync(client):
    r = client.post("/api/blast", json={"sequence": QUERY, "database": "nucldb"})
    assert r.status_code == 200
    assert r.get_json()["results"][0][:2] == ["query", "seq0"]
    assert client.post("/api/blast", json={"sequence": QUERY, "database": "x"}).status_code == 400
    assert client.post("/api/blast", data="not json").status_code == 400


def test_api_jobs(client):
    r = client.post("/api/jobs", json={"sequence": QUERY, "database": "nucldb"})
    assert r.status_code == 202
    body = wait_for_job(client, r.get_json()["job_id"])
    assert body["state"] == "done"
    hits = body["blast_output"]["BlastOutput2"][0]["report"]["results"]["search"]["hits"]
    assert hits[0]["description"][0]["id"] == "seq0"


def test_subject_is_plain_text(client):
    r = client.get("/subject/nucldb/seq1")
    assert r.content_type.startswith("text/plain")
    assert r.data.startswith(b">seq1")
    assert client.get("/subject/evil/seq1").status_code == 404


def test_old_jobs_are_removed(app, client):
    job_id = submit(client).headers["Location"].rsplit("/", 1)[-1]
    wait_for_job(client, job_id)
    path = os.path.join(app.config["job_dir"], job_id)
    old = time.time() - 8 * 86400
    os.utime(path, (old, old))
    submit(client)
    assert not os.path.exists(path)
