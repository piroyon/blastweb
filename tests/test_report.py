"""Unit tests for alignment coordinates and summary values (no BLAST+ needed)."""
from blastweb.report import build_report, _build_hsp, _query_cover, _ticks


def hsp(**kw):
    base = {"bit_score": 100.0, "score": 50, "evalue": 1e-20, "identity": 10, "align_len": 10, "gaps": 0}
    base.update(kw)
    return base


def coords(alignment):
    """[(start, end)] for Query and Sbjct lines."""
    result = []
    for line in alignment.splitlines():
        if line.startswith(("Query", "Sbjct")):
            parts = line.split()
            result.append((parts[0], int(parts[1]), int(parts[3])))
    return result


def test_plus_strand_blocks():
    seq = "A" * 130
    h = _build_hsp(hsp(qseq=seq, hseq=seq, midline="|" * 130, align_len=130, identity=130,
                       query_from=1, query_to=130, hit_from=101, hit_to=230,
                       query_strand="Plus", hit_strand="Plus"))
    assert coords(h["alignment"]) == [
        ("Query", 1, 60), ("Sbjct", 101, 160),
        ("Query", 61, 120), ("Sbjct", 161, 220),
        ("Query", 121, 130), ("Sbjct", 221, 230),
    ]
    assert h["orientation"] == "Strand = Plus/Plus"


def test_minus_strand_hit_counts_down():
    seq = "A" * 70
    h = _build_hsp(hsp(qseq=seq, hseq=seq, midline="|" * 70, align_len=70, identity=70,
                       query_from=1, query_to=70, hit_from=300, hit_to=231,
                       query_strand="Plus", hit_strand="Minus"))
    assert coords(h["alignment"]) == [("Query", 1, 60), ("Sbjct", 300, 241), ("Query", 61, 70), ("Sbjct", 240, 231)]


def test_gaps_do_not_advance_coordinates():
    h = _build_hsp(hsp(qseq="AC--GT", hseq="ACTTGT", midline="||  ||", align_len=6, identity=4, gaps=2,
                       query_from=1, query_to=4, hit_from=1, hit_to=6))
    assert coords(h["alignment"]) == [("Query", 1, 4), ("Sbjct", 1, 6)]


def test_translated_negative_frame():
    seq = "M" * 70
    h = _build_hsp(hsp(qseq=seq, hseq=seq, midline=seq, align_len=70, identity=70, positive=70,
                       query_from=1, query_to=210, query_frame=-1, hit_from=21, hit_to=90))
    assert coords(h["alignment"]) == [("Query", 210, 31), ("Sbjct", 21, 80), ("Query", 30, 1), ("Sbjct", 81, 90)]
    assert h["orientation"] == "Frame = -1"
    assert h["positive_pct"] == 100.0


def test_query_cover_merges_overlaps():
    hsps = [{"query_from": 1, "query_to": 50}, {"query_from": 40, "query_to": 60}, {"query_from": 81, "query_to": 100}]
    assert _query_cover(hsps, 100) == 80.0


def test_ticks():
    assert _ticks(500) == [1, 100, 200, 300, 400, 500]
    assert _ticks(30) == [1, 5, 10, 15, 20, 25, 30]


def test_build_report_no_hits():
    data = {"BlastOutput2": [{"report": {"results": {"search": {
        "query_id": "Query_1", "query_title": "q1", "query_len": 10, "message": "No hits found"}}}}]}
    [q] = build_report(data)
    assert q["title"] == "q1" and q["hits"] == [] and q["message"] == "No hits found"
