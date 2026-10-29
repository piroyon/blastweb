"""Turn BLAST outfmt 15 JSON into data for the result page:
hit table, graphical overview and pairwise alignments (wwwblast style)."""

ALIGNMENT_WIDTH = 60
OVERVIEW_MAX_HITS = 100
OVERVIEW_WIDTH = 700

# Bit score classes used by NCBI BLAST for the graphical overview
SCORE_CLASSES = [
    (40, "#000000", "< 40"),
    (50, "#0000ff", "40 - 50"),
    (80, "#00a000", "50 - 80"),
    (200, "#d000d0", "80 - 200"),
    (float("inf"), "#e00000", ">= 200"),
]


def score_color(bit_score):
    for limit, color, _ in SCORE_CLASSES:
        if bit_score < limit:
            return color
    return SCORE_CLASSES[-1][1]


def score_legend():
    return [(color, label) for _, color, label in SCORE_CLASSES]


def build_report(data):
    """Return a list of queries, each with its hits, overview and alignments."""
    queries = []
    for qi, entry in enumerate(data.get("BlastOutput2", [])):
        report = entry["report"]
        search = report["results"]["search"]
        query_len = search["query_len"]
        hits = []
        for hi, hit in enumerate(search.get("hits", [])):
            hits.append(_build_hit(qi, hi, hit, query_len))
        queries.append({
            "index": qi,
            "anchor": f"q{qi}",
            "title": search.get("query_title") or search.get("query_id", ""),
            "len": query_len,
            "message": search.get("message", ""),
            "hits": hits,
            "overview": _build_overview(hits, query_len),
        })
    return queries


def _build_hit(qi, hi, hit, query_len):
    desc = hit["description"][0]
    hsps = hit["hsps"]
    best = max(hsps, key=lambda h: h["bit_score"])
    return {
        "anchor": f"q{qi}-h{hi}",
        "id": desc.get("id", ""),
        "accession": desc.get("accession", desc.get("id", "")),
        "title": desc.get("title", ""),
        "len": hit["len"],
        "max_score": best["bit_score"],
        "total_score": sum(h["bit_score"] for h in hsps),
        "query_cover": _query_cover(hsps, query_len),
        "evalue": min(h["evalue"] for h in hsps),
        "identity": 100.0 * best["identity"] / best["align_len"],
        "hsps": [_build_hsp(h) for h in hsps],
    }


def _query_cover(hsps, query_len):
    """Percent of the query covered by the union of HSP ranges."""
    ranges = sorted((min(h["query_from"], h["query_to"]), max(h["query_from"], h["query_to"])) for h in hsps)
    covered = 0
    cur_start, cur_end = ranges[0]
    for start, end in ranges[1:]:
        if start > cur_end + 1:
            covered += cur_end - cur_start + 1
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    covered += cur_end - cur_start + 1
    return 100.0 * covered / query_len if query_len else 0.0


def _build_overview(hits, query_len):
    """Bars along the query for each hit, positioned in OVERVIEW_WIDTH units."""
    scale = OVERVIEW_WIDTH / query_len if query_len else 0
    rows = []
    for hit in hits[:OVERVIEW_MAX_HITS]:
        segments = []
        for hsp in hit["hsps"]:
            start = min(hsp["query_from"], hsp["query_to"])
            end = max(hsp["query_from"], hsp["query_to"])
            segments.append({
                "x": (start - 1) * scale,
                "width": max((end - start + 1) * scale, 1),
                "color": score_color(hsp["bit_score"]),
            })
        rows.append({"anchor": hit["anchor"], "label": f"{hit['id']} {hit['title']}  score={hit['max_score']:.1f}  E={hit['evalue']:.2g}",
                     "segments": segments})
    return {
        "width": OVERVIEW_WIDTH,
        "ticks": [{"x": (pos - 1) * scale, "label": pos} for pos in _ticks(query_len)],
        "rows": rows,
        "truncated": len(hits) > OVERVIEW_MAX_HITS,
    }


def _ticks(length):
    """At most 8 evenly spaced 'nice' tick positions from 1 to length."""
    if length <= 1:
        return [1]
    magnitude = 1
    while True:
        for base in (1, 2, 5):
            step = base * magnitude
            if length / step <= 8:
                return [1] + list(range(step, length + 1, step))
        magnitude *= 10


def _coord_walker(hsp, side):
    """Return (start, direction, step) for walking along one sequence of an HSP.

    Minus strand hits come with from > to; negative frames keep from < to,
    so both are normalized here. Translated sequences move 3 bases per residue.
    """
    start, end = hsp[f"{side}_from"], hsp[f"{side}_to"]
    frame = hsp.get(f"{side}_frame")
    minus = (frame is not None and frame < 0) or hsp.get(f"{side}_strand") == "Minus" or start > end
    step = 3 if frame is not None else 1
    if minus:
        return max(start, end), -1, step
    return min(start, end), 1, step


def _build_hsp(hsp):
    qseq, hseq, midline = hsp["qseq"], hsp["hseq"], hsp["midline"]
    align_len = hsp["align_len"]

    stats = {
        "bit_score": hsp["bit_score"],
        "score": hsp["score"],
        "evalue": hsp["evalue"],
        "query_from": hsp["query_from"],
        "query_to": hsp["query_to"],
        "identity": hsp["identity"],
        "align_len": align_len,
        "identity_pct": 100.0 * hsp["identity"] / align_len,
        "gaps": hsp.get("gaps", 0),
        "gaps_pct": 100.0 * hsp.get("gaps", 0) / align_len,
        "positive": hsp.get("positive"),
        "positive_pct": 100.0 * hsp["positive"] / align_len if "positive" in hsp else None,
    }
    if "query_frame" in hsp or "hit_frame" in hsp:
        frames = [f"{hsp[k]:+d}" for k in ("query_frame", "hit_frame") if k in hsp]
        stats["orientation"] = "Frame = " + "/".join(frames)
    elif "hit_strand" in hsp:
        stats["orientation"] = f"Strand = {hsp.get('query_strand', 'Plus')}/{hsp['hit_strand']}"
    else:
        stats["orientation"] = ""

    q_pos, q_dir, q_step = _coord_walker(hsp, "query")
    h_pos, h_dir, h_step = _coord_walker(hsp, "hit")
    blocks = []
    for i in range(0, len(qseq), ALIGNMENT_WIDTH):
        q_chunk, h_chunk = qseq[i:i + ALIGNMENT_WIDTH], hseq[i:i + ALIGNMENT_WIDTH]
        q_start, q_end, q_pos = _advance(q_pos, q_dir, q_step, q_chunk)
        h_start, h_end, h_pos = _advance(h_pos, h_dir, h_step, h_chunk)
        blocks.append((q_start, q_chunk, q_end, midline[i:i + ALIGNMENT_WIDTH], h_start, h_chunk, h_end))

    width = max(len(str(n)) for b in blocks for n in (b[0], b[2], b[4], b[6]))
    lines = []
    for q_start, q_chunk, q_end, mid, h_start, h_chunk, h_end in blocks:
        lines.append(f"Query  {q_start:<{width}}  {q_chunk}  {q_end}")
        lines.append(f"       {'':<{width}}  {mid}")
        lines.append(f"Sbjct  {h_start:<{width}}  {h_chunk}  {h_end}")
        lines.append("")
    stats["alignment"] = "\n".join(lines).rstrip()
    return stats


def _advance(pos, direction, step, chunk):
    """Return (block start, block end, next position) for one alignment row."""
    residues = len(chunk) - chunk.count("-")
    if residues == 0:
        return pos, pos, pos
    end = pos + direction * (residues * step - 1)
    return pos, end, end + direction
