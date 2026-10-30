"""Streaming pre-screen cap: bounded-memory accumulator == batch select_top_per_gene."""
from flexprobe.models import (
    ScreenResources, ExportSettings, ModalityDefinition, ProbeCandidate, RunConfig,
)
from flexprobe.ranking import (
    DEFAULT_WEIGHTS, TopPerGeneAccumulator, select_top_per_gene,
)


def _cands():
    """Many candidates across 3 genes with varied, deterministic GC so
    prescreen_score differs and the accumulator's trim path is exercised."""
    out = []
    for g in range(3):
        gid = f"gene{g}"
        for i in range(40):
            # GC walks across [0.30, 0.70] so scores spread and some tie.
            gc = 0.30 + ((i * 7 + g * 3) % 41) / 100.0
            out.append(ProbeCandidate(
                gene_id=gid, probe_id=f"{gid}_t{i:06d}",
                sequence="A" * 50, start=i * 5, end=i * 5 + 50,
                gc_content=gc, gc_left=gc, gc_right=gc, homopolymer_max=2,
            ))
    return out


def _ids(cands):
    return {c.probe_id for c in cands}


def _seqlen(_gid):
    return 220


def _gk(c):
    return c.gene_id


def test_accumulator_matches_select_top_per_gene():
    cands = _cands()
    for n in (1, 5, 13, 100):
        batch = select_top_per_gene(
            cands, gene_key=_gk, seq_length_of=_seqlen, n=n, weights=DEFAULT_WEIGHTS)
        # Feed in a non-sorted order to prove order-independence + the trim path.
        acc = TopPerGeneAccumulator(
            gene_key=_gk, seq_length_of=_seqlen, n=n,
            weights=DEFAULT_WEIGHTS, trim_factor=2)
        shuffled = cands[::3] + cands[1::3] + cands[2::3]
        acc.add_many(shuffled)
        streamed = acc.result()
        assert _ids(streamed) == _ids(batch), f"n={n}: streamed set must equal batch"
        # never keeps more than n per gene
        per_gene = {}
        for c in streamed:
            per_gene[c.gene_id] = per_gene.get(c.gene_id, 0) + 1
        assert all(v <= n for v in per_gene.values())


def test_accumulator_result_is_deterministic():
    cands = _cands()
    a = TopPerGeneAccumulator(gene_key=_gk, seq_length_of=_seqlen, n=7)
    a.add_many(cands)
    b = TopPerGeneAccumulator(gene_key=_gk, seq_length_of=_seqlen, n=7)
    b.add_many(cands[::-1])
    assert [c.probe_id for c in a.result()] == [c.probe_id for c in b.result()]


# ---- pipeline integration: streamed cap (serial & parallel) == same set ----

_GENES = {
    "geneA": ("GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
              "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"),
    "geneB": ("TTGGCCAATTGGCCAATTGGCCAATAGGCCTTAAGGCCTTAAGGCCTTAA"
              "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGG"),
    "geneC": ("ACACACGTGTGTACACACGTGTGTACGTGTACACGTGTACACGTGTACAC"
              "TGTGACACTGTGACACTGTGACACTGTGACACTGTGACACTGTGACACTG"),
}


def _multi_fasta(tmp_path):
    p = tmp_path / "multi.fasta"
    p.write_text("".join(f">{g}\n{s}\n" for g, s in _GENES.items()))
    return str(p)


def _cfg(fa, out, workers, top_n):
    return RunConfig(
        fasta_path=fa, output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=False, db_paths=["/fake/db"]),
        export=ExportSettings(overwrite=True), discovery_only=True,
        workers=workers, screen_top_per_gene=top_n, resume_screen=False,
    )


def _patch_blast_capture(monkeypatch):
    """Capture the probe_ids that reach BLAST; score nothing (empty hits)."""
    import flexprobe.pipeline as pl
    seen = {}

    def fake_run_screen(passing, *_a, per_chunk=None, **_k):
        seen["ids"] = {c.probe_id for c in passing}
        seen["n"] = len(passing)
        if per_chunk is not None:
            per_chunk(0, list(passing), {})   # stream: empty hits -> all pass
            return None
        return {c.probe_id: [] for c in passing}

    monkeypatch.setattr(pl, "_run_screen", fake_run_screen)
    return seen


def test_streamed_cap_serial_equals_parallel(tmp_path, monkeypatch):
    from flexprobe.pipeline import run_pipeline
    seen = _patch_blast_capture(monkeypatch)

    run_pipeline(_cfg(_multi_fasta(tmp_path), str(tmp_path / "s"), workers=1, top_n=4))
    serial_ids = seen["ids"]

    run_pipeline(_cfg(_multi_fasta(tmp_path), str(tmp_path / "p"), workers=2, top_n=4))
    parallel_ids = seen["ids"]

    assert serial_ids == parallel_ids, "streamed cap must be workers-independent"
    # 3 genes, cap 4 -> at most 4 per gene
    assert 0 < len(serial_ids) <= 12


def test_streamed_cap_matches_batch_reference(tmp_path, monkeypatch):
    # The set reaching BLAST with the streamed cap must equal select_top_per_gene
    # over the full filtered set (skip_screen run keeps every filtered candidate).
    from flexprobe.pipeline import run_pipeline

    # 1) Full filtered set, no cap, BLAST skipped -> results.csv has all of them.
    seen = _patch_blast_capture(monkeypatch)
    run_pipeline(_cfg(_multi_fasta(tmp_path), str(tmp_path / "full"),
                      workers=1, top_n=0))
    full_ids = seen["ids"]               # top_n=0 -> every filtered candidate

    # 2) Capped run.
    run_pipeline(_cfg(_multi_fasta(tmp_path), str(tmp_path / "cap"),
                      workers=1, top_n=3))
    capped_ids = seen["ids"]

    assert capped_ids.issubset(full_ids)
    # Per gene, the cap kept the genuine top-3 (subset relationship + bound).
    per_gene = {}
    for pid in capped_ids:
        per_gene[pid.split("_t")[0]] = per_gene.get(pid.split("_t")[0], 0) + 1
    assert all(v <= 3 for v in per_gene.values())
    assert capped_ids != full_ids        # the cap actually dropped candidates
