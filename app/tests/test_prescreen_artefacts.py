"""The pre-screen candidate set has to survive the run as a real artefact.

``results.csv`` holds only the candidates that SURVIVED the specificity screen,
so a rank-enrichment null computed from it is conditioned on the outcome. Fixing
that needs the set that was CONSIDERED. The pipeline builds that set, checkpoints
it in ``.candidate_cache/`` — and then deletes it on success, which is why the
July 2026 whole-transcriptome runs cannot be corrected without re-running them.

These tests pin the artefacts that close that hole, and in particular the two
things it is easy to get subtly wrong:
  * what "considered" means relative to the ``screen_top_per_gene`` cap, and
  * that the parallel filter path records the same thing as the serial one.
"""
import csv
import gzip
import os

import pytest

from flexprobe.errors import ConfigError
from flexprobe.models import (
    ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
)
from flexprobe.prescreen_export import (
    CANDIDATES_FILE, PRECAP_FILE, SUMMARY_FILE, TRANSCRIPTS_FILE,
)

# Two genes, each ~200 bp of GC-balanced, low-homopolymer sequence, so step-1
# tiling yields far more filter-passing windows per gene than any cap below.
_SEQ_A = ("ACGTACGTGCATGCATACGTACGTGCATGCATTGCATGCAACGTTGCAAC"
          "GTACGTGCATGCATACGTGTCAGCATGCATACGTACGTGCATTGCATGCA"
          "ACGTACGTGCATGCATACGTACGTGCATGCATTGCAACGTTGCATGCATA"
          "CGTACGTGCATGCATACGTGTCAGCATGCATACGTACGTGCATTGCATGC")
_SEQ_B = ("TGCATGCAACGTACGTGCATTGCATGCAACGTACGTGCATGCATACGTAC"
          "GTGCATGCATTGCATGCAACGTTGCAACGTACGTGCATGCATACGTGTCA"
          "GCATGCATACGTACGTGCATTGCATGCAACGTACGTGCATGCATACGTAC"
          "GTGCATGCATTGCAACGTTGCATGCATACGTACGTGCATGCATACGTGTC")


def _cfg(tmp_path, *, top_n=5, workers=1, **kw):
    fa = tmp_path / "g.fasta"
    fa.write_text(f">GENE1 cdna\n{_SEQ_A}\n>GENE2 cdna\n{_SEQ_B}\n")
    return RunConfig(
        fasta_path=str(fa), output_dir=str(tmp_path / "out"),
        modality=ModalityDefinition(name="flex", version="v2"),
        # BLAST is monkeypatched out; db_paths just has to be non-empty for the
        # cap (which only applies when BLAST will actually run) to engage.
        screen=ScreenResources(skip_screen=False, db_paths=["/fake/db"]),
        export=ExportSettings(overwrite=True), discovery_only=True,
        screen_top_per_gene=top_n, workers=workers, **kw,
    )


def _patch_blast(monkeypatch):
    import flexprobe.pipeline as pl
    seen = {}

    def fake_blast(passing, *_a, **_k):
        seen["n"] = len(passing)
        seen["ids"] = {c.probe_id for c in passing}
        return {c.probe_id: [] for c in passing}

    monkeypatch.setattr(pl, "_run_screen", fake_blast)
    return seen


def _read(path):
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def _out(tmp_path, name):
    return os.path.join(str(tmp_path / "out"), name)


# ---- default is off ----

def test_no_artefacts_by_default(tmp_path, monkeypatch):
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    run_pipeline(_cfg(tmp_path))
    for name in (CANDIDATES_FILE, PRECAP_FILE, SUMMARY_FILE):
        assert not os.path.exists(_out(tmp_path, name)), name


# ---- the capped (screened) set ----

def test_emits_exactly_what_blast_screened(tmp_path, monkeypatch):
    """The artefact must be the BLAST denominator — not a superset, not a subset."""
    from flexprobe.pipeline import run_pipeline
    seen = _patch_blast(monkeypatch)
    run_pipeline(_cfg(tmp_path, emit_prescreen_candidates=True))
    rows = _read(_out(tmp_path, CANDIDATES_FILE))
    assert len(rows) == seen["n"] == 10          # 2 genes x cap 5
    assert {r["probe_id"] for r in rows} == seen["ids"]


def test_artefact_columns(tmp_path, monkeypatch):
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    run_pipeline(_cfg(tmp_path, emit_prescreen_candidates=True))
    rows = _read(_out(tmp_path, CANDIDATES_FILE))
    assert set(rows[0]) == {
        "probe_id", "transcript_id", "gene", "start", "end",
        "gc_content", "gc_left", "gc_right", "homopolymer_max",
        "junction_base", "prescreen_score", "rank_in_gene",
    }
    r = rows[0]
    assert r["gene"] and r["transcript_id"]
    assert 0.0 <= float(r["prescreen_score"]) <= 1.0
    assert int(r["end"]) > int(r["start"])


def test_rank_in_gene_reproduces_the_caps_ordering(tmp_path, monkeypatch):
    """rank_in_gene is what makes the cap's cut point legible."""
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    run_pipeline(_cfg(tmp_path, emit_prescreen_candidates=True))
    rows = _read(_out(tmp_path, CANDIDATES_FILE))
    by_gene = {}
    for r in rows:
        by_gene.setdefault(r["gene"], []).append(r)
    for gene, rs in by_gene.items():
        ranks = [int(r["rank_in_gene"]) for r in rs]
        assert sorted(ranks) == list(range(1, len(rs) + 1)), gene
        scores = [float(r["prescreen_score"])
                  for r in sorted(rs, key=lambda x: int(x["rank_in_gene"]))]
        assert scores == sorted(scores, reverse=True), gene


# ---- the pre-cap summary: the cap is itself a selection step ----

def test_summary_records_precap_and_screened_counts(tmp_path, monkeypatch):
    """Post-cap alone hides a selection step. The summary is what exposes it."""
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    run_pipeline(_cfg(tmp_path, top_n=5, emit_prescreen_candidates=True))
    rows = _read(_out(tmp_path, SUMMARY_FILE))
    assert len(rows) == 2                       # one row per gene
    for r in rows:
        passed = int(r["candidates_passed_filters"])
        screened = int(r["candidates_screened"])
        assert screened == 5
        assert passed > screened, (
            "the cap discarded candidates and the summary must say so")
        assert abs(float(r["kept_fraction"]) - screened / passed) < 1e-6


def test_summary_screened_totals_match_the_candidate_file(tmp_path, monkeypatch):
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    run_pipeline(_cfg(tmp_path, emit_prescreen_candidates=True))
    summary = _read(_out(tmp_path, SUMMARY_FILE))
    cands = _read(_out(tmp_path, CANDIDATES_FILE))
    assert sum(int(r["candidates_screened"]) for r in summary) == len(cands)


# ---- the optional full pre-cap dump ----

def _reconstruct_precap(tmp_path, tile_length=50):
    """Rebuild (probe_id, gene) for every pre-cap row using ONLY the documented
    reconstruction rules — which is the real test of the compact schema."""
    tx = {r["t"]: r for r in _read(_out(tmp_path, TRANSCRIPTS_FILE))}
    out = []
    for r in _read(_out(tmp_path, PRECAP_FILE)):
        t = tx[r["t"]]
        start = int(r["start"])
        out.append({
            "probe_id": f"{t['transcript_id']}_t{start:06d}",
            "gene": t["gene"],
            "end": start + tile_length,
            "score": float(r["prescreen_score"]),
            "junction_base": r["junction_base"],
        })
    return out


def test_precap_dump_is_strictly_larger(tmp_path, monkeypatch):
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    run_pipeline(_cfg(tmp_path, emit_prescreen_candidates=True,
                      emit_prescreen_precap_rows=True))
    precap = _reconstruct_precap(tmp_path)
    capped = _read(_out(tmp_path, CANDIDATES_FILE))
    assert len(precap) > len(capped)
    # The reconstruction rule has to actually work: every candidate BLAST saw
    # must be findable in the pre-cap set by its rebuilt probe_id.
    assert {r["probe_id"] for r in capped} <= {r["probe_id"] for r in precap}
    # and the row count must agree with the counts the summary reports
    summary = _read(_out(tmp_path, SUMMARY_FILE))
    assert sum(int(r["candidates_passed_filters"]) for r in summary) == len(precap)


def test_precap_schema_is_compact(tmp_path, monkeypatch):
    """Four columns, and no repeated transcript-id string.

    The transcript id inline would be the single largest cost in a 100M-row
    file; this is the guard on the dictionary encoding not quietly regressing.
    """
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    run_pipeline(_cfg(tmp_path, emit_prescreen_candidates=True,
                      emit_prescreen_precap_rows=True))
    rows = _read(_out(tmp_path, PRECAP_FILE))
    assert set(rows[0]) == {"t", "start", "prescreen_score", "junction_base"}
    assert all(r["t"].isdigit() for r in rows)
    tx = _read(_out(tmp_path, TRANSCRIPTS_FILE))
    assert set(tx[0]) == {"t", "transcript_id", "gene", "length"}
    assert {r["transcript_id"] for r in tx} == {"GENE1", "GENE2"}
    assert all(int(r["length"]) > 0 for r in tx)


def test_precap_scores_match_the_postcap_file(tmp_path, monkeypatch):
    """Both files must report the same prescreen_score for the same candidate.

    They are computed by different call paths (streaming vs batch), so this is
    the check that the compact file is a genuine superset and not a re-derived
    approximation.
    """
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    run_pipeline(_cfg(tmp_path, emit_prescreen_candidates=True,
                      emit_prescreen_precap_rows=True))
    precap = {r["probe_id"]: r for r in _reconstruct_precap(tmp_path)}
    for r in _read(_out(tmp_path, CANDIDATES_FILE)):
        assert abs(precap[r["probe_id"]]["score"]
                   - float(r["prescreen_score"])) < 5e-5
        assert precap[r["probe_id"]]["gene"] == r["gene"]
        assert precap[r["probe_id"]]["end"] == int(r["end"])


# ---- the run must say which set it captured ----

@pytest.mark.parametrize("kw,expected", [
    ({}, "none"),
    ({"emit_prescreen_candidates": True}, "post_cap"),
    ({"emit_prescreen_candidates": True, "emit_prescreen_precap_rows": True},
     "pre_cap"),
])
def test_capture_mode_is_recorded_in_the_run_stats(tmp_path, monkeypatch,
                                                   kw, expected):
    """'pre-cap or post-cap?' must never be inferred from which files exist."""
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    result = run_pipeline(_cfg(tmp_path, **kw))
    assert result["stats"]["prescreen_capture"] == expected


@pytest.mark.parametrize("precap,expected", [(False, "post_cap"),
                                             (True, "pre_cap")])
def test_index_file_states_the_capture_mode(tmp_path, monkeypatch,
                                            precap, expected):
    import json
    from flexprobe.pipeline import run_pipeline
    from flexprobe.prescreen_export import INDEX_FILE
    _patch_blast(monkeypatch)
    run_pipeline(_cfg(tmp_path, emit_prescreen_candidates=True,
                      emit_prescreen_precap_rows=precap))
    idx = json.load(open(_out(tmp_path, INDEX_FILE)))
    assert idx["capture"] == expected
    assert expected.replace("_", "-") in idx["capture_meaning"] \
        or expected.split("_")[0] in idx["capture_meaning"]
    assert idx["tile_length"] == 50
    assert idx["screen_top_per_gene"] == 5
    # the reconstruction rules the compact schema depends on are written down
    assert idx["reconstruction"]["end"] == "start + tile_length"
    assert "probe_id" in idx["reconstruction"]
    assert CANDIDATES_FILE in idx["files"]
    assert (PRECAP_FILE in idx["files"]) is precap


def test_precap_needs_the_main_artefact(tmp_path):
    cfg = _cfg(tmp_path, emit_prescreen_precap_rows=True)
    with pytest.raises(ConfigError, match="emit_prescreen_precap_rows"):
        cfg.validate()


# ---- the parallel filter path must record the same thing ----

def test_parallel_path_produces_the_same_artefacts(tmp_path, monkeypatch):
    """The tally hangs off the filter stream, which has two implementations.

    A per-gene count collected on only one of them would be silently wrong on
    every real (workers > 1) HPC run and right in every test.
    """
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    serial = tmp_path / "s"
    parallel = tmp_path / "p"
    serial.mkdir(); parallel.mkdir()
    run_pipeline(_cfg(serial, workers=1, emit_prescreen_candidates=True))
    run_pipeline(_cfg(parallel, workers=2, emit_prescreen_candidates=True))
    for name in (CANDIDATES_FILE, SUMMARY_FILE):
        a = _read(_out(serial, name))
        b = _read(_out(parallel, name))
        assert a == b, f"{name} differs between the serial and parallel paths"


# ---- the belt-and-braces cache ----

def test_candidate_cache_is_deleted_by_default(tmp_path, monkeypatch):
    from flexprobe import candidate_cache
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    cfg = _cfg(tmp_path)
    run_pipeline(cfg)
    assert not os.path.exists(candidate_cache.cache_dir(cfg.output_dir))


def test_keep_candidate_cache_survives_a_successful_run(tmp_path, monkeypatch):
    from flexprobe import candidate_cache
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    cfg = _cfg(tmp_path, keep_candidate_cache=True)
    run_pipeline(cfg)
    cdir = candidate_cache.cache_dir(cfg.output_dir)
    assert os.path.isdir(cdir)
    assert os.path.exists(os.path.join(cdir, "candidates.jsonl"))


def test_screen_cache_is_still_cleared_when_keeping_candidates(tmp_path, monkeypatch):
    """Only the small cache is kept — the BLAST cache is the huge one."""
    from flexprobe import screen_cache
    from flexprobe.pipeline import run_pipeline
    _patch_blast(monkeypatch)
    cfg = _cfg(tmp_path, keep_candidate_cache=True)
    run_pipeline(cfg)
    assert not os.path.exists(
        os.path.join(cfg.output_dir, screen_cache.CACHE_DIRNAME))


# ---- the new flags must not invalidate a resume checkpoint ----

def test_new_flags_are_not_in_the_cache_fingerprint(tmp_path):
    """Turning the artefact on mid-chain must not force a re-tile.

    A whole-transcriptome tile+filter pass costs hours; invalidating a good
    checkpoint to produce a byte-identical candidate set would be an expensive
    way to change nothing.
    """
    from flexprobe import candidate_cache
    fa = tmp_path / "g.fasta"
    fa.write_text(f">GENE1\n{_SEQ_A}\n")
    base = RunConfig(fasta_path=str(fa))
    for field in ("emit_prescreen_candidates", "emit_prescreen_precap_rows",
                  "keep_candidate_cache"):
        other = RunConfig(fasta_path=str(fa), **{field: True})
        assert candidate_cache.fingerprint(base) == \
            candidate_cache.fingerprint(other), field
