"""Real-blastn resume: verdict cache survives interruption (needs real blastn).

The unit tests in test_screen_cache.py fake blastn. This exercises the whole
streaming path against a REAL blastn — the thing the HPC job chain bets on,
where each link resumes from `.screen_cache` written by the previous link.

A cache of raw hits on a human transcriptome reaches gigabytes per chunk and
can fill a disk partway through a run, so the cache holds verdicts only.
"""
import json
import os
import shutil

import pytest

from flexprobe import screen_cache
from flexprobe.screen import ENGINE_BLAST
from flexprobe.blast import build_db, discover_blast
from flexprobe.models import (
    ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
    SpecificityThresholds,
)
from flexprobe.pipeline import run_pipeline

pytestmark = pytest.mark.skipif(
    shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
    reason="ncbi-blast+ (blastn/makeblastdb) not installed",
)

TARGET = (
    "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGTTACCGATGCATGC"
    "ATGCATGCATGCATGCATGCAACGTACGTACGTACGTACGTTGCATGCATG"
    "CATGCATGCATG"
)


def _keep_cache(monkeypatch):
    """Stop the pipeline clearing the cache on success (pipeline.py:1353), so a
    test can inspect it. Resume in the wild only happens when a run did NOT
    finish (walltime kill), which is exactly when the cache survives."""
    monkeypatch.setattr(screen_cache, "clear", lambda _out: None)


def _setup(tmp_path):
    fa = tmp_path / "target.fasta"
    fa.write_text(f">GENE1\n{TARGET}\n")
    dbfa = tmp_path / "off.fasta"
    # Share sequence with the target so real hits are produced.
    dbfa.write_text(f">OFFTGT\n{TARGET}\n")
    db = str(tmp_path / "offdb")
    build_db(discover_blast()["makeblastdb"], str(dbfa), db, title="off")
    return fa, db


def _cfg(fa, out, db):
    return RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=False, db_paths=[db],
                             engine=ENGINE_BLAST),
        # The junction-duplex gate is k-mer-only, so a BLAST config now has
        # to disable it explicitly (models.RunConfig.validate refuses the
        # combination rather than reporting a gate that never ran).
        specificity=SpecificityThresholds(off_target_junction_duplex_bp=0),
        export=ExportSettings(overwrite=True), discovery_only=True,
    )


def _results(out):
    """(status, margin) per probe from results.csv — the user-visible answer."""
    import csv
    with open(os.path.join(str(out), "results.csv")) as fh:
        return sorted(
            (r["mrna_target_50bp"], r["specificity_status"])
            for r in csv.DictReader(fh)
        )


def test_real_screen_cache_holds_verdicts_not_hits(tmp_path, monkeypatch):
    """The cache must contain scored verdict fields, never raw BLAST hits."""
    _keep_cache(monkeypatch)
    fa, db = _setup(tmp_path)
    out = tmp_path / "out"
    run_pipeline(_cfg(fa, out, db))

    cdir = screen_cache.cache_dir(str(out))
    chunks = [f for f in os.listdir(cdir) if f.startswith("chunk_")]
    assert chunks, "expected at least one cached chunk"
    payload = json.load(open(os.path.join(cdir, chunks[0])))
    probe_id, verdict = next(iter(payload.items()))

    # verdict shape, not hit shape
    assert isinstance(verdict, dict), f"{probe_id} cached a list => raw hits!"
    assert set(screen_cache.VERDICT_FIELDS) <= set(verdict)
    # the raw-hit fingerprint must be absent
    for hit_field in ("qseqid", "sseqid", "bitscore", "mismatch", "gapopen"):
        assert hit_field not in verdict, f"raw hit field {hit_field!r} in cache"


def test_real_blast_resume_reproduces_uninterrupted_results(tmp_path, monkeypatch):
    """Delete nothing, re-run: every chunk resumes from cache, zero blastn
    calls, and results.csv is byte-for-byte equivalent."""
    _keep_cache(monkeypatch)
    fa, db = _setup(tmp_path)
    out = tmp_path / "out"
    run_pipeline(_cfg(fa, out, db))
    first = _results(out)

    calls = {"n": 0}
    # Count at the DEFINITION site: screen.BlastScreen reaches run_blastn
    # through the ``blast`` module, so that is where an invocation is
    # observable. Wrapping a name bound in ``pipeline`` would count zero
    # regardless of whether BLAST ran, and this test's whole assertion is that
    # the count is zero.
    from flexprobe import blast as fp_blast
    real = fp_blast.run_blastn

    def counting(*a, **k):
        calls["n"] += 1
        return real(*a, **k)

    fp_blast.run_blastn = counting
    try:
        run_pipeline(_cfg(fa, out, db))     # same output dir => resume
        second = _results(out)
    finally:
        fp_blast.run_blastn = real

    assert calls["n"] == 0, f"resume re-BLASTed {calls['n']} chunk(s)"
    assert second == first, "resumed results differ from the uninterrupted run"


def test_real_blast_stale_cache_is_not_reused_across_thresholds(tmp_path, monkeypatch):
    """A verdict is only valid for the thresholds that produced it. Changing a
    specificity threshold must invalidate the cache (fingerprint), not silently
    reuse verdicts scored under the old one."""
    from flexprobe.models import SpecificityThresholds
    _keep_cache(monkeypatch)
    fa, db = _setup(tmp_path)
    out = tmp_path / "out"
    run_pipeline(_cfg(fa, out, db))
    cdir = screen_cache.cache_dir(str(out))
    fp_before = json.load(open(os.path.join(cdir, "_meta.json")))["fingerprint"]
    assert [f for f in os.listdir(cdir) if f.startswith("chunk_")]

    # Same probes, same DB -- only a scoring threshold moves. Under schema 1
    # (raw hits) this legitimately reused the cache and re-scored. Under
    # schema 2 the cached VERDICTS were scored at the OLD threshold, so reuse
    # would silently return wrong answers: the fingerprint must change.
    cfg = _cfg(fa, out, db)
    # Junction gate still off — this config is on the BLAST engine, which
    # RunConfig.validate refuses to pair with a live gate.
    cfg.specificity = SpecificityThresholds(off_target_pident_min=99.0,
                                            off_target_junction_duplex_bp=0)
    run_pipeline(cfg)

    fp_after = json.load(open(os.path.join(cdir, "_meta.json")))["fingerprint"]
    assert fp_after != fp_before, (
        "changing a specificity threshold did NOT invalidate the verdict cache "
        "-- verdicts scored at the old threshold would be silently reused"
    )
