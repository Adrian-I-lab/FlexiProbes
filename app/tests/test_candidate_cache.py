"""Candidate checkpoint (resume-skip-tiling) + streamed-specificity parity."""
from flexprobe import blast as fp_blast
from flexprobe.screen import ENGINE_BLAST
import dataclasses
import os
import threading

import pytest

from flexprobe import candidate_cache
from flexprobe.errors import CancellationError
from flexprobe.models import (
    ScreenResources, ExportSettings, FilterThresholds, ModalityDefinition,
    ProbeCandidate, RunConfig, SpecificityThresholds,
)
from flexprobe.pipeline import run_pipeline
from flexprobe.specificity import batch_evaluate

from conftest import SYNTH_200, write_fasta


# ---------------------------------------------------------------------------
# candidate_cache module
# ---------------------------------------------------------------------------

def _filters():
    return FilterThresholds(gc_min=0.30, gc_max=0.80, max_n_per_half=2,
                            max_homopolymer_run=6, require_tn_junction=False)


def _cfg(fa, **over):
    kw = dict(fasta_path=fa, output_dir="/tmp/out", tile_length=50, tile_step=5,
              filters=_filters(), screen_top_per_gene=0, group_by_gene=False,
              ranking_weights=None)
    kw.update(over)
    return RunConfig(**kw)


def _fp(fa, **over):
    return candidate_cache.fingerprint(_cfg(fa, **over))


def test_fingerprint_stable_and_sensitive(tmp_path):
    fa = str(write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)]))
    base = _fp(fa)
    assert base == _fp(fa)                                   # deterministic
    assert base != _fp(fa, tile_step=1)                      # tiling changed
    assert base != _fp(fa, screen_top_per_gene=50)            # cap changed
    assert base != _fp(fa, group_by_gene=True)               # grouping changed
    assert base != _fp(fa, filters=FilterThresholds(gc_min=0.10))
    # editing the input FASTA (size/mtime) invalidates the checkpoint
    write_fasta(tmp_path / "g.fa", [("G", SYNTH_200 + "ACGT")])
    assert base != _fp(fa)


# ---------------------------------------------------------------------------
# fingerprint fails CLOSED: every output-affecting config field invalidates it
# ---------------------------------------------------------------------------
#
# The bug this pins down: the fingerprint was an allowlist of seven hand-picked
# values, so the pre-tiling trim options and `filters.reject_junction_base_n`
# were invisible to it. In the documented HPC workflow (walltime kill ->
# resubmit with a flag added) such a flag produced an IDENTICAL digest: the
# stale checkpoint was reloaded, tiling and filtering were skipped, and probes
# designed into the poly-A tail came back as if the flag had been honoured.
#
# The parametrization is DERIVED from RunConfig rather than hand-listed, so a
# field added to the config in future automatically acquires a test case and
# must either invalidate the cache or be justified in FINGERPRINT_IGNORED.

def _leaf_paths(obj, prefix=""):
    """Every dotted leaf path of a RunConfig instance."""
    for f in dataclasses.fields(obj):
        val = getattr(obj, f.name)
        path = f"{prefix}{f.name}"
        if dataclasses.is_dataclass(val):
            yield from _leaf_paths(val, prefix=f"{path}.")
        else:
            yield path, val


def _mutate(val):
    """A different value of the same shape."""
    if isinstance(val, bool):
        return not val
    if isinstance(val, int):
        return val + 7
    if isinstance(val, float):
        return val + 0.11
    if isinstance(val, str):
        return (val or "x") + "_CHANGED"
    if isinstance(val, list):
        return list(val) + ["CHANGED"]
    if isinstance(val, dict):
        return {**val, "gc_proximity": 0.99}
    if val is None:
        return {"gc_proximity": 0.99}
    raise AssertionError(f"no mutation rule for {type(val)!r}")


def _set_path(cfg, path, value):
    obj = cfg
    parts = path.split(".")
    for p in parts[:-1]:
        obj = getattr(obj, p)
    setattr(obj, parts[-1], value)


def _all_leaf_paths():
    return [p for p, _v in _leaf_paths(RunConfig())]


def _ignored(path):
    """True if this leaf is covered by a FINGERPRINT_IGNORED entry."""
    return path.split(".")[0] in candidate_cache.FINGERPRINT_IGNORED


_RELEVANT = [p for p in _all_leaf_paths() if not _ignored(p)]
_IGNORED = [p for p in _all_leaf_paths() if _ignored(p)]


@pytest.mark.parametrize("path", _RELEVANT)
def test_every_output_affecting_field_invalidates_the_cache(tmp_path, path):
    fa = str(write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)]))
    base_cfg = _cfg(fa)
    base = candidate_cache.fingerprint(base_cfg)

    mutated = _cfg(fa)
    current = mutated
    for p in path.split(".")[:-1]:
        current = getattr(current, p)
    _set_path(mutated, path, _mutate(getattr(current, path.split(".")[-1])))

    assert candidate_cache.fingerprint(mutated) != base, (
        f"changing {path} left the candidate-checkpoint digest unchanged — a "
        f"resumed run would reload candidates built under the OLD value. Either "
        f"fold it into the fingerprint or justify it in FINGERPRINT_IGNORED."
    )


@pytest.mark.parametrize("path", _IGNORED)
def test_denylisted_fields_do_not_invalidate_the_cache(tmp_path, path):
    """The denylist is small and deliberate; this is its executable spec.

    Every entry is a field that provably cannot change which candidates the
    tile+filter+cap stage emits (see FINGERPRINT_IGNORED for the per-field
    justification), so it must NOT force a whole-transcriptome retile.
    """
    fa = str(write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)]))
    base = candidate_cache.fingerprint(_cfg(fa))

    mutated = _cfg(fa)
    current = mutated
    for p in path.split(".")[:-1]:
        current = getattr(current, p)
    _set_path(mutated, path, _mutate(getattr(current, path.split(".")[-1])))

    assert candidate_cache.fingerprint(mutated) == base, \
        f"{path} is in FINGERPRINT_IGNORED but changes the digest"


def test_the_three_fields_the_old_allowlist_missed(tmp_path):
    """Named explicitly so the specific reported defects can't silently return."""
    fa = str(write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)]))
    base = _fp(fa)
    assert base != _fp(fa, filters=FilterThresholds(reject_junction_base_n=False))


def test_schema_bump_invalidates_existing_caches(tmp_path, monkeypatch):
    """A checkpoint written by schema-1 code must never be matched by schema 2."""
    fa = str(write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)]))
    assert candidate_cache.CACHE_SCHEMA >= 2
    now = _fp(fa)
    monkeypatch.setattr(candidate_cache, "CACHE_SCHEMA", 1)
    assert _fp(fa) != now


def _cand(i):
    return ProbeCandidate(gene_id="G", probe_id=f"G_t{i:06d}",
                          sequence="ACGTACGTAC" * 5, start=i, end=i + 50,
                          gc_content=0.5, gc_left=0.48, gc_right=0.52,
                          homopolymer_max=3)


def test_candidate_round_trip(tmp_path):
    out = str(tmp_path)
    cands = [_cand(i) for i in range(5)]
    candidate_cache.save(out, "fpA", cands, tiles_total=999, filter_passed=42)
    hit = candidate_cache.load(out, "fpA")
    assert hit is not None
    loaded, meta = hit
    assert [c.probe_id for c in loaded] == [c.probe_id for c in cands]
    assert loaded[0].sequence == cands[0].sequence
    assert loaded[2].gc_content == cands[2].gc_content
    assert meta["tiles_total"] == 999 and meta["filter_passed"] == 42


def test_load_rejects_wrong_fingerprint(tmp_path):
    out = str(tmp_path)
    candidate_cache.save(out, "fpA", [_cand(0)])
    assert candidate_cache.load(out, "fpB") is None      # stale config -> miss


def test_load_rejects_incomplete_meta(tmp_path):
    import json
    out = str(tmp_path)
    candidate_cache.save(out, "fpA", [_cand(0)])
    # Simulate a save killed before the meta committed.
    meta_path = os.path.join(candidate_cache.cache_dir(out), "_meta.json")
    # Keep the stamp save() wrote: this test is about the `complete` flag, and
    # an unstamped meta is refused before the flag is read (test_cache_stamp).
    with open(meta_path) as fh:
        stamp = json.load(fh)["stamp"]
    with open(meta_path, "w") as fh:
        json.dump({"fingerprint": "fpA", "complete": False, "stamp": stamp}, fh)
    assert candidate_cache.load(out, "fpA") is None


def test_clear_removes_cache(tmp_path):
    out = str(tmp_path)
    candidate_cache.save(out, "fpA", [_cand(0)])
    assert os.path.isdir(candidate_cache.cache_dir(out))
    candidate_cache.clear(out)
    assert not os.path.isdir(candidate_cache.cache_dir(out))


# ---------------------------------------------------------------------------
# streamed-specificity parity (per-chunk == one big batch)
# ---------------------------------------------------------------------------

def _mock_blast(monkeypatch, pl):
    """Fake blastn: one off-target hit for the first probe of every chunk."""
    monkeypatch.setattr(fp_blast, "discover_blast",
                        lambda: {"blastn": "b", "makeblastdb": "m"})
    monkeypatch.setattr(fp_blast, "check_version", lambda *a, **k: "2.13.0")
    monkeypatch.setattr(pl, "is_db_prefix", lambda _p: True)

    def fake_run_blastn(_b, query_fa, _db, **_k):
        tsv = query_fa + ".tsv"
        first = None
        with open(query_fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    first = line[1:].strip()
                    break
        with open(tsv, "w") as fh:
            if first:
                fh.write(f"{first}\tsubjX\t100.0\t25\t0\t0\t1\t25\t1\t25\t1e-5\t50\n")
        return tsv

    def fake_parse(tsv):
        rows = []
        with open(tsv) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 12:
                    rows.append({"qseqid": p[0], "sseqid": p[1],
                                 "pident": p[2], "length": p[3]})
        return rows

    monkeypatch.setattr(fp_blast, "run_blastn", fake_run_blastn)
    monkeypatch.setattr(fp_blast, "parse_hits_tsv", fake_parse)
    monkeypatch.setattr("flexprobe.constants.BLAST_QUERY_CHUNK", 2)


def test_streaming_specificity_matches_batch(tmp_path, monkeypatch):
    import flexprobe.pipeline as pl
    _mock_blast(monkeypatch, pl)

    # Batch path: accumulate every hit, then evaluate the whole set at once.
    batch_cands = [_cand(i) for i in range(5)]
    out_a = tmp_path / "a"; out_a.mkdir()
    hits = pl._run_screen(batch_cands, ["/db"], str(out_a),
                                min_pident=80.0, min_aln=20, resume=False)
    assert hits is not None    # no per_chunk -> full dict returned
    batch_evaluate(batch_cands, hits, set())
    batch = {c.probe_id: c.specificity_pass for c in batch_cands}

    # Streamed path: score each chunk as its hits arrive, holding nothing.
    stream_cands = [_cand(i) for i in range(5)]
    out_b = tmp_path / "b"; out_b.mkdir()

    def sink(_ci, chunk, chunk_hits):
        batch_evaluate(chunk, chunk_hits, set())

    ret = pl._run_screen(stream_cands, ["/db"], str(out_b),
                               min_pident=80.0, min_aln=20, resume=False,
                               per_chunk=sink)
    assert ret is None         # streamed -> nothing accumulated/returned
    stream = {c.probe_id: c.specificity_pass for c in stream_cands}

    assert stream == batch
    # Sanity on what the mock produced: it emits a hit for the FIRST query of
    # each chunk only, and the query order is <probe>__L then <probe>__R — so
    # probe 0 of each chunk has an off-target on its LHS half and nothing on its
    # RHS half. Per CG000621/CG000839 ("at least five mismatches in at least one
    # of the LHS or RHS probes") that pair still passes: the clean RHS half
    # cannot hybridise there, so the nick never forms and no ligation product
    # (hence no off-target signal) is made.
    assert batch["G_t000000"] is True and batch["G_t000001"] is True
    # ...and the legacy "either half rejects" mode still rejects it, so the
    # reproducibility switch really does restore the old verdict.
    legacy = [_cand(i) for i in range(5)]
    out_c = tmp_path / "c"; out_c.mkdir()
    legacy_hits = pl._run_screen(legacy, ["/db"], str(out_c),
                                       min_pident=80.0, min_aln=20, resume=False)
    batch_evaluate(legacy, legacy_hits, set(),
                   off_target_require_both_halves=False)
    assert {c.probe_id: c.specificity_pass for c in legacy}["G_t000000"] is False


# ---------------------------------------------------------------------------
# full-pipeline resume: a killed run leaves a checkpoint; the rerun skips tiling
# ---------------------------------------------------------------------------

def _resumable_cfg(fa, out):
    return RunConfig(
        fasta_path=fa, output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        filters=_filters(),
        screen=ScreenResources(db_paths=["/db/prefix"], skip_screen=False,
                             engine=ENGINE_BLAST),
        # The junction-duplex gate is k-mer-only, so a BLAST config now has
        # to disable it explicitly (models.RunConfig.validate refuses the
        # combination rather than reporting a gate that never ran).
        specificity=SpecificityThresholds(off_target_junction_duplex_bp=0),
        export=ExportSettings(overwrite=True),
        tile_length=50, tile_step=5, resume_screen=True,
    )


def test_pipeline_resumes_candidates_after_kill(tmp_path, monkeypatch):
    import flexprobe.pipeline as pl
    _mock_blast(monkeypatch, pl)

    # Count real tiling work so we can prove the rerun skips it.
    real_tile = pl.tile_sequence
    calls = {"tile": 0}

    def counting_tile(*a, **k):
        calls["tile"] += 1
        return real_tile(*a, **k)

    monkeypatch.setattr(pl, "tile_sequence", counting_tile)

    fa = str(write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)]))
    out = str(tmp_path / "out")
    cfg = _resumable_cfg(fa, out)

    # Run 1 — request cancellation the instant the candidate checkpoint lands
    # (Stage 3b), simulating a walltime kill just before BLAST.
    ev = threading.Event()

    def kill_after_checkpoint(e):
        if "Checkpointed" in e.message:
            ev.set()

    with pytest.raises(CancellationError):
        run_pipeline(cfg, cancel_event=ev, progress_callback=kill_after_checkpoint)

    assert calls["tile"] > 0, "run 1 must have tiled"
    assert os.path.isdir(candidate_cache.cache_dir(out)), \
        "the checkpoint must survive the kill (not cleared on failure)"

    # Run 2 — resume. Must reload candidates, skip tiling, and finish cleanly.
    calls["tile"] = 0
    msgs = []
    result = run_pipeline(cfg, progress_callback=lambda e: msgs.append(e.message))

    assert calls["tile"] == 0, "resume must skip tiling/filtering"
    assert any("Resumed" in m and "checkpoint" in m for m in msgs)
    assert result["stats"]["pairs_output"] > 0
    # a fully successful run drops both checkpoints so they never reach the tarball
    assert not os.path.isdir(candidate_cache.cache_dir(out))


def test_resumed_run_results_identical_to_clean_run(tmp_path, monkeypatch):
    # Audit #38: prove OUTPUT identity, not just the resume mechanics — a run
    # killed just before BLAST and resumed must produce a byte-identical
    # results.csv to a clean single-shot run on the same input.
    import flexprobe.pipeline as pl
    _mock_blast(monkeypatch, pl)
    fa = str(write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)]))

    # Clean single-shot run.
    clean_out = str(tmp_path / "clean")
    run_pipeline(_resumable_cfg(fa, clean_out))
    with open(os.path.join(clean_out, "results.csv")) as fh:
        clean_csv = fh.read()

    # Killed-then-resumed run to a separate output dir.
    resume_out = str(tmp_path / "resume")
    cfg = _resumable_cfg(fa, resume_out)
    ev = threading.Event()

    def kill_after_checkpoint(e):
        if "Checkpointed" in e.message:
            ev.set()

    with pytest.raises(CancellationError):
        run_pipeline(cfg, cancel_event=ev, progress_callback=kill_after_checkpoint)
    run_pipeline(cfg)   # resume to completion
    with open(os.path.join(resume_out, "results.csv")) as fh:
        resume_csv = fh.read()

    assert resume_csv == clean_csv, \
        "resumed results.csv must be row-identical to a clean run"
