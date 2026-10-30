"""Resumable BLAST: checkpoint cache unit tests + pipeline resume behaviour."""
from flexprobe import blast as fp_blast
import math
import os

from flexprobe import screen_cache
from flexprobe.models import ProbeCandidate


# ---- cache module ----

def _ids(n):
    return [f"g1_{i:04d}" for i in range(n)]


def test_fingerprint_is_stable_and_sensitive():
    a = screen_cache.fingerprint(_ids(5), db_arg="db", chunk_size=2,
                                perc_identity=80.0, qcov_hsp_perc=70.0)
    # identical inputs -> identical digest (resume across runs)
    b = screen_cache.fingerprint(_ids(5), db_arg="db", chunk_size=2,
                                perc_identity=80.0, qcov_hsp_perc=70.0)
    assert a == b
    # any change to candidates / DB / thresholds -> different digest
    assert a != screen_cache.fingerprint(_ids(6), db_arg="db", chunk_size=2,
                                        perc_identity=80.0, qcov_hsp_perc=70.0)
    assert a != screen_cache.fingerprint(_ids(5), db_arg="other", chunk_size=2,
                                        perc_identity=80.0, qcov_hsp_perc=70.0)
    assert a != screen_cache.fingerprint(_ids(5), db_arg="db", chunk_size=5,
                                        perc_identity=80.0, qcov_hsp_perc=70.0)


def test_fingerprint_sensitive_to_ref_content_hash():
    """`db_arg` is a PATH; `ref_content_hash` is what catches the reference's
    FILE changing underneath an unchanged path (audit item 21 / E1)."""
    base = dict(db_arg="db", chunk_size=2, perc_identity=80.0, qcov_hsp_perc=70.0)
    a = screen_cache.fingerprint(_ids(5), **base, ref_content_hash="a" * 64)
    b = screen_cache.fingerprint(_ids(5), **base, ref_content_hash="a" * 64)
    assert a == b
    assert a != screen_cache.fingerprint(_ids(5), **base, ref_content_hash="b" * 64)
    # default (omitted) must agree with an explicit "" -- both mean "unknown".
    assert screen_cache.fingerprint(_ids(5), **base) == \
        screen_cache.fingerprint(_ids(5), **base, ref_content_hash="")


def test_blast_screen_ref_content_hash_changes_with_fasta_content(tmp_path, monkeypatch):
    """BlastScreen.fingerprint_parts() must reflect the reference's actual
    bytes, not merely its (unchanged) path. Direct unit coverage of the fix,
    independent of whether makeblastdb ever rebuilds the .nin/.nhr/.nsq at
    that path -- see `screen._blast_reference_content_hash`."""
    from flexprobe.screen import BlastScreen

    monkeypatch.setattr(fp_blast, "discover_blast",
                        lambda: {"blastn": "b", "makeblastdb": "m"})
    monkeypatch.setattr(fp_blast, "check_version", lambda *_a, **_k: "2.13.0")

    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">refA\nACGTACGTACGTACGTACGTACGTACGTACGT\n")

    s = BlastScreen()
    s.prepare([str(ref_fa)])
    hash_a = s.fingerprint_parts()["ref_content_hash"]
    assert hash_a, "a real FASTA on disk must produce a real content hash"

    # Reference B's CONTENT, copied over the SAME path.
    ref_fa.write_text(">refB\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n")
    s2 = BlastScreen()
    s2.prepare([str(ref_fa)])
    hash_b = s2.fingerprint_parts()["ref_content_hash"]

    assert hash_b != hash_a, (
        "reference FASTA content changed at an unchanged path but "
        "ref_content_hash didn't -- db_arg alone can't catch this"
    )
    # Stable, not merely "always different": identical content -> identical hash.
    s3 = BlastScreen()
    s3.prepare([str(ref_fa)])
    assert s3.fingerprint_parts()["ref_content_hash"] == hash_b


def test_kmer_screen_ref_content_hash_changes_with_fasta_content(tmp_path):
    """Same guarantee for the k-mer engine (the default). Its content hash is
    a free reuse of `ref.composition.source_sha256`, already computed by
    `KmerIndex.build_or_load` to key its OWN index cache -- see
    `KmerScreen.fingerprint_parts`."""
    from flexprobe.screen import KmerScreen

    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">refA\n" + "ACGT" * 40 + "\n")
    s = KmerScreen(cache_dir=str(tmp_path / "idxA"), mmap=False)
    s.prepare([str(ref_fa)])
    hash_a = s.fingerprint_parts()["ref_content_hash"]
    assert hash_a

    ref_fa.write_text(">refB\n" + "TTGG" * 40 + "\n")
    s2 = KmerScreen(cache_dir=str(tmp_path / "idxB"), mmap=False)
    s2.prepare([str(ref_fa)])
    hash_b = s2.fingerprint_parts()["ref_content_hash"]
    assert hash_b != hash_a


def test_kmer_reference_content_change_invalidates_resume_cache(tmp_path):
    """End-to-end reproduction of the audit's verified unsafe case, through
    the engine the pipeline actually ships by default (kmer): run a screen
    against reference A, copy reference B's CONTENT over the SAME path
    reference A was screened against, run again, and confirm the resume
    cache is invalidated (different fingerprint, fresh chunks) rather than
    silently reused."""
    import json

    import flexprobe.pipeline as pl
    from flexprobe.screen import ENGINE_KMER

    ref_fa = tmp_path / "reference.fa"
    ref_fa.write_text(">refA\n" + ("ACGTTGCA" * 60) + "\n")

    out = tmp_path / "out"
    out.mkdir()
    passing = [_candidate(0)]

    pl._run_screen(passing, [str(ref_fa)], str(out), min_pident=80.0,
                        min_aln=20, resume=True, engine=ENGINE_KMER)
    cdir = screen_cache.cache_dir(str(out))
    fp_a = json.load(open(os.path.join(cdir, "_meta.json")))["fingerprint"]
    assert [f for f in os.listdir(cdir) if f.startswith("chunk_")]

    # Reference B's CONTENT, copied over the SAME path reference A was
    # screened against -- db_arg (the path) never changes.
    ref_fa.write_text(">refB\n" + ("TTGGCCAA" * 60) + "\n")

    pl._run_screen(passing, [str(ref_fa)], str(out), min_pident=80.0,
                        min_aln=20, resume=True, engine=ENGINE_KMER)
    fp_b = json.load(open(os.path.join(cdir, "_meta.json")))["fingerprint"]

    assert fp_b != fp_a, (
        "reference content changed at an unchanged path but the resume-cache "
        "fingerprint didn't -- stale off-target results would be silently "
        "reused, the exact hazard ref_content_hash exists to close"
    )


def test_chunk_round_trip(tmp_path):
    cdir = screen_cache.open_cache(str(tmp_path), "fp1")
    assert screen_cache.load_chunk(cdir, 3) is None
    hits = {"g1_0001": [{"sseqid": "x", "pident": "100"}]}
    screen_cache.save_chunk(cdir, 3, hits)
    assert screen_cache.load_chunk(cdir, 3) == hits


def test_open_cache_resumes_on_match_resets_on_mismatch(tmp_path):
    out = str(tmp_path)
    cdir = screen_cache.open_cache(out, "fp1")
    screen_cache.save_chunk(cdir, 0, {"g1_0000": []})
    # same fingerprint -> cache kept (resume)
    screen_cache.open_cache(out, "fp1")
    assert screen_cache.load_chunk(cdir, 0) is not None
    # different fingerprint -> stale cache wiped
    screen_cache.open_cache(out, "fp2")
    assert screen_cache.load_chunk(cdir, 0) is None


def test_clear_removes_cache(tmp_path):
    out = str(tmp_path)
    screen_cache.open_cache(out, "fp1")
    assert os.path.isdir(screen_cache.cache_dir(out))
    screen_cache.clear(out)
    assert not os.path.isdir(screen_cache.cache_dir(out))


# ---- pipeline resume behaviour ----

def _candidate(i):
    return ProbeCandidate(gene_id="g1", probe_id=f"g1_{i:04d}",
                          sequence="ACGT" * 12 + "AC", start=i, end=i + 50)


def test_blast_resume_skips_completed_chunks(tmp_path, monkeypatch):
    """A second _run_screen over the same candidates does zero blastn
    calls — every chunk is served from the checkpoint cache."""
    import flexprobe.pipeline as pl

    calls = {"n": 0}

    def fake_run_blastn(_blastn, query_fa, _db, **_kw):
        calls["n"] += 1
        tsv = query_fa + ".tsv"
        # Emit one off-target-looking hit for the first probe in the chunk.
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
                    rows.append({"qseqid": p[0], "sseqid": p[1], "pident": p[2],
                                 "length": p[3]})
        return rows

    monkeypatch.setattr(fp_blast, "discover_blast", lambda: {"blastn": "blastn", "makeblastdb": "mkdb"})
    monkeypatch.setattr(fp_blast, "check_version", lambda *_a, **_k: "2.13.0")
    monkeypatch.setattr(pl, "is_db_prefix", lambda _p: True)   # treat as ready DB
    monkeypatch.setattr(fp_blast, "run_blastn", fake_run_blastn)
    monkeypatch.setattr(fp_blast, "parse_hits_tsv", fake_parse)
    # Small chunk size so a few candidates span multiple chunks.
    monkeypatch.setattr("flexprobe.constants.BLAST_QUERY_CHUNK", 2)

    out = tmp_path / "out"
    out.mkdir()
    passing = [_candidate(i) for i in range(5)]   # 5 candidates -> 3 chunks

    first = pl._run_screen(passing, ["/db/prefix"], str(out),
                                 min_pident=80.0, min_aln=20, resume=True)
    assert calls["n"] == 3, "first run should BLAST all 3 chunks"
    assert os.path.isdir(screen_cache.cache_dir(str(out)))

    calls["n"] = 0
    second = pl._run_screen(passing, ["/db/prefix"], str(out),
                                  min_pident=80.0, min_aln=20, resume=True)
    assert calls["n"] == 0, "resume should reuse every cached chunk"
    assert second == first, "resumed hits must match the original"


def test_blast_no_resume_writes_no_cache(tmp_path, monkeypatch):
    import flexprobe.pipeline as pl
    monkeypatch.setattr(fp_blast, "discover_blast", lambda: {"blastn": "b", "makeblastdb": "m"})
    monkeypatch.setattr(fp_blast, "check_version", lambda *_a, **_k: "2.13.0")
    monkeypatch.setattr(pl, "is_db_prefix", lambda _p: True)
    monkeypatch.setattr(fp_blast, "run_blastn", lambda _b, q, _d, **_k: open(q + ".tsv", "w").close() or (q + ".tsv"))
    monkeypatch.setattr(fp_blast, "parse_hits_tsv", lambda _t: [])
    out = tmp_path / "out"; out.mkdir()
    pl._run_screen([_candidate(0)], ["/db/prefix"], str(out), resume=False)
    assert not os.path.isdir(screen_cache.cache_dir(str(out)))


# ---- verdict cache (schema 2) ----
# The streaming path must cache verdicts, not the raw hit dict. At
# whole-transcriptome scale the raw hits run to gigabytes per chunk and can
# fill a disk over a full run.

def _scored(i, *, passed=True):
    c = _candidate(i)
    c.specificity_pass = passed
    c.specificity_status = "passed" if passed else "failed_off_target"
    c.specificity_margin = 12.5
    c.offtarget_nearmiss_count = 3
    c.best_offtarget_pident = 67.5
    c.best_offtarget_bitscore = 40.0
    c.best_offtarget_evalue = 1e-4
    c.isoform_total = 4
    c.isoform_hits = 3
    c.isoform_coverage = 0.75
    c.isoform_hit_ids = ["tx1", "tx2", "tx3"]
    return c


def test_verdict_round_trip_restores_every_field(tmp_path):
    """Every field batch_evaluate() writes must survive save -> load -> apply.
    A field we forget to persist silently reverts to its model default."""
    cdir = screen_cache.open_cache(str(tmp_path), "fp1")
    assert screen_cache.load_verdicts(cdir, 0) is None
    saved = [_scored(0), _scored(1, passed=False)]
    screen_cache.save_verdicts(cdir, 0, saved)

    fresh = [_candidate(0), _candidate(1)]
    n_pass = screen_cache.apply_verdicts(fresh, screen_cache.load_verdicts(cdir, 0))
    assert n_pass == 1
    for orig, restored in zip(saved, fresh):
        for f in screen_cache.VERDICT_FIELDS:
            got, want = getattr(restored, f), getattr(orig, f)
            # NaN is a real value here, not a missing one: offtarget_min_qcov
            # uses it for "there were no non-target hits at all". `nan == nan`
            # is False, so compare NaN-ness explicitly rather than letting a
            # correctly restored sentinel read as a round-trip failure. Both
            # must be NaN — a NaN restored over a real number still fails.
            if isinstance(want, float) and math.isnan(want):
                assert isinstance(got, float) and math.isnan(got), \
                    f"{f} not restored (expected NaN, got {got!r})"
            else:
                assert got == want, f"{f} not restored"


def test_apply_verdicts_rejects_partial_chunk(tmp_path):
    """A chunk missing a probe must re-BLAST, not half-apply — otherwise the
    absent probes keep specificity_pass=False and vanish from the final set."""
    cdir = screen_cache.open_cache(str(tmp_path), "fp1")
    screen_cache.save_verdicts(cdir, 0, [_scored(0)])          # only probe 0
    fresh = [_candidate(0), _candidate(1)]                     # chunk has 2
    assert screen_cache.apply_verdicts(fresh, screen_cache.load_verdicts(cdir, 0)) is None
    # nothing applied — probe 0 must NOT have been mutated
    assert fresh[0].specificity_status == "not_screened"


def test_apply_verdicts_rejects_missing_field(tmp_path):
    cdir = screen_cache.open_cache(str(tmp_path), "fp1")
    screen_cache.save_verdicts(cdir, 0, [_scored(0)])
    v = screen_cache.load_verdicts(cdir, 0)
    del v["g1_0000"]["isoform_hit_ids"]
    assert screen_cache.apply_verdicts([_candidate(0)], v) is None


def test_fingerprint_includes_schema_and_spec_params():
    """Schema in the hash stops v2 code reading a v1 raw-hit cache. spec_params
    in the hash stops verdicts being reused across different thresholds."""
    base = dict(db_arg="db", chunk_size=2, perc_identity=80.0, qcov_hsp_perc=68.0)
    a = screen_cache.fingerprint(_ids(5), **base, spec_params={"off_target_pident_min": 80.0})
    assert a == screen_cache.fingerprint(_ids(5), **base,
                                        spec_params={"off_target_pident_min": 80.0})
    # a different scoring threshold MUST invalidate verdicts
    assert a != screen_cache.fingerprint(_ids(5), **base,
                                        spec_params={"off_target_pident_min": 85.0})
    # verdict-cached vs hit-cached runs must never share a digest
    assert a != screen_cache.fingerprint(_ids(5), **base, spec_params=None)
    # set ordering must not change the digest (stability across runs)
    x = screen_cache.fingerprint(_ids(5), **base, spec_params={"ids": {"a", "b"}})
    y = screen_cache.fingerprint(_ids(5), **base, spec_params={"ids": {"b", "a"}})
    assert x == y


def test_schema_bump_invalidates_old_cache(tmp_path, monkeypatch):
    """The real trap: v1 wrote raw hits under the same digest inputs. If the
    schema weren't hashed, v2 would json.load a v1 chunk and hand raw hits to
    code expecting verdicts."""
    base = dict(db_arg="db", chunk_size=2, perc_identity=80.0, qcov_hsp_perc=68.0)
    v2 = screen_cache.fingerprint(_ids(3), **base)
    monkeypatch.setattr(screen_cache, "CACHE_SCHEMA", 1)
    v1 = screen_cache.fingerprint(_ids(3), **base)
    assert v1 != v2, "schema must participate in the fingerprint"


# ---- streaming: scale invariance + interruption/resume ----

def _streaming_harness(monkeypatch, hits_per_probe, fail_on_chunk=None):
    """Fake blastn emitting `hits_per_probe` off-target hits per query half."""
    import flexprobe.pipeline as pl
    calls = {"n": 0}

    def fake_run_blastn(_b, query_fa, _db, **_kw):
        calls["n"] += 1
        if fail_on_chunk is not None and calls["n"] > fail_on_chunk:
            raise RuntimeError("simulated walltime kill / crash")
        qids = []
        with open(query_fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    qids.append(line[1:].strip())
        tsv = query_fa + ".tsv"
        with open(tsv, "w") as fh:
            for q in qids:
                for j in range(hits_per_probe):
                    fh.write(f"{q}\tsubj{j}\t95.0\t25\t0\t0\t1\t25\t1\t25\t1e-5\t50\n")
        return tsv

    def fake_parse(tsv):
        rows = []
        with open(tsv) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 12:
                    rows.append({"qseqid": p[0], "sseqid": p[1],
                                 "pident": p[2], "length": p[3],
                                 "bitscore": p[11], "evalue": p[10]})
        return rows

    monkeypatch.setattr(fp_blast, "discover_blast",
                        lambda: {"blastn": "blastn", "makeblastdb": "mkdb"})
    monkeypatch.setattr(fp_blast, "check_version", lambda *_a, **_k: "2.13.0")
    monkeypatch.setattr(pl, "is_db_prefix", lambda _p: True)
    monkeypatch.setattr(fp_blast, "run_blastn", fake_run_blastn)
    monkeypatch.setattr(fp_blast, "parse_hits_tsv", fake_parse)
    monkeypatch.setattr("flexprobe.constants.BLAST_QUERY_CHUNK", 2)
    return calls


def _run_streaming(out, passing, monkeypatch, hits_per_probe, fail_on_chunk=None,
                   cache_mode="verdicts"):
    """Drive the real streaming path; returns (blastn call count, pass tally)."""
    import flexprobe.pipeline as pl
    from flexprobe.specificity import batch_evaluate
    calls = _streaming_harness(monkeypatch, hits_per_probe, fail_on_chunk)
    tally = {"passed": 0}
    spec_kwargs = dict(off_target_pident_min=80.0, off_target_aln_min=20,
                       on_target_pident_min=90.0, on_target_aln_min=21)

    def sink(ci, cands, hits):
        if hits is None:                       # resumed -> tally only
            tally["passed"] += sum(1 for c in cands if c.specificity_pass)
            return
        tally["passed"] += batch_evaluate(cands, hits, set(), **spec_kwargs)

    pl._run_screen(
        passing, ["/db/prefix"], str(out),
        min_pident=80.0, min_aln=20, resume=True,
        per_chunk=sink, spec_params={**spec_kwargs, "target_ids": set()},
        cache_mode=cache_mode,
    )
    return calls["n"], tally["passed"]


def _cache_bytes(out):
    cdir = screen_cache.cache_dir(str(out))
    return sum(os.path.getsize(os.path.join(cdir, f))
               for f in os.listdir(cdir) if f.startswith("chunk_"))


def test_verdict_cache_size_is_independent_of_hit_volume(tmp_path, monkeypatch):
    """Cache size must not grow with the number of hits per probe.

    A cache of the raw hit dict scales linearly with hits/probe, which at
    whole-transcriptome scale (hundreds of hits per query) means gigabytes per
    chunk. Verdicts are a fixed ~11 fields per candidate, so 200x the hits must
    NOT mean a bigger cache.
    """
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir(); b.mkdir()
    _run_streaming(a, [_candidate(i) for i in range(4)], monkeypatch, hits_per_probe=1)
    small = _cache_bytes(a)
    _run_streaming(b, [_candidate(i) for i in range(4)], monkeypatch, hits_per_probe=200)
    big = _cache_bytes(b)

    # 200x the hits, same verdicts -> essentially identical cache footprint.
    assert big == small, (
        f"verdict cache grew with hit volume ({small} -> {big} bytes): the cache "
        "is storing hits again, which is exactly what filled a 460 GB disk"
    )


def test_interrupted_streaming_run_resumes_to_identical_result(tmp_path, monkeypatch):
    """Interrupt mid-run, resume, and land on the same answer as an
    uninterrupted run — with zero re-BLAST of completed chunks."""
    ref, out = tmp_path / "ref", tmp_path / "out"
    ref.mkdir(); out.mkdir()

    # baseline: clean run, 5 candidates -> 3 chunks
    ref_calls, ref_passed = _run_streaming(
        ref, [_candidate(i) for i in range(5)], monkeypatch, hits_per_probe=3)
    assert ref_calls == 3

    # interrupted: blow up after chunk 2
    import pytest
    with pytest.raises(RuntimeError):
        _run_streaming(out, [_candidate(i) for i in range(5)], monkeypatch,
                       hits_per_probe=3, fail_on_chunk=2)
    assert len(os.listdir(screen_cache.cache_dir(str(out)))) >= 2  # 2 chunks + meta

    # resume: only the unfinished chunk should BLAST
    calls, passed = _run_streaming(
        out, [_candidate(i) for i in range(5)], monkeypatch, hits_per_probe=3)
    assert calls == 1, f"resume re-BLASTed {calls} chunks; expected only the lost one"
    assert passed == ref_passed, "resumed tally must match the uninterrupted run"


# ---- screen_cache_mode knob ----

def test_hits_mode_cache_scales_with_hit_volume(tmp_path, monkeypatch):
    """The knob's whole point: `hits` trades disk for re-scorability. Assert the
    cost is real, so nobody reaches for it at transcriptome scale by accident."""
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir(); b.mkdir()
    _run_streaming(a, [_candidate(i) for i in range(4)], monkeypatch,
                   hits_per_probe=1, cache_mode="hits")
    small = _cache_bytes(a)
    _run_streaming(b, [_candidate(i) for i in range(4)], monkeypatch,
                   hits_per_probe=200, cache_mode="hits")
    big = _cache_bytes(b)
    assert big > small * 50, (
        "hits mode should grow with hit volume; if it doesn't, the modes are "
        "wired the same way and the knob is a lie"
    )


def test_hits_mode_rescoring_is_free_verdicts_mode_is_not(tmp_path):
    """`hits` keeps spec_params OUT of the fingerprint so a threshold change
    re-scores from cache; `verdicts` puts it IN so the same change re-BLASTs.
    This asymmetry IS the trade-off the knob exposes."""
    base = dict(db_arg="db", chunk_size=2, perc_identity=80.0, qcov_hsp_perc=68.0)
    # verdicts: threshold change must change the digest (=> re-BLAST)
    v1 = screen_cache.fingerprint(_ids(3), **base, cache_mode="verdicts",
                                 spec_params={"off_target_pident_min": 80.0})
    v2 = screen_cache.fingerprint(_ids(3), **base, cache_mode="verdicts",
                                 spec_params={"off_target_pident_min": 85.0})
    assert v1 != v2
    # hits: spec_params is omitted by the caller, so the digest is unchanged
    # across thresholds => the cache is reused and re-scored for free.
    h1 = screen_cache.fingerprint(_ids(3), **base, cache_mode="hits")
    h2 = screen_cache.fingerprint(_ids(3), **base, cache_mode="hits")
    assert h1 == h2
    # and the two modes must never share a cache
    assert h1 != v1


def test_switching_mode_invalidates_cache(tmp_path, monkeypatch):
    """Flipping the knob changes the on-disk format. The fingerprint must
    invalidate the old cache rather than misread it as the other shape."""
    out = tmp_path / "out"; out.mkdir()
    cands = [_candidate(i) for i in range(4)]
    _run_streaming(out, cands, monkeypatch, hits_per_probe=3, cache_mode="verdicts")
    import json
    cdir = screen_cache.cache_dir(str(out))
    first = json.load(open(os.path.join(cdir, "_meta.json")))["fingerprint"]

    calls, _ = _run_streaming(out, [_candidate(i) for i in range(4)], monkeypatch,
                              hits_per_probe=3, cache_mode="hits")
    second = json.load(open(os.path.join(cdir, "_meta.json")))["fingerprint"]
    assert second != first, "mode must participate in the fingerprint"
    assert calls > 0, "switching mode must re-BLAST, not reuse the old format"


def test_invalid_mode_is_rejected(tmp_path, monkeypatch):
    import pytest
    out = tmp_path / "out"; out.mkdir()
    with pytest.raises(ValueError, match="screen_cache_mode"):
        _run_streaming(out, [_candidate(0)], monkeypatch, hits_per_probe=1,
                       cache_mode="raw")


# ---- HPC core detection ----

def test_usable_cores_respects_scheduler_allocation(monkeypatch):
    """os.cpu_count() reports the NODE, ignoring cgroup/cpuset limits: a PBS job
    with ncpus=16 on a 64-core node would hand blastn 63 threads on 16 cores."""
    import flexprobe.pipeline as pl

    # Python 3.13+: process_cpu_count honours affinity -> must win.
    monkeypatch.setattr(pl.os, "cpu_count", lambda: 64)
    monkeypatch.setattr(pl.os, "process_cpu_count", lambda: 16, raising=False)
    assert pl._usable_cores() == 16

    # Older Python on Linux: fall back to sched_getaffinity (cpuset-aware).
    monkeypatch.delattr(pl.os, "process_cpu_count", raising=False)
    monkeypatch.setattr(pl.os, "sched_getaffinity", lambda _pid: set(range(16)),
                        raising=False)
    assert pl._usable_cores() == 16

    # macOS/Windows: neither exists -> cpu_count is all we have.
    monkeypatch.delattr(pl.os, "sched_getaffinity", raising=False)
    assert pl._usable_cores() == 64


def test_usable_cores_never_returns_zero(monkeypatch):
    import flexprobe.pipeline as pl
    monkeypatch.delattr(pl.os, "process_cpu_count", raising=False)
    monkeypatch.delattr(pl.os, "sched_getaffinity", raising=False)
    monkeypatch.setattr(pl.os, "cpu_count", lambda: None)
    assert pl._usable_cores() >= 1


def test_fingerprint_records_the_EFFECTIVE_mode(tmp_path, monkeypatch):
    """Non-streaming always caches hits, even when `verdicts` is requested.
    The fingerprint must describe what's actually on disk, or a later streaming
    run could match the digest and try to read hits as verdicts."""
    import flexprobe.pipeline as pl
    import json
    _streaming_harness(monkeypatch, hits_per_probe=2)
    out = tmp_path / "out"; out.mkdir()
    # non-streaming (no per_chunk) but asking for verdicts
    pl._run_screen([_candidate(i) for i in range(4)], ["/db/prefix"],
                         str(out), min_pident=80.0, min_aln=20, resume=True,
                         cache_mode="verdicts")
    cdir = screen_cache.cache_dir(str(out))
    non_streaming_fp = json.load(open(os.path.join(cdir, "_meta.json")))["fingerprint"]

    # a genuine `hits` streaming run over the same candidates must agree, since
    # both actually hold hits; a `verdicts` run must NOT.
    out2 = tmp_path / "out2"; out2.mkdir()
    _run_streaming(out2, [_candidate(i) for i in range(4)], monkeypatch,
                   hits_per_probe=2, cache_mode="verdicts")
    verdict_fp = json.load(open(os.path.join(
        screen_cache.cache_dir(str(out2)), "_meta.json")))["fingerprint"]
    assert non_streaming_fp != verdict_fp, (
        "a hits-on-disk cache shares a digest with a verdicts cache — a "
        "streaming run could match it and misread hits as verdicts"
    )


def test_VERDICT_FIELDS_covers_everything_the_scorer_writes():
    """Drift guard, checked against the SOURCE rather than a hand-kept list.

    A field that specificity.py sets but screen_cache doesn't persist would
    silently revert to its ProbeCandidate default on resume — probes quietly
    losing specificity_pass, isoform coverage, etc. That's the worst failure
    mode here: no crash, just wrong results in a 127-hour run. If you add a
    `candidate.<x> = ...` to specificity.py, add it to VERDICT_FIELDS too.
    """
    import ast
    import pathlib
    src_file = pathlib.Path(screen_cache.__file__).parent / "specificity.py"
    tree = ast.parse(src_file.read_text())
    written = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        for tgt in node.targets:
            if (isinstance(tgt, ast.Attribute)
                    and isinstance(tgt.value, ast.Name)
                    and tgt.value.id in ("candidate", "c")):
                written.add(tgt.attr)
    missing = written - set(screen_cache.VERDICT_FIELDS)
    assert not missing, (
        f"specificity.py writes {sorted(missing)} but screen_cache.VERDICT_FIELDS "
        "does not persist them — they will silently revert to model defaults on "
        "resume. Add them to VERDICT_FIELDS."
    )


def test_verdict_cache_requires_spec_params(tmp_path, monkeypatch):
    """Caching verdicts without keying them to the scoring config would yield a
    digest identical to a non-streaming (hits-on-disk) run at the same defaults:
    two on-disk formats under one fingerprint."""
    import flexprobe.pipeline as pl
    import pytest
    _streaming_harness(monkeypatch, hits_per_probe=1)
    out = tmp_path / "out"; out.mkdir()
    with pytest.raises(ValueError, match="spec_params is required"):
        pl._run_screen([_candidate(0)], ["/db/prefix"], str(out),
                             min_pident=80.0, min_aln=20, resume=True,
                             per_chunk=lambda ci, c, h: None,
                             cache_mode="verdicts")   # no spec_params
