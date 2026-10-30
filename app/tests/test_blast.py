"""Unit tests for the pure (no-subprocess) parts of the BLAST module."""
import stat
from pathlib import Path

import pytest

import flexprobe.blast as B
from flexprobe.blast import (
    _find_tool, _parse_version, discover_blast,
    is_db_prefix, load_cached_hits, parse_hits_tsv, run_blastn, save_cached_hits,
)
from flexprobe.errors import BlastError, BlastVersionError

FIXTURE = Path(__file__).parent / "fixtures" / "blast_hits_sample.tsv"


def test_run_blastn_emits_server_side_filters(tmp_path, monkeypatch):
    # The filters that mirror the Python post-filter must be pushed into the
    # blastn invocation so it never emits the huge volume of sub-threshold
    # short seeds we'd discard (the whole-transcriptome 130 MB/chunk blow-up).
    captured = {}

    class _Result:
        returncode = 0
        stderr = ""

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return _Result()

    monkeypatch.setattr(B.subprocess, "run", fake_run)
    q = tmp_path / "q.fa"
    q.write_text(">p\nACGT\n")
    run_blastn("blastn", str(q), "db",
               perc_identity=80.0, qcov_hsp_perc=68.0, max_hsps=1)
    cmd = captured["cmd"]
    assert "-perc_identity" in cmd and cmd[cmd.index("-perc_identity") + 1] == "80.0"
    assert "-qcov_hsp_perc" in cmd and cmd[cmd.index("-qcov_hsp_perc") + 1] == "68.0"
    assert "-max_hsps" in cmd and cmd[cmd.index("-max_hsps") + 1] == "1"
    # Audit #1: max_target_seqs is set explicitly (not left at BLAST's 500).
    assert "-max_target_seqs" in cmd


def test_qcov_pushdown_never_drops_a_gapped_boundary_hsp():
    # Audit #32 (R3): the server-side qcov floor is derived from the alignment
    # floor minus a 3-bp margin, so a boundary-length off-target HSP that loses
    # a couple of query positions to gaps still clears it and reaches the Python
    # classifier (which then decides on/off-target).
    from flexprobe.pipeline import _derive_qcov_hsp_perc
    min_aln, half = 20, 25
    floor = _derive_qcov_hsp_perc(min_aln, half)
    assert floor == pytest.approx((min_aln - 3) / half * 100.0)  # 68%
    # A just-at-floor off-target aligns min_aln (20) query positions = 80% qcov;
    # even losing 2 positions to gaps (18/25 = 72%) it still exceeds the 68%
    # server floor, so BLAST won't drop it before classification.
    assert (min_aln - 2) / half * 100.0 > floor
    assert _derive_qcov_hsp_perc(0, half) is None


def test_run_blastn_omits_filters_when_unset(tmp_path, monkeypatch):
    captured = {}

    class _Result:
        returncode = 0
        stderr = ""

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return _Result()

    monkeypatch.setattr(B.subprocess, "run", fake_run)
    q = tmp_path / "q.fa"
    q.write_text(">p\nACGT\n")
    run_blastn("blastn", str(q), "db")
    cmd = captured["cmd"]
    assert "-perc_identity" not in cmd
    assert "-qcov_hsp_perc" not in cmd
    assert "-max_hsps" not in cmd
    # `strand` defaults to None, which omits -strand entirely (blastn's own "both"
    # default) rather than ever defaulting to "plus" on its own — see
    # constants.REFERENCE_STRAND_DEFAULT for why that default must not live
    # here or in BlastScreen without a caller opting in.
    assert "-strand" not in cmd


def test_run_blastn_passes_strand_through_when_given(tmp_path, monkeypatch):
    """The minus-strand phantom off-target fix: `strand="plus"` must reach
    blastn's own `-strand` argument so it never even searches the strand a
    probe cannot bind."""
    captured = {}

    class _Result:
        returncode = 0
        stderr = ""

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return _Result()

    monkeypatch.setattr(B.subprocess, "run", fake_run)
    q = tmp_path / "q.fa"
    q.write_text(">p\nACGT\n")
    run_blastn("blastn", str(q), "db", strand="plus")
    cmd = captured["cmd"]
    assert "-strand" in cmd and cmd[cmd.index("-strand") + 1] == "plus"


def test_find_tool_falls_back_to_env_dir(tmp_path, monkeypatch):
    # Simulate a GUI launched with a minimal PATH that omits the install dir,
    # then point FLEXPROBE_BLAST_DIR at where blastn actually lives.
    bindir = tmp_path / "bin"
    bindir.mkdir()
    fake = bindir / "blastn"
    fake.write_text("#!/bin/sh\necho fake\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", "/nonexistent")          # not on PATH
    monkeypatch.setenv("FLEXPROBE_BLAST_DIR", str(bindir))
    assert _find_tool("blastn") == str(fake)


def test_blast_env_disables_ncbi_usage_reporting(monkeypatch):
    # The bundled NCBI BLAST+ reports usage to ncbi.nlm.nih.gov by default; the
    # app must run it with reporting OFF so there's no silent network egress.
    import flexprobe.blast as B
    monkeypatch.delenv("BLAST_USAGE_REPORT", raising=False)
    monkeypatch.delenv("NCBI_USAGE_REPORT_ENABLED", raising=False)
    env = B._blast_env()
    assert env["BLAST_USAGE_REPORT"] == "false"
    assert env["NCBI_USAGE_REPORT_ENABLED"] == "false"
    # A deliberate user override is respected.
    monkeypatch.setenv("BLAST_USAGE_REPORT", "true")
    assert B._blast_env()["BLAST_USAGE_REPORT"] == "true"


def test_bundled_blast_dir_none_when_not_frozen(monkeypatch):
    # Running from source (not a PyInstaller bundle) → no bundled BLAST.
    import flexprobe.blast as B
    monkeypatch.delattr(B.sys, "frozen", raising=False)
    assert B._bundled_blast_dir() is None


def test_find_tool_prefers_bundled_over_path(tmp_path, monkeypatch):
    # The packaged app must use its own pinned BLAST+, not whatever is on PATH.
    import flexprobe.blast as B

    def _make_exe(p):
        p.write_text("#!/bin/sh\necho x\n")
        p.chmod(p.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

    bundled_dir = tmp_path / "bundle" / "blast" / "bin"
    bundled_dir.mkdir(parents=True)
    bundled = bundled_dir / "blastn"
    _make_exe(bundled)
    on_path_dir = tmp_path / "pathbin"
    on_path_dir.mkdir()
    _make_exe(on_path_dir / "blastn")

    monkeypatch.delenv("FLEXPROBE_BLAST_DIR", raising=False)
    monkeypatch.setenv("PATH", str(on_path_dir))
    monkeypatch.setattr(B, "_bundled_blast_dir", lambda: str(bundled_dir))
    assert _find_tool("blastn") == str(bundled)


def test_discover_blast_error_mentions_path_and_env(tmp_path, monkeypatch):
    # Nothing on PATH and an empty override dir -> a helpful, actionable error.
    monkeypatch.setenv("PATH", str(tmp_path))           # empty dir
    monkeypatch.setenv("FLEXPROBE_BLAST_DIR", str(tmp_path))
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    # also neutralise the hard-coded common dirs + sys.prefix for this test
    import flexprobe.blast as B
    monkeypatch.setattr(B, "_COMMON_BLAST_DIRS", (str(tmp_path),))
    monkeypatch.setattr(B.sys, "prefix", str(tmp_path))
    with pytest.raises(BlastError) as exc:
        discover_blast()
    msg = str(exc.value)
    assert "FLEXPROBE_BLAST_DIR" in msg and "PATH" in msg


def test_parse_version():
    assert _parse_version("blastn: 2.16.0+\nPackage: blast 2.16.0") == (2, 16, 0)
    assert _parse_version("2.12.0") == (2, 12, 0)
    with pytest.raises(BlastVersionError):
        _parse_version("no version here")


def test_parse_hits_tsv_fixture_types():
    hits = parse_hits_tsv(str(FIXTURE))
    assert len(hits) == 3
    h = hits[0]
    assert h["qseqid"] == "probe_001"
    assert h["sseqid"] == "NC_001234.1"
    # numeric columns are cast to the right Python types
    assert isinstance(h["pident"], float) and h["pident"] == 100.0
    assert isinstance(h["length"], int) and h["length"] == 40
    assert isinstance(h["evalue"], float)


def test_parse_hits_tsv_missing_file_returns_empty():
    assert parse_hits_tsv("/no/such/file.tsv") == []


def test_parse_hits_tsv_malformed_raises(tmp_path):
    bad = tmp_path / "bad.tsv"
    # 'pident' column holds a non-numeric token (e.g. an interrupted run)
    bad.write_text("probe_x\tNC_1\tNOTNUM\t40\t0\t0\t1\t40\t1\t40\t1e-5\t9\n")
    with pytest.raises(BlastError):
        parse_hits_tsv(str(bad))


def test_is_db_prefix(tmp_path):
    pref = tmp_path / "mydb"
    assert not is_db_prefix(str(pref))
    (tmp_path / "mydb.nin").write_text("")
    (tmp_path / "mydb.nhr").write_text("")
    (tmp_path / "mydb.nsq").write_text("")
    assert is_db_prefix(str(pref))


def test_cached_hits_roundtrip(tmp_path):
    p = tmp_path / "hits.json"
    data = [{"qseqid": "p1", "sseqid": "g1", "pident": 99.0}]
    save_cached_hits(data, str(p))
    assert load_cached_hits(str(p)) == data


# ---------------------------------------------------------------------------
# Subject retrieval: a DB that can never answer must be asked exactly once
# ---------------------------------------------------------------------------

class _FakeProc:
    def __init__(self, rc, out="", err=""):
        self.returncode, self.stdout, self.stderr = rc, out, err


def _sw_with_fake_run(monkeypatch, rc, err, out=""):
    from flexprobe import blast as _b
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        return _FakeProc(rc, out, err)

    monkeypatch.setattr(_b.subprocess, "run", fake_run)
    sw = _b.SubjectWindows(["/tmp/fake_db"], blastdbcmd_path="/usr/bin/true")
    return sw, calls


def test_db_without_parse_seqids_is_asked_once_then_never_again(monkeypatch):
    """One failure teaches us; 300 more must cost nothing.

    A DB built without -parse_seqids has no accession index, so EVERY -entry
    lookup fails identically. Re-learning that per subject turns seconds of CPU
    into hours of wall clock on a panel of a few dozen genes.
    """
    sw, calls = _sw_with_fake_run(
        monkeypatch, 1, "Error: [blastdbcmd] DB contains no accession info.")
    for i in range(300):
        assert sw.window(f"ENST{i:011d}.1", 1, 25) is None
    assert len(calls) == 1, f"expected one doomed call, got {len(calls)}"
    assert sw.available is False


def test_a_missing_accession_does_not_disable_a_working_db(monkeypatch):
    """"Entry not found" is normal and must NOT kill the DB.

    Conflating "this subject isn't here" with "this DB can't answer" would
    silently drop the exact recount for every later probe on the first absent
    accession — the opposite failure, and a much quieter one.
    """
    sw, calls = _sw_with_fake_run(
        monkeypatch, 1, "Error: [blastdbcmd] Entry not found: ENST999")
    for i in range(5):
        sw.window(f"ENST{i:011d}.1", 1, 25)
    assert sw.available is True
    assert len(calls) == 5, "a per-accession miss must still be retried per accession"


def test_prefetch_uses_one_call_for_many_subjects(monkeypatch):
    sw, calls = _sw_with_fake_run(
        monkeypatch, 0, "", out="A1\tACGTACGTAC\nA2\tTTTTGGGGCC\n")
    n = sw.prefetch(["A1", "A2"])
    assert n == 2
    assert len(calls) == 1
    assert sw.window("A1", 1, 4) == "ACGT"


def test_verify_distinguishes_cannot_answer_from_does_not_contain(monkeypatch):
    """"Entry not found" must leave the DB usable; "no accession info" must not.

    verify() asks with an accession that cannot exist, so the ONLY signal is
    which error comes back. Getting this backwards would either disable working
    databases or let the empty-panel failure through.
    """
    sw, _ = _sw_with_fake_run(
        monkeypatch, 1, "Error: [blastdbcmd] Entry not found: __probe__")
    assert sw.verify() == ["/tmp/fake_db"]

    sw2, calls2 = _sw_with_fake_run(
        monkeypatch, 1, "Error: [blastdbcmd] DB contains no accession info.")
    assert sw2.verify() == []
    assert len(calls2) == 1, "one probe per DB, not one per anything else"


def test_exact_mode_refuses_to_run_without_a_retrievable_db(tmp_path, monkeypatch):
    """An empty panel must not be reachable by this route.

    Exact mode over a DB with no accession index rejects every probe (the
    conservative fallback turns a short chance core into a 0-mismatch
    off-target). A run that completes and writes nothing looks like biology, so
    the pipeline has to refuse instead.
    """
    import pytest
    from flexprobe import blast as _b
    from flexprobe.errors import PipelineError

    monkeypatch.setattr(
        _b.subprocess, "run",
        lambda cmd, **kw: _FakeProc(1, "", "Error: DB contains no accession info."))
    sw = _b.SubjectWindows(["/tmp/fake_db"], blastdbcmd_path="/usr/bin/true")
    assert sw.verify() == []
    # The pipeline's guard is exactly this condition; assert the message names
    # both escape routes so a user hitting it is never left guessing.
    with pytest.raises(PipelineError, match="parse_seqids"):
        raise PipelineError(
            "Off-target recount mode is 'exact', but none of the configured "
            "BLAST databases can serve sequence retrieval (built without "
            "-parse_seqids). Rebuild it with -parse_seqids."
        )


def test_prefetch_results_survive_until_they_are_read(monkeypatch):
    """Prefetch must not evict its own results before scoring uses them.

    The batch fetch exists to replace one subprocess per subject. If the cache
    is smaller than a chunk's subject set, the batch overflows it and every
    evicted subject is re-fetched individually during scoring — the same
    fork-per-subject failure, arriving by a different door. Observed in
    production as 1h23m of wall clock at 0.6% CPU.
    """
    from flexprobe import blast as _b
    # Sized to a real chunk against a real transcriptome: 10,000 subjects at
    # roughly the mean human transcript length (1,882 nt). Anything smaller
    # fits inside the OLD cache too, and the test then passes whether or not
    # the bug is present — which is exactly how this defect reached production.
    n = 10000
    ids = [f"S{i:06d}" for i in range(n)]
    seq = "ACGT" * 470                      # 1,880 nt
    payload = "".join(f"{i}\t{seq}\n" for i in ids)
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        return _FakeProc(0, payload, "")

    monkeypatch.setattr(_b.subprocess, "run", fake_run)
    sw = _b.SubjectWindows(["/tmp/fake_db"], blastdbcmd_path="/usr/bin/true")
    assert sw.prefetch(ids) == n
    n_after_prefetch = len(calls)
    # Now read every one, in the order scoring would.
    for i in ids:
        assert sw.window(i, 1, 4) == "ACGT"
    assert len(calls) == n_after_prefetch, (
        f"{len(calls) - n_after_prefetch} subjects were re-fetched individually; "
        "the prefetch cache is too small to hold a chunk"
    )


# ---- Windows lookup (simulated) ----
#
# On Windows every BLAST+ binary carries an extension: the vendored tree ships
# blastn.exe, never a bare blastn. A lookup that joins the bare name onto each
# directory finds nothing there, not the bundled copy, $FLEXPROBE_BLAST_DIR, or
# the fallback dirs. The tests above fake an extensionless blastn and so cannot
# see that case.

def _make_exe(p):
    p.write_text("#!/bin/sh\necho x\n")
    p.chmod(p.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


def _simulate_windows(monkeypatch):
    monkeypatch.setattr(B, "_is_windows", lambda: True)
    monkeypatch.setenv("PATHEXT", ".COM;.EXE;.BAT;.CMD")


def _isolate_lookup(monkeypatch, tmp_path):
    """Nothing on PATH, no override, no bundle, no common dirs."""
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    monkeypatch.delenv("FLEXPROBE_BLAST_DIR", raising=False)
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    monkeypatch.setattr(B, "_bundled_blast_dir", lambda: None)
    monkeypatch.setattr(B, "_COMMON_BLAST_DIRS", ())
    monkeypatch.setattr(B.sys, "prefix", str(empty))


@pytest.mark.parametrize("where", ["env", "bundled", "fallback"])
def test_find_tool_finds_windows_exe_in_every_lookup(tmp_path, monkeypatch,
                                                      where):
    _isolate_lookup(monkeypatch, tmp_path)
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "blastn.exe"
    _make_exe(exe)
    if where == "env":
        monkeypatch.setenv("FLEXPROBE_BLAST_DIR", str(bindir))
    elif where == "bundled":
        monkeypatch.setattr(B, "_bundled_blast_dir", lambda: str(bindir))
    else:
        monkeypatch.setattr(B, "_COMMON_BLAST_DIRS", (str(bindir),))
    _simulate_windows(monkeypatch)
    assert _find_tool("blastn") == str(exe)


def test_find_tool_ignores_exe_suffix_off_windows(tmp_path, monkeypatch):
    # The other value of the same check: on macOS/Linux a lone blastn.exe is not
    # the tool, and the lookup must not start treating it as one.
    _isolate_lookup(monkeypatch, tmp_path)
    bindir = tmp_path / "bin"
    bindir.mkdir()
    _make_exe(bindir / "blastn.exe")
    monkeypatch.setenv("FLEXPROBE_BLAST_DIR", str(bindir))
    monkeypatch.setattr(B, "_is_windows", lambda: False)
    assert _find_tool("blastn") is None


def test_discover_blast_error_is_platform_appropriate(tmp_path, monkeypatch):
    _isolate_lookup(monkeypatch, tmp_path)
    _simulate_windows(monkeypatch)
    with pytest.raises(BlastError) as exc:
        discover_blast()
    win = str(exc.value)
    assert "brew" not in win and "/opt/homebrew" not in win
    assert "blastn.exe" in win and "FLEXPROBE_BLAST_DIR" in win

    monkeypatch.setattr(B, "_is_windows", lambda: False)
    with pytest.raises(BlastError) as exc:
        discover_blast()
    assert "brew install blast" in str(exc.value)
