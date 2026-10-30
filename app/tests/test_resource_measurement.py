"""An absent measurement must be written as null, never as 0.0.

Every Windows run in the lab-meeting dry launch recorded ``cpu_s``,
``cores_mean``, ``peak_rss_mb`` and the tree peak as 0.0. Two causes: the
``resource`` module does not exist there, and the tree sampler shells out to
``ps``, which is not a Windows command. Both failures were then written as the
number zero, and the tree report's ``method`` still said it had summed ps RSS.
A zero reads as a measurement; it was the absence of one.

The "nothing measurable" tests reproduce that state on any platform by removing
``resource`` and emptying PATH, with the Windows back end switched off. The
Windows back end itself cannot be exercised off Windows, so its tree walk is
tested against a fake process table and its loader is tested for failing
quietly.
"""
import json
import os
import time

import pytest

import flexprobe.pipeline as P
from flexprobe.pipeline import _StageTimer, _TreeRssSampler


def _no_measurement(monkeypatch, tmp_path):
    """No resource module, no ps on PATH, no Windows back end."""
    monkeypatch.setattr(P, "resource", None)
    monkeypatch.setattr(P, "_is_windows", lambda: False, raising=False)
    empty = tmp_path / "empty-path"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))


def test_stage_cpu_and_memory_are_null_when_unmeasurable(monkeypatch,
                                                         tmp_path):
    _no_measurement(monkeypatch, tmp_path)
    t = _StageTimer()
    t.enter("tile")
    t.enter("filter")
    res = t.resources()
    for name in ("tile", "filter", "total"):
        assert res[name]["cpu_s"] is None, name
        assert res[name]["cores_mean"] is None, name
        assert res[name]["peak_rss_mb"] is None, name
        assert res[name]["wall_s"] >= 0.0          # wall clock always works


def test_timeline_step_cpu_is_null_when_unmeasurable(monkeypatch, tmp_path):
    _no_measurement(monkeypatch, tmp_path)
    sink = tmp_path / "run_timeline.jsonl"
    t = _StageTimer()
    t.open_sink(str(sink))
    try:
        t.enter("tile")
        t.enter("filter")
        t.enter("rank")
    finally:
        t._tree_rss.stop()
    recs = [json.loads(line) for line in sink.read_text().splitlines()]
    assert recs
    for r in recs:
        assert r["step_cpu_s"] is None
        assert r["peak_rss_mb"] is None


def test_tree_sampler_without_ps_reports_null_and_does_not_claim_ps(
        monkeypatch, tmp_path):
    _no_measurement(monkeypatch, tmp_path)
    s = _TreeRssSampler(interval_s=0.05)
    s.start()
    time.sleep(0.3)
    s.stop()
    rep = s.report()
    assert rep["samples"] == 0 and rep["failed_samples"] > 0
    assert rep["peak_tree_rss_mib"] is None
    assert rep["peak_tree_rss_gib"] is None
    assert rep["processes_at_peak"] is None
    assert rep.get("available") is False and rep.get("reason")
    assert "sum of ps RSS" not in (rep["method"] or "")


def test_run_timings_json_writes_null_with_availability_flags(monkeypatch,
                                                              tmp_path):
    _no_measurement(monkeypatch, tmp_path)
    t = _StageTimer()
    t.open_sink(str(tmp_path / "run_timeline.jsonl"))
    t.enter("tile")
    t.enter("filter")
    path = tmp_path / "run_timings.json"
    t.write_report(str(path))
    data = json.loads(path.read_text())

    assert data["stages"]["total"]["peak_rss_mb"] is None
    assert data["stages"]["total"]["cpu_s"] is None
    m = data.get("measurement", {})
    for what in ("cpu", "peak_rss"):
        assert m.get(what, {}).get("available") is False, what
        assert m[what]["reason"], what
        assert m[what]["method"] is None, what
    tree = data["peak_tree_rss"]
    assert tree["available"] is False
    assert tree["parent_exact_rss_mib"] is None
    assert tree["recommended_floor_mib"] is None
    assert tree["recommended_floor_gib"] is None
    assert tree["which_is_higher"] is None


def test_measured_platform_still_reports_numbers(tmp_path):
    # The other value of every check above: where resource and ps exist, the
    # same report carries real numbers and says how it got them.
    if P.resource is None:
        pytest.skip("needs the POSIX resource module and ps")
    t = _StageTimer()
    t.open_sink(str(tmp_path / "run_timeline.jsonl"))
    t.enter("tile")
    time.sleep(0.5)
    t.enter("filter")
    path = tmp_path / "run_timings.json"
    t.write_report(str(path))
    data = json.loads(path.read_text())

    assert isinstance(data["stages"]["tile"]["cpu_s"], float)
    assert isinstance(data["stages"]["tile"]["cores_mean"], float)
    assert data["stages"]["total"]["peak_rss_mb"] > 0.0
    assert data["measurement"]["cpu"]["available"] is True
    assert data["measurement"]["peak_rss"]["available"] is True
    tree = data["peak_tree_rss"]
    assert tree["available"] is True and tree["samples"] > 0
    assert tree["peak_tree_rss_mib"] > 0.0
    assert "sum of ps RSS" in tree["method"]
    assert tree["recommended_floor_mib"] >= tree["parent_exact_rss_mib"]


# ---- Windows back end ----

@pytest.mark.skipif(os.name == "nt", reason="exercises the non-Windows path")
def test_windows_probe_loads_only_on_windows_and_never_raises(monkeypatch):
    monkeypatch.setattr(P, "_is_windows", lambda: False)
    assert P._windows_probe() is None

    # Claimed Windows, but there is no kernel32 here: the loader must swallow
    # that, and the timer must degrade to null rather than crash.
    monkeypatch.setattr(P, "_is_windows", lambda: True)
    P._load_windows_probe.cache_clear()
    try:
        assert P._windows_probe() is None
        monkeypatch.setattr(P, "resource", None)
        assert _StageTimer._cpu_now() is None
        assert _StageTimer._peak_rss_mb() is None
    finally:
        P._load_windows_probe.cache_clear()


def test_windows_tree_walk_sums_descendants_and_skips_reused_parent_pids():
    mib = 1024 * 1024
    procs = [(10, 1), (11, 10), (12, 11), (13, 10), (99, 10), (50, 2)]
    # pid -> (creation time, working set bytes)
    info = {
        10: (100, 5 * mib),
        11: (110, 3 * mib),
        12: (120, 2 * mib),                 # grandchild, reached through 11
        13: (130, 1 * mib),
        # Lists 10 as its parent but was created BEFORE 10 existed: a Windows
        # parent pid is never cleared, so this is an orphan of an earlier
        # holder of pid 10, not a descendant.
        99: (50, 7 * mib),
        50: (10, 9 * mib),                  # unrelated process
    }
    probe = object.__new__(P._WinProcessProbe)
    probe._processes = lambda: procs
    probe._probe = lambda pid: info.get(pid)
    assert probe.tree_rss(10) == (11.0, 4)

    probe._probe = lambda pid: None         # cannot open even the root
    assert probe.tree_rss(10) is None
    probe._processes = lambda: None         # snapshot failed
    assert probe.tree_rss(10) is None
