"""results.csv loading: off the UI thread, bounded in memory, cancellable.

A whole-transcriptome results.csv can exceed a gigabyte. Parsed inline on the
UI thread into an unbounded cache, it would freeze the app or run the machine
out of memory, and every run opened in a session would stay resident until
the process exited.

These tests assert both halves of the design: the parse really happens on a
worker QThread (and its signals really arrive), and the cache really evicts.
"""
import csv
import importlib.util
import json
import os
import pathlib
import threading

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


# ---------------------------------------------------------------------------
# The cache — no Qt needed
# ---------------------------------------------------------------------------

def _tiny_run(out_dir, n_rows=3, pad=""):
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "results.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["probe_id", "gene_id", "rank", "start", "end", "pad"])
        for i in range(n_rows):
            w.writerow([f"p{i}", "G1", i + 1, i * 400, i * 400 + 50, pad])
    with open(out_dir / "manifest.json", "w") as fh:
        json.dump({"params": {}}, fh)
    return str(out_dir)


@pytest.fixture
def clean_cache():
    from flexprobe import select_export as se
    se.clear_results_cache()
    yield se
    se.clear_results_cache()


def test_cache_hits_are_returned_without_rereading(clean_cache, tmp_path):
    se = clean_cache
    run = _tiny_run(tmp_path / "r")
    first = se.load_results_rows(run)
    assert se.load_results_rows(run) is first      # same object: a real hit
    assert se.results_cache_paths() == [os.path.join(run, "results.csv")]


def test_cache_evicts_past_its_entry_cap(clean_cache, tmp_path, monkeypatch):
    """The cap is the whole point: an unbounded cache is the memory bug."""
    se = clean_cache
    monkeypatch.setattr(se, "_RESULTS_CACHE_MAX_ENTRIES", 3)
    runs = [_tiny_run(tmp_path / f"r{i}") for i in range(5)]
    for r in runs:
        se.load_results_rows(r)
    resident = se.results_cache_paths()
    assert len(resident) == 3
    # Least-recently-used first, so the two oldest are the ones that went.
    assert resident == [os.path.join(r, "results.csv") for r in runs[2:]]


def test_cache_eviction_is_least_recently_USED_not_least_recently_added(
        clean_cache, tmp_path, monkeypatch):
    se = clean_cache
    monkeypatch.setattr(se, "_RESULTS_CACHE_MAX_ENTRIES", 2)
    a, b, c = (_tiny_run(tmp_path / n) for n in "abc")
    se.load_results_rows(a)
    se.load_results_rows(b)
    se.load_results_rows(a)          # re-touch a: b is now the LRU
    se.load_results_rows(c)
    assert se.results_cache_paths() == [
        os.path.join(a, "results.csv"), os.path.join(c, "results.csv")]


def test_cache_evicts_on_the_byte_budget_too(clean_cache, tmp_path, monkeypatch):
    """Four small runs are within the entry cap and still must not all stay
    resident if they blow the byte budget between them."""
    se = clean_cache
    monkeypatch.setattr(se, "_RESULTS_CACHE_MAX_ENTRIES", 100)
    monkeypatch.setattr(se, "_RESULTS_CACHE_MAX_BYTES", 4000)
    runs = [_tiny_run(tmp_path / f"b{i}", n_rows=20, pad="x" * 60)
            for i in range(4)]
    sizes = [os.path.getsize(os.path.join(r, "results.csv")) for r in runs]
    assert all(500 < sz < 4000 for sz in sizes), sizes
    for r in runs:
        se.load_results_rows(r)
    resident = se.results_cache_paths()
    assert 0 < len(resident) < 4
    total = sum(os.path.getsize(p) for p in resident)
    assert total <= 4000


def test_a_file_larger_than_the_budget_is_never_admitted(
        clean_cache, tmp_path, monkeypatch):
    """The 1.17 GB case. Admitting it would evict everything else and still
    blow the budget, so it is re-read instead of cached."""
    se = clean_cache
    monkeypatch.setattr(se, "_RESULTS_CACHE_MAX_BYTES", 200)
    small = _tiny_run(tmp_path / "small", n_rows=1)
    big = _tiny_run(tmp_path / "big", n_rows=50, pad="y" * 80)
    se.load_results_rows(small)
    rows_1 = se.load_results_rows(big)
    rows_2 = se.load_results_rows(big)
    assert rows_1 is not rows_2                       # re-read, not cached
    assert os.path.join(big, "results.csv") not in se.results_cache_paths()
    assert os.path.join(small, "results.csv") in se.results_cache_paths()


def test_a_changed_file_is_a_miss_even_when_the_path_is_resident(
        clean_cache, tmp_path):
    se = clean_cache
    run = _tiny_run(tmp_path / "r")
    first = se.load_results_rows(run)
    assert len(first) == 3
    _tiny_run(tmp_path / "r", n_rows=7)
    os.utime(os.path.join(run, "results.csv"), (0, 0))
    assert len(se.load_results_rows(run)) == 7


def test_progress_and_cancel_reach_the_streaming_reader(clean_cache, tmp_path):
    """Cancel must stop the READ, not just hide a dialog."""
    from flexprobe.select_export import ResultsLoadCancelled
    se = clean_cache
    # Comfortably over the 1 MiB poll interval so the hooks fire mid-file.
    run = _tiny_run(tmp_path / "big", n_rows=4000, pad="z" * 900)
    assert os.path.getsize(os.path.join(run, "results.csv")) > 2 * 1024 * 1024

    seen = []
    se.load_results_rows(run, progress=seen.append)
    assert seen and seen[-1] == 100
    assert any(0 < p < 100 for p in seen), seen

    se.clear_results_cache()
    with pytest.raises(ResultsLoadCancelled):
        se.load_results_rows(run, should_cancel=lambda: True)
    # A cancelled read must not leave a partial row list in the cache.
    assert se.results_cache_paths() == []


# ---------------------------------------------------------------------------
# The worker — Qt
# ---------------------------------------------------------------------------

qt = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)


@qt
def test_run_data_worker_parses_off_the_main_thread(qapp, clean_cache, tmp_path):
    """The load must not happen on the UI thread. Asserted by capturing the
    thread identity INSIDE the parse, not by trusting that a QThread was
    constructed."""
    from PySide6.QtCore import QEventLoop
    from flexprobe import select_export as se
    from flexprobe_gui.review_dialog import _RunDataWorker

    run = _tiny_run(tmp_path / "r", n_rows=5)
    main_thread = threading.get_ident()
    parse_threads = []

    real_load = se.load_results_rows

    def spy(*a, **kw):
        parse_threads.append(threading.get_ident())
        return real_load(*a, **kw)

    import flexprobe_gui.review_dialog as rd
    orig = rd.load_results_rows
    rd.load_results_rows = spy
    try:
        worker = _RunDataWorker(run)
        loop = QEventLoop()
        got = {}
        worker.done.connect(lambda d: (got.setdefault("data", d), loop.quit()))
        worker.failed.connect(lambda m: (got.setdefault("err", m), loop.quit()))
        worker.start()
        loop.exec()
        worker.wait()
    finally:
        rd.load_results_rows = orig

    assert "err" not in got, got.get("err")
    assert parse_threads and all(t != main_thread for t in parse_threads)
    assert list(got["data"].rows_by_gene) == ["G1"]
    assert len(got["data"].rows_by_gene["G1"]) == 5


@qt
def test_review_dialog_loads_through_the_worker(qapp, clean_cache, tmp_path):
    """End to end: the dialog's own data arrives via the worker's done signal."""
    from flexprobe_gui.review_dialog import ReviewDialog, _RunDataWorker
    from flexprobe_gui.run_history import RunRecord

    run = _tiny_run(tmp_path / "r", n_rows=4)
    started = []
    real_start = _RunDataWorker.start

    def spy_start(self, *a, **kw):
        started.append(self)
        return real_start(self, *a, **kw)

    _RunDataWorker.start = spy_start
    try:
        dlg = ReviewDialog(RunRecord(name="t", output_dir=run,
                                     status="ready_for_review"))
    finally:
        _RunDataWorker.start = real_start

    assert len(started) == 1
    assert dlg._load_aborted is False
    assert len(dlg._rows_by_gene["G1"]) == 4
    assert dlg.gene_combo.count() == 1
    dlg.deleteLater()


@qt
def test_cancelling_the_load_refuses_to_show_an_empty_dialog(
        qapp, clean_cache, tmp_path):
    """A cancelled load must not be presented as 'this run has no candidates'."""
    from PySide6.QtWidgets import QDialog
    from flexprobe_gui.review_dialog import ReviewDialog, _RunDataWorker
    from flexprobe_gui.run_history import RunRecord

    run = _tiny_run(tmp_path / "r", n_rows=4)
    # Cancel the moment the worker starts, before it can read anything.
    real_start = _RunDataWorker.start

    def cancel_then_start(self, *a, **kw):
        self.request_cancel()
        return real_start(self, *a, **kw)

    _RunDataWorker.start = cancel_then_start
    try:
        dlg = ReviewDialog(RunRecord(name="t", output_dir=run,
                                     status="ready_for_review"))
    finally:
        _RunDataWorker.start = real_start

    if dlg._load_aborted:                 # the cancel won the race
        assert dlg.exec() == QDialog.Rejected
        assert dlg._rows_by_gene == {}
    else:                                 # tiny file finished first — fine
        assert dlg.gene_combo.count() == 1
    dlg.deleteLater()


@qt
def test_hpc_import_runs_off_the_main_thread(qapp, tmp_path):
    """The import extracts a cluster archive; that cannot block the UI either."""
    from PySide6.QtCore import QEventLoop
    from flexprobe_gui.hpc_dialog import _ImportWorker

    import flexprobe_gui.hpc_dialog as hd

    run = _tiny_run(tmp_path / "cluster_out", n_rows=2)
    main_thread = threading.get_ident()
    seen = {}
    work_threads = []

    real_import = hd.import_hpc_output

    def spy(*a, **kw):
        work_threads.append(threading.get_ident())
        return real_import(*a, **kw)

    hd.import_hpc_output = spy
    try:
        worker = _ImportWorker(run)
        loop = QEventLoop()
        worker.finished_ok.connect(
            lambda info: (seen.update(info=info, tid=threading.get_ident()),
                          loop.quit()))
        worker.failed.connect(lambda m, u: (seen.update(err=m), loop.quit()))
        worker.start()
        loop.exec()
        worker.wait()
    finally:
        hd.import_hpc_output = real_import

    assert "err" not in seen, seen.get("err")
    assert seen["info"]["output_dir"] == run
    assert seen["info"]["status"] == "ready_for_review"
    # The extraction/scan ran on the worker; only the signal came back to the UI.
    assert work_threads and all(t != main_thread for t in work_threads)
    assert seen["tid"] == main_thread
    assert worker.isFinished()
