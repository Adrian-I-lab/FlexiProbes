"""Importing run folders, and finding them again after they move.

What this file guards, all on one row of the home table:

* a run folder that is moved or deleted must not leave a row whose Open does
  nothing, whose Review fails with "results.csv not found", and which a
  project bulk export skips without a word;
* a desktop run folder (copied from another machine, or listed in a history
  that was lost) can be brought back into the list without going through the
  HPC importer, so it is not badged HPC, re-runs with its own FASTA and
  settings, fills its Pairs column, and is not imported twice.

Every run folder here is made by the REAL pipeline, on a tiny FASTA with the
off-target screen skipped, so the files the importer reads (manifest.json,
run_config.json, design_summary.json, results.csv) have their true shape rather
than a hand-written imitation of it. Imports are done inside each test so that,
against code without the feature, the tests that drive an existing entry point
fail on their assertions rather than all at once on an ImportError.
"""
import datetime
import importlib.util
import json
import os
import shutil

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)

STATUS_COL = 5


def _actions_col() -> int:
    # Read from the app rather than restated, so a column move cannot break it.
    from flexprobe_gui._home_completed import COL_ACTIONS
    return COL_ACTIONS

# 1,500 bp. Enough to leave candidates after the filters; the pipeline takes a
# fraction of a second on it with the screen skipped.
_SEQ = ("AGCTGTAGCTAGCTAGCATGCATGCATGCAGTCAGTCGATGCATCGTAGCATCGATCGTACTAGCTAGCAT"
        * 30)[:1500]
_NON_DEFAULT_TILE_STEP = 7


@pytest.fixture(scope="module")
def run_templates(tmp_path_factory):
    """Three finished desktop runs, ``{"alpha", "beta", "gamma"}`` -> folder.

    Built once per module; tests COPY them (``_place``) so each test owns its
    folders. The runs differ in their FASTA header, so each has its own FASTA
    SHA-256, and in their manifest timestamp.
    """
    from flexprobe.models import ExportSettings, RunConfig, ScreenResources
    from flexprobe.pipeline import run_pipeline
    from flexprobe.summary import write_ok_summary
    from flexprobe_gui.modality import MODALITY_PRESETS

    base = tmp_path_factory.mktemp("run_templates")
    out = {}
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("FLEXPROBE_HOME", str(base / ".flexprobe_home"))
        for i, name in enumerate(("alpha", "beta", "gamma")):
            fasta = base / f"{name}.fa"
            fasta.write_text(f">G{i + 1}\n{_SEQ}\n")
            run_dir = base / name
            cfg = RunConfig(
                fasta_path=str(fasta), output_dir=str(run_dir),
                modality=MODALITY_PRESETS["Flex v2 (10X GEM-X)"],
                discovery_only=True, workers=1,
                tile_step=_NON_DEFAULT_TILE_STEP,
                screen=ScreenResources(skip_screen=True),
                export=ExportSettings(overwrite=True),
            )
            result = run_pipeline(cfg)
            # What the desktop worker writes after every run (run_manager).
            write_ok_summary(str(run_dir), result["stats"])
            out[name] = str(run_dir)
    return out


def _place(templates, key, dest):
    dest = str(dest)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    shutil.copytree(templates[key], dest)
    return dest


def _json(path):
    with open(path) as fh:
        return json.load(fh)


def _buttons(win, row):
    from PySide6.QtWidgets import QPushButton
    return [b.text() for b in
            win.table.cellWidget(row, _actions_col()).findChildren(QPushButton)]


def _status_texts(win, row):
    from PySide6.QtWidgets import QLabel
    widget = win.table.cellWidget(row, STATUS_COL)
    if widget is None:
        return [win.table.item(row, STATUS_COL).text()]
    return [lb.text() for lb in widget.findChildren(QLabel)]


def _silence_boxes(monkeypatch, answer=None):
    """Auto-answer every modal box; returns the list warnings are recorded in."""
    from PySide6.QtWidgets import QMessageBox
    warnings = []
    monkeypatch.setattr(QMessageBox, "exec",
                        lambda self: answer if answer is not None else QMessageBox.Ok)
    monkeypatch.setattr(QMessageBox, "warning",
                        staticmethod(lambda *a, **k: warnings.append(a)))
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda *a, **k: QMessageBox.Ok))
    return warnings


# ---------------------------------------------------------------------------
# The shared folder reader
# ---------------------------------------------------------------------------

def test_describing_a_run_folder_reads_its_counts_and_identity(run_templates,
                                                                tmp_path):
    from flexprobe.hpc_export import _describe_run
    run = _place(run_templates, "alpha", tmp_path / "alpha")
    info = _describe_run(run)
    summary = _json(os.path.join(run, "design_summary.json"))
    manifest = _json(os.path.join(run, "manifest.json"))

    assert summary["stats"]["pairs_output"] > 0
    assert info["stats"]["pairs_output"] == summary["stats"]["pairs_output"]
    assert info["stats"]["imported"] is True
    assert info["manifest_timestamp"] == manifest["timestamp"]
    assert info["fasta_sha256"] == manifest["inputs"]["fasta_sha256"]
    assert info["app_version"] == manifest["app_version"]
    assert info["has_run_config"] is True
    assert info["cluster_evidence"] == ""


# ---------------------------------------------------------------------------
# I1. Import runs…
# ---------------------------------------------------------------------------

def test_importing_one_run_folder_fills_the_record_from_its_own_files(
        run_templates, tmp_path):
    from flexprobe.config_io import load_config
    from flexprobe_gui.run_history import RunHistory
    from flexprobe_gui.run_import import import_runs, import_source

    run = _place(run_templates, "alpha", tmp_path / "runs" / "alpha")
    history = RunHistory(str(tmp_path / "runs.json"))
    summary = import_runs(run, history)

    assert [r.name for r in summary.imported] == ["alpha"]
    assert summary.already_present == [] and summary.skipped == []
    rec = history.list_completed()[0]
    cfg = load_config(os.path.join(run, "run_config.json"))
    manifest = _json(os.path.join(run, "manifest.json"))
    design = _json(os.path.join(run, "design_summary.json"))

    assert rec.output_dir == os.path.abspath(run)
    assert rec.fasta_path == cfg.fasta_path and rec.fasta_path
    assert rec.modality_label == "Flex v2 (10X GEM-X)"
    assert rec.status == "ready_for_review"
    assert rec.stats["discovery_only"] is True
    assert rec.stats["pairs_output"] == design["stats"]["pairs_output"]
    assert rec.stats["app_version"] == manifest["app_version"]
    assert rec.manifest_timestamp == manifest["timestamp"]
    assert rec.fasta_sha256 == manifest["inputs"]["fasta_sha256"]
    assert import_source(rec) == "desktop"
    assert rec.config_summary["tile_step"] == _NON_DEFAULT_TILE_STEP
    started = datetime.datetime.fromisoformat(rec.started_at)
    completed = datetime.datetime.fromisoformat(rec.completed_at)
    assert started <= completed

    # Persisted, identity included, through the ordinary atomic save.
    again = RunHistory(str(tmp_path / "runs.json")).list_completed()
    assert [(r.run_id, r.manifest_timestamp) for r in again] == \
        [(rec.run_id, manifest["timestamp"])]


def test_a_run_without_run_config_resolves_the_manifest_fasta_path(
        run_templates, tmp_path):
    # The manifest names its input FASTA relative to the run folder.
    # With no run_config.json the importer falls back to the manifest, and a
    # relative path copied as is would point nowhere from the GUI's own cwd.
    from flexprobe_gui.run_history import RunHistory
    from flexprobe_gui.run_import import import_runs

    run = _place(run_templates, "alpha", tmp_path / "alpha")
    src_fasta = os.path.join(os.path.dirname(run_templates["alpha"]), "alpha.fa")
    shutil.copy(src_fasta, tmp_path / "alpha.fa")    # same layout as the run
    os.remove(os.path.join(run, "run_config.json"))
    fasta_in_manifest = _json(os.path.join(run, "manifest.json"))["inputs"]["fasta"]
    assert not os.path.isabs(fasta_in_manifest)    # the relative-path shape this covers

    history = RunHistory(str(tmp_path / "runs.json"))
    import_runs(run, history)
    rec = history.list_completed()[0]

    assert rec.fasta_path == str(tmp_path / "alpha.fa")
    assert os.path.isfile(rec.fasta_path)


def test_an_imported_run_re_runs_with_its_own_fasta_and_settings(
        qapp, run_templates, tmp_path):
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunHistory
    from flexprobe_gui.run_import import import_runs

    run = _place(run_templates, "alpha", tmp_path / "alpha")
    history = RunHistory(str(tmp_path / "runs.json"))
    rec = import_runs(run, history).imported[0]

    win = HomeWindow()
    try:
        draft = win._record_to_draft(rec)
    finally:
        win.close()
    assert draft["fasta_path"] == rec.fasta_path and draft["fasta_path"]
    assert draft["tile_step"] == _NON_DEFAULT_TILE_STEP
    assert draft["modality"] == "Flex v2 (10X GEM-X)"


def test_every_setting_re_run_reads_is_rebuilt_from_run_config():
    """Drift guard: the draft map gains a key, the import must produce it."""
    from flexprobe.models import RunConfig
    from flexprobe_gui._nq_build import CONFIG_SUMMARY_TO_DRAFT
    from flexprobe_gui.run_import import (
        SUMMARY_KEYS_NOT_FROM_CONFIG, summary_from_config,
    )
    produced = set(summary_from_config(RunConfig(fasta_path="x.fa")))
    wanted = {f.summary for f in CONFIG_SUMMARY_TO_DRAFT}
    assert wanted - produced == set(SUMMARY_KEYS_NOT_FROM_CONFIG)


def test_importing_a_parent_finds_nested_runs_and_says_what_it_skipped(
        run_templates, tmp_path):
    from flexprobe_gui.run_history import RunHistory
    from flexprobe_gui.run_import import import_runs

    base = tmp_path / "base"
    _place(run_templates, "alpha", base / "projA" / "alpha")
    _place(run_templates, "beta", base / "projA" / "beta")
    _place(run_templates, "gamma", base / "projB" / "gamma")
    bare = base / "projB" / "bare"
    bare.mkdir(parents=True)
    (bare / "results.csv").write_text("rank,probe_id,gene_id\n1,p1,G1\n")
    # Four levels down: past the scan depth, so it must appear NOWHERE. A copy
    # of a real run, so a scan that went too deep would report it rather than
    # pass it by in silence.
    deep = _place(run_templates, "alpha", base / "d1" / "d2" / "d3" / "deep")

    history = RunHistory(str(tmp_path / "runs.json"))
    summary = import_runs(str(base), history)

    assert sorted(r.name for r in summary.imported) == ["alpha", "beta", "gamma"]
    assert [p for p, _why in summary.skipped] == [str(bare)]
    assert "no manifest.json" in summary.skipped[0][1]
    everything = ([r.output_dir for r in summary.imported]
                  + [p for p, _ in summary.already_present + summary.skipped
                     + summary.relinked])
    assert deep not in everything
    assert len(history.list_completed()) == 3
    assert "Imported 3 runs" in summary.message()


def test_a_second_import_and_a_copied_run_are_not_added_again(run_templates,
                                                              tmp_path):
    from flexprobe_gui.run_history import RunHistory
    from flexprobe_gui.run_import import import_runs

    base = tmp_path / "base"
    _place(run_templates, "alpha", base / "p" / "alpha")
    _place(run_templates, "beta", base / "p" / "beta")
    history = RunHistory(str(tmp_path / "runs.json"))
    assert len(import_runs(str(base), history).imported) == 2

    again = import_runs(str(base), history)
    assert again.imported == []
    assert sorted(n for _p, n in again.already_present) == ["alpha", "beta"]

    # The same run copied to a second place: a different path, the same
    # manifest timestamp and FASTA hash.
    copy = _place(run_templates, "alpha", tmp_path / "elsewhere" / "alpha_copy")
    third = import_runs(copy, history)
    assert third.imported == []
    assert third.already_present == [(copy, "alpha")]
    assert len(history.list_completed()) == 2


def test_a_cluster_output_folder_keeps_the_hpc_badge(run_templates, tmp_path):
    from flexprobe_gui.run_history import RunHistory
    from flexprobe_gui.run_import import import_runs, import_source

    selected = _place(run_templates, "alpha", tmp_path / "a" / "output")
    with open(os.path.join(selected, "selected_probes_per_gene.csv"), "w") as fh:
        fh.write("gene_id,probe_id\n")
    in_bundle = _place(run_templates, "beta", tmp_path / "bundle" / "output")
    (tmp_path / "bundle" / ".flexiprobes_DONE").write_text("")
    desktop = _place(run_templates, "gamma", tmp_path / "c" / "gamma")

    history = RunHistory(str(tmp_path / "runs.json"))
    for folder in (selected, in_bundle, desktop):
        import_runs(folder, history)
    sources = {r.output_dir: import_source(r) for r in history.list_completed()}
    assert sources == {selected: "hpc", in_bundle: "hpc", desktop: "desktop"}


def test_the_hpc_dialog_refuses_a_run_it_already_imported(qapp, run_templates,
                                                          tmp_path, monkeypatch):
    from flexprobe_gui.hpc_dialog import ImportHpcDialog
    from flexprobe_gui.run_history import RunHistory
    _silence_boxes(monkeypatch)

    run = _place(run_templates, "alpha", tmp_path / "alpha")
    history = RunHistory(str(tmp_path / "runs.json"))
    for _ in range(2):
        dlg = ImportHpcDialog(history)
        dlg.src_edit.setText(run)
        dlg._on_import()
        dlg.deleteLater()

    records = history.list_completed()
    assert len(records) == 1
    # And the HPC path now fills the same fields the desktop import does.
    assert records[0].fasta_path
    assert records[0].stats.get("pairs_output")
    assert records[0].manifest_timestamp


def test_the_import_badge_says_imported_or_hpc(qapp, run_templates, tmp_path,
                                               monkeypatch):
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunRecord
    _silence_boxes(monkeypatch)

    win = HomeWindow()
    try:
        win._on_import_runs(_place(run_templates, "alpha", tmp_path / "alpha"))
        win.manager.history.append(RunRecord(
            name="cluster", status="ready_for_review", stats={"imported": True}))
        win._refresh_completed()
        rows = {win.table.item(r, 2).text(): r
                for r in range(win.table.rowCount())}
        assert "Imported" in _status_texts(win, rows["alpha"])
        assert "HPC" not in _status_texts(win, rows["alpha"])
        assert "HPC" in _status_texts(win, rows["cluster"])
    finally:
        win.close()


# ---------------------------------------------------------------------------
# I2. Missing folders and Locate…
# ---------------------------------------------------------------------------

def test_a_moved_run_shows_missing_then_locate_repairs_it_and_review_opens(
        qapp, run_templates, tmp_path, monkeypatch):
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.review_dialog import ReviewDialog
    from flexprobe_gui.run_history import RunHistory
    _silence_boxes(monkeypatch)

    old = _place(run_templates, "alpha", tmp_path / "old" / "alpha")
    win = HomeWindow()
    try:
        win._on_import_runs(old)
        rec = win.manager.history.list_completed()[0]
        present = _buttons(win, 0)
        assert {"Review", "Open", "Filters"} <= set(present)
        assert "Locate…" not in present

        new = tmp_path / "new" / "alpha"
        new.parent.mkdir()
        shutil.move(old, str(new))
        win._refresh()

        gone = _buttons(win, 0)
        assert "Locate…" in gone
        assert not {"Review", "Open", "Filters"} & set(gone)
        assert "Folder missing" in _status_texts(win, 0)

        assert win._on_locate_run(rec, folder=str(new)) == 1
        assert rec.output_dir == str(new)
        assert "Review" in _buttons(win, 0)
        assert "Folder missing" not in _status_texts(win, 0)
        # Persisted: a fresh read of the store sees the new place.
        reloaded = RunHistory().get(rec.run_id)
        assert reloaded is not None and reloaded.output_dir == str(new)

        review = ReviewDialog(rec)
        assert review.gene_combo.count() == 1
        review.deleteLater()
    finally:
        win.close()


def test_locate_refuses_a_folder_holding_a_different_run(
        qapp, run_templates, tmp_path, monkeypatch):
    from flexprobe_gui.home_window import HomeWindow
    warnings = _silence_boxes(monkeypatch)

    old = _place(run_templates, "alpha", tmp_path / "old" / "alpha")
    win = HomeWindow()
    try:
        win._on_import_runs(old)
        rec = win.manager.history.list_completed()[0]
        shutil.rmtree(old)
        win._refresh()
        # Same folder NAME, different run inside: the name alone must not pass.
        impostor = _place(run_templates, "beta", tmp_path / "other" / "alpha")

        assert win._on_locate_run(rec, folder=impostor) == 0
        assert rec.output_dir == os.path.abspath(old)
        assert warnings and "different run" in warnings[0][2]
        assert "Locate…" in _buttons(win, 0)
    finally:
        win.close()


def test_a_record_without_a_stored_identity_is_matched_by_name_and_confirmed(
        qapp, run_templates, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunRecord
    from flexprobe_gui.run_import import CONFIRM, REJECT, match_folder

    found = _place(run_templates, "alpha", tmp_path / "found" / "alpha")
    renamed = _place(run_templates, "alpha", tmp_path / "found" / "zeta")
    rec = RunRecord(name="old run", status="complete",
                    output_dir=str(tmp_path / "gone" / "alpha"),
                    stats={"discovery_only": True})
    assert not rec.manifest_timestamp
    assert match_folder(rec, found).verdict == CONFIRM
    assert match_folder(rec, renamed).verdict == REJECT

    win = HomeWindow()
    try:
        win.manager.history.append(rec)
        _silence_boxes(monkeypatch, answer=QMessageBox.No)
        assert win._on_locate_run(rec, folder=found) == 0
        assert rec.output_dir == str(tmp_path / "gone" / "alpha")

        _silence_boxes(monkeypatch, answer=QMessageBox.Yes)
        assert win._on_locate_run(rec, folder=found) == 1
        assert rec.output_dir == found
        # The identity is taken now, so the next move is checked properly.
        assert rec.manifest_timestamp == \
            _json(os.path.join(found, "manifest.json"))["timestamp"]
    finally:
        win.close()


def test_locating_one_run_offers_to_repair_the_others_that_moved_with_it(
        qapp, run_templates, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from flexprobe_gui.home_window import HomeWindow
    _silence_boxes(monkeypatch, answer=QMessageBox.Yes)

    base = tmp_path / "base"
    _place(run_templates, "alpha", base / "projA" / "alpha")
    _place(run_templates, "beta", base / "projA" / "beta")
    _place(run_templates, "gamma", base / "projB" / "gamma")
    win = HomeWindow()
    try:
        win._on_import_runs(str(base))
        by_name = {r.name: r for r in win.manager.history.list_completed()}
        moved = tmp_path / "moved"
        shutil.move(str(base), str(moved))
        win._refresh()

        # beta sits beside alpha; gamma is under a different project folder
        # and is found through the moved root.
        assert win._on_locate_run(
            by_name["alpha"], folder=str(moved / "projA" / "alpha")) == 3
        assert {n: r.output_dir for n, r in by_name.items()} == {
            "alpha": str(moved / "projA" / "alpha"),
            "beta": str(moved / "projA" / "beta"),
            "gamma": str(moved / "projB" / "gamma"),
        }
        assert all("Locate…" not in _buttons(win, row)
                   for row in range(win.table.rowCount()))
    finally:
        win.close()


def test_an_import_relinks_a_missing_run_instead_of_adding_it_twice(
        run_templates, tmp_path):
    from flexprobe_gui.run_history import RunHistory
    from flexprobe_gui.run_import import import_runs

    old = _place(run_templates, "alpha", tmp_path / "old" / "alpha")
    history = RunHistory(str(tmp_path / "runs.json"))
    rec = import_runs(old, history).imported[0]
    new = str(tmp_path / "new" / "alpha")
    os.makedirs(os.path.dirname(new))
    shutil.move(old, new)

    summary = import_runs(new, history)
    assert summary.imported == [] and summary.relinked == [(new, "alpha")]
    assert [r.output_dir for r in history.list_completed()] == [new]
    assert history.get(rec.run_id).output_dir == new


def test_a_finished_run_records_its_identity_once_it_completes(
        qapp, run_templates, tmp_path):
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunHistory, RunRecord

    run = _place(run_templates, "alpha", tmp_path / "alpha")
    win = HomeWindow()
    try:
        running = RunRecord(name="r", status="running", output_dir=run)
        win.manager.history.append(running)
        win._on_record_changed(running.run_id)
        assert running.manifest_timestamp == ""      # control: not finished

        running.status = "ready_for_review"
        win.manager.history.update(running)
        win._on_record_changed(running.run_id)
        manifest = _json(os.path.join(run, "manifest.json"))
        assert running.manifest_timestamp == manifest["timestamp"]
        assert running.fasta_sha256 == manifest["inputs"]["fasta_sha256"]
        assert RunHistory().get(running.run_id).manifest_timestamp == \
            manifest["timestamp"]
    finally:
        win.close()


def test_open_and_bulk_export_say_when_a_folder_is_missing(
        qapp, tmp_path, monkeypatch):
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.reveal import reveal_in_finder
    from flexprobe_gui.run_history import RunRecord
    warnings = _silence_boxes(monkeypatch)

    assert reveal_in_finder(str(tmp_path / "nope")) is False
    win = HomeWindow()
    try:
        proj = win.projects.create("P")
        rec = RunRecord(name="moved away", status="complete",
                        project_id=proj.project_id,
                        output_dir=str(tmp_path / "nope"))
        win.manager.history.append(rec)
        win._on_open_folder(rec)
        assert warnings and "no longer at" in warnings[-1][2]

        win._refresh_project_filter()
        win.project_filter_combo.setCurrentIndex(
            win.project_filter_combo.findData(proj.project_id))
        win._on_bulk_export()
        assert "moved away" in warnings[-1][2]
    finally:
        win.close()


# ---------------------------------------------------------------------------
# I3. Manual only
# ---------------------------------------------------------------------------

def test_nothing_scans_for_runs_or_reads_manifests_at_start_up(
        qapp, run_templates, tmp_path, monkeypatch):
    """Start-up and Refresh may stat each record's folder, and nothing more."""
    from flexprobe import hpc_export
    from flexprobe_gui import home_window, run_import
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunHistory, RunRecord
    _silence_boxes(monkeypatch)

    present = _place(run_templates, "alpha", tmp_path / "alpha")
    history = RunHistory()
    history.append(RunRecord(name="here", status="complete", output_dir=present))
    history.append(RunRecord(name="gone", status="complete",
                             output_dir=str(tmp_path / "gone")))

    def tripwire(*_a, **_k):
        raise AssertionError("scanned or read a run folder at start-up")

    for mod, attr in ((hpc_export, "find_run_folders"),
                      (hpc_export, "describe_run_folder"),
                      (run_import, "find_run_folders"),
                      (run_import, "describe_run_folder"),
                      (run_import, "read_identity"),
                      (home_window, "import_runs")):
        monkeypatch.setattr(mod, attr, tripwire)

    win = HomeWindow()
    try:
        win._refresh()
        rows = {win.table.item(r, 2).text(): r
                for r in range(win.table.rowCount())}
        assert "Folder missing" in _status_texts(win, rows["gone"])
        assert "Folder missing" not in _status_texts(win, rows["here"])
        # Positive control: the tripwire is armed, so its silence above means
        # something.
        with pytest.raises(AssertionError, match="start-up"):
            win._on_import_runs(str(tmp_path))
    finally:
        win.close()
