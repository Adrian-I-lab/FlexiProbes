"""GUI: New query form layout, run naming, fan-out and the output folder.

Covers the New query dialog: inputs first, a run name
that follows the target FASTA, the strain-unique "Every sequence" fan-out, and
an output folder that sticks as the default.

Skipped when PySide6 is unavailable. Headless via offscreen Qt.
"""
import importlib.util
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)


def _dialog(qapp, tmp_path, **kwargs):
    # Point all stores at a temp dir so the test never touches ~/.flexprobe.
    os.environ["FLEXPROBE_HOME"] = str(tmp_path)
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    return NewQueryDialog(**kwargs)


def _top_level_groups(dlg):
    """Group boxes in the order the scroll content lays them out."""
    from PySide6.QtWidgets import QGroupBox, QScrollArea
    layout = dlg.findChild(QScrollArea).widget().layout()
    out = []
    for i in range(layout.count()):
        w = layout.itemAt(i).widget()
        if isinstance(w, QGroupBox):
            out.append(w)
    return out


# ---------------------------------------------------------------------------
# Inputs first
# ---------------------------------------------------------------------------

def test_inputs_group_leads_the_form(qapp, tmp_path):
    """The target FASTA and the off-target references come before any setting."""
    dlg = _dialog(qapp, tmp_path)
    groups = _top_level_groups(dlg)
    titles = [g.title() for g in groups]
    assert titles == [
        "Inputs", "Output", "Chemistry", "Homology mode", "Workflow",
        "Off-target screening",
    ]
    inputs, output = groups[0], groups[1]
    assert inputs.isAncestorOf(dlg.fasta_edit)
    assert inputs.isAncestorOf(dlg.db_edit)
    assert output.isAncestorOf(dlg.out_edit)
    assert output.isAncestorOf(dlg.overwrite)
    assert output.isAncestorOf(dlg.run_name)
    # The remaining screening options stay in their own group.
    assert groups[-1].isAncestorOf(dlg.screen_engine)
    assert not groups[-1].isAncestorOf(dlg.db_edit)


def test_bulk_fasta_list_sits_in_the_inputs_group(qapp, tmp_path):
    a = tmp_path / "a.fasta"
    b = tmp_path / "b.fasta"
    for p in (a, b):
        p.write_text(">g\n" + "ACGT" * 25 + "\n")
    dlg = _dialog(qapp, tmp_path, bulk_fasta_paths=[str(a), str(b)])
    inputs = _top_level_groups(dlg)[0]
    assert inputs.title() == "Inputs"
    assert inputs.isAncestorOf(dlg.bulk_list)
    assert dlg.bulk_list.count() == 2


# ---------------------------------------------------------------------------
# The run name follows the target FASTA
# ---------------------------------------------------------------------------

_SEQ_A = ("GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
          "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG")
_SEQ_B = _SEQ_A[:-1] + "A"


def _segment_fasta(tmp_path, segment):
    """A two-strain FASTA for one segment, named like the lab's inputs."""
    fa = tmp_path / f"OROV_{segment}_2strain.fasta"
    fa.write_text(
        f">STRAIN_A|{segment}|ACC1\n{_SEQ_A}\n"
        f">STRAIN_B|{segment}|ACC2\n{_SEQ_B}\n"
    )
    return fa


def _select_strategy(dlg, prefix):
    dlg.strain_mode_combo.setCurrentIndex(1)            # Homology consensus
    idx = next(i for i in range(dlg.strategy_combo.count())
               if dlg.strategy_combo.itemText(i).startswith(prefix))
    dlg.strategy_combo.setCurrentIndex(idx)


def _name_part(output_dir):
    """The run-name part of '<YYYYmmdd>_<HHMMSS>_<name>'."""
    return os.path.basename(output_dir).split("_", 2)[2]


def test_rerun_of_a_small_run_follows_a_new_large_fasta(qapp, tmp_path):
    """Re-run a Small-segment run, point it at the Large FASTA: no "Small".

    The dry launch produced a Large run named
    "..._Small_2strain__strainunique-discriminating__..._Large_...": the draft
    stored the generated name, Re-run restored it into the Run name box, and
    picking a new FASTA never touched it.
    """
    small = _segment_fasta(tmp_path, "Small")
    large = _segment_fasta(tmp_path, "Large")
    dlg = _dialog(qapp, tmp_path)
    _select_strategy(dlg, "Sequence-unique")
    dlg.fasta_edit.setText(str(small))
    dlg.out_edit.setText(str(tmp_path / "out"))
    rec, _cfg, _c = dlg.build_runs()[0]
    assert rec.name.startswith("OROV_Small_2strain")

    # Re-run: home_window applies the record's stored draft to a new dialog.
    rerun = _dialog(qapp, tmp_path)
    rerun.apply_draft(rec.query_draft)
    rerun.fasta_edit.setText(str(large))
    rec2, cfg2, _c = rerun.build_runs()[0]
    assert cfg2.fasta_path == str(large)
    assert _name_part(cfg2.output_dir).startswith("OROV_Large_2strain")
    assert "Small" not in os.path.basename(cfg2.output_dir)
    assert "Small" not in rec2.name


def test_a_legacy_suffixed_draft_name_does_not_leak(qapp, tmp_path):
    """Older records can hold the full generated name, suffix included.

    Both routes into the form carry it: a stored query_draft, and the
    config_summary fallback, which copies the record's name.
    """
    small = _segment_fasta(tmp_path, "Small")
    large = _segment_fasta(tmp_path, "Large")
    base = _dialog(qapp, tmp_path)
    _select_strategy(base, "Sequence-unique")
    draft = base.to_draft()
    draft.update({
        "fasta_path": str(small),
        "run_name": ("OROV_Small_2strain__strainunique-discriminating"
                     "__STRAIN_A_Small_ACC1"),
    })

    rerun = _dialog(qapp, tmp_path)
    rerun.apply_draft(draft)
    # The stem was derived, not typed, so the box falls back to its placeholder.
    assert rerun.run_name.text() == ""
    rerun.fasta_edit.setText(str(large))
    rerun.out_edit.setText(str(tmp_path / "out"))
    _rec, cfg, _c = rerun.build_runs()[0]
    assert _name_part(cfg.output_dir).startswith("OROV_Large_2strain")
    assert "Small" not in os.path.basename(cfg.output_dir)


def test_a_legacy_typed_name_keeps_only_what_was_typed(qapp, tmp_path):
    from flexprobe_gui._nq_build import strip_generated_run_name
    small = _segment_fasta(tmp_path, "Small")
    dlg = _dialog(qapp, tmp_path)
    draft = dlg.to_draft()
    draft.update({"fasta_path": str(small),
                  "run_name": "orov_panel__conserved-panstrain__STRAIN_B"})
    dlg.apply_draft(draft)
    assert dlg.run_name.text() == "orov_panel"
    # Every generated tag is cut at its earliest occurrence, constitutive too.
    assert strip_generated_run_name("x__constitutive") == "x"
    assert strip_generated_run_name(
        "x__strainunique-discriminating__A__conserved-panstrain__B") == "x"
    assert strip_generated_run_name("plain_name") == "plain_name"


def test_a_restored_name_is_cleared_when_the_fasta_changes(qapp, tmp_path):
    small = _segment_fasta(tmp_path, "Small")
    large = _segment_fasta(tmp_path, "Large")
    dlg = _dialog(qapp, tmp_path)
    draft = dlg.to_draft()
    draft.update({"fasta_path": str(small), "run_name": "orov_panel"})
    dlg.apply_draft(draft)
    assert dlg.run_name.text() == "orov_panel"      # same FASTA: kept
    dlg.fasta_edit.setText(str(large))
    assert dlg.run_name.text() == ""                # new FASTA: cleared


def test_a_name_typed_in_this_dialog_survives_a_fasta_change(qapp, tmp_path):
    from PySide6.QtTest import QTest
    small = _segment_fasta(tmp_path, "Small")
    large = _segment_fasta(tmp_path, "Large")
    dlg = _dialog(qapp, tmp_path)
    draft = dlg.to_draft()
    draft.update({"fasta_path": str(small), "run_name": "orov_panel"})
    dlg.apply_draft(draft)
    dlg.run_name.clear()
    QTest.keyClicks(dlg.run_name, "my_large_panel")
    dlg.fasta_edit.setText(str(large))
    assert dlg.run_name.text() == "my_large_panel"
    dlg.out_edit.setText(str(tmp_path / "out"))
    _rec, cfg, _c = dlg.build_runs()[0]
    assert _name_part(cfg.output_dir) == "my_large_panel"


def test_the_stored_draft_holds_only_the_typed_name(qapp, tmp_path):
    small = _segment_fasta(tmp_path, "Small")
    dlg = _dialog(qapp, tmp_path)
    _select_strategy(dlg, "Sequence-unique")
    dlg.fasta_edit.setText(str(small))
    dlg.out_edit.setText(str(tmp_path / "out"))
    rec, _cfg, _c = dlg.build_runs()[0]
    assert "__" in rec.name                         # the record is suffixed...
    assert rec.query_draft["run_name"] == ""        # ...the draft is not

    dlg.run_name.setText("orov_panel")
    rec, _cfg, _c = dlg.build_runs()[0]
    assert rec.name.startswith("orov_panel__strainunique-discriminating")
    assert rec.query_draft["run_name"] == "orov_panel"


# ---------------------------------------------------------------------------
# Strain-unique "Every sequence" fan-out
# ---------------------------------------------------------------------------

def test_strain_unique_every_sequence_fans_out_one_run_per_record(
        qapp, tmp_path):
    """build_runs can fan a strain-unique run out, and the picker must offer
    the sentinel so the GUI can reach it."""
    from flexprobe_gui._nq_build import ALL_TEMPLATES
    fa = _segment_fasta(tmp_path, "Large")
    dlg = _dialog(qapp, tmp_path)
    _select_strategy(dlg, "Sequence-unique")
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(str(tmp_path / "out"))
    idx = dlg.unique_target_combo.findData(ALL_TEMPLATES)
    assert idx >= 0, "the Every sequence entry is not offered"
    dlg.unique_target_combo.setCurrentIndex(idx)

    runs = dlg.build_runs()
    assert len(runs) == 2, "expected one run per input record"
    by_target = {cfg.strain_unique_target: cfg for _r, cfg, _c in runs}
    assert sorted(by_target) == [0, 1]
    for _rec, cfg, _c in runs:
        assert cfg.consensus_strategy == "strain_unique"
        assert cfg.consensus_reference == 0
        cfg.validate()
    outs = [os.path.basename(cfg.output_dir) for _r, cfg, _c in runs]
    assert len(set(outs)) == 2
    assert _name_part(by_target[0].output_dir).endswith(
        "OROV_Large_2strain__strainunique-discriminating__STRAIN_A_Large_ACC1")
    assert _name_part(by_target[1].output_dir).endswith(
        "OROV_Large_2strain__strainunique-discriminating__STRAIN_B_Large_ACC2")


def test_every_sequence_is_offered_last_and_is_not_the_default(qapp, tmp_path):
    from flexprobe_gui._nq_build import ALL_TEMPLATES
    fa = _segment_fasta(tmp_path, "Large")
    dlg = _dialog(qapp, tmp_path)
    _select_strategy(dlg, "Sequence-unique")
    dlg.fasta_edit.setText(str(fa))
    combo = dlg.unique_target_combo
    assert [combo.itemData(i) for i in range(combo.count())] == [
        0, 1, ALL_TEMPLATES]
    assert combo.currentData() == 0
    # One record: nothing to fan out to.
    solo = tmp_path / "solo.fasta"
    solo.write_text(f">SOLO|Large|ACC9\n{_SEQ_A}\n")
    dlg.fasta_edit.setText(str(solo))
    assert [combo.itemData(i) for i in range(combo.count())] == [0]


@pytest.mark.parametrize("choice", ["every", 1])
def test_unique_target_round_trips_through_a_draft(qapp, tmp_path, choice):
    from flexprobe_gui._nq_build import ALL_TEMPLATES
    want = ALL_TEMPLATES if choice == "every" else choice
    fa = _segment_fasta(tmp_path, "Large")
    dlg = _dialog(qapp, tmp_path)
    _select_strategy(dlg, "Sequence-unique")
    dlg.fasta_edit.setText(str(fa))
    dlg.unique_target_combo.setCurrentIndex(
        dlg.unique_target_combo.findData(want))
    draft = dlg.to_draft()
    assert draft["unique_target"] == want

    # A fresh dialog: the picker is empty until the draft's FASTA fills it.
    other = _dialog(qapp, tmp_path)
    other.apply_draft(draft)
    assert other.unique_target_combo.currentData() == want
    other.out_edit.setText(str(tmp_path / "out"))
    n = len(other.build_runs())
    assert n == (2 if choice == "every" else 1)


def test_a_single_template_rerun_is_not_widened_to_every_sequence(
        qapp, tmp_path):
    """The template choice was selected into a still-empty combo and lost, so
    the refill defaulted to "Every sequence" and the Re-run queued N runs."""
    fa = _segment_fasta(tmp_path, "Large")
    dlg = _dialog(qapp, tmp_path)
    _select_strategy(dlg, "Conserved")
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(str(tmp_path / "out"))
    dlg.template_combo.setCurrentIndex(dlg.template_combo.findData(1))
    rec, _cfg, _c = dlg.build_runs()[0]

    rerun = _dialog(qapp, tmp_path)
    rerun.apply_draft(rec.query_draft)
    assert rerun.template_combo.currentData() == 1
    runs = rerun.build_runs()
    assert len(runs) == 1
    assert runs[0][1].consensus_reference == 1


# ---------------------------------------------------------------------------
# The output folder sticks, and a draft does not override it
# ---------------------------------------------------------------------------

def _settings(tmp_path):
    from flexprobe_gui.settings import AppSettings
    return AppSettings(str(tmp_path / "settings.json"))


def test_a_browsed_output_folder_becomes_the_default(
        qapp, tmp_path, monkeypatch):
    from flexprobe_gui import new_query_dialog as nqd
    from flexprobe_gui.tools_dialog import PreferencesDialog
    chosen = str(tmp_path / "lab_outputs")

    class _Picker:
        @staticmethod
        def getExistingDirectory(*_a, **_k):
            return chosen

    monkeypatch.setattr(nqd, "QFileDialog", _Picker)
    dlg = _dialog(qapp, tmp_path, settings=_settings(tmp_path))
    assert dlg.out_edit.text() != chosen
    dlg._browse_out()
    assert dlg.out_edit.text() == chosen

    # A new query, and Preferences, each read the folder from disk afresh.
    assert _dialog(qapp, tmp_path,
                   settings=_settings(tmp_path)).out_edit.text() == chosen
    assert PreferencesDialog(_settings(tmp_path)).out_edit.text() == chosen


def test_a_typed_output_folder_is_saved_when_the_run_starts(qapp, tmp_path):
    from PySide6.QtWidgets import QDialog
    fa = _segment_fasta(tmp_path, "Large")
    typed = str(tmp_path / "typed_outputs")
    dlg = _dialog(qapp, tmp_path, settings=_settings(tmp_path))
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(typed)
    # Typing alone does not save it: only starting the run does.
    assert _settings(tmp_path).default_output_dir() != typed
    dlg._validate_and_accept()
    assert dlg.result() == QDialog.DialogCode.Accepted.value
    assert _settings(tmp_path).default_output_dir() == typed
    assert _dialog(qapp, tmp_path,
                   settings=_settings(tmp_path)).out_edit.text() == typed


@pytest.mark.parametrize("route", ["rerun_draft", "loaded_draft"])
def test_a_draft_does_not_override_the_current_default_folder(
        qapp, tmp_path, route):
    """Re-run and Load draft both go through apply_draft (decision Q1)."""
    small = _segment_fasta(tmp_path, "Small")
    old = str(tmp_path / "old_run_folder")
    current = str(tmp_path / "current_default")
    first = _dialog(qapp, tmp_path, settings=_settings(tmp_path))
    first.fasta_edit.setText(str(small))
    first.out_edit.setText(old)
    if route == "rerun_draft":
        draft = first.build_runs()[0][0].query_draft
    else:
        draft = first.to_draft()
    assert draft["output_dir"] == old

    _settings(tmp_path).set_default_output_dir(current)
    dlg = _dialog(qapp, tmp_path, settings=_settings(tmp_path))
    dlg.apply_draft(draft)
    assert dlg.out_edit.text() == current
    _rec, cfg, _c = dlg.build_runs()[0]
    assert cfg.output_dir.startswith(current + os.sep)
