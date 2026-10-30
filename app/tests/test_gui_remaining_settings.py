"""Design settings the command line has that the GUI must also expose.

Probe spacing, the Tm ceiling,
a homopolymer limit above 25, FASTA header parsing, the seventh ranking weight,
null-model calibration and the strain-unique junction duplex. Plus the tile
size box, which offered 15-200 bp while every chemistry drops any tile that is
not 50 bp (pipeline.py, "Skip tiles that aren't exactly the chemistry's
expected length"), so any other value would produce no probes.
"""
import importlib.util
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture()
def dlg(qapp, tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / "home"))
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    d = NewQueryDialog()
    d.fasta_edit.setText(str(tmp_path / "in.fa"))
    d.out_edit.setText(str(tmp_path / "out"))
    yield d
    d.close()


def _cfg(d):
    return d.build_runs()[0][1]


def test_an_untouched_form_still_builds_the_engine_defaults(dlg):
    from flexprobe.models import FilterThresholds, RunConfig
    cfg = _cfg(dlg)
    r = RunConfig(fasta_path="", output_dir="")
    for f in ("min_probe_gap_bp", "relaxed_probe_gap_bp", "floor_probe_gap_bp",
              "split_header_on_pipe", "null_model_shuffles",
              "null_model_max_windows", "null_model_seed",
              "strain_unique_junction_duplex_bp", "ranking_weights"):
        assert getattr(cfg, f) == getattr(r, f), f
    assert cfg.filters.tm_max == FilterThresholds().tm_max
    cfg.validate()


def test_each_new_control_reaches_the_config(dlg):
    dlg.min_gap.setValue(300)
    dlg.relaxed_gap.setValue(80)
    dlg.floor_gap.setValue(20)
    dlg.tm_max.setValue(70.0)
    dlg.max_hp.setValue(50)
    dlg.split_header.setChecked(True)
    dlg.null_shuffles.setValue(20)
    dlg.null_windows.setValue(150)
    dlg.null_seed.setValue(7)
    cfg = _cfg(dlg)
    assert (cfg.min_probe_gap_bp, cfg.relaxed_probe_gap_bp,
            cfg.floor_probe_gap_bp) == (300, 80, 20)
    assert cfg.filters.tm_max == 70.0
    assert cfg.filters.max_homopolymer_run == 50
    assert cfg.split_header_on_pipe is True
    assert (cfg.null_model_shuffles, cfg.null_model_max_windows,
            cfg.null_model_seed) == (20, 150, 7)
    cfg.validate()


def test_a_tm_ceiling_of_zero_means_off(dlg):
    dlg.tm_max.setValue(0.0)
    assert _cfg(dlg).filters.tm_max is None


def test_the_spacing_rungs_cannot_be_ordered_wrongly(dlg):
    """floor <= relaxed <= target. An equal value removes that rung, which is
    exactly what the tooltip says (the old CLI help said it 'disabled the
    relaxation', which was false)."""
    dlg.min_gap.setValue(100)
    dlg.relaxed_gap.setValue(500)
    assert dlg.relaxed_gap.value() <= dlg.min_gap.value()
    dlg.floor_gap.setValue(500)
    assert dlg.floor_gap.value() <= dlg.relaxed_gap.value()


def test_the_seventh_ranking_weight_is_exposed(dlg):
    from flexprobe.ranking import DEFAULT_WEIGHTS
    sp = dlg._ranking_spins["mismatch_placement"]
    assert sp.value() == pytest.approx(DEFAULT_WEIGHTS.mismatch_placement)
    dlg.ranking_box.setChecked(True)
    sp.setValue(0.4)
    assert _cfg(dlg).ranking_weights["mismatch_placement"] == pytest.approx(0.4)


def test_tile_size_is_locked_to_the_chemistry(dlg):
    from flexprobe.chemistry import get_chemistry
    from flexprobe_gui.modality import MODALITY_PRESETS
    for label in MODALITY_PRESETS:
        dlg.modality.setCurrentText(label)
        tl = get_chemistry(MODALITY_PRESETS[label]).tile_length
        assert (dlg.tile_size.minimum(), dlg.tile_size.maximum()) == (tl, tl)
        assert _cfg(dlg).tile_length == tl


def test_strain_unique_junction_duplex_is_settable(dlg):
    from flexprobe_gui._nq_build import STRAIN_MODE_HOMOLOGY, STRATEGY_UNIQUE
    dlg.strain_mode_combo.setCurrentText(STRAIN_MODE_HOMOLOGY)
    dlg.strategy_combo.setCurrentText(STRATEGY_UNIQUE)
    assert not dlg.unique_jd_row.isHidden()
    assert (dlg.unique_jd_bp.minimum(), dlg.unique_jd_bp.maximum()) == (6, 25)
    dlg.unique_jd_bp.setValue(11)
    assert _cfg(dlg).strain_unique_junction_duplex_bp == 11
    dlg.unique_jd.setChecked(False)
    assert _cfg(dlg).strain_unique_junction_duplex_bp == 0


def test_the_strain_unique_duplex_row_is_hidden_outside_that_mode(dlg):
    assert dlg.unique_jd_row.isHidden()


def test_the_form_refuses_settings_the_engine_refuses(dlg, monkeypatch):
    """Tm ceiling below the floor: RunConfig.validate refuses it, so the form
    must say so before accepting, not leave the run to fail after Start."""
    from PySide6.QtWidgets import QMessageBox
    shown = []
    monkeypatch.setattr(QMessageBox, "warning",
                        lambda *a, **k: shown.append(a[2] if len(a) > 2 else ""))
    dlg.tm_min.setValue(60.0)
    dlg.tm_max.setValue(55.0)
    assert dlg._validate_inputs() is False
    assert shown and "tm_max" in shown[0], shown


def test_a_draft_and_a_re_run_restore_every_new_setting(dlg, qapp, tmp_path):
    from flexprobe_gui._nq_build import summary_to_draft
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    dlg.min_gap.setValue(250)
    dlg.relaxed_gap.setValue(60)
    dlg.floor_gap.setValue(25)
    dlg.tm_max.setValue(68.0)
    dlg.split_header.setChecked(True)
    dlg.null_shuffles.setValue(10)
    for source in ("draft", "record"):
        d = (dlg.to_draft() if source == "draft"
             else {**dlg.to_draft(),
                   **summary_to_draft(dlg.build_runs()[0][0].config_summary)})
        fresh = NewQueryDialog()
        try:
            fresh.apply_draft(d)
            fresh.fasta_edit.setText(str(tmp_path / "in.fa"))
            fresh.out_edit.setText(str(tmp_path / "out"))
            c = _cfg(fresh)
            assert (c.min_probe_gap_bp, c.relaxed_probe_gap_bp,
                    c.floor_probe_gap_bp) == (250, 60, 25), source
            assert c.filters.tm_max == 68.0, source
            assert c.split_header_on_pipe is True, source
            assert c.null_model_shuffles == 10, source
        finally:
            fresh.close()
