"""The New-query form exposes the off-target specificity rules.

Cross-subject pooling, either-half rejection and the junction-duplex
threshold are settable from the command line, and the form must be able to
set them too rather than always building SpecificityThresholds() and the
engine's default reference handling. These tests pin that the form can
express every one of those settings, that an untouched form still builds
exactly the engine defaults, and that a re-run restores them.
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


def _pick(combo, data):
    i = combo.findData(data)
    assert i >= 0, data
    combo.setCurrentIndex(i)


def test_an_untouched_form_builds_the_engine_defaults(dlg):
    from flexprobe.models import ScreenResources, SpecificityThresholds
    cfg = _cfg(dlg)
    assert cfg.specificity == SpecificityThresholds()
    dflt = ScreenResources()
    assert cfg.screen.exclude_alt_scaffolds == dflt.exclude_alt_scaffolds
    assert cfg.screen.reference_strand == dflt.reference_strand
    cfg.validate()


def test_either_half_rule_builds_the_either_half_config(dlg):
    """Either half: off_target_require_both_halves False, same_subject left
    at its default (it has no effect then)."""
    from flexprobe_gui._nq_build import OT_RULE_EITHER_HALF
    _pick(dlg.ot_rule, OT_RULE_EITHER_HALF)
    sp = _cfg(dlg).specificity
    assert sp.off_target_require_both_halves is False
    assert sp.off_target_same_subject is True


def test_cross_subject_rule_builds_the_pooled_config(dlg):
    """Whole probe: both halves, pooled across subjects."""
    from flexprobe_gui._nq_build import OT_RULE_CROSS_SUBJECT
    _pick(dlg.ot_rule, OT_RULE_CROSS_SUBJECT)
    sp = _cfg(dlg).specificity
    assert sp.off_target_require_both_halves is True
    assert sp.off_target_same_subject is False


def test_every_other_control_reaches_the_config(dlg):
    dlg.ot_min_mm.setValue(3)
    dlg.junction_duplex_bp.setValue(12)
    dlg.exclude_alt.setChecked(False)
    _pick(dlg.ref_strand, "plus")
    cfg = _cfg(dlg)
    assert cfg.specificity.off_target_min_mismatch == 3
    assert cfg.specificity.off_target_junction_duplex_bp == 12
    assert cfg.screen.exclude_alt_scaffolds is False
    assert cfg.screen.reference_strand == "plus"
    cfg.validate()
    dlg.junction_duplex.setChecked(False)
    assert _cfg(dlg).specificity.off_target_junction_duplex_bp == 0


def test_the_controls_cannot_build_a_config_the_engine_refuses(dlg):
    """Ranges follow the engine: duplex 0 or 6-25 bp (models.validate refuses
    1-5, and a 25 nt half cannot hold more), mismatches 1-6 (the k-mer screen
    reports at most 5 mismatches, screen.py MAX_MISMATCH_REPORTED)."""
    assert (dlg.junction_duplex_bp.minimum(), dlg.junction_duplex_bp.maximum()) == (6, 25)
    assert (dlg.ot_min_mm.minimum(), dlg.ot_min_mm.maximum()) == (1, 6)


def test_blast_switches_the_duplex_filter_off_and_disables_it(dlg, qapp):
    from flexprobe.screen import ENGINE_BLAST, ENGINE_KMER
    _pick(dlg.screen_engine, ENGINE_BLAST)
    qapp.processEvents()
    assert not dlg.junction_duplex.isEnabled()
    cfg = _cfg(dlg)
    assert cfg.specificity.off_target_junction_duplex_bp == 0
    cfg.validate()
    _pick(dlg.screen_engine, ENGINE_KMER)
    qapp.processEvents()
    assert dlg.junction_duplex.isEnabled()


def test_local_alignment_counting_is_offered_only_under_blast(dlg, qapp):
    """Under k-mer every hit is a full 25 nt match with <= 5 mismatches, so
    the 80%/20 bp local rule would not be a local alignment at all."""
    from flexprobe.screen import ENGINE_BLAST, ENGINE_KMER
    from flexprobe_gui._nq_build import OT_COUNT_LOCAL
    idx = dlg.ot_counting.findData(OT_COUNT_LOCAL)
    model = dlg.ot_counting.model()
    assert not model.item(idx).isEnabled()
    _pick(dlg.screen_engine, ENGINE_BLAST)
    qapp.processEvents()
    assert model.item(idx).isEnabled()
    dlg.ot_counting.setCurrentIndex(idx)
    assert _cfg(dlg).specificity.off_target_full_half is False
    _pick(dlg.screen_engine, ENGINE_KMER)
    qapp.processEvents()
    assert dlg.ot_counting.currentData() != OT_COUNT_LOCAL
    assert _cfg(dlg).specificity.off_target_full_half is True


def test_a_draft_round_trips_every_specificity_control(dlg, qapp, tmp_path):
    from flexprobe_gui._nq_build import OT_RULE_CROSS_SUBJECT
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    _pick(dlg.ot_rule, OT_RULE_CROSS_SUBJECT)
    dlg.ot_min_mm.setValue(4)
    dlg.junction_duplex_bp.setValue(14)
    dlg.exclude_alt.setChecked(False)
    _pick(dlg.ref_strand, "plus")
    draft = dlg.to_draft()
    fresh = NewQueryDialog()
    try:
        fresh.apply_draft(draft)
        a, b = _cfg(dlg), None
        fresh.fasta_edit.setText(str(tmp_path / "in.fa"))
        fresh.out_edit.setText(str(tmp_path / "out"))
        b = _cfg(fresh)
        assert b.specificity == a.specificity
        assert b.screen.exclude_alt_scaffolds is False
        assert b.screen.reference_strand == "plus"
    finally:
        fresh.close()


def test_a_re_run_restores_the_recorded_specificity(dlg):
    from flexprobe_gui._nq_build import OT_RULE_EITHER_HALF, summary_to_draft
    from flexprobe_gui._nq_build import OT_RULE_SAME_SUBJECT
    _pick(dlg.ot_rule, OT_RULE_EITHER_HALF)
    dlg.junction_duplex.setChecked(False)
    rec = dlg.build_runs()[0][0]
    draft = summary_to_draft(rec.config_summary)
    assert draft["ot_rule"] == OT_RULE_EITHER_HALF
    assert draft["junction_duplex"] is False
    # A record made before these settings existed was built from the engine
    # defaults, so that is what a re-run restores.
    old = summary_to_draft({})
    assert old["ot_rule"] == OT_RULE_SAME_SUBJECT
    assert old["junction_duplex"] is True


def test_an_imported_either_half_run_restores_on_re_run():
    """An HPC run imported from its run_config.json records the same keys a
    run made here does (one shared writer), so Re-run restores it."""
    from dataclasses import replace
    from flexprobe.models import RunConfig, SpecificityThresholds
    from flexprobe_gui._nq_build import OT_RULE_EITHER_HALF, summary_to_draft
    from flexprobe_gui.run_import import summary_from_config
    cfg = RunConfig(fasta_path="x.fa", output_dir="o")
    cfg = replace(cfg, specificity=SpecificityThresholds(
        off_target_require_both_halves=False, off_target_junction_duplex_bp=0))
    draft = summary_to_draft(summary_from_config(cfg))
    assert draft["ot_rule"] == OT_RULE_EITHER_HALF
    assert draft["junction_duplex"] is False
