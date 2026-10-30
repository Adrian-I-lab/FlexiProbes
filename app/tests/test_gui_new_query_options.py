"""GUI: the new RunConfig options are wired into the New-query dialog.

Skipped when PySide6 is unavailable. Headless via offscreen Qt.
"""
import importlib.util
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)


def _dialog(qapp, tmp_path):
    # Point all stores at a temp dir so the test never touches ~/.flexprobe.
    os.environ["FLEXPROBE_HOME"] = str(tmp_path)
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    return NewQueryDialog()


def test_constitutive_toggle_shows_isoform_min(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    dlg.strain_mode_combo.setCurrentIndex(1)            # multi-strain
    # find + select the constitutive strategy
    idx = dlg.strategy_combo.findText(
        next(dlg.strategy_combo.itemText(i)
             for i in range(dlg.strategy_combo.count())
             if dlg.strategy_combo.itemText(i).startswith("Constitutive"))
    )
    dlg.strategy_combo.setCurrentIndex(idx)
    # The dialog is never shown headless, so isVisible() is always False;
    # isHidden() reflects the explicit setVisible() toggle we care about.
    # Checked on the ROWS: the whole row (spin box plus its unit text) is
    # hidden, so a hidden row cannot leave its unit label floating in the
    # form. See test_gui_new_query_layout.py.
    assert not dlg.isoform_min_row.isHidden()
    assert dlg.unique_min_mm_row.isHidden()
    # the hint warns about the alignment caveat for constitutive mode
    hint = dlg.strain_hint.text()
    assert "isoform-collapse" in hint and "MAFFT" in hint


def test_new_options_draft_round_trip(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    dlg.workers.setValue(6)
    dlg.extra_on_target.setPlainText("ISO_A\nISO_B")
    dlg.ranking_box.setChecked(True)
    dlg._ranking_spins["specificity_margin"].setValue(0.5)
    dlg.strain_mode_combo.setCurrentIndex(1)
    sidx = next(i for i in range(dlg.strategy_combo.count())
               if dlg.strategy_combo.itemText(i).startswith("Constitutive"))
    dlg.strategy_combo.setCurrentIndex(sidx)
    dlg.isoform_min.setValue(3)

    draft = dlg.to_draft()
    assert draft["workers"] == 6
    assert draft["extra_on_target"] == "ISO_A\nISO_B"
    assert draft["ranking_enabled"] is True
    assert draft["ranking_weights"]["specificity_margin"] == 0.5
    assert draft["isoform_min"] == 3
    assert draft["strategy"].startswith("Constitutive")

    # A fresh dialog restores every new field from the draft.
    dlg2 = _dialog(qapp, tmp_path)
    dlg2.apply_draft(draft)
    assert dlg2.workers.value() == 6
    assert dlg2.extra_on_target.toPlainText() == "ISO_A\nISO_B"
    assert dlg2.ranking_box.isChecked()
    assert dlg2._ranking_spins["specificity_margin"].value() == 0.5
    assert dlg2.isoform_min.value() == 3


def test_scalar_fields_round_trip_via_shared_spec(qapp, tmp_path):
    # Every scalar spin/check field flows through one shared spec, so
    # to_draft/apply_draft can't drift. Set non-default values, round-trip, and
    # confirm each restores.
    dlg = _dialog(qapp, tmp_path)
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    dlg.tile_step.setValue(7)
    dlg.gc_min.setValue(0.40)
    dlg.max_hp.setValue(6)
    dlg.workers.setValue(4)
    dlg.screen_top_per_gene.setValue(50)
    dlg.group_by_gene.setChecked(True)
    dlg.skip_screen.setChecked(True)
    draft = dlg.to_draft()
    # Every spec field is present in the draft.
    for key, _attr, _kind, _default in NewQueryDialog._SCALAR_FIELDS:
        assert key in draft
    dlg2 = _dialog(qapp, tmp_path)
    dlg2.apply_draft(draft)
    assert dlg2.tile_step.value() == 7
    assert dlg2.gc_min.value() == 0.40
    assert dlg2.max_hp.value() == 6
    assert dlg2.workers.value() == 4
    assert dlg2.screen_top_per_gene.value() == 50
    assert dlg2.group_by_gene.isChecked()
    assert dlg2.skip_screen.isChecked()


def test_primer3_default_off_when_key_missing_in_draft(qapp, tmp_path):
    # A legacy draft must not silently change the primer3 setting on
    # round-trip. The two defaults are DELIBERATELY different, which makes
    # keeping them straight more important, not less:
    #
    #   fresh dialog      -> ON  (modality preset; matches the engine default)
    #   draft w/o the key -> OFF (it describes a run made WITHOUT primer3, and
    #                             reproducing it must not change what it did)
    dlg = _dialog(qapp, tmp_path)
    assert dlg.primer3.isChecked() is True         # fresh dialog: engine default
    dlg.apply_draft({})                            # draft predating the key
    assert dlg.primer3.isChecked() is False, (
        "a legacy draft must not inherit the new ON default — it describes a "
        "run that had primer3 off"
    )
    # Both explicit values still round-trip.
    dlg.apply_draft({"primer3": True})
    assert dlg.primer3.isChecked() is True
    dlg.apply_draft({"primer3": False})
    assert dlg.primer3.isChecked() is False


def test_group_by_gene_options_round_trip_and_reach_config(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    # Grouping is ON by default, so the gate starts enabled.
    assert dlg.group_by_gene.isChecked()
    assert dlg.min_isoform_coverage.isEnabled()
    # The coverage gate is disabled while grouping is off.
    dlg.group_by_gene.setChecked(False)
    assert not dlg.min_isoform_coverage.isEnabled()
    dlg.group_by_gene.setChecked(True)
    assert dlg.min_isoform_coverage.isEnabled()
    dlg.min_isoform_coverage.setValue(0.8)

    draft = dlg.to_draft()
    assert draft["group_by_gene"] is True
    assert draft["min_isoform_coverage"] == 0.8

    # Round-trips into a fresh dialog.
    dlg2 = _dialog(qapp, tmp_path)
    dlg2.apply_draft(draft)
    assert dlg2.group_by_gene.isChecked()
    assert dlg2.min_isoform_coverage.value() == 0.8

    # And reaches the built RunConfig.
    fa = tmp_path / "g.fa"
    fa.write_text(">G1\n" + "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCATGCATG")
    dlg2.fasta_edit.setText(str(fa))
    dlg2.out_edit.setText(str(tmp_path / "out"))
    cfg = dlg2.build_runs()[0][1]
    assert cfg.group_by_gene is True
    assert cfg.min_isoform_coverage == 0.8


def test_isoform_coverage_is_a_tunable_ranking_weight(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    assert "isoform_coverage" in dlg._ranking_spins
    iso = dlg._ranking_spins["isoform_coverage"]
    assert iso.value() == 0.20
    # Greyed out while gene grouping is off (it is on by default), then
    # editable again when it is turned back on.
    dlg.group_by_gene.setChecked(False)
    assert not iso.isEnabled()
    dlg.group_by_gene.setChecked(True)
    assert iso.isEnabled()

    # With the ranking panel on, an explicit value flows to the RunConfig.
    dlg.ranking_box.setChecked(True)
    iso.setValue(0.40)
    draft = dlg.to_draft()
    assert draft["ranking_weights"]["isoform_coverage"] == 0.40

    fa = tmp_path / "g.fa"
    fa.write_text(">G1\n" + "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCATGCATG")
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(str(tmp_path / "out"))
    cfg = dlg.build_runs()[0][1]
    assert cfg.ranking_weights["isoform_coverage"] == 0.40
    cfg.validate()   # isoform_coverage is an accepted ranking-weight key


def test_flex_v1_multiplex_defaults_to_four_barcodes(qapp, tmp_path):
    # Regression: selecting Flex v1 + Multiplex must produce a genuinely
    # multiplexed (multi-barcode) config from the sample-count selector, which
    # defaults to 4 samples (BC001-BC004) — never silently a single barcode.
    dlg = _dialog(qapp, tmp_path)
    for i in range(dlg.modality.count()):
        if dlg.modality.itemText(i).startswith("Flex v1"):
            dlg.modality.setCurrentIndex(i)
            break
    dlg.mode.setCurrentText("Multiplex")
    assert dlg._is_flex_v1_multiplex()
    # defaults to 4 samples -> BC001..BC004, no manual entry required
    assert dlg._selected_multi_barcodes() == ["BC001", "BC002", "BC003", "BC004"]

    fa = tmp_path / "g.fa"
    fa.write_text(">G1\n" + "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCATGCATG")
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(str(tmp_path / "out"))
    cfg = dlg.build_runs()[0][1]
    assert cfg.multiplex is True
    assert cfg.probe_barcodes == ["BC001", "BC002", "BC003", "BC004"]


def test_flex_v1_multiplex_sample_count(qapp, tmp_path):
    # The sample-count selector maps N samples -> BC001..BC00N.
    dlg = _dialog(qapp, tmp_path)
    for i in range(dlg.modality.count()):
        if dlg.modality.itemText(i).startswith("Flex v1"):
            dlg.modality.setCurrentIndex(i)
            break
    dlg.mode.setCurrentText("Multiplex")
    dlg.multi_count.setValue(16)
    bcs = dlg._selected_multi_barcodes()
    assert len(bcs) == 16
    assert bcs[0] == "BC001" and bcs[-1] == "BC016"
    dlg.multi_count.setValue(8)
    assert dlg._selected_multi_barcodes() == [f"BC{i:03d}" for i in range(1, 9)]


def test_build_runs_threads_new_options_into_config(qapp, tmp_path):
    fa = tmp_path / "toy.fasta"
    fa.write_text(
        ">toy|gene|ACC1\n"
        "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
        "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG\n"
    )
    dlg = _dialog(qapp, tmp_path)
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(str(tmp_path / "out"))
    dlg.workers.setValue(2)
    dlg.extra_on_target.setPlainText("ISO_X")
    dlg.ranking_box.setChecked(True)
    dlg._ranking_spins["gc_balance"].setValue(0.4)

    runs = dlg.build_runs()
    assert runs, "expected at least one (record, config) run"
    cfg = runs[0][1]
    assert cfg.workers == 2
    assert cfg.extra_on_target_ids == ["ISO_X"]
    assert cfg.ranking_weights is not None
    assert cfg.ranking_weights["gc_balance"] == 0.4
    # the assembled config must pass validation
    cfg.validate()


def test_strain_unique_per_half_defaults_on_and_reaches_config(qapp, tmp_path):
    """The GUI must offer the same per-half strain-discrimination choice as
    --strain-unique-whole-window, and default to the per-half reading."""
    dlg = _dialog(qapp, tmp_path)
    dlg.strain_mode_combo.setCurrentIndex(1)            # multi-strain
    idx = next(i for i in range(dlg.strategy_combo.count())
               if dlg.strategy_combo.itemText(i).startswith("Sequence-unique"))
    dlg.strategy_combo.setCurrentIndex(idx)

    # Shown only for the strain-unique strategy, and ticked by default.
    assert not dlg.unique_per_half.isHidden()
    assert dlg.unique_per_half.isChecked()

    # Hidden again once the strategy no longer uses it.
    idx_c = next(i for i in range(dlg.strategy_combo.count())
                 if dlg.strategy_combo.itemText(i).startswith("Constitutive"))
    dlg.strategy_combo.setCurrentIndex(idx_c)
    assert dlg.unique_per_half.isHidden()


def test_strain_unique_per_half_round_trips_through_a_draft(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    dlg.unique_per_half.setChecked(False)
    draft = dlg.to_draft()
    assert draft["unique_per_half"] is False

    dlg2 = _dialog(qapp, tmp_path)
    dlg2.apply_draft(draft)
    assert dlg2.unique_per_half.isChecked() is False

    # A draft saved before this option existed must restore the safe default,
    # not silently fall back to the legacy whole-window count.
    dlg3 = _dialog(qapp, tmp_path)
    draft.pop("unique_per_half")
    dlg3.apply_draft(draft)
    assert dlg3.unique_per_half.isChecked() is True


def _two_strain_fasta(tmp_path):
    """Two records of one fragment — the shape sequence-unique mode is for."""
    fa = tmp_path / "two_strain.fasta"
    a = ("GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
         "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG")
    b = ("GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
         "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATA")
    fa.write_text(f">STRAIN_A|Large|ACC1\n{a}\n>STRAIN_B|Large|ACC2\n{b}\n")
    return fa


def _select_sequence_unique(dlg):
    dlg.strain_mode_combo.setCurrentIndex(1)            # Homology consensus
    idx = next(i for i in range(dlg.strategy_combo.count())
               if dlg.strategy_combo.itemText(i).startswith("Sequence-unique"))
    dlg.strategy_combo.setCurrentIndex(idx)


def test_strain_unique_target_survives_a_closed_dialog(qapp, tmp_path):
    """The picked target sequence must reach RunConfig.

    build_runs() is called by home_window AFTER exec() returns — i.e. on a
    dialog that has already been closed — and QWidget.isVisible() is False for
    every child of a hidden window. Reading the picker through an isVisible()
    guard therefore forced strain_unique_target=0 on every run regardless of
    the user's choice, so the GUI could only ever design probes unique to the
    FIRST record. The pipeline itself was fine: the same config at target 0 vs
    1 yields different probe sets. Guard against the widget-visibility form
    coming back.
    """
    fa = _two_strain_fasta(tmp_path)
    dlg = _dialog(qapp, tmp_path)
    _select_sequence_unique(dlg)
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(str(tmp_path / "out"))

    # The picker is populated from the FASTA even though the dialog was never
    # shown (the second route to a silent 0 was an unpopulated combo).
    # Two records plus the trailing "Every sequence" fan-out entry.
    assert dlg.unique_target_combo.count() == 3
    assert dlg.unique_target_combo.itemText(1) == "STRAIN_B|Large|ACC2"

    dlg.unique_target_combo.setCurrentIndex(1)          # pick the SECOND record
    assert dlg.unique_target_combo.currentData() == 1

    # Headless, the dialog is never shown, so isVisible() is already False —
    # exactly the state build_runs() sees in the real app after exec() returns.
    assert not dlg.unique_target_combo.isVisible()

    runs = dlg.build_runs()
    assert runs
    cfg = runs[0][1]
    assert cfg.consensus_strategy == "strain_unique"
    assert cfg.strain_unique_target == 1, (
        "the picked target sequence was dropped between the dialog and RunConfig"
    )
    cfg.validate()


def test_strain_unique_target_is_zero_for_other_strategies(qapp, tmp_path):
    """Only sequence-unique mode consults the picker; others pin target 0."""
    fa = _two_strain_fasta(tmp_path)
    dlg = _dialog(qapp, tmp_path)
    _select_sequence_unique(dlg)
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(str(tmp_path / "out"))
    dlg.unique_target_combo.setCurrentIndex(1)

    idx_c = next(i for i in range(dlg.strategy_combo.count())
                 if dlg.strategy_combo.itemText(i).startswith("Constitutive"))
    dlg.strategy_combo.setCurrentIndex(idx_c)

    cfg = dlg.build_runs()[0][1]
    assert cfg.consensus_strategy == "constitutive"
    assert cfg.strain_unique_target == 0


def _select_strategy(dlg, prefix):
    dlg.strain_mode_combo.setCurrentIndex(1)            # Homology consensus
    idx = next(i for i in range(dlg.strategy_combo.count())
               if dlg.strategy_combo.itemText(i).startswith(prefix))
    dlg.strategy_combo.setCurrentIndex(idx)


def test_consensus_template_defaults_to_fan_out_and_lists_records(qapp, tmp_path):
    """The template picker offers every record plus an "all" fan-out entry."""
    from flexprobe_gui._nq_build import ALL_TEMPLATES
    fa = _two_strain_fasta(tmp_path)
    dlg = _dialog(qapp, tmp_path)
    _select_strategy(dlg, "Conserved")
    dlg.fasta_edit.setText(str(fa))

    assert not dlg.template_combo.isHidden()
    # fan-out entry + one per record
    assert dlg.template_combo.count() == 3
    assert dlg.template_combo.itemData(0) == ALL_TEMPLATES
    assert dlg.template_combo.itemData(1) == 0
    assert dlg.template_combo.itemData(2) == 1
    assert dlg.template_combo.itemText(1) == "STRAIN_A|Large|ACC1"

    # Hidden for sequence-unique, which picks its sequence via "Target sequence".
    _select_strategy(dlg, "Sequence-unique")
    assert dlg.template_combo.isHidden()


def test_single_template_reaches_config(qapp, tmp_path):
    fa = _two_strain_fasta(tmp_path)
    dlg = _dialog(qapp, tmp_path)
    _select_strategy(dlg, "Constitutive")
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(str(tmp_path / "out"))
    dlg.template_combo.setCurrentIndex(2)               # the SECOND record
    assert dlg.template_combo.currentData() == 1

    runs = dlg.build_runs()
    assert len(runs) == 1
    cfg = runs[0][1]
    assert cfg.consensus_reference == 1
    cfg.validate()


def test_every_template_fans_out_one_run_per_record(qapp, tmp_path):
    """A single consensus run only ever emits the template's flavour of each
    variant window, so the other records' equally valid probes are never
    generated. The fan-out is what makes the complete set reachable."""
    from flexprobe_gui._nq_build import ALL_TEMPLATES
    fa = _two_strain_fasta(tmp_path)
    dlg = _dialog(qapp, tmp_path)
    _select_strategy(dlg, "Conserved")
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(str(tmp_path / "out"))
    dlg.template_combo.setCurrentIndex(
        dlg.template_combo.findData(ALL_TEMPLATES))

    runs = dlg.build_runs()
    assert len(runs) == 2, "expected one run per input record"
    refs = sorted(cfg.consensus_reference for _rec, cfg, _c in runs)
    assert refs == [0, 1]
    # Output dirs must not collide, and each run is named for its template.
    outs = {cfg.output_dir for _rec, cfg, _c in runs}
    assert len(outs) == 2
    assert any("STRAIN_A" in o for o in outs)
    assert any("STRAIN_B" in o for o in outs)
    for _rec, cfg, _c in runs:
        cfg.validate()


def test_fan_out_is_per_fasta_in_bulk(qapp, tmp_path):
    """Record counts differ between files, so the fan-out must read each FASTA
    rather than assume the first one's record count."""
    from flexprobe_gui._nq_build import ALL_TEMPLATES
    two = _two_strain_fasta(tmp_path)
    one = tmp_path / "single.fasta"
    one.write_text(">SOLO|Large|ACC9\n" + "ACGT" * 25 + "\n")
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    os.environ["FLEXPROBE_HOME"] = str(tmp_path)
    dlg = NewQueryDialog(bulk_fasta_paths=[str(two), str(one)])
    _select_strategy(dlg, "Conserved")
    dlg.out_edit.setText(str(tmp_path / "out"))
    dlg.template_combo.setCurrentIndex(
        dlg.template_combo.findData(ALL_TEMPLATES))

    runs = dlg.build_runs()
    # 2 templates for the two-record file + 1 for the single-record file.
    assert len(runs) == 3
    by_fasta = {}
    for _rec, cfg, _c in runs:
        by_fasta.setdefault(os.path.basename(cfg.fasta_path), []).append(
            cfg.consensus_reference)
    assert sorted(by_fasta["two_strain.fasta"]) == [0, 1]
    assert by_fasta["single.fasta"] == [0]


def test_consensus_template_round_trips_through_a_draft(qapp, tmp_path):
    from flexprobe_gui._nq_build import ALL_TEMPLATES
    fa = _two_strain_fasta(tmp_path)
    dlg = _dialog(qapp, tmp_path)
    _select_strategy(dlg, "Conserved")
    dlg.fasta_edit.setText(str(fa))
    dlg.template_combo.setCurrentIndex(
        dlg.template_combo.findData(ALL_TEMPLATES))
    draft = dlg.to_draft()
    assert draft["consensus_template"] == ALL_TEMPLATES

    dlg2 = _dialog(qapp, tmp_path)
    dlg2.apply_draft(draft)
    dlg2.fasta_edit.setText(str(fa))
    assert dlg2.template_combo.currentData() == ALL_TEMPLATES


# ---------------------------------------------------------------------------
# Off-target engine selector (k-mer is the default)
# ---------------------------------------------------------------------------

def test_engine_selector_offers_both_and_defaults_to_the_engine_default(
        qapp, tmp_path):
    """Desktop users must be able to reach the fallback engine.

    The k-mer default refuses a reference over ~2.1 Gbases and its error names
    `--screen-engine blast` — a CLI flag a GUI user cannot reach. Without this
    control the people most likely to point at a genome have no way out.
    """
    from flexprobe.constants import SCREEN_ENGINE_DEFAULT
    from flexprobe.screen import ENGINE_BLAST, ENGINE_KMER
    dlg = _dialog(qapp, tmp_path)
    offered = {dlg.screen_engine.itemData(i)
               for i in range(dlg.screen_engine.count())}
    assert offered == {ENGINE_KMER, ENGINE_BLAST}
    assert dlg.screen_engine.currentData() == SCREEN_ENGINE_DEFAULT


def test_engine_choice_round_trips_through_a_draft(qapp, tmp_path):
    from flexprobe.screen import ENGINE_BLAST
    dlg = _dialog(qapp, tmp_path)
    idx = dlg.screen_engine.findData(ENGINE_BLAST)
    dlg.screen_engine.setCurrentIndex(idx)
    draft = dlg.to_draft()
    assert draft["screen_engine"] == ENGINE_BLAST

    dlg2 = _dialog(qapp, tmp_path)
    dlg2.apply_draft(draft)
    assert dlg2.screen_engine.currentData() == ENGINE_BLAST


def test_draft_without_an_engine_key_restores_blast(qapp, tmp_path):
    """A draft saved before the setting existed describes a BLAST run.

    Same reasoning as the `primer3` default in _nq_draft._SCALAR_FIELDS:
    restoring it as today's default would change the run the draft claims to
    reproduce.
    """
    from flexprobe.screen import ENGINE_BLAST
    dlg = _dialog(qapp, tmp_path)
    draft = dlg.to_draft()
    draft.pop("screen_engine", None)
    dlg2 = _dialog(qapp, tmp_path)
    dlg2.apply_draft(draft)
    assert dlg2.screen_engine.currentData() == ENGINE_BLAST


def test_selected_engine_reaches_the_built_config(qapp, tmp_path):
    """The widget must actually change the run, not just the form."""
    from flexprobe.screen import ENGINE_BLAST
    fa = tmp_path / "t.fa"
    fa.write_text(">g\n" + "ACGT" * 40 + "\n")
    dlg = _dialog(qapp, tmp_path)
    dlg.fasta_edit.setText(str(fa))
    dlg.out_edit.setText(str(tmp_path / "out"))
    dlg.db_edit.setPlainText(str(fa))
    dlg.screen_engine.setCurrentIndex(dlg.screen_engine.findData(ENGINE_BLAST))
    runs = dlg.build_runs()
    assert runs, "expected the dialog to build a run"
    cfg = runs[0][1] if isinstance(runs[0], tuple) else runs[0]
    cfg = getattr(cfg, "config", cfg)
    assert cfg.screen.engine == ENGINE_BLAST
