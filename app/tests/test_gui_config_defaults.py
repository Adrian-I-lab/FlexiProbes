# flexiprobes — probe design for 10x Genomics Flex v1 / Flex v2 / Visium.
# Copyright (C) 2026 Adrian Ilich and contributors
#
# This file is part of flexiprobes.
#
# flexiprobes is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# flexiprobes is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""The GUI must not keep its own copy of the engine's defaults.

Two ways the app has drifted from the engine while saying nothing:

  * ``MODALITY_DEFAULTS`` restated all 21 engine-owned values (7 filter
    thresholds x 3 modalities) as literals. Two of them went stale.
  * ``HomeWindow._record_to_draft`` — the "Re-run" path for a record saved
    before drafts existed — carried a hand-written subset of the keys the
    dialog consumes, so eight design-relevant settings recorded in
    ``config_summary`` were silently dropped on the way back into the form.

Both are pinned REFLECTIVELY here, over ``dataclasses.fields``, because a
hand-listed test drifts exactly the way the hand-listed code did.
"""
from __future__ import annotations

import dataclasses

import pytest

pytest.importorskip("PySide6")

from flexprobe.models import FilterThresholds

from flexprobe_gui.modality import MODALITY_DEFAULTS


# ---------------------------------------------------------------------------
# MODALITY_DEFAULTS is the engine's defaults, not a copy of them
# ---------------------------------------------------------------------------

# A modality entry may legitimately differ from the engine on a NON-filter key:
# Flex v1 embeds a per-sample Probe Barcode, so its preset opens singleplex.
# Nothing in CG000621/CG000839 asks the three chemistries to FILTER differently,
# so no FilterThresholds field belongs in here. Adding one is a decision that
# needs a tech-note citation, not a quiet literal.
_CHEMISTRY_OVERRIDES_ALLOWED: dict = {}


@pytest.mark.parametrize("label", sorted(MODALITY_DEFAULTS))
@pytest.mark.parametrize("field", [f.name
                                   for f in dataclasses.fields(FilterThresholds)])
def test_modality_preset_carries_the_engine_filter_default(label, field):
    """Derived from FilterThresholds, not hand-listed.

    A field added to the engine gets its own case automatically, and a preset
    that restates a value the engine has since changed fails here rather than
    designing differently from the CLI in silence.
    """
    assert field not in _CHEMISTRY_OVERRIDES_ALLOWED, (
        f"{field} is listed as a chemistry override; assert its reason instead"
    )
    assert field in MODALITY_DEFAULTS[label], (
        f"{label} has no {field!r}: the preset no longer covers every engine "
        f"filter threshold, so the dialog cannot open on the engine's value."
    )
    assert MODALITY_DEFAULTS[label][field] == getattr(FilterThresholds(), field), (
        f"{label}[{field!r}] = {MODALITY_DEFAULTS[label][field]!r} but the "
        f"engine defaults to {getattr(FilterThresholds(), field)!r} — the same "
        f"untouched form would design differently in the app than on the CLI."
    )


def test_the_only_cross_modality_differences_are_chemistry_facts():
    """Whatever DOES differ between the three presets must be a chemistry fact.

    Guards the emptiness that is the point of ``_modality_defaults``: today the
    presets differ only on the Probe Barcode and the multiplex flag that follows
    from it, both of which the engine has no opinion about.
    """
    keys = set().union(*(set(d) for d in MODALITY_DEFAULTS.values()))
    differing = {
        k for k in keys
        if len({repr(d.get(k)) for d in MODALITY_DEFAULTS.values()}) > 1
    }
    assert differing == {"probe_barcode", "multiplex"}, differing


# ---------------------------------------------------------------------------
# config_summary -> draft: "Re-run" must rebuild the design it shows
# ---------------------------------------------------------------------------

@pytest.fixture()
def dialog(qapp, tmp_path):
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    dlg = NewQueryDialog()
    dlg.fasta_edit.setText(str(tmp_path / "in.fa"))
    dlg.out_edit.setText(str(tmp_path / "out"))
    return dlg


@pytest.fixture()
def default_summary(dialog):
    """What build_runs records for an untouched form — the live key list."""
    rec, _cfg, _cached = dialog.build_runs()[0]
    return dict(rec.config_summary)


def _record_to_draft(summary, modality="Flex v2 (10X GEM-X)", fasta="in.fa"):
    """Call HomeWindow._record_to_draft unbound, off a stub host.

    The method needs only ``self.settings.default_output_dir()``, so this avoids
    standing up a whole main window (and its stores) for a pure mapping test.
    """
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunRecord

    class _Settings:
        def default_output_dir(self):
            return "/tmp/out"

    class _Host:
        settings = _Settings()

    rec = RunRecord(name="r", modality_label=modality, fasta_path=str(fasta),
                    project_id="", config_summary=dict(summary))
    return HomeWindow._record_to_draft(_Host(), rec)


# --- completeness: the map cannot fall behind either side -------------------

def test_every_recorded_config_summary_key_is_mapped_or_excluded(default_summary):
    """Derived from a REAL built record, not a hand-listed set.

    A key added to build_runs' config_summary shows up here automatically and
    must be either mapped into the draft or listed as deliberately excluded —
    the alternative is what this finding was: recorded, then silently dropped.
    """
    from flexprobe_gui._nq_build import (
        CONFIG_SUMMARY_NOT_IN_DRAFT, CONFIG_SUMMARY_TO_DRAFT,
    )
    mapped = {f.summary for f in CONFIG_SUMMARY_TO_DRAFT}
    unaccounted = set(default_summary) - mapped - CONFIG_SUMMARY_NOT_IN_DRAFT
    assert not unaccounted, (
        f"config_summary records {sorted(unaccounted)} and 'Re-run' throws it "
        f"away, so the rebuilt run is a different design from the recorded one."
    )


def test_every_filter_thresholds_field_is_mapped_or_excluded():
    """The reflective guard the hand-written map needed and did not have.

    Adding a field to FilterThresholds fails this until it is either carried
    through config_summary into the draft or listed as GUI-invisible.
    """
    from flexprobe_gui._nq_build import (
        CONFIG_SUMMARY_TO_DRAFT, FILTER_FIELDS_NOT_IN_SUMMARY,
    )
    mapped = {f.summary for f in CONFIG_SUMMARY_TO_DRAFT}
    unaccounted = {
        f.name for f in dataclasses.fields(FilterThresholds)
        if f.name not in mapped and f.name not in FILTER_FIELDS_NOT_IN_SUMMARY
    }
    assert not unaccounted, (
        f"FilterThresholds grew {sorted(unaccounted)}; add it to "
        f"CONFIG_SUMMARY_TO_DRAFT (with the engine's own default) or to "
        f"FILTER_FIELDS_NOT_IN_SUMMARY with the reason it has no form field."
    )


def test_the_exclusion_lists_never_hide_a_key_that_is_mapped():
    """An exclusion must mean 'considered and excluded', never 'both'."""
    from flexprobe_gui._nq_build import (
        CONFIG_SUMMARY_NOT_IN_DRAFT, CONFIG_SUMMARY_TO_DRAFT,
        FILTER_FIELDS_NOT_IN_SUMMARY,
    )
    mapped = {f.summary for f in CONFIG_SUMMARY_TO_DRAFT}
    assert not (mapped & CONFIG_SUMMARY_NOT_IN_DRAFT)
    assert not (mapped & FILTER_FIELDS_NOT_IN_SUMMARY)


def test_every_mapped_draft_key_is_one_apply_draft_consumes(dialog):
    """A draft key the dialog does not read is a drop with extra steps."""
    from flexprobe_gui._nq_build import CONFIG_SUMMARY_TO_DRAFT
    known = set(dialog.to_draft()) | {"multi_barcodes"}  # legacy-only key
    unknown = {f.draft for f in CONFIG_SUMMARY_TO_DRAFT} - known
    assert not unknown, f"draft keys nothing restores: {sorted(unknown)}"


def test_combo_texts_match_the_dialog(dialog):
    """apply_draft restores these combos by EXACT text match.

    The constants are mirrored from new_query_dialog, so a reworded item there
    must fail here rather than silently restoring nothing.
    """
    from flexprobe_gui import _nq_build as B

    def _items(combo):
        return [combo.itemText(i) for i in range(combo.count())]

    assert _items(dialog.run_mode) == [B.RUN_MODE_MANUAL, B.RUN_MODE_AUTO]
    assert _items(dialog.strain_mode_combo) == [B.STRAIN_MODE_SINGLE,
                                                B.STRAIN_MODE_HOMOLOGY]
    assert _items(dialog.strategy_combo) == [B.STRATEGY_CONSERVED,
                                             B.STRATEGY_UNIQUE,
                                             B.STRATEGY_CONSTITUTIVE]
    # ...and the code the dialog stores for each of those texts round-trips.
    for code, text in B.STRATEGY_TEXT_BY_CODE.items():
        assert B.strategy_code(text) == code


# --- round trip: a recorded non-default must come back unchanged ------------

# One non-default, in-range value per mapped config_summary key. Non-default is
# the point: a map that merely echoed the defaults would pass a test written
# with default values while still dropping every key.
_NON_DEFAULT = {
    "tile_size": 45,
    "tile_step": 3,
    "gc_min": 0.30,
    "gc_max": 0.80,
    "enforce_gc_range": True,
    "tm_min": 61.5,
    "max_n_per_half": 3,
    "max_homopolymer_run": 6,
    "require_tn_junction": False,
    "tn_junction_mode": "off",
    "screen_workers": 6,
    "screen_engine": "blast",
    "skip_screen": True,
    "screen_top_per_gene": 7,
    "primer3": False,
    "group_by_gene": False,          # True is the default
    "min_isoform_coverage": 0.35,
    "multiplex": False,
    "probe_barcode": "BC003",
    "probe_barcodes": ["BC001", "BC002", "BC003"],
    "auto_select_per_gene": 4,
    "multi_strain": True,
    "consensus_strategy": "strain_unique",
    "conserved_max_degenerate_per_half": 2,
    "strain_unique_min_mismatches": 8,
    "strain_unique_per_half": False,
    "consensus_reference": 2,
    "workers": 5,
    "isoform_min_count": 3,
    "multi_count": 3,
    "run_mode": "Auto-pick top N per gene (whole-transcriptome)",
    "ranking_weights": {"gc_proximity": 0.9, "specificity_margin": 0.05},
    "extra_on_target_ids": ["NM_000001.1", "NM_000002.2"],
    "off_target_rule": "either_half",
    "off_target_counting": "local",
    "off_target_min_mismatch": 3,
    "off_target_junction_duplex_bp": 12,
    "exclude_alt_scaffolds": False,
    "reference_strand": "plus",
    "tm_max": 68.0,
    "min_probe_gap_bp": 250,
    "relaxed_probe_gap_bp": 60,
    "floor_probe_gap_bp": 25,
    "split_header_on_pipe": True,
    "null_model_shuffles": 20,
    "null_model_max_windows": 150,
    "null_model_seed": 7,
    "strain_unique_junction_duplex_bp": 11,
}


def test_the_non_default_table_covers_every_mapped_key():
    from flexprobe_gui._nq_build import CONFIG_SUMMARY_TO_DRAFT
    mapped = {f.summary for f in CONFIG_SUMMARY_TO_DRAFT}
    assert mapped == set(_NON_DEFAULT), mapped ^ set(_NON_DEFAULT)


@pytest.mark.parametrize("key", sorted(_NON_DEFAULT))
def test_a_recorded_non_default_changes_the_draft(default_summary, key):
    """Every recorded key must MOVE something in the rebuilt draft.

    This is the finding stated as a test: eight keys were recorded and then had
    no effect at all on the draft, so 'Re-run' silently rebuilt a different
    design. Written generically so a ninth cannot be added quietly.
    """
    baseline = _record_to_draft(default_summary)
    changed = _record_to_draft({**default_summary, key: _NON_DEFAULT[key]})
    assert changed != baseline, (
        f"config_summary[{key!r}] = {_NON_DEFAULT[key]!r} leaves the re-run "
        f"draft untouched, so re-running rebuilds a different design."
    )


def test_the_whole_non_default_summary_reaches_the_form(qapp, default_summary):
    """End to end: recorded run -> draft -> a fresh dialog's widgets.

    The eight keys the old map dropped are asserted by name, because 'it
    reaches the draft' is not the claim a user cares about — 'the form I am
    about to re-run shows the run I clicked' is.
    """
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    from flexprobe_gui import _nq_build as B

    summary = {**default_summary, **_NON_DEFAULT}
    dlg = NewQueryDialog()
    dlg.apply_draft(_record_to_draft(summary, modality="Flex v1 (Chromium)"))

    # --- the eight that were dropped ---
    assert dlg.enforce_gc.isChecked() is True
    assert dlg.tm_min.value() == pytest.approx(61.5)
    assert dlg.screen_engine.currentData() == "blast"
    assert dlg.screen_top_per_gene.value() == 7
    assert dlg.strain_mode_combo.currentText() == B.STRAIN_MODE_HOMOLOGY
    assert dlg.strategy_combo.currentText() == B.STRATEGY_UNIQUE
    assert dlg.unique_min_mm.value() == 8
    assert dlg.unique_per_half.isChecked() is False
    # --- and the ninth: auto-pick came back as a manual-review run ---
    assert dlg.run_mode.currentText() == B.RUN_MODE_AUTO
    assert dlg.auto_n.value() == 4
    # --- the seven config_summary did not record at all ---
    #     (consensus_reference needs a populated template combo; see below)
    assert dlg.workers.value() == 5
    assert dlg.isoform_min.value() == 3
    assert dlg.multi_count.value() == 3
    assert dlg.ranking_box.isChecked() is True
    assert dlg._ranking_spins["gc_proximity"].value() == pytest.approx(0.9)
    assert dlg._ranking_spins["specificity_margin"].value() == pytest.approx(0.05)
    assert dlg.extra_on_target.toPlainText() == "NM_000001.1\nNM_000002.2"
    # --- what already worked must keep working ---
    # The tile box is locked to the chemistry's 50 bp. A recorded 45 would
    # produce no probes (every non-50 tile is dropped at assembly), so the
    # form shows the only size that works.
    assert dlg.tile_size.value() == 50
    assert dlg.tile_step.value() == 3
    assert dlg.gc_min.value() == pytest.approx(0.30)
    assert dlg.gc_max.value() == pytest.approx(0.80)
    assert dlg.max_n.value() == 3
    assert dlg.max_hp.value() == 6
    assert dlg.screen_workers.value() == 6
    assert dlg.skip_screen.isChecked() is True
    assert dlg.primer3.isChecked() is False
    assert dlg.group_by_gene.isChecked() is False     # non-default
    assert dlg.min_isoform_coverage.value() == pytest.approx(0.35)
    assert dlg.conserved_degen.value() == 2
    assert dlg.mode.currentText() == "Singleplex"
    assert dlg._tn_mode() == "off"


# --- the None-vs-falsy rules, one test each ---------------------------------

def test_a_deliberately_recorded_zero_survives(default_summary):
    """``or`` would promote every one of these to the engine default.

    conserved_max_degenerate_per_half=0 is exact-identity conserved windows,
    which on the OROV L segment is the difference between 44 orderable probes
    and none; max_n_per_half=0 rejects every degenerate consensus window.
    """
    draft = _record_to_draft({
        **default_summary,
        "conserved_max_degenerate_per_half": 0,
        "max_n_per_half": 0,
        "max_homopolymer_run": 0,
        "min_isoform_coverage": 0.0,
        "screen_top_per_gene": 0,
    })
    assert draft["conserved_degen"] == 0
    assert draft["max_n"] == 0
    assert draft["max_hp"] == 0
    assert draft["min_isoform_coverage"] == 0.0
    assert draft["screen_top_per_gene"] == 0


def test_a_recorded_false_survives(default_summary):
    draft = _record_to_draft({
        **default_summary,
        "strain_unique_per_half": False,
        "enforce_gc_range": False,
        "require_tn_junction": False,
        "multiplex": False,
    })
    assert draft["unique_per_half"] is False
    assert draft["enforce_gc"] is False
    assert draft["tn_junction"] is False
    assert draft["mode"] == "Singleplex"


def test_a_recorded_none_means_not_applicable_and_takes_the_engine_default(
        default_summary):
    """A non-conserved run records None, which is not 'the user chose 0'."""
    from flexprobe.models import RunConfig
    engine = RunConfig(fasta_path="", output_dir="")
    draft = _record_to_draft({
        **default_summary,
        "conserved_max_degenerate_per_half": None,
        "strain_unique_min_mismatches": None,
        "strain_unique_per_half": None,
    })
    assert draft["conserved_degen"] == engine.conserved_max_degenerate_per_half == 1
    assert draft["unique_min_mm"] == engine.strain_unique_min_mismatches
    assert draft["unique_per_half"] is engine.strain_unique_per_half


def test_a_recorded_tm_min_of_none_means_the_gate_was_off(default_summary):
    """The one None that is a user choice, not a 'not applicable'.

    build_runs writes tm_min=None when the spin box reads 0 (the documented
    'off'). Restoring the engine's 56.0 would silently re-arm a gate the user
    turned off, the mirror image of the None-means-not-applicable case above.
    """
    assert _record_to_draft({**default_summary, "tm_min": None})["tm_min"] == 0.0


def test_an_absent_key_takes_the_engine_default_not_a_literal():
    """A record predating a key gets the engine's value, not a stale copy.

    max_n_per_half named explicitly: this path defaulted it to 0 while the
    engine defaults to 1, so re-running an old record changed what the screen
    threw away.
    """
    draft = _record_to_draft({})
    engine = FilterThresholds()
    assert draft["max_n"] == engine.max_n_per_half == 1
    assert draft["max_hp"] == engine.max_homopolymer_run
    assert draft["gc_min"] == pytest.approx(engine.gc_min)
    assert draft["gc_max"] == pytest.approx(engine.gc_max)
    assert draft["tm_min"] == pytest.approx(engine.tm_min)
    assert draft["enforce_gc"] is engine.enforce_gc_range
    assert draft["tn_junction"] is engine.require_tn_junction


def test_a_record_predating_a_feature_omits_it_so_the_draft_default_applies():
    """The engine default would be WRONG for these two.

    A record with no primer3 key describes a run made before structure scoring
    existed, and one with no screen_engine key ran under BLAST. The engine now
    defaults both the other way, so the draft must carry no key at all and let
    apply_draft's own legacy reasoning stand. Otherwise restoring the record
    would silently change the run it describes.
    """
    draft = _record_to_draft({})
    assert "primer3" not in draft
    assert "screen_engine" not in draft

    from flexprobe_gui.new_query_dialog import NewQueryDialog
    dlg = NewQueryDialog()
    dlg.apply_draft(draft)
    assert dlg.primer3.isChecked() is False
    assert dlg.screen_engine.currentData() == "blast"


# --- the seven settings config_summary did not record at all ---------------

def test_an_absent_new_key_takes_the_engine_default_too():
    """The seven added later get the same treatment: engine value, no literal."""
    from flexprobe.models import RunConfig
    engine = RunConfig(fasta_path="", output_dir="")
    draft = _record_to_draft({})
    assert draft["workers"] == engine.workers
    assert draft["isoform_min"] == engine.isoform_min_count
    assert draft["consensus_template"] == engine.consensus_reference
    assert draft["extra_on_target"] == "\n".join(engine.extra_on_target_ids) == ""
    assert draft["ranking_enabled"] is (engine.ranking_weights is not None) is False


def test_ranking_off_records_none_and_restores_as_off(default_summary):
    """None is 'no custom ranking', which is not an empty weight set.

    Flattening it to {} on the way out would read back as 'ranking on, every
    weight zero' — a different design, and one the form cannot express.
    """
    draft = _record_to_draft({**default_summary, "ranking_weights": None})
    assert draft["ranking_enabled"] is False
    # No weights to restore, so the dialog keeps its own spin-box values.
    assert "ranking_weights" not in draft


def test_ranking_weights_recorded_are_json_and_survive(default_summary):
    """config_summary is persisted as JSON in the run history."""
    import json
    weights = {"gc_proximity": 0.9, "specificity_margin": 0.05}
    summary = {**default_summary, "ranking_weights": weights}
    assert json.loads(json.dumps(summary))["ranking_weights"] == weights
    draft = _record_to_draft(summary)
    assert draft["ranking_enabled"] is True
    assert draft["ranking_weights"] == weights


def test_the_whole_recorded_summary_is_json_serializable(dialog):
    """Every value build_runs records has to survive the run history round trip."""
    import json
    rec, _cfg, _cached = dialog.build_runs()[0]
    assert json.loads(json.dumps(rec.config_summary)) == rec.config_summary


def test_the_record_does_not_alias_the_config(dialog):
    """A recorded list/dict must not share an object with the live RunConfig."""
    rec, cfg, _cached = dialog.build_runs()[0]
    assert rec.config_summary["extra_on_target_ids"] is not cfg.extra_on_target_ids
    assert rec.config_summary["probe_barcodes"] is not cfg.probe_barcodes


# --- legacy records must keep their existing fallbacks ---------------------

def test_a_record_without_run_mode_still_infers_it(default_summary):
    """The explicit key wins, but a record predating it keeps the inference."""
    from flexprobe_gui import _nq_build as B
    legacy = {k: v for k, v in default_summary.items() if k != "run_mode"}
    assert _record_to_draft({**legacy, "auto_select_per_gene": 4})["run_mode"] \
        == B.RUN_MODE_AUTO
    assert _record_to_draft({**legacy, "auto_select_per_gene": 0})["run_mode"] \
        == B.RUN_MODE_MANUAL
    # ...and an explicitly recorded mode overrides what the count implies.
    assert _record_to_draft({**default_summary, "run_mode": B.RUN_MODE_AUTO,
                             "auto_select_per_gene": 0})["run_mode"] \
        == B.RUN_MODE_AUTO


def test_a_record_without_multi_count_still_falls_back_to_the_barcodes(
        qapp, default_summary):
    """apply_draft counts multi_barcodes only while no multi_count is present.

    Writing a default multi_count for a legacy record would disable that
    fallback and silently change its sample count.
    """
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    legacy = {k: v for k, v in default_summary.items() if k != "multi_count"}
    draft = _record_to_draft(
        {**legacy, "multiplex": True,
         "probe_barcodes": ["BC001", "BC002", "BC003"]},
        modality="Flex v1 (Chromium)")
    assert "multi_count" not in draft
    dlg = NewQueryDialog()
    dlg.apply_draft(draft)
    assert dlg.multi_count.value() == 3


# --- consensus template: recorded, and restored into the form ---------------

def test_the_consensus_template_this_run_used_reaches_the_draft(default_summary):
    """The RESOLVED per-run index, not the picker's "every sequence" sentinel.

    A fan-out turns one form into N runs; this record is one of them, so
    re-running it must repeat that run rather than re-queue the whole batch.
    """
    draft = _record_to_draft({**default_summary, "consensus_reference": 2})
    assert draft["consensus_template"] == 2


def test_the_recorded_consensus_template_reaches_the_form(qapp, tmp_path,
                                                          default_summary):
    """Was a strict xfail. apply_draft selected the template at step 4, before
    step 6 sets the FASTA that fills the combo, so the choice found an empty
    combo and the refill defaulted to the ALL_TEMPLATES fan-out. The choice is
    now deferred until the FASTA fills the picker (_select_after_fill)."""
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    from flexprobe_gui import _nq_build as B
    from conftest import SYNTH_200, write_fasta

    fa = write_fasta(tmp_path / "g.fa", [("S1", SYNTH_200),
                                         ("S2", SYNTH_200[::-1]),
                                         ("S3", SYNTH_200)])
    draft = _record_to_draft(
        {**default_summary, "multi_strain": True,
         "consensus_strategy": "conserved", "consensus_reference": 2},
        fasta=fa)
    assert draft["consensus_template"] == 2
    dlg = NewQueryDialog()
    dlg.apply_draft(draft)
    assert dlg.strain_mode_combo.currentText() == B.STRAIN_MODE_HOMOLOGY
    assert dlg.template_combo.currentData() == 2
