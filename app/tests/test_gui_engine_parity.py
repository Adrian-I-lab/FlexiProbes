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

"""The untouched New-query form must build the config an unflagged CLI run does.

If the dialog opens on one value while the engine defaults to another, the
same untouched form produces a different design in the app than on the
command line. ``conserved_max_degenerate_per_half`` and ``max_n_per_half``
are two examples: a form opening on 0 while the engine defaults to 1.

Such drift is silent: nothing in the output says which value was used, and on
real strain alignments the difference is between orderable probes and none. This
module pins the whole comparison rather than those two fields, so any other
mismatch is caught by a failing test instead of by a user.
"""
from __future__ import annotations

import dataclasses

import pytest

pytest.importorskip("PySide6")

from flexprobe.chemistry import get_chemistry, resolve_tn_mode
from flexprobe.constants import TN_MODE_OFF, TN_MODE_SOFT, TN_MODE_STRICT_A
from flexprobe.models import FilterThresholds, RunConfig


@pytest.fixture()
def dialog(qapp, tmp_path):
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    dlg = NewQueryDialog()
    dlg.fasta_edit.setText(str(tmp_path / "in.fa"))
    dlg.out_edit.setText(str(tmp_path / "out"))
    return dlg


def _built_filters(dlg) -> FilterThresholds:
    _rec, cfg, _cached = dlg.build_runs()[0]
    return cfg.filters


# --- the parity that keeps breaking ----------------------------------------

# Fields the dialog is ALLOWED to differ from FilterThresholds on, each with the
# reason. Anything not listed here must match, so a new engine default that the
# dialog does not follow fails this test rather than shipping.
_DIVERGENCE_ALLOWED = {
    # The dialog always states the gate explicitly, so it never emits "auto".
    # strict_a is what auto resolves to under the default require_tn_junction,
    # so the effective gate is identical.
    "tn_junction_mode",
}


def test_untouched_form_matches_engine_filter_defaults(dialog):
    built = _built_filters(dialog)
    engine = FilterThresholds()
    mismatched = {
        f.name: (getattr(built, f.name), getattr(engine, f.name))
        for f in dataclasses.fields(FilterThresholds)
        if f.name not in _DIVERGENCE_ALLOWED
        and getattr(built, f.name) != getattr(engine, f.name)
    }
    assert not mismatched, (
        "The New query dialog opens on values the engine does not default to, "
        "so the same untouched form designs differently in the app and on the "
        f"CLI. field: (gui, engine) = {mismatched}"
    )


def test_untouched_tn_gate_resolves_to_the_engine_default(dialog):
    """The one allowed divergence must still resolve to the same gate."""
    _rec, cfg, _cached = dialog.build_runs()[0]
    chem = get_chemistry(cfg.modality)
    assert resolve_tn_mode(cfg.filters, chem) == \
        resolve_tn_mode(FilterThresholds(), chem) == TN_MODE_STRICT_A


def test_max_n_per_half_opens_at_the_engine_default(dialog):
    """Named explicitly: it opened at 0 for three modalities.

    At 0, homology mode discards every consensus window carrying a two-way
    IUPAC code — which is most of them on a real strain alignment — while the
    same run on the CLI keeps them.
    """
    assert dialog.max_n.value() == FilterThresholds().max_n_per_half == 1


def test_every_modality_preset_matches_the_engine_default(dialog):
    from flexprobe_gui.modality import MODALITY_PRESETS
    engine = FilterThresholds()
    for label in MODALITY_PRESETS:
        dialog.modality.setCurrentText(label)
        assert dialog.max_n.value() == engine.max_n_per_half, label
        assert dialog.max_hp.value() == engine.max_homopolymer_run, label
        assert dialog.gc_min.value() == pytest.approx(engine.gc_min), label
        assert dialog.gc_max.value() == pytest.approx(engine.gc_max), label
        assert dialog.tm_min.value() == pytest.approx(engine.tm_min), label


# --- TN junction mode is reachable from the dialog --------------------------

def test_tn_mode_picker_offers_all_three_gates(dialog):
    offered = {dialog.tn_junction_mode.itemData(i)
               for i in range(dialog.tn_junction_mode.count())}
    assert offered == {TN_MODE_STRICT_A, TN_MODE_SOFT, TN_MODE_OFF}


@pytest.mark.parametrize("mode,expect_require", [
    (TN_MODE_STRICT_A, True),
    (TN_MODE_SOFT, True),      # soft IS a gate; only "off" means none
    (TN_MODE_OFF, False),
])
def test_tn_mode_reaches_the_config(dialog, mode, expect_require):
    dialog._set_tn_mode(mode)
    built = _built_filters(dialog)
    assert built.tn_junction_mode == mode
    assert built.require_tn_junction is expect_require


def test_soft_gate_is_flex_only_and_never_builds_for_visium(dialog):
    """resolve_tn_mode raises ConfigError for soft on Visium, so the dialog
    must not be able to produce that config at all."""
    dialog.modality.setCurrentText("Flex v2 (10X GEM-X)")
    dialog._set_tn_mode(TN_MODE_SOFT)
    assert dialog._tn_mode() == TN_MODE_SOFT
    dialog.modality.setCurrentText("Visium spatial")
    assert dialog._tn_mode() == TN_MODE_STRICT_A
    _rec, cfg, _cached = dialog.build_runs()[0]
    # Would raise if the soft mode had survived the modality switch.
    assert resolve_tn_mode(cfg.filters, get_chemistry(cfg.modality)) \
        == TN_MODE_STRICT_A


def test_soft_draft_cannot_smuggle_itself_into_visium(dialog):
    dialog.apply_draft({"modality": "Visium spatial",
                        "tn_junction_mode": TN_MODE_SOFT})
    assert dialog._tn_mode() == TN_MODE_STRICT_A
    # The note explains the substitution. Asserted on its text, not on
    # isVisibleTo: it lives inside the "Advanced options" section, which opens
    # collapsed, so every child of it reads as not-visible until the user
    # expands it.
    assert "Soft is not available" in dialog.tn_junction_note.text()


# --- draft round-trip, including the pre-mode boolean -----------------------

def test_tn_mode_round_trips_through_a_draft(dialog):
    dialog._set_tn_mode(TN_MODE_SOFT)
    draft = dialog.to_draft()
    assert draft["tn_junction_mode"] == TN_MODE_SOFT
    # Written alongside so a draft saved now still loads in an older build.
    assert draft["tn_junction"] is True

    from flexprobe_gui.new_query_dialog import NewQueryDialog
    other = NewQueryDialog()
    other.apply_draft(draft)
    assert other._tn_mode() == TN_MODE_SOFT


@pytest.mark.parametrize("legacy,expect", [
    (True, TN_MODE_STRICT_A),
    (False, TN_MODE_OFF),
])
def test_pre_mode_draft_maps_the_boolean_not_the_current_default(
        qapp, legacy, expect):
    """A draft saved before the mode existed can only have meant hard or off.

    Restoring it as anything else would change the run the draft claims to
    reproduce.
    """
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    dlg = NewQueryDialog()
    dlg.apply_draft({"tn_junction": legacy})
    assert dlg._tn_mode() == expect


# --- screen parallelism -----------------------------------------------------

def test_screen_workers_reaches_the_config(dialog):
    assert dialog.screen_workers.value() == RunConfig(
        fasta_path="x", output_dir="y").screen_workers == 1
    dialog.screen_workers.setValue(6)
    _rec, cfg, _cached = dialog.build_runs()[0]
    assert cfg.screen_workers == 6
    # ...and is recorded, so a re-run of this run restores it.
    assert _rec.config_summary["screen_workers"] == 6


def test_screen_workers_round_trips_through_a_draft(dialog):
    dialog.screen_workers.setValue(4)
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    other = NewQueryDialog()
    other.apply_draft(dialog.to_draft())
    assert other.screen_workers.value() == 4


# --- engine warning ---------------------------------------------------------

def test_engine_warning_appears_only_off_the_default_engine(dialog):
    from flexprobe.screen import ENGINE_BLAST, ENGINE_KMER
    dialog.screen_engine.setCurrentIndex(
        dialog.screen_engine.findData(ENGINE_KMER))
    assert not dialog.engine_warning.isVisibleTo(dialog)
    dialog.screen_engine.setCurrentIndex(
        dialog.screen_engine.findData(ENGINE_BLAST))
    assert dialog.engine_warning.isVisibleTo(dialog)
    # The ranking consequence is the part a user cannot see in the output.
    assert "rank-comparable" in dialog.engine_warning.text()


# --- HPC export: the dialog and the engine must agree on the defaults -------

def test_hpc_dialog_opens_on_the_engine_chain_defaults(qapp, tmp_path):
    """The dialog's defaults come from flexprobe.hpc_export.

    A separate copy in the dialog could drift below what a
    whole-transcriptome job needs and get it killed. One source means one
    place to change and nothing to keep in step by hand.
    """
    from flexprobe.hpc_export import (
        DEFAULT_CHAIN_LINKS, DEFAULT_CORES, DEFAULT_MEM_GB, DEFAULT_WALLTIME,
    )
    from flexprobe_gui.hpc_dialog import HpcExportDialog
    from flexprobe_gui.run_history import RunRecord

    class _Settings:
        """Empty store, so the dialog shows defaults rather than a saved run's."""

        def __init__(self):
            self._d = {}

        def get(self, key, default=None):
            return self._d.get(key, default)

        def set(self, key, value):
            self._d[key] = value

    cfg = RunConfig(fasta_path="x.fa", output_dir=str(tmp_path / "out"))
    dlg = HpcExportDialog([(RunRecord(name="r"), cfg, None)], _Settings())
    assert dlg.walltime_edit.text() == DEFAULT_WALLTIME == "24:00:00"
    assert dlg.mem_spin.value() == DEFAULT_MEM_GB == 28
    assert dlg.cores_spin.value() == DEFAULT_CORES == 8
    assert dlg.links_spin.value() == DEFAULT_CHAIN_LINKS == 12


@pytest.mark.parametrize("walltime,hours", [
    ("24:00:00", 24.0),
    ("100:00:00", 100.0),   # >24 h is a walltime, not an error
    ("12:30", 12.5),        # HH:MM is accepted by both schedulers
    ("0:30:00", 0.5),
    ("", None),
    ("later", None),
    ("-1:00:00", None),
])
def test_walltime_parser(qapp, walltime, hours):
    from flexprobe_gui.hpc_dialog import HpcExportDialog
    assert HpcExportDialog._walltime_hours(walltime) == hours


def test_budget_line_multiplies_walltime_by_links(qapp, tmp_path):
    """walltime x links is what decides whether the run finishes, and it is the
    one number neither field shows on its own."""
    from flexprobe_gui.hpc_dialog import HpcExportDialog
    from flexprobe_gui.run_history import RunRecord

    class _Settings:
        def __init__(self):
            self._d = {}

        def get(self, key, default=None):
            return self._d.get(key, default)

        def set(self, key, value):
            self._d[key] = value

    cfg = RunConfig(fasta_path="x.fa", output_dir=str(tmp_path / "out"))
    dlg = HpcExportDialog([(RunRecord(name="r"), cfg, None)], _Settings())
    assert "288 h" in dlg.budget_label.text()          # 12 x 24
    dlg.links_spin.setValue(30)
    assert "720 h" in dlg.budget_label.text()          # 30 x 24
    dlg.walltime_edit.setText("nope")
    assert "HH:MM:SS" in dlg.budget_label.text()       # says so, doesn't crash
