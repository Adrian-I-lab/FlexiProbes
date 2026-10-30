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

"""RunConfig assembly for NewQueryDialog.

Split out of new_query_dialog.py (audit #19): ``BuildMixin`` turns the dialog's
form state into ``(RunRecord, RunConfig, cached_path)`` tuples. It reads the
host dialog's widgets via ``self`` and holds no state; mixed into NewQueryDialog
ahead of QDialog in the MRO."""
from __future__ import annotations

import logging
import os
from dataclasses import replace
from datetime import datetime  # noqa: F401  (used when stamping run records)
from typing import Any, List, NamedTuple, Optional, Tuple  # noqa: F401

from flexprobe.constants import TN_MODE_OFF
from flexprobe.models import (
    ScreenResources, ExportSettings, FilterThresholds, RunConfig,
    SpecificityThresholds,
)
from .modality import MODALITY_PRESETS
from .run_history import RunRecord

logger = logging.getLogger(__name__)

# userData sentinel on the consensus-template combo meaning "design once per
# input record and queue them all", as opposed to an int index meaning "use just
# this record". Not an index, so it can never collide with one. Defined here
# rather than in new_query_dialog because that module imports this one.
ALL_TEMPLATES = "__all__"

#: Short, unambiguous tag per multi-strain strategy, for output directory names.
#: Spelled out rather than abbreviated: "cons" would be read as consensus by
#: half the lab and conserved by the other half, and conserved IS the consensus
#: design, which is precisely the confusion this is here to stop.
_STRATEGY_TAG = {
    "conserved": "conserved-panstrain",
    "strain_unique": "strainunique-discriminating",
}

#: Every tag build_runs appends to a run name after "__". Constitutive has no
#: entry above, so it is appended under its own code and listed here by hand.
_GENERATED_NAME_TAGS = tuple(_STRATEGY_TAG.values()) + ("constitutive",)


def strip_generated_run_name(name: str) -> str:
    """Cut a run name back to what the user typed, or to the FASTA stem.

    build_runs appends ``__<strategy tag>`` and then ``__<record id>``, so
    everything from the EARLIEST generated tag onwards is machine-made. That
    suffix names the FASTA and record it was built from, and it must never be
    carried to another input: a Re-run of a Small-segment run that was pointed
    at the Large FASTA came out named
    "..._Small_2strain__strainunique-discriminating__..._Large_...", because the
    suffixed name was stored in the draft, restored into the Run name box, and
    only the strategy tag was ever stripped from it.
    """
    name = (name or "").strip()
    cuts = [c for c in (name.find(f"__{t}") for t in _GENERATED_NAME_TAGS)
            if c != -1]
    return name[:min(cuts)] if cuts else name


# ---------------------------------------------------------------------------
# config_summary -> draft: the ONE map
# ---------------------------------------------------------------------------
#
# ``build_runs`` writes ``RunRecord.config_summary``; ``HomeWindow._record_to_draft``
# reads it back when the user hits Re-run on a record saved before drafts
# existed. Those two lists were maintained by hand in two files and had drifted:
# eight design-relevant settings were recorded and then silently dropped on the
# way back into the form — screen_engine, tm_min, enforce_gc_range,
# screen_top_per_gene, multi_strain, consensus_strategy,
# strain_unique_min_mismatches, strain_unique_per_half — so "Re-run" rebuilt a
# DIFFERENT design from the one the user was looking at, and said nothing.
#
# One table now defines the correspondence, and the reflective tests in
# tests/test_gui_config_defaults.py fail if a key is added to config_summary,
# or a field to FilterThresholds, without an entry here.
#
# Every default below is READ from the engine dataclasses, never restated, for
# the same reason the modality presets are (see modality.py).

# Combo item texts the draft is written in. apply_draft restores these combos by
# findText, i.e. by EXACT match, so a near-miss silently restores nothing —
# which is the failure mode this map exists to remove. They are mirrored here
# because the dialog builds them in new_query_dialog; the pinning test
# test_combo_texts_match_the_dialog fails the moment the two diverge.
RUN_MODE_MANUAL = "Manual review (default)"
RUN_MODE_AUTO = "Auto-pick top N per gene (whole-transcriptome)"
STRAIN_MODE_SINGLE = "Single sequence per gene (default)"
STRAIN_MODE_HOMOLOGY = (
    "Homology consensus (one FASTA = N related sequences of one fragment)"
)
STRATEGY_CONSERVED = "Conserved — probes hit every sequence"
STRATEGY_UNIQUE = "Sequence-unique — probes distinguish one sequence"
STRATEGY_CONSTITUTIVE = (
    "Constitutive (isoform-collapse) — probes shared by ≥N records"
)

# consensus_strategy code <-> strategy combo text. build_runs derives the code
# from the text; the re-run path needs the inverse, so both live here.
STRATEGY_TEXT_BY_CODE = {
    "conserved": STRATEGY_CONSERVED,
    "strain_unique": STRATEGY_UNIQUE,
    "constitutive": STRATEGY_CONSTITUTIVE,
}


# Off-target specificity rule, as ONE choice. The engine has two booleans
# (off_target_require_both_halves, off_target_same_subject), but same_subject
# has no effect once either half is enough (models.SpecificityThresholds), so
# three of the four combinations are distinct and the fourth is a duplicate.
# These codes are what config_summary records and a draft restores.
OT_RULE_SAME_SUBJECT = "same_subject"      # both halves, one off-target gene
OT_RULE_CROSS_SUBJECT = "cross_subject"    # both halves, pooled across genes
OT_RULE_EITHER_HALF = "either_half"        # one disqualifying half rejects

# Mismatch counting: the 10x rule over the full 25 nt half, or the legacy
# >= 80% / >= 20 bp local alignment (off_target_full_half False).
OT_COUNT_FULL_HALF = "full_half"
OT_COUNT_LOCAL = "local"


def ot_rule_code(require_both_halves: bool, same_subject: bool) -> str:
    """The two engine booleans -> the one rule code the form offers."""
    if not require_both_halves:
        return OT_RULE_EITHER_HALF
    return OT_RULE_SAME_SUBJECT if same_subject else OT_RULE_CROSS_SUBJECT


def specificity_summary(specificity: Any, screen: Any) -> dict:
    """The config_summary keys for the specificity controls.

    ONE writer, used by build_runs and by run_import.summary_from_config, so a
    run made here and a run imported from the cluster record these settings
    identically and a Re-run restores both the same way.
    """
    return {
        "off_target_rule": ot_rule_code(
            specificity.off_target_require_both_halves,
            specificity.off_target_same_subject),
        "off_target_counting": (OT_COUNT_FULL_HALF
                                if specificity.off_target_full_half
                                else OT_COUNT_LOCAL),
        "off_target_min_mismatch": int(specificity.off_target_min_mismatch),
        # What actually ran: 0 under BLAST, whatever the box said.
        "off_target_junction_duplex_bp":
            int(specificity.off_target_junction_duplex_bp),
        "exclude_alt_scaffolds": bool(screen.exclude_alt_scaffolds),
        "reference_strand": screen.reference_strand,
    }


def advanced_summary(cfg: Any) -> dict:
    """config_summary keys for the settings added on 2026-09-22 (spacing, Tm
    ceiling, header parsing, null model, strain-unique duplex). ONE writer,
    used by build_runs and by run_import.summary_from_config, like
    specificity_summary."""
    return {
        "tm_max": cfg.filters.tm_max,
        "min_probe_gap_bp": int(cfg.min_probe_gap_bp),
        "relaxed_probe_gap_bp": int(cfg.relaxed_probe_gap_bp),
        "floor_probe_gap_bp": int(cfg.floor_probe_gap_bp),
        "split_header_on_pipe": bool(cfg.split_header_on_pipe),
        "null_model_shuffles": int(cfg.null_model_shuffles),
        "null_model_max_windows": int(cfg.null_model_max_windows),
        "null_model_seed": int(cfg.null_model_seed),
        "strain_unique_junction_duplex_bp":
            int(cfg.strain_unique_junction_duplex_bp),
    }


def strategy_code(text: str) -> str:
    """Strategy combo text -> the code stored in RunConfig/config_summary."""
    if text.startswith("Sequence-unique"):
        return "strain_unique"
    if text.startswith("Constitutive"):
        return "constitutive"
    return "conserved"


# Returned by a converter to mean "write no draft key at all", which hands the
# decision to apply_draft's own default. That is deliberate wherever the draft's
# default already encodes what a record MISSING the key must have meant, so the
# reasoning lives in exactly one place (see primer3 / screen_engine below).
_OMIT: Any = object()

# Engine defaults, read live rather than restated. A run_config's two path
# arguments are required but irrelevant to its defaults.
_F = FilterThresholds()
_R = RunConfig(fasta_path="", output_dir="")
_S = ScreenResources()
_SP = SpecificityThresholds()
_OT_RULE_DEFAULT = ot_rule_code(_SP.off_target_require_both_halves,
                                _SP.off_target_same_subject)
_OT_COUNT_DEFAULT = (OT_COUNT_FULL_HALF if _SP.off_target_full_half
                     else OT_COUNT_LOCAL)


def _int_or(default: Any) -> Any:
    """None (the key was recorded as 'not applicable') -> the engine default.

    Tested against None, never with ``or``: a recorded 0 is a real, deliberate
    setting — max_n_per_half=0, conserved_max_degenerate_per_half=0 — and ``or``
    would silently promote it to the default. On a short or variable segment
    that one substitution can be the difference between usable probes and
    none.
    """
    return lambda v: default if v is None else int(v)


def _float_or(default: Any) -> Any:
    return lambda v: default if v is None else float(v)


def _bool_or(default: Any) -> Any:
    """Same None-not-falsy rule: a recorded False must survive as False."""
    return lambda v: default if v is None else bool(v)


def _tm_min_to_draft(v: Any) -> Any:
    """A RECORDED None means the Tm gate was switched off, not 'unknown'.

    build_runs writes ``tm_min=None`` when the spin box reads 0, which is the
    documented "off" value (a real per-half Tm floor is never 0 C). Mapping that
    None to the engine's 56.0 would silently re-arm a gate the user turned off,
    so it maps back to the spin box's 0. A record that never recorded the key at
    all is a different case and takes the engine default (see the table).
    """
    return 0.0 if v is None else float(v)


def _mode_text(v: Any) -> str:
    return "Singleplex" if v is False else "Multiplex"


def _multi_barcodes(v: Any) -> Any:
    """Only a genuinely multiplexed run (>1 barcode) restores a barcode list.

    apply_draft consults this only for drafts predating the sample-count spin,
    where the list length IS the sample count; one barcode is singleplex and
    carries no count.
    """
    bcs = [b for b in (v or []) if b]
    return bcs if len(bcs) > 1 else []


def _strain_mode_text(v: Any) -> str:
    return STRAIN_MODE_HOMOLOGY if v else STRAIN_MODE_SINGLE


def _strategy_text(v: Any) -> Any:
    # "" is what build_runs records for a non-homology run: no strategy was in
    # play, so claim none and let the dialog keep its own default.
    return STRATEGY_TEXT_BY_CODE.get(v or "", _OMIT)


def _auto_n(v: Any) -> Any:
    # 0 means the run was in manual review, and the spin box's minimum is 1, so
    # there is no honest value to restore — run_mode below carries the fact.
    return int(v) if v else _OMIT


def _lines(v: Any) -> str:
    """A recorded id list -> the free-text box the dialog splits back on."""
    return "\n".join(v or [])


def _ranking_enabled(v: Any) -> bool:
    """None is the engine's "no custom ranking", not "an empty weight set"."""
    return v is not None


def _ranking_weights(v: Any) -> Any:
    # Nothing to restore when ranking was off: apply_draft iterates the dict, so
    # an empty one is a no-op anyway, and omitting says so honestly.
    return dict(v) if v else _OMIT


def _run_mode_text(v: Any) -> str:
    """auto_select_per_gene drives the Run-mode combo as well as the spin box.

    Restoring only the spin box (what this path used to do) put an auto-pick run
    back into the form as a manual-review run with a decorative "probes per
    gene" number: build_runs reads auto_n only when the combo says Auto, so the
    re-run silently became a different workflow.
    """
    return RUN_MODE_AUTO if v else RUN_MODE_MANUAL


class _SummaryField(NamedTuple):
    """One config_summary key and where it lands in a draft."""

    summary: str                 # key in RunRecord.config_summary
    draft: str                   # key apply_draft consumes
    convert: Any                 # recorded value -> draft value (or _OMIT)
    missing: Any                 # used when the key is ABSENT; _OMIT = no key


CONFIG_SUMMARY_TO_DRAFT: tuple = (
    # --- tiling ---
    _SummaryField("tile_size", "tile_size", _int_or(_R.tile_length), _R.tile_length),
    _SummaryField("tile_step", "tile_step", _int_or(_R.tile_step), _R.tile_step),
    # --- filter funnel (defaults straight off FilterThresholds) ---
    _SummaryField("gc_min", "gc_min", _float_or(_F.gc_min), _F.gc_min),
    _SummaryField("gc_max", "gc_max", _float_or(_F.gc_max), _F.gc_max),
    _SummaryField("enforce_gc_range", "enforce_gc",
                  _bool_or(_F.enforce_gc_range), _F.enforce_gc_range),
    _SummaryField("tm_min", "tm_min", _tm_min_to_draft, _F.tm_min),
    _SummaryField("max_n_per_half", "max_n",
                  _int_or(_F.max_n_per_half), _F.max_n_per_half),
    _SummaryField("max_homopolymer_run", "max_hp",
                  _int_or(_F.max_homopolymer_run), _F.max_homopolymer_run),
    _SummaryField("require_tn_junction", "tn_junction",
                  _bool_or(_F.require_tn_junction), _F.require_tn_junction),
    # Absent means the run predates the three-way mode, so it can only have been
    # the hard gate or none at all — never soft. apply_draft derives exactly
    # that from `tn_junction`, so pass a mode through only when one was recorded.
    _SummaryField("tn_junction_mode", "tn_junction_mode",
                  lambda v: v or _OMIT, _OMIT),
    # --- off-target screen ---
    _SummaryField("screen_workers", "screen_workers",
                  _int_or(_R.screen_workers), _R.screen_workers),
    _SummaryField("workers", "workers", _int_or(_R.workers), _R.workers),
    # Absent means the run predates the setting, i.e. it ran under BLAST.
    # apply_draft already defaults a key-less draft to ENGINE_BLAST for that
    # exact reason, so omit rather than restate it here.
    _SummaryField("screen_engine", "screen_engine", lambda v: v or _OMIT, _OMIT),
    _SummaryField("skip_screen", "skip_screen",
                  _bool_or(_S.skip_screen), _S.skip_screen),
    _SummaryField("screen_top_per_gene", "screen_top_per_gene",
                  _int_or(_R.screen_top_per_gene), _R.screen_top_per_gene),
    # --- off-target specificity rules (added 2026-09-22) ---
    # A record without these keys predates the controls, and was therefore
    # built from SpecificityThresholds() and ScreenResources() defaults, so
    # the engine default is what a re-run restores.
    _SummaryField("off_target_rule", "ot_rule",
                  lambda v: v or _OT_RULE_DEFAULT, _OT_RULE_DEFAULT),
    _SummaryField("off_target_counting", "ot_counting",
                  lambda v: v or _OT_COUNT_DEFAULT, _OT_COUNT_DEFAULT),
    _SummaryField("off_target_min_mismatch", "ot_min_mm",
                  _int_or(_SP.off_target_min_mismatch),
                  _SP.off_target_min_mismatch),
    # One recorded number, two form fields: 0 is "off", anything else is on
    # at that many bp. A recorded 0 leaves the bp box at its default.
    _SummaryField("off_target_junction_duplex_bp", "junction_duplex",
                  lambda v: (_SP.off_target_junction_duplex_bp > 0
                             if v is None else int(v) > 0),
                  _SP.off_target_junction_duplex_bp > 0),
    _SummaryField("off_target_junction_duplex_bp", "junction_duplex_bp",
                  lambda v: int(v) if v else _OMIT,
                  _SP.off_target_junction_duplex_bp),
    _SummaryField("exclude_alt_scaffolds", "exclude_alt",
                  _bool_or(_S.exclude_alt_scaffolds), _S.exclude_alt_scaffolds),
    _SummaryField("reference_strand", "ref_strand",
                  lambda v: v or _S.reference_strand, _S.reference_strand),
    # --- the rest of the CLI-only design settings (added 2026-09-22) ---
    # tm_max None is "no ceiling", which the spin box shows as 0.
    _SummaryField("tm_max", "tm_max",
                  lambda v: 0.0 if v is None else float(v),
                  0.0 if _F.tm_max is None else float(_F.tm_max)),
    _SummaryField("min_probe_gap_bp", "min_gap",
                  _int_or(_R.min_probe_gap_bp), _R.min_probe_gap_bp),
    _SummaryField("relaxed_probe_gap_bp", "relaxed_gap",
                  _int_or(_R.relaxed_probe_gap_bp), _R.relaxed_probe_gap_bp),
    _SummaryField("floor_probe_gap_bp", "floor_gap",
                  _int_or(_R.floor_probe_gap_bp), _R.floor_probe_gap_bp),
    _SummaryField("split_header_on_pipe", "split_header",
                  _bool_or(_R.split_header_on_pipe), _R.split_header_on_pipe),
    _SummaryField("null_model_shuffles", "null_shuffles",
                  _int_or(_R.null_model_shuffles), _R.null_model_shuffles),
    _SummaryField("null_model_max_windows", "null_windows",
                  _int_or(_R.null_model_max_windows), _R.null_model_max_windows),
    _SummaryField("null_model_seed", "null_seed",
                  _int_or(_R.null_model_seed), _R.null_model_seed),
    # One number, two fields, as for the off-target duplex: 0 is off.
    _SummaryField("strain_unique_junction_duplex_bp", "unique_jd",
                  lambda v: (_R.strain_unique_junction_duplex_bp > 0
                             if v is None else int(v) > 0),
                  _R.strain_unique_junction_duplex_bp > 0),
    _SummaryField("strain_unique_junction_duplex_bp", "unique_jd_bp",
                  lambda v: int(v) if v else _OMIT,
                  _R.strain_unique_junction_duplex_bp),
    # --- scoring / grouping ---
    # DELIBERATELY omitted when absent, so the draft's own False applies. The
    # engine defaults primer3_scoring to True, but a record with no primer3 key
    # describes a run made BEFORE structure scoring existed, i.e. without it.
    # Restoring True would change the run the re-run claims to repeat (audit #12).
    _SummaryField("primer3", "primer3", _bool_or(False), _OMIT),
    _SummaryField("group_by_gene", "group_by_gene",
                  _bool_or(_R.group_by_gene), _R.group_by_gene),
    _SummaryField("min_isoform_coverage", "min_isoform_coverage",
                  _float_or(_R.min_isoform_coverage), _R.min_isoform_coverage),
    _SummaryField("isoform_min_count", "isoform_min",
                  _int_or(_R.isoform_min_count), _R.isoform_min_count),
    # Which extra records count as on-target. Changes what the screen rejects,
    # so a re-run without it screens against a different definition of "self".
    _SummaryField("extra_on_target_ids", "extra_on_target",
                  _lines, _lines(_R.extra_on_target_ids)),
    # One recorded value, two form fields: the checkbox says custom ranking was
    # on at all, the dict says with what weights.
    _SummaryField("ranking_weights", "ranking_enabled",
                  _ranking_enabled, _ranking_enabled(_R.ranking_weights)),
    _SummaryField("ranking_weights", "ranking_weights",
                  _ranking_weights, _OMIT),
    # --- chemistry ---
    _SummaryField("multiplex", "mode", _mode_text, "Multiplex"),
    _SummaryField("probe_barcode", "probe_barcode", lambda v: v or "", ""),
    _SummaryField("probe_barcodes", "multi_barcodes", _multi_barcodes, []),
    # OMITTED when absent, deliberately: apply_draft falls back to counting
    # `multi_barcodes` only while no multi_count is present, and that fallback
    # exists for exactly the legacy records that do not record this key. A
    # default here would write the key for them too and disable it.
    _SummaryField("multi_count", "multi_count", _int_or(_OMIT), _OMIT),
    # --- workflow: one recorded number drives two form fields ---
    _SummaryField("auto_select_per_gene", "auto_n", _auto_n, _OMIT),
    _SummaryField("auto_select_per_gene", "run_mode",
                  _run_mode_text, RUN_MODE_MANUAL),
    # AFTER the inference above, so an explicitly recorded mode wins and a
    # legacy record with no run_mode key still gets the inferred one.
    _SummaryField("run_mode", "run_mode", lambda v: v or _OMIT, _OMIT),
    # --- homology mode ---
    _SummaryField("multi_strain", "strain_mode",
                  _strain_mode_text, STRAIN_MODE_SINGLE),
    _SummaryField("consensus_strategy", "strategy", _strategy_text, _OMIT),
    _SummaryField("consensus_reference", "consensus_template",
                  _int_or(_R.consensus_reference), _R.consensus_reference),
    _SummaryField("conserved_max_degenerate_per_half", "conserved_degen",
                  _int_or(_R.conserved_max_degenerate_per_half),
                  _R.conserved_max_degenerate_per_half),
    _SummaryField("strain_unique_min_mismatches", "unique_min_mm",
                  _int_or(_R.strain_unique_min_mismatches),
                  _R.strain_unique_min_mismatches),
    _SummaryField("strain_unique_per_half", "unique_per_half",
                  _bool_or(_R.strain_unique_per_half),
                  _R.strain_unique_per_half),
)

# Recorded, and deliberately NOT restored. These describe how THIS record was
# batched, not how its probes were designed: the re-run dialog derives its own
# batching from the FASTAs it is handed. Listed rather than merely absent so the
# completeness test can tell "considered and excluded" from "forgotten".
CONFIG_SUMMARY_NOT_IN_DRAFT = frozenset({
    "bulk_batch_index",
    "bulk_batch_size",
})

# FilterThresholds fields with no corresponding form field, so nothing to record
# or restore. Same purpose: the completeness test reads this, so adding a field
# to FilterThresholds fails a test until it is either mapped or listed here.
FILTER_FIELDS_NOT_IN_SUMMARY = frozenset({
    "reject_junction_base_n",        # junction rules are not user-tunable here
    "degenerate_junction_clearance",
    "reject_junction_mismatch",
    "junction_window",
})


def summary_to_draft(summary: dict) -> dict:
    """Turn a ``RunRecord.config_summary`` back into a draft for apply_draft.

    Covers only what the summary records; the caller adds the record-level
    fields (modality, project, run name, paths) that live on the RunRecord
    itself rather than in the summary.
    """
    draft: dict = {}
    for field in CONFIG_SUMMARY_TO_DRAFT:
        if field.summary in summary:
            value = field.convert(summary[field.summary])
        else:
            value = field.missing
        if value is not _OMIT:
            draft[field.draft] = value
    return draft


from typing import TYPE_CHECKING

if TYPE_CHECKING:  # give the mixin the host widget's type for mypy
    from PySide6.QtWidgets import QDialog as _HostBase
else:
    _HostBase = object


def _fasta_record_ids(path: str) -> list:
    """Record ids of ``path``, or [] if it cannot be read.

    Silent on failure by design: the caller falls back to a single run, and the
    user gets a precise error from the pipeline rather than a dialog-time one.
    """
    try:
        from flexprobe.errors import FlexprobeError
        from flexprobe.io_fasta import parse_fasta
        return [gid for gid, _desc, _seq in parse_fasta((path or "").strip())]
    except (OSError, ValueError, FlexprobeError):
        return []


def _safe_name(sid: str) -> str:
    """Record id -> a path-safe run-name fragment.

    Ids routinely carry pipes and dots (``OROV_BeAn19991|Large|KP052850.1``),
    which must not reach a directory name.
    """
    return "".join(c if c.isalnum() or c in "_-" else "_" for c in sid)[:40]


class BuildMixin(_HostBase):
    def build_runs(self) -> list:
        """Return a list of ``(record, config, cached_path)`` tuples.

        In single mode the list has one entry. In bulk mode every dropped /
        picked FASTA gets its own (record, config) sharing the dialog's
        chemistry / filters / BLAST settings. The output directory for
        each run is auto-derived from the FASTA filename so the bulk runs
        don't collide.
        """
        out_base = self.out_edit.text().strip()
        # Common config — everything that doesn't vary per FASTA.
        modality_label = self.modality.currentText()
        modality = MODALITY_PRESETS[modality_label]
        _tm_min = float(self.tm_min.value())
        _tn_mode = self._tn_mode()
        filters = FilterThresholds(
            gc_min=float(self.gc_min.value()),
            gc_max=float(self.gc_max.value()),
            enforce_gc_range=self.enforce_gc.isChecked(),
            # The spinbox has no null state; 0 is the documented "off" value
            # (a real per-half Tm floor is never 0 C).
            tm_min=(_tm_min if _tm_min > 0 else None),
            # Same "0 = off" convention as tm_min.
            tm_max=(float(self.tm_max.value())
                    if self.tm_max.value() > 0 else None),
            max_n_per_half=int(self.max_n.value()),
            max_homopolymer_run=int(self.max_hp.value()),
            # Both are recorded. resolve_tn_mode lets the explicit mode win,
            # and require_tn_junction stays the honest "is a junction gate
            # engaged at all" flag that config_summary, older drafts and
            # every committed run_config.json are written in terms of.
            require_tn_junction=(_tn_mode != TN_MODE_OFF),
            tn_junction_mode=_tn_mode,
        )
        db_paths = [
            line.strip() for line in self.db_edit.toPlainText().splitlines()
            if line.strip()
        ]
        cached_path = self.cached_hits.text().strip()
        export = ExportSettings(
            write_idt_opools=True,
            write_cellranger_probe_set=True,
            write_results_table=True,
            write_manifest=True,
            overwrite=self.overwrite.isChecked(),
        )
        multiplex = self.mode.currentText() == "Multiplex"
        # Flex v1 multiplex -> Probe Barcodes BC001..BC00N for the chosen sample
        # count; otherwise fall back to the single combo.
        multi_labels = self._selected_multi_barcodes()
        if multi_labels:
            probe_barcodes = list(multi_labels)
        else:
            single = self.barcode_combo.currentData() or ""
            probe_barcodes = [single] if single else []
        probe_bc = probe_barcodes[0] if probe_barcodes else ""
        auto_mode = self.run_mode.currentText().startswith("Auto")
        auto_n = int(self.auto_n.value()) if auto_mode else 0
        workers = int(self.workers.value())
        screen_workers = int(self.screen_workers.value())
        extra_on_target_ids = [
            line.strip() for line in self.extra_on_target.toPlainText().splitlines()
            if line.strip()
        ]
        ranking_weights = (
            {k: float(sp.value()) for k, sp in self._ranking_spins.items()}
            if self.ranking_box.isChecked() else None
        )
        project_id = self.project_combo.currentData() or ""
        project_safe = ""
        if project_id:
            proj = self._projects.get(project_id)
            if proj is not None:
                project_safe = "".join(
                    c if c.isalnum() or c in "_-." else "_" for c in proj.name
                )

        # Build the per-FASTA list.
        if self._is_bulk:
            fastas = list(self._bulk_paths)
        else:
            fastas = [self.fasta_edit.text().strip()]

        # Strategy is fasta-independent, so resolve it once — the template
        # fan-out below needs it before the per-run loop starts.
        multi_strain = self.strain_mode_combo.currentText().startswith("Homology")
        # Via the shared helper, so the code written into config_summary and the
        # combo text a re-run restores from it can never disagree.
        consensus_strategy = strategy_code(self.strategy_combo.currentText())

        # Expand each FASTA into one job per consensus template. "Every
        # sequence" is the only way to obtain the complete pan-strain set: a
        # single template emits only its own flavour of each variant window.
        # Record ids are read per FASTA, so a bulk batch of files with different
        # record counts fans out correctly instead of assuming the first file's.
        jobs: List[Tuple[str, int, Optional[str]]] = []
        # ONE selector drives the fan-out, whichever multi-strain mode is on.
        # The two modes pick their sequence with DIFFERENT parameters --
        # consensus_reference for conserved (a coordinate frame) and
        # strain_unique_target for strain_unique (the strain being
        # discriminated FOR) -- but the fan-out question is identical: "this
        # record, or every record in turn?". Until 2026-09-11 only the
        # conserved path could answer "every", so a strain_unique matrix had to
        # be launched by hand once per strain, and on 2026-09-11 that produced
        # two byte-identical runs at target 0 instead of one per strain.
        _sel = (
            self.unique_target_combo.currentData()
            if multi_strain and consensus_strategy == "strain_unique"
            else self.template_combo.currentData()
            if multi_strain
            else 0
        )
        _tpl = _sel
        for fasta in fastas:
            if _tpl == ALL_TEMPLATES:
                ids = _fasta_record_ids(fasta)
                if len(ids) > 1:
                    jobs.extend((fasta, i, sid) for i, sid in enumerate(ids))
                    continue
            jobs.append((fasta, int(_tpl or 0) if _tpl != ALL_TEMPLATES else 0, None))
        # More than one job per FASTA means the batch behaves like a bulk run
        # (unique output dirs, per-run names) even if a single file was chosen.
        _fanned = len(jobs) > len(fastas)

        # Base draft shared by every run in this batch; per-FASTA fields
        # (fasta_path, run_name) are overridden below so each run re-opens
        # pre-filled with its own input.
        base_draft = self.to_draft()

        stamp_base = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_runs: list = []
        for i, (fasta, template_idx, template_id) in enumerate(jobs):
            name = (
                self.run_name.text().strip()
                if not self._is_bulk and not _fanned and self.run_name.text().strip()
                else os.path.splitext(os.path.basename(fasta))[0]
            )
            # NAME THE STRATEGY, NOT JUST THE TEMPLATE. Until 2026-09-11 the
            # folder name carried a `__<record>` suffix for a CONSERVED fan-out
            # and nothing at all for STRAIN_UNIQUE, so the presence of a suffix
            # read as "this is the per-strain, discriminating one" when it meant
            # the exact opposite. A real run set on 2026-09-11 was labelled
            # backwards because of it, and a deck nearly went out that way.
            #
            # The two modes also differ in which parameter picks the sequence:
            # `consensus_reference` for conserved (a coordinate frame, and only
            # the template combo fans out over it), `strain_unique_target` for
            # strain_unique. Naming the target for strain_unique too means two
            # runs at different targets get DISTINCT directories instead of
            # relying on the timestamp to keep them apart.
            # STRIP ANY SUFFIX A PREVIOUS RUN LEFT BEHIND BEFORE ADDING ONE.
            # The dialog keeps the Run name field between runs, so appending
            # compounds: on 2026-09-11 a second run came out named
            # "..._strainunique-discriminating__PEIAM..._strainunique-discriminating__BeAn...",
            # carrying a STALE strain name ahead of the real one, which is
            # exactly the misreading the tag was added to prevent.
            name = strip_generated_run_name(name)
            if multi_strain:
                name = f"{name}__{_STRATEGY_TAG.get(consensus_strategy, consensus_strategy)}"
            if template_id is not None:
                # Name the run after its template so the queue is readable and
                # the output dirs of a fan-out cannot collide.
                name = f"{name}__{_safe_name(template_id)}"
            elif multi_strain and consensus_strategy == "strain_unique":
                _ids = _fasta_record_ids(fasta)
                _ti = int(template_idx)   # already resolved by the fan-out
                if 0 <= _ti < len(_ids):
                    name = f"{name}__{_safe_name(_ids[_ti])}"
            # In bulk mode add a sub-second nonce so two FASTAs created in
            # the same wall-clock second still get unique output dirs.
            stamp = (f"{stamp_base}_{i:03d}"
                     if (self._is_bulk or _fanned) else stamp_base)
            sub = f"{stamp}_{name}"
            out_dir = (
                os.path.join(out_base, project_safe, sub) if project_safe
                else os.path.join(out_base, sub)
            )
            screen = ScreenResources(
                db_paths=list(db_paths),
                # Recorded explicitly, never left to the dataclass default: a
                # config file is a reproducibility artefact, and one that omits
                # the engine silently adopts whatever the default is when it is
                # next loaded (see config_io._pin_legacy_screen_engine).
                engine=self.screen_engine.currentData(),
                skip_screen=self.skip_screen.isChecked(),
                exclude_alt_scaffolds=self.exclude_alt.isChecked(),
                reference_strand=self.ref_strand.currentData(),
                use_cached_hits=bool(cached_path),
                cached_hits_path=cached_path or None,
            )
            # The junction-duplex gate is measurable only by the k-mer engine:
            # the sites it exists to catch carry 6+ mismatches per half, so
            # blastn reports no alignment there at all and the records the gate
            # needs provably cannot be recovered from a hit list. RunConfig
            # .validate() therefore REFUSES engine=blast while the gate is on.
            #
            # The dialog has an engine picker but no junction-duplex control, so
            # left alone every BLAST run the GUI builds would be rejected with an
            # error the user has no widget to act on. Turn the gate off for the
            # BLAST path explicitly -- never silently: `_blast_gate_note` is
            # surfaced to the user, and 0 is what the CLI error tells you to pass.
            _rule = self.ot_rule.currentData()
            _both = _rule != OT_RULE_EITHER_HALF
            specificity = SpecificityThresholds(
                off_target_require_both_halves=_both,
                # No effect under either-half, so it keeps the engine default
                # there, which is what the CLI writes for the same run.
                off_target_same_subject=(
                    _rule == OT_RULE_SAME_SUBJECT if _both
                    else _SP.off_target_same_subject),
                off_target_full_half=(
                    self.ot_counting.currentData() != OT_COUNT_LOCAL),
                off_target_min_mismatch=int(self.ot_min_mm.value()),
                off_target_junction_duplex_bp=(
                    int(self.junction_duplex_bp.value())
                    if self.junction_duplex.isChecked() else 0),
            )
            blast_gate_note = ""
            if screen.engine == "blast" and (
                    specificity.off_target_junction_duplex_bp > 0):
                blast_gate_note = (
                    "Screening engine is BLAST, so the junction-duplex filter "
                    "is switched off for this run (BLAST cannot measure it). "
                    "The off-target screen is weaker than the default. Use the "
                    "k-mer engine to keep the filter."
                )
                specificity = replace(
                    specificity, off_target_junction_duplex_bp=0)
                logger.warning("%s (run %r)", blast_gate_note, name)
                # Read by home_window after build_runs() returns, so the user is
                # told rather than only the log file.
                self.blast_gate_note = blast_gate_note
            # Gate on the strategy (robust) rather than widget.isVisible().
            # build_runs() runs AFTER exec() returns, i.e. on a closed dialog,
            # and isVisible() is False for every child of a hidden window — so
            # the visibility form silently forced target 0 on every run no
            # matter which sequence the user picked. Same reasoning as
            # NewQueryDialog._selected_multi_barcodes.
            # `template_idx` carries the FANNED selection for whichever mode is
            # active (see the _sel comment above), so a fan-out has already
            # resolved it to a concrete record index. Route it to the parameter
            # that mode actually reads, and never to both.
            if consensus_strategy == "strain_unique":
                unique_target_idx = template_idx
                consensus_reference_idx = 0
            else:
                unique_target_idx = 0
                consensus_reference_idx = template_idx
            cfg = RunConfig(
                fasta_path=fasta,
                output_dir=out_dir,
                modality=modality,
                filters=filters,
                screen=screen,
                specificity=specificity,
                export=export,
                tile_length=int(self.tile_size.value()),
                tile_step=int(self.tile_step.value()),
                multiplex=multiplex,
                probe_barcodes=list(probe_barcodes),
                discovery_only=True,
                auto_select_per_gene=auto_n,
                multi_strain=multi_strain,
                consensus_strategy=consensus_strategy,
                isoform_min_count=int(self.isoform_min.value()),
                # Only the conserved rule has a degenerate tolerance; the
                # other two strategies design from a concrete sequence, so
                # send 0 rather than a stale spin-box value. Resolved here,
                # once, for the whole batch, like every other strategy field.
                conserved_max_degenerate_per_half=(
                    int(self.conserved_degen.value())
                    if multi_strain and consensus_strategy == "conserved" else 0
                ),
                # Routed above, never both: strain_unique IGNORES
                # consensus_reference (models.py says so), so writing the
                # target index into it would record a parameter that had no
                # effect and invite exactly the mode confusion this session
                # spent an afternoon untangling.
                consensus_reference=int(consensus_reference_idx),
                strain_unique_target=int(unique_target_idx or 0),
                strain_unique_min_mismatches=int(self.unique_min_mm.value()),
                strain_unique_per_half=bool(self.unique_per_half.isChecked()),
                primer3_scoring=self.primer3.isChecked(),
                workers=workers,
                screen_workers=screen_workers,
                extra_on_target_ids=list(extra_on_target_ids),
                group_by_gene=self.group_by_gene.isChecked(),
                min_isoform_coverage=float(self.min_isoform_coverage.value()),
                screen_top_per_gene=int(self.screen_top_per_gene.value()),
                ranking_weights=ranking_weights,
                min_probe_gap_bp=int(self.min_gap.value()),
                relaxed_probe_gap_bp=int(self.relaxed_gap.value()),
                floor_probe_gap_bp=int(self.floor_gap.value()),
                split_header_on_pipe=self.split_header.isChecked(),
                null_model_shuffles=int(self.null_shuffles.value()),
                null_model_max_windows=int(self.null_windows.value()),
                null_model_seed=int(self.null_seed.value()),
                strain_unique_junction_duplex_bp=(
                    int(self.unique_jd_bp.value())
                    if self.unique_jd.isChecked() else 0),
            )
            record = RunRecord(
                name=name,
                modality_label=modality_label,
                fasta_path=fasta,
                output_dir=out_dir,
                project_id=project_id,
                config_summary={
                    "tile_size": cfg.tile_length,
                    "tile_step": cfg.tile_step,
                    "gc_min": filters.gc_min,
                    "gc_max": filters.gc_max,
                    "enforce_gc_range": filters.enforce_gc_range,
                    "tm_min": filters.tm_min,
                    "max_n_per_half": filters.max_n_per_half,
                    "max_homopolymer_run": filters.max_homopolymer_run,
                    "require_tn_junction": filters.require_tn_junction,
                    "tn_junction_mode": filters.tn_junction_mode,
                    "screen_workers": cfg.screen_workers,
                    # The filter/tiling worker count, distinct from the screen's
                    # own pool above. Not a design setting on its own, but a run
                    # re-run at 1 worker on a whole transcriptome is a different
                    # afternoon, so it is restored rather than reset.
                    "workers": cfg.workers,
                    "screen_engine": screen.engine,
                    "skip_screen": screen.skip_screen,
                    "primer3": cfg.primer3_scoring,
                    "group_by_gene": cfg.group_by_gene,
                    "min_isoform_coverage": cfg.min_isoform_coverage,
                    "isoform_min_count": cfg.isoform_min_count,
                    "screen_top_per_gene": cfg.screen_top_per_gene,
                    **specificity_summary(specificity, screen),
                    **advanced_summary(cfg),
                    # Copied, not aliased: the record is serialized to the run
                    # history and must not share a mutable object with the
                    # config the pipeline is about to run.
                    "extra_on_target_ids": list(cfg.extra_on_target_ids),
                    # None is meaningful — it is what the engine stores when
                    # custom ranking is off — so it is recorded as None rather
                    # than flattened to {}, which would mean "on, all zero".
                    "ranking_weights": (dict(cfg.ranking_weights)
                                        if cfg.ranking_weights is not None
                                        else None),
                    "multiplex": multiplex,
                    "probe_barcode": probe_bc,
                    "probe_barcodes": list(probe_barcodes),
                    # The sample count itself, not just the barcodes it
                    # produced: it is the field the form actually has, and the
                    # barcode list only implies it for a Flex v1 multiplex run.
                    "multi_count": int(self.multi_count.value()),
                    "auto_select_per_gene": auto_n,
                    # Recorded rather than left to be inferred from
                    # auto_select_per_gene above. The inference is exact today
                    # (the auto_n spin box cannot go below 1, so auto mode
                    # always writes a non-zero count) but it is exact only by
                    # accident of that range; widening the spin box to 0 would
                    # silently turn every auto run into a manual one on re-run.
                    "run_mode": self.run_mode.currentText(),
                    "multi_strain": multi_strain,
                    "consensus_strategy": consensus_strategy if multi_strain else "",
                    # The template this record was actually designed against —
                    # the RESOLVED per-run index, not the picker's selection.
                    # "Every sequence" fans one form out into N runs, and this
                    # record is one of them; restoring the sentinel would re-run
                    # the whole fan-out instead of repeating the run clicked.
                    "consensus_reference": cfg.consensus_reference,
                    "conserved_max_degenerate_per_half": (
                        int(self.conserved_degen.value())
                        if multi_strain and consensus_strategy == "conserved"
                        else None
                    ),
                    "strain_unique_min_mismatches": (
                        int(self.unique_min_mm.value()) if multi_strain else None
                    ),
                    "strain_unique_per_half": (
                        bool(self.unique_per_half.isChecked())
                        if multi_strain else None
                    ),
                    "bulk_batch_index": i if self._is_bulk else None,
                    "bulk_batch_size": len(fastas) if self._is_bulk else 1,
                },
                # The draft keeps only the name the user TYPED, never `name`:
                # that carries this FASTA's stem and record id, and a Re-run
                # pointed at a different FASTA would inherit both. A bulk or
                # fanned run never used the box, so its name re-derives.
                query_draft={
                    **base_draft, "fasta_path": fasta,
                    "run_name": ("" if self._is_bulk or _fanned
                                 else base_draft.get("run_name", "")),
                },
            )
            out_runs.append((record, cfg, cached_path or None))
        return out_runs

    # Backwards-compat shim for callers that still expect a single tuple.
    def build(self) -> Tuple[RunRecord, RunConfig, Optional[str]]:
        runs = self.build_runs()
        if not runs:
            raise RuntimeError("build_runs() returned an empty list")
        return runs[0]
