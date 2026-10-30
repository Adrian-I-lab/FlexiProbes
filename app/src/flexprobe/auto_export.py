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

"""Auto-select and export: one implementation, two front ends.

``flexiprobes design --auto-select-per-gene N`` and the GUI's whole-
transcriptome mode are the same operation: pick the top N spaced probes per
gene from a finished run's ``results.csv``, save the selection, and write the
order CSVs. Both front ends call this module, so they cannot disagree on:

1. **A raised exporter.** Reported as ``EXPORT_FAILED``, never as a finished,
   exported run.
2. **An empty pick.** Reported as ``NOTHING_SELECTED``, never silently.
3. **Spacing reporting.** The same green/orange/red tier lines and the same
   red-tier warning, whichever front end prints them.
4. **Where the gaps come from.** The caller passes the run's actual
   ``--min-probe-gap`` / ``--relaxed-probe-gap`` / ``--floor-probe-gap``.
   This module never substitutes the defaults (150/40/15) on its own, because
   that would silently pick at a spacing the user did not ask for.

The selection rule itself is :func:`flexprobe.auto_select.auto_pick_with_report`.
The risk this module removes is not in which probes are picked but in what each
front end then does and says about the result, for example reporting a failed
export as a success.

This module therefore owns the operation and its *facts*; the two callers own
only their output channel (stdout/stderr for the CLI, the logger and the run
record for the GUI). Reporting text lives here so the two cannot say different
things about the same run.
"""
from __future__ import annotations

import enum
from dataclasses import dataclass
from typing import List, Mapping, Optional

from .auto_select import AutoPick, auto_pick_with_report
from .constants import (
    GAP_TIER_RED, MIN_PROBE_GAP_BP, MIN_PROBE_GAP_FLOOR_BP,
    PROBE_GAP_RELAX_MIN_PROBES, RELAXED_PROBE_GAP_BP,
)
from . import select_export


class AutoExportStatus(enum.Enum):
    """What the auto-export actually managed to do.

    A bool cannot carry this: "did not export" splits into three cases that
    want three different things said, and a caller deciding a run's final
    status has to be able to tell them apart. Only ``EXPORTED`` means order
    CSVs exist on disk; every other value leaves the candidates in
    ``results.csv`` for a manual pick.
    """

    EXPORTED = "exported"                    # order CSVs written
    NO_RESULTS = "no_results"                # results.csv missing entirely
    NOTHING_SELECTED = "nothing_selected"    # the picker returned no probes
    EXPORT_FAILED = "export_failed"          # the writer raised


def spacing_tier_lines(counts: Mapping[str, int], min_gap: int,
                       relaxed_gap: int, floor_gap: int) -> List[str]:
    """One line per spacing tier, always all three, always in the same order.

    Emitted even when a tier is empty: "orange 0" is the reassurance a user is
    looking for, and a line that vanishes when the news is good makes them go
    and check whether the tool simply forgot to say.

    Returned rather than printed so the CLI can put them on stdout and the GUI
    can put them in its log without the wording drifting between the two.
    """
    total = sum(counts.values()) or 1
    rows = (("green", f">= {min_gap} bp", "as designed"),
            ("orange", f"{relaxed_gap}-{min_gap} bp", "tight; usable"),
            ("red", f"{floor_gap}-{relaxed_gap} bp", "crowded; inspect these"))
    lines = ["probe spacing:"]
    for name, band, note in rows:
        n = counts.get(name, 0)
        lines.append(f"  {name:<7}{band:>12}  {n:>8,}  "
                     f"({100 * n / total:4.1f}%)  {note}")
    return lines


@dataclass
class AutoExportResult:
    """The outcome, plus everything either front end needs to report it."""

    status: AutoExportStatus
    n_per_gene: int = 0
    min_gap: int = MIN_PROBE_GAP_BP
    relaxed_gap: int = RELAXED_PROBE_GAP_BP
    floor_gap: int = MIN_PROBE_GAP_FLOOR_BP
    pick: Optional[AutoPick] = None
    error: Optional[BaseException] = None

    @property
    def exported(self) -> bool:
        return self.status is AutoExportStatus.EXPORTED

    @property
    def n_probes(self) -> int:
        if self.pick is None:
            return 0
        return sum(len(v) for v in self.pick.selections.values())

    @property
    def n_genes(self) -> int:
        return 0 if self.pick is None else len(self.pick.selections)

    def report_lines(self) -> List[str]:
        """The good news: what was picked, and how it is spaced. Empty unless
        the order CSVs actually got written."""
        if not self.exported or self.pick is None:
            return []
        lines = [f"auto-selected {self.n_probes} probes across "
                 f"{self.n_genes} genes (top {self.n_per_gene}/gene)."]
        lines.extend(spacing_tier_lines(
            self.pick.tier_counts(), self.min_gap, self.relaxed_gap,
            self.floor_gap))
        return lines

    def warning_lines(self) -> List[str]:
        """Everything the user has to be told, whether or not it went well.

        The red-tier warning names genes because a crowded gene is something a
        user can act on (order fewer probes for it); the orange count in the
        banner is left as a count because "tight but usable" needs no action.
        """
        if self.status is AutoExportStatus.NO_RESULTS:
            return ["no results.csv to auto-select from; nothing was exported."]
        if self.status is AutoExportStatus.NOTHING_SELECTED:
            return ["auto-select found no probes to pick."]
        if self.status is AutoExportStatus.EXPORT_FAILED:
            return [f"WARNING: auto-export failed ({self.error}); "
                    f"results.csv left for manual review."]
        if self.pick is None:
            return []
        counts = self.pick.tier_counts()
        if not counts.get("red"):
            return []
        red = self.pick.genes_by_tier().get(GAP_TIER_RED, [])
        shown = ", ".join(red[:5]) + (" ..." if len(red) > 5 else "")
        return [f"WARNING: {counts['red']:,} probe(s) across {len(red):,} "
                f"gene(s) sit closer than {self.relaxed_gap} bp — the last "
                f"rung of the spacing ladder. Their transcripts cannot hold "
                f"{PROBE_GAP_RELAX_MIN_PROBES} separated probes; consider "
                f"fewer probes for these genes: {shown}"]


def auto_export(outdir: str, n_per_gene: int, *,
                min_gap: int = MIN_PROBE_GAP_BP,
                relaxed_gap: int = RELAXED_PROBE_GAP_BP,
                floor_gap: int = MIN_PROBE_GAP_FLOOR_BP,
                write_idt: bool = True,
                write_cellranger: bool = True) -> AutoExportResult:
    """Pick the top ``n_per_gene`` spaced probes per gene and write the orders.

    Reads ``results.csv`` from ``outdir``, saves ``selections.json`` and writes
    the IDT / Cell Ranger order CSVs. Never raises for an expected failure —
    the caller gets an :class:`AutoExportResult` and decides what to say and
    what status to record.

    The gaps are parameters with no clever defaulting on purpose: this function
    cannot see whether its caller *has* the run's real spacing flags, so a
    caller that silently substitutes the module defaults must do so in the
    open, at its own call site, where it can say it is doing it.
    """
    def _res(status: AutoExportStatus, pick: Optional[AutoPick] = None,
             error: Optional[BaseException] = None) -> AutoExportResult:
        return AutoExportResult(
            status, n_per_gene=n_per_gene, min_gap=min_gap,
            relaxed_gap=relaxed_gap, floor_gap=floor_gap, pick=pick,
            error=error)

    try:
        rows = select_export.load_results_rows(outdir)
    except FileNotFoundError:
        return _res(AutoExportStatus.NO_RESULTS)
    # Real transcript lengths where the run recorded them: the spread target
    # scales with the gene, and inferring the length from the furthest passing
    # probe understates a transcript whose 3' end yielded nothing.
    try:
        lengths = select_export.load_gene_lengths(outdir)
    except Exception:
        lengths = {}
    pick = auto_pick_with_report(
        rows, n_per_gene=n_per_gene,
        min_gap=min_gap, relaxed_gap=relaxed_gap, floor_gap=floor_gap,
        gene_length_of=(lambda g: lengths.get(g, 0)) if lengths else None,
    )
    if not pick.selections:
        return _res(AutoExportStatus.NOTHING_SELECTED, pick)
    select_export.save_selections(outdir, pick.selections, note="auto-picked")
    try:
        # overwrite=True on purpose: these are the order files THIS run just
        # wrote in full, replaced by the selected subset. Whether pre-existing
        # outputs may be overwritten was already decided when the run wrote
        # them (--force). The two switches are the run's own (--no-idt /
        # --no-cellranger), so a run that declined an order format does not
        # get it written here.
        select_export.export_selected_probes(
            outdir, pick.selections, overwrite=True,
            write_idt=write_idt, write_cellranger=write_cellranger)
    except Exception as exc:  # noqa: BLE001 — handed back to the caller verbatim
        return _res(AutoExportStatus.EXPORT_FAILED, pick, exc)
    return _res(AutoExportStatus.EXPORTED, pick)
