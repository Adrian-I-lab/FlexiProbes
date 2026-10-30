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

"""Filter explorer dialog (U8).

Shows why tiles were rejected for a run: the rejection counts per category
plus concrete example rejected tiles, loaded from ``filter_report.json``.
"""
from __future__ import annotations

import json
import os
from typing import Optional

from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QHBoxLayout, QHeaderView, QLabel,
    QTableWidget, QTableWidgetItem, QVBoxLayout,
)

from .theme import muted_hex


def _examples_for(examples: dict, reason: str) -> list:
    """Every example tile filed under ``reason``, including its per-half keys.

    ``filters.py`` COUNTS under a bare key ('gc_range', 'tm_range', 'n_count')
    but FILES the example under a side-suffixed one ('gc_range_left',
    'tm_range_right', 'n_count_left'). Looking the bare key up directly — which
    is what this dialog used to do — therefore returned nothing for the three
    commonest rejection reasons, so the example browser sat permanently empty
    on exactly the runs a user opens it to understand.
    """
    out = list(examples.get(reason) or [])
    prefix = reason + "_"
    for key in sorted(examples):
        if key.startswith(prefix):
            out.extend(examples[key])
    return out


_REASON_HELP = {
    "gap_or_invalid": "The window contains an alignment gap ('-' or '.').",
    "junction_base": ("A degenerate base (N or another IUPAC code) lies "
                      "within the junction clearance band, 4 bp either side "
                      "of the ligation nick by default."),
    "junction_mismatch": ("On a homology run, the input records disagree within "
                          "the junction clearance. The ligase will not seal "
                          "across a mispair there, so the probe's behaviour on "
                          "the other records is unpredictable."),
    "degenerate_unresolved": ("A degenerate consensus parent arrived without a "
                              "usable per-record variant set. This is a "
                              "pipeline fault, not a chemistry rejection — "
                              "report it if the count is non-zero."),
    "tn_junction": ("The ligation-junction base failed the active gate — "
                    "'strict_a' rejects anything but A at mRNA position 25, "
                    "'soft' rejects only C."),
    "n_count": ("Too many degenerate bases (N or other IUPAC codes) in "
                "one half."),
    "gc_range": "GC content of a half is outside the allowed window.",
    "tm_range": ("A half's melting temperature is outside the allowed window. "
                 "Relative Tm (primer3, PCR-primer salt), not an in-assay value."),
    "homopolymer": "A homopolymer run exceeds the limit.",
    "passed": "Tiles that passed every filter.",
    "total": "Total tiles considered.",
}


class FilterReportDialog(QDialog):
    """Read-only viewer for one run's filter_report.json."""

    def __init__(self, output_dir: str, run_name: str = "", parent=None):
        super().__init__(parent)
        self._output_dir = output_dir
        self.setWindowTitle(f"Filter report: {run_name}" if run_name else "Filter report")
        self.setModal(True)
        self.resize(880, 600)
        self.setMinimumSize(560, 400)
        self.setSizeGripEnabled(True)

        report = self._load()
        outer = QVBoxLayout(self)

        if report is None:
            outer.addWidget(QLabel(
                "No readable filter_report.json for this run (missing or "
                "damaged). Runs made by older versions do not write one."
            ))
            bb = QDialogButtonBox(QDialogButtonBox.Close)
            bb.rejected.connect(self.reject); bb.accepted.connect(self.accept)
            outer.addWidget(bb)
            return

        counts = report.get("counts", {})
        thr = report.get("thresholds", {})
        total = counts.get("total", 0)
        passed = counts.get("passed", 0)

        # Headline summary.
        #
        # The per-rule counts OVERLAP. filters.py checks every rule on every
        # tile and counts a failing tile once under EACH rule it fails, so the
        # rule counts normally sum to more than the rejected tiles. Rejections
        # are therefore total - passed, never a subtraction of the rule sum.
        # (Until 2026-09-22 this printed total - passed - sum as "rejected for
        # a reason not listed below", which was negative on real runs: -127 on
        # the OROV S conserved panel run, -3937 on the ProbeST runs.)
        rejected = total - passed
        rule_failures = sum(
            v for k, v in counts.items() if k not in ("total", "passed")
        )
        # Every rejected tile fails at least one counted rule, so the rule sum
        # can never honestly fall BELOW the rejections. If it does, some
        # rejected tiles carry no recorded reason: an accounting gap in the
        # engine, and the one case worth a line of its own.
        unexplained = max(0, rejected - rule_failures)
        summary = QLabel(
            (f"<b>{total}</b> tiles considered · <b>{passed}</b> passed "
             f"({(100*passed/total):.1f}%) · <b>{rejected}</b> rejected"
             + (f" · <b>{unexplained}</b> rejected tiles have no recorded "
                "reason (please report this)"
                if unexplained else "")
             + "<br>")
            if total else "No tiles."
        )
        self.summary_label = summary
        outer.addWidget(summary)

        # The GC band is a ranking term unless enforced, so say which, in
        # words: this line used to print Python's own "gate: False".
        _gc_gate = thr.get("enforce_gc_range")
        _gc_gate_text = ("enforced" if _gc_gate is True
                         else "not enforced" if _gc_gate is False else "?")
        # tm_min None means the Tm gate was switched off, and printed as
        # "Tm/half >= None" until 2026-09-22.
        if "tm_min" not in thr:
            _tm_text = ">= ?"
        elif thr.get("tm_min") is None:
            _tm_text = "off"
        else:
            _tm_text = f">= {thr['tm_min']}"
        thr_lbl = QLabel(
            "Thresholds: GC {gc_min} to {gc_max} ({gc_gate}), "
            "Tm/half {tm_text}, max degenerate bases/half {max_n_per_half}, "
            "max homopolymer {max_homopolymer_run}, "
            "TN junction {tn_junction_mode}".format(
                **{**{"gc_min": "?", "gc_max": "?", "max_n_per_half": "?",
                      "max_homopolymer_run": "?",
                      "tm_min": "?", "tn_junction_mode": "?"}, **thr,
                   "gc_gate": _gc_gate_text, "tm_text": _tm_text}
            )
        )
        thr_lbl.setStyleSheet(f"color:{muted_hex()};")
        outer.addWidget(thr_lbl)

        # Counts table.
        # Must cover every key filters.new_filter_counters() can raise, or a
        # whole rejection class disappears from the report and the arithmetic
        # (total - passed) stops reconciling.
        order = ["gc_range", "tm_range", "homopolymer", "tn_junction",
                 "n_count", "junction_base", "junction_mismatch",
                 "degenerate_unresolved", "gap_or_invalid"]
        present = [k for k in order if counts.get(k)]
        ct = QTableWidget(len(present), 3)
        ct.setHorizontalHeaderLabels(
            ["Reason", "Tiles failing this rule", "What it means"])
        ct.verticalHeader().setVisible(False)
        ct.setEditTriggers(QTableWidget.NoEditTriggers)
        for i, k in enumerate(present):
            ct.setItem(i, 0, QTableWidgetItem(k))
            ct.setItem(i, 1, QTableWidgetItem(str(counts.get(k, 0))))
            ct.setItem(i, 2, QTableWidgetItem(_REASON_HELP.get(k, "")))
        ct.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        ct.resizeColumnsToContents()
        ct.setMaximumHeight(190)
        outer.addWidget(ct)
        overlap_note = QLabel(
            "A tile can fail several rules at once, and is counted under each "
            "one, so this column does not add up to the rejected total.")
        overlap_note.setWordWrap(True)
        overlap_note.setStyleSheet(f"color:{muted_hex()};")
        outer.addWidget(overlap_note)

        # Example explorer.
        row = QHBoxLayout()
        row.addWidget(QLabel("Examples for:"))
        self.reason_combo = QComboBox()
        examples = report.get("examples", {})
        for k in present:
            n_ex = len(_examples_for(examples, k))
            if n_ex:
                self.reason_combo.addItem(f"{k}  ({n_ex} shown)", userData=k)
        row.addWidget(self.reason_combo, 1)
        outer.addLayout(row)

        self.ex_table = QTableWidget(0, 6)
        self.ex_table.setHorizontalHeaderLabels(
            ["probe_id", "start", "GC", "GC L", "GC R", "sequence"]
        )
        self.ex_table.verticalHeader().setVisible(False)
        self.ex_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.ex_table.horizontalHeader().setSectionResizeMode(5, QHeaderView.Stretch)
        outer.addWidget(self.ex_table, 1)

        self._examples = examples
        self.reason_combo.currentIndexChanged.connect(self._show_examples)
        if self.reason_combo.count():
            self._show_examples(0)

        bb = QDialogButtonBox(QDialogButtonBox.Close)
        bb.rejected.connect(self.reject); bb.accepted.connect(self.accept)
        outer.addWidget(bb)

    def _load(self) -> Optional[dict]:
        path = os.path.join(self._output_dir, "filter_report.json")
        if not os.path.exists(path):
            return None
        try:
            with open(path) as fh:
                return json.load(fh)
        except (OSError, json.JSONDecodeError):
            return None

    def _show_examples(self, _idx) -> None:
        key = self.reason_combo.currentData()
        rows = _examples_for(self._examples, key) if key else []
        self.ex_table.setRowCount(0)
        for r in rows:
            i = self.ex_table.rowCount()
            self.ex_table.insertRow(i)
            self.ex_table.setItem(i, 0, QTableWidgetItem(str(r.get("probe_id", ""))))
            self.ex_table.setItem(i, 1, QTableWidgetItem(str(r.get("start", ""))))
            self.ex_table.setItem(i, 2, QTableWidgetItem(f"{r.get('gc_content', 0):.2f}"))
            self.ex_table.setItem(i, 3, QTableWidgetItem(f"{r.get('gc_left', 0):.2f}"))
            self.ex_table.setItem(i, 4, QTableWidgetItem(f"{r.get('gc_right', 0):.2f}"))
            self.ex_table.setItem(i, 5, QTableWidgetItem(str(r.get("sequence", ""))))
        self.ex_table.resizeColumnsToContents()
