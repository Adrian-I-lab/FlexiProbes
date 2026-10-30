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

"""Completed-runs table for HomeWindow.

Split out of home_window.py (audit #19): the paginated history table, its cell
renderers and page navigation. ``CompletedTableMixin`` reads the host window's
widgets/state via ``self``; the module-level ``_item``/``_status_*``/``_duration``
helpers and ``PAGE_SIZE`` are shared with home_window."""
from __future__ import annotations

import os
from datetime import datetime

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QHBoxLayout, QLabel, QPushButton, QTableWidgetItem, QWidget,
)

from .run_history import RunRecord
from .run_import import SOURCE_HPC, import_source, is_folder_missing
from .run_manager import is_resumable
from .theme import link_hex

PAGE_SIZE = 20


def _section_label(text: str) -> QLabel:
    lbl = QLabel(text)
    f = lbl.font(); f.setBold(True); f.setPointSize(11); lbl.setFont(f)
    lbl.setStyleSheet("margin-top: 4px;")
    return lbl


def _item(text: str) -> QTableWidgetItem:
    item = QTableWidgetItem(text or "")
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    return item


def _status_label(status: str) -> str:
    """Friendlier text for the home-page Status column."""
    return {
        "ready_for_review": "ready for review",
        "interrupted": "interrupted",
    }.get(status, status)


def _status_cell_widget(status_text: str, accent: str,
                        tooltip: str = "", badge: str = "HPC",
                        badge_tip: str = "Imported from an HPC run") -> QWidget:
    """Status text followed by a small pill (imported-run badge).

    ``badge`` is "HPC" for a cluster run and "Imported" for a desktop run folder
    brought in by "Import runs…"; the rule is in ``run_import``.

    This is an OVERLAY. ``setCellWidget`` does not replace the cell's item, it
    draws on top of it, and this widget has no background of its own, so the
    item underneath must be blank or both strings paint (see the caller). The
    widget owns the status text for imported runs; nothing else may also carry
    it.

    ``tooltip`` goes on the widget rather than the item because the widget
    swallows the hover, which would otherwise make a failed imported run's
    error message unreachable.
    """
    w = QWidget()
    lay = QHBoxLayout(w)
    lay.setContentsMargins(4, 0, 4, 0)
    lay.setSpacing(6)
    lay.addWidget(QLabel(status_text))
    pill = QLabel(badge)
    pill.setToolTip(badge_tip)
    pill.setStyleSheet(
        f"background: {accent}; color: white; border-radius: 6px; "
        f"padding: 0px 6px; font-size: 10px; font-weight: bold;"
    )
    lay.addWidget(pill)
    lay.addStretch(1)
    if tooltip:
        w.setToolTip(tooltip)
    return w


# A traceback is thousands of characters; a tooltip that long is unreadable
# and on some platforms simply does not render. Show the tail, which is where
# the actual exception message is.
_ERROR_TOOLTIP_CHARS = 1200


def _error_tooltip(error: str) -> str:
    text = (error or "").strip()
    if len(text) > _ERROR_TOOLTIP_CHARS:
        text = "…\n" + text[-_ERROR_TOOLTIP_CHARS:]
    return "Why this run failed:\n\n" + text


#: Column positions in the home table. The per-run buttons sit BEFORE Output:
#: Output holds a full folder path, and as the last-but-one column it pushed
#: every button off the right edge of a 1280 px window.
COL_ACTIONS = 8
COL_OUTPUT = 9


def _started_text(stamp: str) -> str:
    """ISO timestamp -> 'YYYY-MM-DD HH:MM'. Unparseable input is shown as is."""
    try:
        return datetime.fromisoformat(stamp).strftime("%Y-%m-%d %H:%M")
    except (ValueError, TypeError):
        return stamp or ""


def _duration(rec: RunRecord) -> str:
    if not rec.started_at or not rec.completed_at:
        return ""
    try:
        s = datetime.fromisoformat(rec.started_at)
        c = datetime.fromisoformat(rec.completed_at)
        secs = int((c - s).total_seconds())
    except (ValueError, TypeError):
        return ""
    if secs < 60:
        return f"{secs}s"
    mins, rem = divmod(secs, 60)
    if mins < 60:
        return f"{mins}m {rem}s"
    hrs, mins = divmod(mins, 60)
    return f"{hrs}h {mins}m"


from typing import TYPE_CHECKING

if TYPE_CHECKING:  # give the mixin the host widget's type for mypy
    from PySide6.QtWidgets import QMainWindow as _HostBase
else:
    _HostBase = object


class CompletedTableMixin(_HostBase):
    _page: int          # set on the host HomeWindow; declared for mypy
    def _refresh_completed(self) -> None:
        completed = self.manager.history.list_completed()
        if self._project_filter:
            completed = [r for r in completed if r.project_id == self._project_filter]
        # Project-name lookup keeps the column readable.
        names_by_id = {p.project_id: p.name for p in self.projects.list()}

        n_pages = max(1, (len(completed) + PAGE_SIZE - 1) // PAGE_SIZE)
        if self._page >= n_pages:
            self._page = n_pages - 1
        start = self._page * PAGE_SIZE
        page_rows = completed[start:start + PAGE_SIZE]

        self.table.setRowCount(0)
        # Suppress the itemChanged handler while we populate check cells.
        self._populating_runs = True
        for rec in page_rows:
            row = self.table.rowCount()
            self.table.insertRow(row)
            chk = QTableWidgetItem()
            chk.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            chk.setCheckState(Qt.Unchecked)
            chk.setData(Qt.UserRole, rec.run_id)   # survives sort/page order
            self.table.setItem(row, 0, chk)
            _stamp = rec.started_at or rec.created_at
            started = _item(_started_text(_stamp))
            started.setToolTip(_stamp or "")
            self.table.setItem(row, 1, started)
            self.table.setItem(row, 2, _item(rec.name))
            self.table.setItem(row, 3, _item(names_by_id.get(rec.project_id, "")))
            self.table.setItem(row, 4, _item(rec.modality_label))
            # One isdir per displayed row, on every rebuild (start-up and
            # Refresh included). A moved or deleted folder used to leave a row
            # whose Open did nothing and whose Review failed with "results.csv
            # not found". Nothing else is scanned here, by decision: runs
            # arrive only through Import runs… and Locate….
            missing = is_folder_missing(rec)
            status_text = _status_label(rec.status)
            # run_manager stores the failure on the record but nothing ever read
            # it back, so a failed run showed the word "failed" and no reason
            # anywhere in the app — including the k-mer engine's "this looks
            # like a genome, switch engines" message, which is actionable.
            tip = (_error_tooltip(rec.error)
                   if rec.status == "failed" and rec.error else "")
            if missing:
                status_text = "Folder missing"
                tip = (f"The output folder is no longer at:\n{rec.output_dir}\n\n"
                       f"Use Locate… to point FlexiProbes at where it went."
                       + (f"\n\n{tip}" if tip else ""))
            source = import_source(rec)
            # Badge runs that were imported from an HPC job so they stand out
            # from locally-run designs.
            #
            # The badge is a cell WIDGET, and a cell widget is an overlay: the
            # view goes on painting the item's own text underneath it. With the
            # status set on both, the two strings rendered a few pixels apart
            # and read as bold, doubled, smudged text. Local runs never took
            # this branch, which is why only imported runs looked wrong. Keep
            # the item — it carries the cell's flags and selection behaviour —
            # but leave it blank and let the widget own the text.
            if source:
                placeholder = _item("")
                self.table.setItem(row, 5, placeholder)
                if source == SOURCE_HPC:
                    badge, badge_tip = "HPC", "Imported from an HPC run"
                else:
                    badge, badge_tip = ("Imported",
                                        "Imported from a run folder on disk")
                self.table.setCellWidget(
                    row, 5, _status_cell_widget(status_text, link_hex(), tip,
                                                badge, badge_tip)
                )
            else:
                status_item = _item(status_text)
                if tip:
                    status_item.setToolTip(tip)
                self.table.setItem(row, 5, status_item)
            pairs = (rec.stats or {}).get("pairs_output", "")
            self.table.setItem(row, 6, _item(str(pairs) if pairs != "" else ""))
            self.table.setItem(row, 7, _item(_duration(rec)))
            out_item = _item(rec.output_dir)
            out_item.setToolTip(rec.output_dir or "")
            if missing:
                out_item.setToolTip(
                    f"{rec.output_dir}\n\nThis folder no longer exists.")
            self.table.setItem(row, COL_OUTPUT, out_item)

            actions = QWidget()
            ah = QHBoxLayout(actions)
            ah.setContentsMargins(2, 0, 2, 0)
            ah.setSpacing(4)
            if missing:
                # Review, Filters and Open all need the folder, so none of them
                # is offered; Locate… is the way back. Re-run and Delete still
                # work from the record alone.
                locate_btn = QPushButton("Locate…")
                locate_btn.setStyleSheet("font-weight: bold;")
                locate_btn.setToolTip(
                    "Point FlexiProbes at where this run's folder is now. The "
                    "folder is checked against this run's manifest before it "
                    "is accepted."
                )
                locate_btn.clicked.connect(
                    lambda _=False, r=rec: self._on_locate_run(r)
                )
                ah.addWidget(locate_btn)
            if not missing and rec.status in ("ready_for_review", "complete") and \
                    (rec.stats or {}).get("discovery_only"):
                review_btn = QPushButton("Review")
                if rec.status == "ready_for_review":
                    review_btn.setStyleSheet("font-weight: bold;")
                review_btn.clicked.connect(
                    lambda _=False, r=rec: self._open_review(r)
                )
                ah.addWidget(review_btn)
            if not missing and \
                    os.path.exists(os.path.join(rec.output_dir, "filter_report.json")):
                filt_btn = QPushButton("Filters")
                filt_btn.setToolTip("Why were tiles rejected? (filter explorer)")
                filt_btn.clicked.connect(
                    lambda _=False, r=rec: self._open_filter_report(r)
                )
                ah.addWidget(filt_btn)
            if rec.status == "failed" and rec.error:
                why_btn = QPushButton("Why?")
                why_btn.setToolTip("Show the error this run failed with.")
                why_btn.clicked.connect(
                    lambda _=False, r=rec: self._show_run_error(r)
                )
                ah.addWidget(why_btn)
            if not missing:
                open_btn = QPushButton("Open")
                open_btn.clicked.connect(
                    lambda _=False, r=rec: self._on_open_folder(r))
                ah.addWidget(open_btn)
            # Resume: only for interrupted/cancelled/failed runs that still have
            # saved progress on disk (candidate/screen checkpoints + config).
            if rec.status in ("interrupted", "cancelled", "failed") and \
                    is_resumable(rec):
                resume_btn = QPushButton("Resume")
                resume_btn.setStyleSheet("font-weight: bold;")
                resume_btn.setToolTip(
                    "Continue this run from where it stopped, reusing the saved "
                    "tiling + screen checkpoints in its output folder."
                )
                resume_btn.clicked.connect(lambda _=False, r=rec: self._on_resume(r))
                ah.addWidget(resume_btn)
            rerun_btn = QPushButton("Re-run…")
            rerun_btn.setToolTip(
                "Re-open New query pre-filled with this run's settings so you "
                "can tweak and launch it again."
            )
            rerun_btn.clicked.connect(lambda _=False, r=rec: self._on_rerun(r))
            ah.addWidget(rerun_btn)
            del_btn = QPushButton("Delete")
            del_btn.setToolTip("Remove this run from the list.")
            del_btn.clicked.connect(lambda _=False, r=rec: self._on_delete_run(r))
            ah.addWidget(del_btn)
            self.table.setCellWidget(row, COL_ACTIONS, actions)
        self._populating_runs = False
        # Auto-size columns once, on the first populated refresh; afterwards
        # leave them alone so a user's manual column resizing survives reloads.
        if page_rows and not getattr(self, "_columns_sized", False):
            self.table.resizeColumnsToContents()
            self._columns_sized = True

        # Rebuilt rows start unchecked; reset the bulk controls to match.
        self._sync_run_bulk_buttons()
        self.page_label.setText(f"Page {self._page + 1} / {n_pages}")
        self.prev_btn.setEnabled(self._page > 0)
        self.next_btn.setEnabled(self._page < n_pages - 1)

    def _show_run_error(self, rec: RunRecord) -> None:
        """Full failure text for a failed run, selectable so it can be pasted."""
        from PySide6.QtWidgets import QMessageBox
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Critical)
        box.setWindowTitle(f"Run failed: {rec.name}")
        box.setText(f"'{rec.name}' did not finish.")
        box.setDetailedText(rec.error or "(no error text was recorded)")
        box.exec()

    def _prev_page(self) -> None:
        if self._page > 0:
            self._page -= 1
            self._refresh_completed()

    def _next_page(self) -> None:
        self._page += 1
        self._refresh_completed()
