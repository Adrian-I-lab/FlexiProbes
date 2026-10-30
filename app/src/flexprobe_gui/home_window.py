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

"""Home window: New-query launcher, active-runs cards, paginated history table."""
from __future__ import annotations

import logging
import os
from typing import Dict

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QFileDialog, QFrame, QHBoxLayout, QHeaderView,
    QInputDialog, QLabel, QMainWindow, QMessageBox, QProgressBar, QPushButton,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from flexprobe.blast import load_cached_hits
from flexprobe.select_export import export_project_bulk, load_selections

from .manage_projects_dialog import ManageProjectsDialog
from .theme import muted_hex
from ._nq_build import summary_to_draft
from .new_query_dialog import NewQueryDialog

logger = logging.getLogger("flexprobe_gui.home_window")
from .project_history import ProjectHistory
from .query_drafts import QueryDraftStore
from .reveal import reveal_in_finder
from ._home_dnd import DragDropMixin
from ._home_completed import CompletedTableMixin, _section_label
from .review_dialog import ReviewDialog
from .run_history import RunHistory, RunRecord, now_iso
from .run_import import (
    CONFIRM, REJECT, ImportSummary, import_runs, is_folder_missing,
    match_folder, nearest_existing_parent, relink, sibling_candidates,
    stamp_identity,
)
from .run_manager import RunManager
from .settings import AppSettings
from .tools_dialog import (
    BuildDBDialog, PreferencesDialog, ReExportDialog, ValidateBlastDialog,
    show_about,
)



class _RunCard(QFrame):
    """A row showing one active run: name, stage, progress, cancel button."""

    def __init__(self, record: RunRecord, on_cancel, parent=None):
        super().__init__(parent)
        self.setFrameShape(QFrame.StyledPanel)
        self.setFrameShadow(QFrame.Raised)
        self._run_id = record.run_id

        lay = QHBoxLayout(self)
        lay.setContentsMargins(10, 6, 10, 6)

        info = QVBoxLayout()
        self.title = QLabel()
        f = self.title.font(); f.setBold(True); self.title.setFont(f)
        info.addWidget(self.title)
        self.stage = QLabel()
        self.stage.setStyleSheet(f"color: {muted_hex()};")
        info.addWidget(self.stage)
        lay.addLayout(info, 3)

        self.bar = QProgressBar()
        self.bar.setRange(0, 100)
        self.bar.setMinimumWidth(220)
        lay.addWidget(self.bar, 3)

        self.cancel_btn = QPushButton("Pause / Cancel")
        self.cancel_btn.setToolTip(
            "Stop this run now. Progress so far (tiling + screen checkpoints) is "
            "kept, so you can pick it back up later with “Resume” — handy for "
            "freeing the machine to work on something else."
        )
        self.cancel_btn.clicked.connect(lambda: on_cancel(self._run_id))
        lay.addWidget(self.cancel_btn)

        self.update_from(record)

    def update_from(self, record: RunRecord) -> None:
        self.title.setText(f"{record.name}  —  {record.modality_label}")
        if record.status == "queued":
            self.stage.setText("queued (waiting for current run to finish)")
            self.bar.setValue(0)
            self.bar.setFormat("queued")
        else:
            self.stage.setText(record.stage or record.status)
            self.bar.setValue(record.progress)
            self.bar.setFormat("%p%")


class HomeWindow(DragDropMixin, CompletedTableMixin, QMainWindow):
    # FASTA extensions accepted by drag-and-drop. Case-insensitive.
    FASTA_EXTS = (".fa", ".fasta", ".fna", ".ffn", ".frn")

    def __init__(self):
        super().__init__()
        self.setWindowTitle("FlexiProbes — Probe Designer")
        self.resize(1100, 800)
        self.setAcceptDrops(True)   # enable drag-and-drop of FASTAs / folders

        history = RunHistory()
        self.manager = RunManager(history, parent=self)
        self.manager.record_changed.connect(self._on_record_changed)
        self.projects = ProjectHistory()
        self.settings = AppSettings()
        self.drafts = QueryDraftStore()

        self._cards: Dict[str, _RunCard] = {}
        # Debounce timer: coalesces bursts of record-changed signals into one
        # completed-table rebuild (audit #14).
        self._completed_refresh_timer = QTimer(self)
        self._completed_refresh_timer.setSingleShot(True)
        self._completed_refresh_timer.setInterval(250)
        self._completed_refresh_timer.timeout.connect(self._refresh_completed)
        self._page = 0
        self._project_filter: str = ""   # "" = All projects
        # Guards the runs-table ``itemChanged`` handler while we rebuild rows,
        # so programmatic check-state changes don't re-enter the selection sync.
        self._populating_runs = False

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(16, 12, 16, 12)
        root.setSpacing(10)

        # --- Top bar: title + New query ---
        top = QHBoxLayout()
        title = QLabel("FlexiProbes")
        tf = title.font(); tf.setPointSize(18); tf.setBold(True); title.setFont(tf)
        top.addWidget(title)
        top.addStretch(1)
        self.new_btn = QPushButton("+ New query")
        self.new_btn.setMinimumHeight(34)
        self.new_btn.setMinimumWidth(140)
        self.new_btn.clicked.connect(self._on_new_query)
        top.addWidget(self.new_btn)
        self.bulk_btn = QPushButton("+ Bulk new query")
        self.bulk_btn.setMinimumHeight(34)
        self.bulk_btn.setMinimumWidth(160)
        self.bulk_btn.setToolTip(
            "Pick multiple FASTA files (or a folder) and run them all with "
            "the same settings. You can also drag-and-drop FASTAs onto the "
            "home window for the same effect."
        )
        self.bulk_btn.clicked.connect(self._on_bulk_new_query_click)
        top.addWidget(self.bulk_btn)
        self.tools_btn = QPushButton("Tools ▾")
        self.tools_btn.setMinimumHeight(34)
        self.tools_btn.setToolTip(
            "BLAST DB management, importing run folders and HPC runs, re-export, "
            "guides and About. The same jobs as the FlexiProbes CLI's "
            "housekeeping commands."
        )
        self.tools_btn.clicked.connect(self._on_tools_clicked)
        top.addWidget(self.tools_btn)
        root.addLayout(top)

        # --- Active runs ---
        root.addWidget(_section_label("Active runs"))
        self.active_wrap = QWidget()
        self.active_layout = QVBoxLayout(self.active_wrap)
        self.active_layout.setContentsMargins(0, 0, 0, 0)
        self.active_layout.setSpacing(6)
        self.active_empty = QLabel("No active runs.")
        self.active_empty.setStyleSheet(f"color: {muted_hex()}; padding: 4px 2px;")
        self.active_layout.addWidget(self.active_empty)
        root.addWidget(self.active_wrap)

        # --- Completed runs (with project filter + bulk export) ---
        root.addWidget(_section_label("Completed runs"))
        filt_row = QHBoxLayout()
        filt_row.addWidget(QLabel("Project filter:"))
        self.project_filter_combo = QComboBox()
        self.project_filter_combo.setMinimumWidth(220)
        filt_row.addWidget(self.project_filter_combo)
        self.bulk_export_btn = QPushButton("Bulk export project")
        self.bulk_export_btn.setEnabled(False)
        self.bulk_export_btn.clicked.connect(self._on_bulk_export)
        filt_row.addWidget(self.bulk_export_btn)
        filt_row.addStretch(1)
        self.manage_projects_btn = QPushButton("Manage projects…")
        self.manage_projects_btn.clicked.connect(self._on_manage_projects)
        filt_row.addWidget(self.manage_projects_btn)
        root.addLayout(filt_row)

        # --- Bulk selection actions (operate on the checked rows) ---
        bulk_row = QHBoxLayout()
        self.select_all_chk = QCheckBox("Select all")
        self.select_all_chk.setToolTip("Check or clear every run shown on this page.")
        self.select_all_chk.clicked.connect(self._on_select_all_runs)
        bulk_row.addWidget(self.select_all_chk)
        self.bulk_move_btn = QPushButton("Move to project…")
        self.bulk_move_btn.setEnabled(False)
        self.bulk_move_btn.setToolTip("Reassign the checked runs to a project.")
        self.bulk_move_btn.clicked.connect(self._on_bulk_move_runs)
        bulk_row.addWidget(self.bulk_move_btn)
        self.bulk_delete_btn = QPushButton("Delete selected")
        self.bulk_delete_btn.setEnabled(False)
        self.bulk_delete_btn.setToolTip("Remove the checked runs from the list.")
        self.bulk_delete_btn.clicked.connect(self._on_bulk_delete_runs)
        bulk_row.addWidget(self.bulk_delete_btn)
        self.bulk_count_lbl = QLabel("")
        self.bulk_count_lbl.setStyleSheet(f"color: {muted_hex()};")
        bulk_row.addWidget(self.bulk_count_lbl)
        bulk_row.addStretch(1)
        root.addLayout(bulk_row)

        # IMPORTANT: build the runs table BEFORE connecting the project-filter
        # signal. The previous order fired ``currentIndexChanged`` on the
        # first ``addItem`` and crashed because ``self.table`` didn't exist
        # yet.

        self.table = QTableWidget(0, 10)
        self.table.setHorizontalHeaderLabels([
            "", "Started", "Name", "Project", "Modality", "Status",
            "Pairs", "Duration", "", "Output",
        ])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.Interactive)
        # Output is last and STRETCHES: Qt gives it exactly the width that is
        # left, so a long path can never widen the table past the window (a
        # 360 px cap still overran). Long paths elide in the middle; the full
        # path is in the cell's tooltip.
        from ._home_completed import COL_OUTPUT
        hh.setSectionResizeMode(COL_OUTPUT, QHeaderView.Stretch)
        self.table.setTextElideMode(Qt.ElideMiddle)
        # Toggling a row's checkbox keeps the bulk-action buttons in sync.
        self.table.itemChanged.connect(self._on_run_item_changed)
        root.addWidget(self.table, 1)

        # --- Pagination ---
        pager = QHBoxLayout()
        pager.addStretch(1)
        self.prev_btn = QPushButton("← Prev")
        self.prev_btn.clicked.connect(self._prev_page)
        pager.addWidget(self.prev_btn)
        self.page_label = QLabel("Page 1 / 1")
        self.page_label.setMinimumWidth(110)
        self.page_label.setAlignment(Qt.AlignCenter)
        pager.addWidget(self.page_label)
        self.next_btn = QPushButton("Next →")
        self.next_btn.clicked.connect(self._next_page)
        pager.addWidget(self.next_btn)
        root.addLayout(pager)

        # Now that the runs table + page widgets exist, populate the project
        # filter and connect its signal. Doing this last avoids the startup
        # AttributeError when the combo's first ``addItem`` fired
        # ``currentIndexChanged`` while ``self.table`` was still None.
        self._refresh_project_filter()
        self.project_filter_combo.currentIndexChanged.connect(
            self._on_project_filter_changed,
        )

        self._install_shortcuts()
        self._refresh()

    def _install_shortcuts(self) -> None:
        """Wire standard keyboard shortcuts. Tooltips advertise them too."""
        def add(seq, handler) -> None:
            QShortcut(QKeySequence(seq), self).activated.connect(handler)

        add(QKeySequence.StandardKey.New, self._on_new_query)
        add("Ctrl+Shift+N", self._on_bulk_new_query_click)
        add("Ctrl+E", self._on_bulk_export)
        add(QKeySequence.StandardKey.Quit, self.close)
        add("Ctrl+W", self.close)
        add(QKeySequence.StandardKey.Refresh, self._refresh)
        # Surface the shortcuts in the buttons' tooltips for discoverability.
        self.new_btn.setToolTip("Start a new probe-design query (Ctrl+N)")
        self.bulk_btn.setToolTip("Queue many FASTAs / a folder (Ctrl+Shift+N)")
        self.bulk_export_btn.setToolTip(
            "Combine a project's saved selections (drafts included) into one "
            "order set (Ctrl+E)"
        )

    # ---- Actions ----

    def _on_new_query(self) -> None:
        self._open_new_query_dialog(bulk_paths=None)

    def _on_bulk_new_query_click(self) -> None:
        # Click-based alternative to drag-and-drop: pick many files at once.
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Pick FASTA files for a bulk run",
            "",
            "FASTA files (*.fa *.fasta *.fna *.ffn *.frn);;All files (*)",
        )
        if not paths:
            return
        valid, bad = self._filter_valid_fastas(paths)
        if bad:
            detail = "\n".join(f"  • {os.path.basename(p)} — {r}" for p, r in bad)
            if not valid:
                QMessageBox.warning(self, "Unreadable FASTA(s)",
                                    "No valid FASTA files selected:\n\n" + detail)
                return
            QMessageBox.warning(self, "Some FASTAs skipped",
                                f"Skipping invalid file(s):\n\n{detail}")
        self._open_new_query_dialog(bulk_paths=valid)

    def _open_new_query_dialog(self, bulk_paths) -> None:
        dlg = NewQueryDialog(
            self, projects=self.projects,
            bulk_fasta_paths=bulk_paths,
            settings=self.settings, drafts=self.drafts,
        )
        # If the home page is currently filtered to a project, default the
        # dialog's project picker to that project — convenient for drag-
        # and-dropping into the right bucket.
        if self._project_filter:
            idx = dlg.project_combo.findData(self._project_filter)
            if idx >= 0:
                dlg.project_combo.setCurrentIndex(idx)
        if dlg.exec() != NewQueryDialog.Accepted:
            return
        runs = dlg.build_runs()
        self._submit_built_runs(runs)
        self._warn_if_blast_gate_disabled(dlg)
        # The user may have created a new project from inside the dialog.
        self._refresh_project_filter()
        if len(runs) > 1:
            QMessageBox.information(
                self, "Bulk runs queued",
                f"Queued {len(runs)} runs. They run one at a time; watch "
                f"the Active runs section for progress.",
            )

    def _submit_built_runs(self, runs) -> None:
        """Submit a list of (record, config, cached_path) tuples to the queue.

        Shared by New query, drag-and-drop and Re-run so cached-hit loading
        and de-duplication stay in one place.
        """
        cached_hits_cache: Dict[str, object] = {}
        for record, config, cached_path in runs:
            cached_hits = None
            if cached_path:
                if cached_path in cached_hits_cache:
                    cached_hits = cached_hits_cache[cached_path]
                else:
                    try:
                        cached_hits = load_cached_hits(cached_path)
                        # Cache only on success. A failed load must NOT be
                        # memoised as None, or an identical path later in this
                        # batch would silently skip re-loading (and re-warning)
                        # and quietly run un-cached.
                        cached_hits_cache[cached_path] = cached_hits
                    except Exception as exc:
                        QMessageBox.warning(
                            self, "Cached hits",
                            f"Could not load cached hits: {exc}",
                        )
            self.manager.submit(record, config, cached_hits=cached_hits)

    def _on_cancel(self, run_id: str) -> None:
        self.manager.cancel(run_id)

    def _on_record_changed(self, run_id: str) -> None:
        # Progress ticks fire frequently and only affect ACTIVE runs, so refresh
        # the (cheap) active cards immediately but debounce the completed-runs
        # table — its rebuild does per-row disk stats (filter_report existence +
        # is_resumable), too heavy to run on every tick (audit #14).
        self._stamp_finished_run(run_id)
        self._refresh_active()
        self._completed_refresh_timer.start()

    def _stamp_finished_run(self, run_id: str) -> None:
        """Keep a finished run's manifest identity on its record.

        That identity is what lets "Locate…" prove a folder found elsewhere is
        this run. Taken once: only while the record has none, and only in the
        two statuses a run reaches after manifest.json is written, so a
        progress tick costs a dict lookup and no file access.
        """
        rec = self.manager.history.get(run_id)
        if rec is None or rec.status not in ("complete", "ready_for_review"):
            return
        if stamp_identity(rec):
            self.manager.history.update(rec)

    # ---- Tools menu (CLI equivalents) ----

    def _on_tools_clicked(self):
        from PySide6.QtWidgets import QMenu
        menu = QMenu(self)
        menu.addAction("Build BLAST database…",
                       lambda: BuildDBDialog(self).exec())
        menu.addAction("Validate BLAST install…",
                       lambda: ValidateBlastDialog(self).exec())
        menu.addAction("Re-export results.csv…",
                       lambda: ReExportDialog(self).exec())
        menu.addAction("Import runs…", self._on_import_runs)
        menu.addAction("Import HPC run…", self._on_import_hpc)
        menu.addSeparator()
        menu.addAction("Preferences…", self._on_preferences)
        menu.addAction("Input data guide…", self._show_input_guide)
        menu.addAction("Off-target reference guide…", self._show_blast_guide)
        menu.addAction("About FlexiProbes", lambda: show_about(self))
        # Anchor the menu just under the button.
        global_pos = self.tools_btn.mapToGlobal(
            self.tools_btn.rect().bottomLeft()
        )
        menu.exec(global_pos)

    def _on_import_hpc(self) -> None:
        from .hpc_dialog import ImportHpcDialog
        ImportHpcDialog(self.manager.history, on_imported=self._refresh, parent=self).exec()

    def _on_import_runs(self, root: str = "") -> "ImportSummary | None":
        """Tools > Import runs…: register the finished run folders under a folder.

        Manual by decision: nothing scans for runs at start-up or on Refresh.
        The badge and duplicate rules are in ``run_import``. ``root`` skips the
        folder picker (tests, and any future drop target).
        """
        if not root:
            root = QFileDialog.getExistingDirectory(
                self, "Folder holding FlexiProbes run folders",
                self.settings.default_output_dir(),
            )
            if not root:
                return None
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            summary = import_runs(root, self.manager.history)
        finally:
            QApplication.restoreOverrideCursor()
        self._refresh()
        box = QMessageBox(self)
        box.setWindowTitle("Import runs")
        box.setIcon(QMessageBox.Information)
        box.setText(summary.message())
        details = summary.details()
        if details:
            box.setDetailedText(details)
        box.exec()
        return summary

    def _on_preferences(self) -> None:
        PreferencesDialog(self.settings, self).exec()

    def _show_blast_guide(self) -> None:
        from .new_query_dialog import _BLAST_TARGETING_HELP_HTML
        box = QMessageBox(self)
        box.setWindowTitle("Off-target reference guide")
        box.setTextFormat(Qt.RichText)
        box.setText(_BLAST_TARGETING_HELP_HTML)
        box.exec()

    def _show_input_guide(self) -> None:
        """Show how to prepare input FASTAs + where the templates live."""
        here = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.normpath(os.path.join(here, "..", "..", "data"))
        templates = os.path.join(data_dir, "templates")
        text = (
            "<b>Preparing input</b><br>"
            "FlexiProbes takes FASTA files (.fa/.fasta/.fna). What each file "
            "holds depends on the mode:<br><br>"
            "<b>Standard (Flex v1/v2, Visium)</b> — one file, one record per "
            "gene; sequence = sense-strand mRNA/cDNA.<br>"
            "<b>Homology mode</b> — one file per fragment; each record is the "
            "same fragment from a different related sequence (strain/isoform). "
            "Put the actual individual sequences in. Every ordered probe is a "
            "concrete sequence taken from those records, so the order files "
            "contain only A, C, G and T.<br><br>"
        )
        # packaging/flexprobe.spec deliberately bundles no sample data, so in a
        # built .app these paths do not exist. Printing them anyway sent users
        # to a folder that was never there, and the "Open templates folder"
        # button silently did nothing. Say which case they are in instead.
        readme = os.path.join(data_dir, "README.md")
        if os.path.isdir(templates):
            text += (
                "Copy-and-edit templates are in:<br>"
                f"<code>{templates}</code><br><br>"
                f"Full guide: <code>{readme}</code>"
            )
        else:
            text += (
                "The example templates are not bundled with the packaged "
                "application. They live in <code>data/templates/</code> in the "
                "FlexiProbes source tree, alongside the full guide in "
                "<code>data/README.md</code>.<br><br>"
                "You do not need them: any FASTA matching the description "
                "above works."
            )
        box = QMessageBox(self)
        box.setWindowTitle("Input data guide")
        box.setTextFormat(Qt.RichText)
        box.setText(text)
        open_btn = None
        if os.path.isdir(templates):
            open_btn = box.addButton("Open templates folder",
                                     QMessageBox.ActionRole)
        box.addButton(QMessageBox.Close)
        box.exec()
        if open_btn is not None and box.clickedButton() is open_btn:
            reveal_in_finder(templates)

    # ---- Open review / filter dialogs for a completed run ----
    # (drag-and-drop FASTA intake lives in DragDropMixin — see _home_dnd.py)

    def _open_review(self, record: RunRecord) -> None:
        # The row hides Review for a missing folder, but the folder can go
        # between the table being drawn and the click.
        if is_folder_missing(record):
            self._warn_folder_missing(record)
            return
        dlg = ReviewDialog(record, self)
        dlg.confirmed.connect(self._on_review_confirmed)
        dlg.exec()

    def _open_filter_report(self, record: RunRecord) -> None:
        if is_folder_missing(record):
            self._warn_folder_missing(record)
            return
        from .filter_report_dialog import FilterReportDialog
        FilterReportDialog(record.output_dir, record.name, self).exec()

    def _on_open_folder(self, record: RunRecord) -> None:
        if not reveal_in_finder(record.output_dir):
            self._warn_folder_missing(record)

    def _warn_folder_missing(self, record: RunRecord) -> None:
        QMessageBox.warning(
            self, "Folder missing",
            f"The output folder for “{record.name}” is no longer at:\n"
            f"{record.output_dir}\n\nUse Locate… on its row to point "
            f"FlexiProbes at where it went.",
        )
        self._refresh_completed()

    # ---- Locate a moved run folder ----

    def _on_locate_run(self, rec: RunRecord, folder: str = "") -> int:
        """Relink a run whose folder moved; returns how many runs were relinked.

        The folder is accepted only if it matches the record (see
        ``run_import.match_folder``). Afterwards any other missing runs found
        in the same new place are offered for repair in one go, since folders
        are usually moved together.
        """
        old_dir = rec.output_dir
        if not folder:
            folder = QFileDialog.getExistingDirectory(
                self, f"Where is the folder for “{rec.name}” now?",
                nearest_existing_parent(old_dir),
            )
            if not folder:
                return 0
        match = match_folder(rec, folder)
        if match.verdict == REJECT:
            QMessageBox.warning(self, "Not this run's folder", match.reason)
            return 0
        if match.verdict == CONFIRM:
            box = QMessageBox(self)
            box.setIcon(QMessageBox.Question)
            box.setWindowTitle("Confirm folder")
            box.setText(f"{match.reason}\n\nUse this folder for “{rec.name}”?")
            box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            box.setDefaultButton(QMessageBox.No)
            if box.exec() != QMessageBox.Yes:
                return 0
        relink(self.manager.history, rec, folder)
        relinked = 1

        siblings = sibling_candidates(self.manager.history, old_dir, folder,
                                      rec.run_id)
        if siblings:
            lines = [
                f"{r.name}  ->  {path}"
                + ("   (matched by folder name only)"
                   if m.verdict == CONFIRM else "")
                for r, path, m in siblings
            ]
            box = QMessageBox(self)
            box.setIcon(QMessageBox.Question)
            box.setWindowTitle("Relink other runs")
            box.setText(
                f"{len(siblings)} other run(s) with a missing folder were found "
                f"in the same new place. Relink them too?")
            box.setDetailedText("\n".join(lines))
            box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            box.setDefaultButton(QMessageBox.Yes)
            if box.exec() == QMessageBox.Yes:
                for r, path, _m in siblings:
                    relink(self.manager.history, r, path)
                    relinked += 1
        self._refresh()
        return relinked

    # ---- Re-run / delete completed runs ----

    def _on_resume(self, rec: RunRecord) -> None:
        if not self.manager.resume(rec.run_id):
            QMessageBox.warning(
                self, "Can't resume",
                "This run can no longer be resumed — its saved progress or "
                "configuration is missing from the output folder. Use “Re-run…” "
                "to start it again.",
            )
            return
        self._refresh()

    def _on_rerun(self, rec: RunRecord) -> None:
        draft = rec.query_draft or self._record_to_draft(rec)
        dlg = NewQueryDialog(
            self, projects=self.projects,
            settings=self.settings, drafts=self.drafts,
        )
        dlg.apply_draft(draft)
        dlg.setWindowTitle(f"Re-run: {rec.name}")
        if dlg.exec() != NewQueryDialog.Accepted:
            return
        runs = dlg.build_runs()
        self._submit_built_runs(runs)
        self._warn_if_blast_gate_disabled(dlg)
        self._refresh_project_filter()
        if len(runs) > 1:
            QMessageBox.information(
                self, "Bulk runs queued",
                f"Queued {len(runs)} runs. They run one at a time; watch "
                f"the Active runs section for progress.",
            )

    def _warn_if_blast_gate_disabled(self, dlg) -> None:
        """Tell the user when picking BLAST silently weakened their screen.

        The junction-duplex gate cannot be measured by BLAST, and RunConfig
        .validate() refuses the combination outright. The dialog has an engine
        picker but no junction-duplex control, so ``build_runs`` turns the gate
        off for the BLAST path to keep the run buildable. That is a real change
        to how strictly off-targets are screened, so it is never left silent.
        """
        note = getattr(dlg, "blast_gate_note", "")
        if note:
            QMessageBox.warning(self, "Junction-duplex filter disabled", note)

    def _record_to_draft(self, rec: RunRecord) -> dict:
        """Best-effort draft for runs created before re-run support.

        The design settings come from ``config_summary`` via the ONE shared map
        in ``_nq_build`` — this method used to keep its own hand-written subset
        of them and had fallen eight keys behind what build_runs records, so
        "Re-run" rebuilt a different design from the one on screen (different
        off-target engine, no Tm floor, homology mode silently dropped) without
        saying so. Every default in that map is read from the engine dataclasses,
        so the two cannot drift apart again.

        Everything below is record-level rather than design-level: it lives on
        the RunRecord itself, or (BLAST DBs) is not recorded at all, so the user
        re-adds it — and can apply the BLAST targeting guidance while doing so.
        """
        draft = summary_to_draft(rec.config_summary or {})
        draft.update({
            "modality": rec.modality_label,
            "project_id": rec.project_id,
            "run_name": rec.name,
            "fasta_path": rec.fasta_path,
            "output_dir": self.settings.default_output_dir(),
            "overwrite": True,
            "blast_dbs": "",
            "cached_hits": "",
        })
        return draft

    def _on_delete_run(self, rec: RunRecord) -> None:
        from PySide6.QtWidgets import QCheckBox
        box = QMessageBox(self)
        box.setWindowTitle("Delete run")
        box.setIcon(QMessageBox.Question)
        box.setText(f"Remove run '{rec.name}' from the list?")
        box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        box.setDefaultButton(QMessageBox.No)
        chk = None
        if rec.output_dir and os.path.isdir(rec.output_dir):
            chk = QCheckBox("Also delete its output folder on disk")
            box.setCheckBox(chk)
        if box.exec() != QMessageBox.Yes:
            return
        if chk is not None and chk.isChecked():
            import shutil
            try:
                shutil.rmtree(rec.output_dir)
            except OSError as exc:
                QMessageBox.warning(
                    self, "Delete run",
                    f"Removed from the list, but could not delete the folder:\n{exc}",
                )
        self.manager.history.delete(rec.run_id)
        self._refresh_completed()

    # ---- Bulk run selection (move / delete checked rows) ----

    def _checked_run_ids(self) -> list:
        """Run-ids of every checked row on the current page."""
        ids = []
        for row in range(self.table.rowCount()):
            it = self.table.item(row, 0)
            if it is not None and it.checkState() == Qt.Checked:
                rid = it.data(Qt.UserRole)
                if rid:
                    ids.append(rid)
        return ids

    def _on_run_item_changed(self, item: QTableWidgetItem) -> None:
        # Only the checkbox column (0) affects the selection state.
        if self._populating_runs or item.column() != 0:
            return
        self._sync_run_bulk_buttons()

    def _on_select_all_runs(self, checked: bool) -> None:
        self._populating_runs = True
        for row in range(self.table.rowCount()):
            it = self.table.item(row, 0)
            if it is not None:
                it.setCheckState(Qt.Checked if checked else Qt.Unchecked)
        self._populating_runs = False
        self._sync_run_bulk_buttons()

    def _sync_run_bulk_buttons(self) -> None:
        n_rows = self.table.rowCount()
        n = len(self._checked_run_ids())
        self.bulk_move_btn.setEnabled(n > 0)
        self.bulk_delete_btn.setEnabled(n > 0)
        self.bulk_count_lbl.setText(f"{n} selected" if n else "")
        # Keep the header checkbox in sync without re-triggering its slot.
        self.select_all_chk.blockSignals(True)
        self.select_all_chk.setChecked(n_rows > 0 and n == n_rows)
        self.select_all_chk.blockSignals(False)

    def _on_bulk_move_runs(self) -> None:
        run_ids = self._checked_run_ids()
        if not run_ids:
            return
        projects = self.projects.list()
        labels = ["(No project)"] + [p.name for p in projects]
        choice, ok = QInputDialog.getItem(
            self, "Move runs to project",
            f"Move {len(run_ids)} run(s) to:", labels, 0, False,
        )
        if not ok:
            return
        target_id = "" if choice == "(No project)" else next(
            (p.project_id for p in projects if p.name == choice), "",
        )
        self.manager.history.set_project_many(run_ids, target_id)
        self._refresh_completed()

    def _on_bulk_delete_runs(self) -> None:
        run_ids = self._checked_run_ids()
        if not run_ids:
            return
        recs = [r for r in (self.manager.history.get(i) for i in run_ids)
                if r is not None]
        if not recs:
            return
        box = QMessageBox(self)
        box.setWindowTitle("Delete runs")
        box.setIcon(QMessageBox.Question)
        box.setText(f"Remove {len(recs)} selected run(s) from the list?")
        box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        box.setDefaultButton(QMessageBox.No)
        have_dirs = any(r.output_dir and os.path.isdir(r.output_dir) for r in recs)
        chk = QCheckBox("Also delete their output folders on disk") if have_dirs else None
        if chk is not None:
            box.setCheckBox(chk)
        if box.exec() != QMessageBox.Yes:
            return
        if chk is not None and chk.isChecked():
            import shutil
            failed = []
            for r in recs:
                if r.output_dir and os.path.isdir(r.output_dir):
                    try:
                        shutil.rmtree(r.output_dir)
                    except OSError as exc:
                        failed.append(f"{r.name}: {exc}")
            if failed:
                QMessageBox.warning(
                    self, "Delete runs",
                    "Removed from the list, but some folders could not be "
                    "deleted:\n" + "\n".join(failed),
                )
        self.manager.history.delete_many([r.run_id for r in recs])
        self._refresh_completed()

    # ---- Project filter / bulk export ----

    def _refresh_project_filter(self) -> None:
        # Block signals while we rebuild the combo so currentIndexChanged
        # doesn't fire mid-population. The blocker must stay alive until
        # the rebuild is done — never ``del`` it.
        self.project_filter_combo.blockSignals(True)
        try:
            current = self._project_filter
            self.project_filter_combo.clear()
            self.project_filter_combo.addItem("All projects", userData="")
            for proj in self.projects.list():
                self.project_filter_combo.addItem(proj.name, userData=proj.project_id)
            idx = self.project_filter_combo.findData(current)
            if idx >= 0:
                self.project_filter_combo.setCurrentIndex(idx)
        finally:
            self.project_filter_combo.blockSignals(False)
        self.bulk_export_btn.setEnabled(bool(self._project_filter))

    def _on_project_filter_changed(self, _idx: int) -> None:
        self._project_filter = self.project_filter_combo.currentData() or ""
        self.bulk_export_btn.setEnabled(bool(self._project_filter))
        self._refresh_completed()

    def _on_manage_projects(self) -> None:
        # Full project manager: create / rename / describe / delete, with
        # per-project run counts. Runs the dialog modally, then re-syncs the
        # home page since the active filter's project may have been renamed
        # or deleted.
        dlg = ManageProjectsDialog(
            projects=self.projects,
            run_history=self.manager.history,
            parent=self,
        )
        dlg.exec()
        # The currently-filtered project may no longer exist after a delete.
        if self._project_filter and self.projects.get(self._project_filter) is None:
            self._project_filter = ""
        self._refresh_project_filter()
        self._refresh_completed()

    def _on_bulk_export(self) -> None:
        if not self._project_filter:
            return
        proj = self.projects.get(self._project_filter)
        if proj is None:
            return
        # Gather every confirmed selection for this project's runs.
        runs_payload = []
        # A run whose folder has moved has no readable selection. It used to be
        # dropped from the export without a word, which reads as "that run had
        # no probes selected" when the truth is "that run was not looked at".
        missing = []
        for rec in self.manager.history.list_completed():
            if rec.project_id != self._project_filter:
                continue
            if is_folder_missing(rec):
                missing.append(rec.name)
                continue
            sel = load_selections(rec.output_dir)
            if not sel:
                continue
            runs_payload.append({
                "name": rec.name,
                "output_dir": rec.output_dir,
                "selections": sel,
            })
        missing_note = (
            "\n\nLeft out because their folder is missing (use Locate… on "
            "their rows): " + ", ".join(missing)) if missing else ""
        if not runs_payload:
            QMessageBox.warning(
                self, "Bulk export",
                f"No runs in '{proj.name}' have a saved selection yet. "
                f"Open each run, save or confirm a selection, then try again."
                + missing_note,
            )
            return
        bulk_dir = QFileDialog.getExistingDirectory(
            self, f"Folder to write '{proj.name}' bulk CSVs into",
        )
        if not bulk_dir:
            return
        try:
            result = export_project_bulk(
                runs_payload, bulk_dir,
                project_name="".join(
                    c if c.isalnum() or c in "_-." else "_" for c in proj.name
                ),
            )
        except Exception as exc:
            QMessageBox.critical(self, "Bulk export failed", str(exc))
            return
        n = sum(result.get("summary", {}).values())
        msg = (
            f"Wrote {n} order rows (probe pairs × barcodes) across "
            f"{len(result['summary'])} runs:\n"
            + "\n".join(result["paths"].values())
        )
        if result.get("skipped"):
            msg += ("\n\nSkipped (no selections, or no results.csv): "
                    + ", ".join(result["skipped"]))
        msg += missing_note
        QMessageBox.information(self, "Bulk export", msg)

    def _on_review_confirmed(self, run_id: str) -> None:
        rec = self.manager.history.get(run_id)
        if rec is None:
            return
        # Once the user confirms a selection the run is done — flip status
        # to "complete" and persist. Re-running review is still allowed
        # because the Review button stays visible for discovery-mode runs.
        rec.status = "complete"
        rec.stage = "exported"
        rec.completed_at = now_iso()
        self.manager.history.update(rec)
        self._refresh()

    # ---- Refresh ----

    def _refresh(self) -> None:
        self._refresh_active()
        self._refresh_completed()

    def _refresh_active(self) -> None:
        active = self.manager.history.list_active()
        active_ids = {r.run_id for r in active}

        # Remove cards no longer active
        for rid in list(self._cards.keys()):
            if rid not in active_ids:
                card = self._cards.pop(rid)
                self.active_layout.removeWidget(card)
                card.deleteLater()

        # Add / update
        for rec in active:
            if rec.run_id in self._cards:
                self._cards[rec.run_id].update_from(rec)
            else:
                card = _RunCard(rec, on_cancel=self._on_cancel)
                self._cards[rec.run_id] = card
                self.active_layout.insertWidget(
                    self.active_layout.count() - 1, card,
                )
        self.active_empty.setVisible(not active)

    # Completed-runs table (refresh/paginate) lives in CompletedTableMixin
    # — see _home_completed.py.
