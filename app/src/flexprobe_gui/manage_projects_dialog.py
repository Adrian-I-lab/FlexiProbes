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

"""Manage-projects dialog.

A small project manager opened from the home page. Lists every project with
its description, creation date and the number of runs attached to it, and
lets the user create, rename, re-describe and delete projects from one place
instead of the old single-shot add-or-delete prompt.

Each row carries a checkbox so several projects can be deleted in one go; a
"Select all" box toggles the whole list. Rename / Edit-description act on a
single checked project (they need exactly one target); Delete acts on every
checked project at once.

Deleting a project never touches the runs themselves — they stay on disk and
simply lose their project filter, matching the long-standing behaviour.
"""
from __future__ import annotations

from collections import Counter
from typing import List, Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QHBoxLayout, QHeaderView,
    QInputDialog, QLabel, QLineEdit, QMessageBox, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout,
)

from .project_history import ProjectHistory
from .run_history import RunHistory


class ManageProjectsDialog(QDialog):
    """Create / rename / describe / delete projects in one window."""

    def __init__(
        self,
        projects: ProjectHistory,
        run_history: Optional[RunHistory] = None,
        parent=None,
    ):
        super().__init__(parent)
        self._projects = projects
        self._run_history = run_history
        # Guards the ``itemChanged`` handler while we rebuild the table so
        # programmatic check-state changes don't re-enter the selection sync.
        self._populating = False
        self.setWindowTitle("Manage projects")
        self.setModal(True)
        self.resize(720, 460)
        self.setMinimumSize(520, 320)
        self.setSizeGripEnabled(True)

        outer = QVBoxLayout(self)
        outer.addWidget(QLabel(
            "Projects group runs together for the home-page filter and bulk "
            "export. Check projects to delete several at once. Deleting a "
            "project leaves its runs on disk; they just stop being filterable "
            "by project."
        ))

        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["", "Name", "Description", "Created", "Runs"]
        )
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.table.doubleClicked.connect(self._on_row_double_clicked)
        self.table.itemChanged.connect(self._on_item_changed)
        outer.addWidget(self.table, 1)

        # --- Action buttons ---
        row = QHBoxLayout()
        self.select_all_chk = QCheckBox("Select all")
        self.select_all_chk.setToolTip("Check or clear every project.")
        self.select_all_chk.clicked.connect(self._on_select_all)
        row.addWidget(self.select_all_chk)

        self.new_btn = QPushButton("New…")
        self.new_btn.clicked.connect(self._on_new)
        row.addWidget(self.new_btn)

        self.rename_btn = QPushButton("Rename…")
        self.rename_btn.clicked.connect(self._on_rename)
        row.addWidget(self.rename_btn)

        self.desc_btn = QPushButton("Edit description…")
        self.desc_btn.clicked.connect(self._on_edit_description)
        row.addWidget(self.desc_btn)

        self.delete_btn = QPushButton("Delete selected")
        self.delete_btn.clicked.connect(self._on_delete)
        row.addWidget(self.delete_btn)

        row.addStretch(1)
        outer.addLayout(row)

        bb = QDialogButtonBox(QDialogButtonBox.Close)
        bb.rejected.connect(self.reject)
        bb.accepted.connect(self.accept)
        outer.addWidget(bb)

        self._reload()

    # ---- Helpers ----

    def _run_counts(self) -> Counter:
        """How many runs (active + completed) reference each project_id."""
        counts: Counter = Counter()
        if self._run_history is None:
            return counts
        for rec in self._run_history.list_active() + self._run_history.list_completed():
            if rec.project_id:
                counts[rec.project_id] += 1
        return counts

    def _reload(self) -> None:
        counts = self._run_counts()
        projects = self._projects.list()
        self._populating = True
        self.table.setRowCount(0)
        for proj in projects:
            i = self.table.rowCount()
            self.table.insertRow(i)
            chk = QTableWidgetItem()
            chk.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            chk.setCheckState(Qt.Unchecked)
            # Stash the id on the checkbox so actions don't depend on sort order.
            chk.setData(Qt.UserRole, proj.project_id)
            self.table.setItem(i, 0, chk)
            self.table.setItem(i, 1, QTableWidgetItem(proj.name))
            self.table.setItem(i, 2, QTableWidgetItem(proj.description))
            self.table.setItem(i, 3, QTableWidgetItem(proj.created_at))
            self.table.setItem(i, 4, QTableWidgetItem(str(counts.get(proj.project_id, 0))))
        self._populating = False
        self.table.resizeColumnToContents(0)
        self.table.resizeColumnToContents(1)
        self.table.resizeColumnToContents(3)
        self.table.resizeColumnToContents(4)
        self._sync_buttons()

    def _checked_ids(self) -> List[str]:
        ids: List[str] = []
        for row in range(self.table.rowCount()):
            it = self.table.item(row, 0)
            if it is not None and it.checkState() == Qt.Checked:
                pid = it.data(Qt.UserRole)
                if pid:
                    ids.append(pid)
        return ids

    def _selected_id(self) -> Optional[str]:
        """The single checked project, or None if zero / more than one."""
        ids = self._checked_ids()
        return ids[0] if len(ids) == 1 else None

    def _pid_at_row(self, row: int) -> Optional[str]:
        it = self.table.item(row, 0)
        return it.data(Qt.UserRole) if it is not None else None

    def _on_item_changed(self, _item: QTableWidgetItem) -> None:
        if self._populating:
            return
        self._sync_buttons()

    def _on_select_all(self, checked: bool) -> None:
        self._populating = True
        for row in range(self.table.rowCount()):
            it = self.table.item(row, 0)
            if it is not None:
                it.setCheckState(Qt.Checked if checked else Qt.Unchecked)
        self._populating = False
        self._sync_buttons()

    def _sync_buttons(self) -> None:
        n_rows = self.table.rowCount()
        n_checked = len(self._checked_ids())
        single = n_checked == 1
        # Rename / describe need exactly one target; delete takes any number.
        self.rename_btn.setEnabled(single)
        self.desc_btn.setEnabled(single)
        self.delete_btn.setEnabled(n_checked > 0)
        self.select_all_chk.blockSignals(True)
        self.select_all_chk.setChecked(n_rows > 0 and n_checked == n_rows)
        self.select_all_chk.blockSignals(False)

    # ---- Actions ----

    def _on_new(self) -> None:
        name, ok = QInputDialog.getText(self, "New project", "Project name:")
        if not ok or not name.strip():
            return
        try:
            self._projects.create(name.strip())
        except Exception as exc:  # ValueError on blank, but stay defensive.
            QMessageBox.warning(self, "Project", str(exc))
            return
        self._reload()

    def _on_row_double_clicked(self, index) -> None:
        # Double-click renames the clicked row regardless of its checkbox.
        pid = self._pid_at_row(index.row())
        if pid is not None:
            self._rename_pid(pid)

    def _on_rename(self) -> None:
        pid = self._selected_id()
        if pid is not None:
            self._rename_pid(pid)

    def _rename_pid(self, pid: str) -> None:
        proj = self._projects.get(pid)
        if proj is None:
            return
        name, ok = QInputDialog.getText(
            self, "Rename project", "New name:",
            QLineEdit.Normal, proj.name,
        )
        if not ok or not name.strip() or name.strip() == proj.name:
            return
        try:
            self._projects.update(pid, name=name.strip())
        except ValueError as exc:
            QMessageBox.warning(self, "Rename project", str(exc))
            return
        self._reload()

    def _on_edit_description(self) -> None:
        pid = self._selected_id()
        if pid is None:
            return
        proj = self._projects.get(pid)
        if proj is None:
            return
        text, ok = QInputDialog.getMultiLineText(
            self, "Edit description",
            f"Description for '{proj.name}':", proj.description,
        )
        if not ok:
            return
        self._projects.update(pid, description=text)
        self._reload()

    def _on_delete(self) -> None:
        pids = self._checked_ids()
        if not pids:
            return
        counts = self._run_counts()
        names = [p.name for p in self._projects.list() if p.project_id in pids]
        n_runs = sum(counts.get(pid, 0) for pid in pids)
        if len(pids) == 1:
            head = f"Delete project '{names[0]}'?"
        else:
            head = f"Delete {len(pids)} projects?\n\n• " + "\n• ".join(names)
        detail = (
            f"\n\n{n_runs} run(s) are attached. They stay on disk; they'll just "
            f"no longer be filterable by project." if n_runs else ""
        )
        ans = QMessageBox.question(
            self, "Delete project(s)",
            f"{head}{detail}",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No,
        )
        if ans == QMessageBox.Yes:
            self._projects.delete_many(pids)
            self._reload()
