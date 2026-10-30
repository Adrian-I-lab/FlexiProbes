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

"""Drag-and-drop FASTA intake for HomeWindow.

Split out of home_window.py (audit #19). ``DragDropMixin`` handles dropped
files/folders and validates FASTAs; it uses only the host window's own state
via ``self`` (e.g. FASTA_EXTS, _open_new_query_dialog) and holds nothing."""
from __future__ import annotations

import os
from typing import TYPE_CHECKING

from PySide6.QtWidgets import QMessageBox

from .new_query_dialog import NewQueryDialog

if TYPE_CHECKING:  # give the mixin the host widget's type for mypy
    from PySide6.QtWidgets import QMainWindow as _HostBase
else:
    _HostBase = object


class DragDropMixin(_HostBase):
    def dragEnterEvent(self, event):  # noqa: N802 (Qt API)
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dragMoveEvent(self, event):   # noqa: N802 (Qt API)
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):       # noqa: N802 (Qt API)
        if not event.mimeData().hasUrls():
            return
        paths = [
            u.toLocalFile() for u in event.mimeData().urls() if u.isLocalFile()
        ]
        fastas = self._gather_fastas(paths)
        if not fastas:
            QMessageBox.warning(
                self, "No FASTA files",
                "Dropped items didn't contain any "
                f"{', '.join(self.FASTA_EXTS)} files.",
            )
            return
        # Fail-fast validation so corrupt/empty files surface now, not three
        # pipeline stages later (U4).
        fastas, bad = self._filter_valid_fastas(fastas)
        if bad:
            detail = "\n".join(f"  • {os.path.basename(p)} — {r}" for p, r in bad)
            if not fastas:
                QMessageBox.warning(
                    self, "Unreadable FASTA(s)",
                    "None of the dropped files look like valid FASTA:\n\n" + detail,
                )
                return
            QMessageBox.warning(
                self, "Some FASTAs skipped",
                f"Skipping {len(bad)} file(s) that don't look like valid "
                f"FASTA:\n\n{detail}\n\nContinuing with {len(fastas)} valid file(s).",
            )
        event.acceptProposedAction()
        # Single file → use the normal single-FASTA dialog (better UX).
        # Multiple files → bulk dialog.
        if len(fastas) == 1:
            dlg = NewQueryDialog(
                self, projects=self.projects,
                settings=self.settings, drafts=self.drafts,
            )
            dlg.fasta_edit.setText(fastas[0])
            if self._project_filter:
                idx = dlg.project_combo.findData(self._project_filter)
                if idx >= 0:
                    dlg.project_combo.setCurrentIndex(idx)
            if dlg.exec() == NewQueryDialog.Accepted:
                # Use the shared submit helper so cached-hit loading + its
                # error reporting stay identical to New query / Re-run (audit M4).
                self._submit_built_runs(dlg.build_runs())
                self._refresh_project_filter()
        else:
            self._open_new_query_dialog(bulk_paths=fastas)

    def _gather_fastas(self, paths):
        """Expand a mix of files and directories into a sorted FASTA list."""
        exts = self.FASTA_EXTS
        seen: list = []
        for p in paths:
            if not p:
                continue
            if os.path.isdir(p):
                for root, _dirs, files in os.walk(p):
                    for f in files:
                        if f.lower().endswith(exts):
                            seen.append(os.path.join(root, f))
            elif p.lower().endswith(exts):
                seen.append(p)
        # Dedupe while preserving order.
        out, marked = [], set()
        for p in sorted(seen):
            if p not in marked:
                marked.add(p)
                out.append(p)
        return out

    @staticmethod
    def _quick_fasta_check(path: str):
        """Cheap fail-fast FASTA sanity check (U4).

        Reads only the first few KB rather than parsing the whole file
        (a human-transcriptome FASTA can be hundreds of MB). Returns
        ``(ok, reason)`` — ``reason`` is empty when ok.
        """
        try:
            with open(path, "r", errors="replace") as fh:
                head = fh.read(8192)
        except OSError as exc:
            return False, f"cannot read ({exc.strerror or exc})"
        stripped = head.lstrip()
        if not stripped:
            return False, "file is empty"
        if not stripped.startswith(">"):
            return False, "does not start with a '>' FASTA header"
        # Must have at least one sequence character after the first header.
        nl = stripped.find("\n")
        if nl == -1 or not stripped[nl:].strip():
            return False, "header has no sequence data"
        return True, ""

    def _filter_valid_fastas(self, fastas):
        """Split a FASTA list into (valid, [(path, reason), ...])."""
        valid, bad = [], []
        for p in fastas:
            ok, reason = self._quick_fasta_check(p)
            (valid if ok else bad).append(p if ok else (p, reason))
        return valid, bad

