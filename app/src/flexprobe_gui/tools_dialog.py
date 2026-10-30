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

"""GUI wrappers around the CLI's housekeeping commands.

Mirrors:
  • flexiprobes build-db
  • flexiprobes validate-blast
  • flexiprobes export
  • flexiprobes version

All four open small focused dialogs; long-running work (BLAST DB builds)
happens on a worker thread so the GUI stays responsive.
"""
from __future__ import annotations

import json
import os
import traceback
from typing import Optional

from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPlainTextEdit, QPushButton,
    QVBoxLayout, QWidget,
)

from flexprobe import __version__
from flexprobe.blast import build_db, check_version, discover_blast
from flexprobe.errors import BlastError

from .theme import muted_hex


# ---------------------------------------------------------------------------
# Build BLAST DB
# ---------------------------------------------------------------------------

class _BuildDBWorker(QThread):
    finished_ok = Signal(str)   # final DB prefix
    failed = Signal(str)        # error message

    def __init__(self, fasta_path: str, db_prefix: str, title: str, parent=None):
        super().__init__(parent)
        self._fasta = fasta_path
        self._prefix = db_prefix
        self._title = title

    def run(self):
        try:
            tools = discover_blast()
            check_version(tools["blastn"])
            prefix = build_db(
                tools["makeblastdb"], self._fasta, self._prefix,
                title=self._title or os.path.basename(self._prefix),
            )
            self.finished_ok.emit(prefix)
        except BlastError as exc:
            self.failed.emit(str(exc))
        except Exception:
            self.failed.emit(traceback.format_exc())


class BuildDBDialog(QDialog):
    """`flexiprobes build-db` — build a BLAST DB from a FASTA."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Build BLAST database")
        self.setModal(True)
        self.resize(560, 240)

        outer = QVBoxLayout(self)
        form = QFormLayout()
        outer.addLayout(form)

        self.fasta_edit = QLineEdit()
        fh = QHBoxLayout(); fh.setContentsMargins(0, 0, 0, 0)
        fh.addWidget(self.fasta_edit)
        b1 = QPushButton("Browse…"); b1.clicked.connect(self._browse_fasta)
        fh.addWidget(b1)
        form.addRow("Source FASTA:", _wrap(fh))

        self.prefix_edit = QLineEdit()
        self.prefix_edit.setPlaceholderText(
            "Output DB prefix, e.g. /Users/me/blast_dbs/human_tx"
        )
        ph = QHBoxLayout(); ph.setContentsMargins(0, 0, 0, 0)
        ph.addWidget(self.prefix_edit)
        b2 = QPushButton("Browse…"); b2.clicked.connect(self._browse_prefix)
        ph.addWidget(b2)
        form.addRow("DB prefix:", _wrap(ph))

        self.title_edit = QLineEdit()
        self.title_edit.setPlaceholderText("(optional; defaults to prefix basename)")
        form.addRow("DB title:", self.title_edit)

        self.status = QLabel(
            "BLAST+ runs makeblastdb -parse_seqids under the hood. SHA-256 of "
            "the FASTA is cached next to the DB so re-builds are skipped "
            "automatically.\n\n"
            "-parse_seqids is required, not preferred: the off-target screen "
            "reads subject bases through blastdbcmd. Duplicate sequence ids "
            "make makeblastdb refuse it, and the build then fails with the fix "
            "rather than quietly producing a DB that cannot be screened against."
            "\n\nThe desktop app bundles blastdbcmd. When running from source "
            "without it, the BLAST engine falls back to a conservative count "
            "that can only over-reject, so install NCBI BLAST+ for exact "
            "counts under BLAST. "
            "The default k-mer engine needs none of this."
        )
        self.status.setWordWrap(True)
        self.status.setStyleSheet(f"color: {muted_hex()};")
        outer.addWidget(self.status)

        btns = QDialogButtonBox(QDialogButtonBox.Cancel)
        self.build_btn = btns.addButton("Build", QDialogButtonBox.AcceptRole)
        btns.accepted.connect(self._on_build)
        btns.rejected.connect(self.reject)
        outer.addWidget(btns)

        # Open at least as tall as the contents need. resize() above is only
        # a starting size, and word-wrapped help text is allowed to squash
        # below its natural height: this dialog's text overlapped or went
        # undrawn at the fixed size (2026-09-22 screenshot pass).
        self.resize(self.width(), max(self.height(), self.sizeHint().height()))

        self._worker: Optional[_BuildDBWorker] = None

    def _browse_fasta(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select FASTA", "",
            "FASTA files (*.fa *.fasta *.fna *.ffn *.frn);;All files (*)",
        )
        if path:
            self.fasta_edit.setText(path)
            if not self.prefix_edit.text():
                # Sensible default: strip the FASTA extension.
                base = path
                for ext in (".fasta", ".fa", ".fna", ".ffn", ".frn"):
                    if base.lower().endswith(ext):
                        base = base[: -len(ext)]
                        break
                self.prefix_edit.setText(base)

    def _browse_prefix(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "DB prefix (any path; sidecars get .nhr/.nin/.nsq)",
            self.prefix_edit.text() or "",
        )
        if path:
            self.prefix_edit.setText(path)

    def _on_build(self):
        fasta = self.fasta_edit.text().strip()
        prefix = self.prefix_edit.text().strip()
        if not fasta or not os.path.exists(fasta):
            QMessageBox.warning(self, "Missing FASTA", "Pick an existing FASTA file.")
            return
        if not prefix:
            QMessageBox.warning(self, "Missing prefix", "Choose an output DB prefix.")
            return
        self.build_btn.setEnabled(False)
        self.status.setText(f"Building DB at {prefix} (this can take a while)…")
        self._worker = _BuildDBWorker(fasta, prefix, self.title_edit.text().strip())
        self._worker.finished_ok.connect(self._on_finished_ok)
        self._worker.failed.connect(self._on_failed)
        self._worker.finished.connect(self._worker.deleteLater)
        self._worker.start()

    def _on_finished_ok(self, prefix: str):
        QMessageBox.information(
            self, "DB built", f"DB ready at:\n{prefix}\n\n"
            "Paste this prefix into 'Off-target references' in New query, "
            "and set Engine to BLAST. The default k-mer engine does not need "
            "a database — point it at the reference FASTA instead.",
        )
        self.accept()

    def _on_failed(self, msg: str):
        self.build_btn.setEnabled(True)
        QMessageBox.critical(self, "Build failed", msg)
        self.status.setText("Build failed — see the message above.")


# ---------------------------------------------------------------------------
# Validate BLAST install + DB
# ---------------------------------------------------------------------------

class ValidateBlastDialog(QDialog):
    """`flexiprobes validate-blast` — check BLAST+ install and (optionally) a DB."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Validate BLAST install")
        self.setModal(True)
        self.resize(560, 360)

        outer = QVBoxLayout(self)
        form = QFormLayout()
        outer.addLayout(form)

        self.db_edit = QLineEdit()
        self.db_edit.setPlaceholderText("(optional)")
        bh = QHBoxLayout(); bh.setContentsMargins(0, 0, 0, 0)
        bh.addWidget(self.db_edit)
        b = QPushButton("Browse…"); b.clicked.connect(self._browse_db)
        bh.addWidget(b)
        form.addRow("Inspect DB prefix:", _wrap(bh))

        self.output = QPlainTextEdit()
        self.output.setReadOnly(True)
        self.output.setStyleSheet(
            "font-family: monospace; font-size: 11px; "
            "background: rgba(0,0,0,0.04);"
        )
        outer.addWidget(self.output, 1)

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        run = btns.addButton("Run check", QDialogButtonBox.AcceptRole)
        run.clicked.connect(self._on_run)
        btns.rejected.connect(self.reject)
        outer.addWidget(btns)

    def _browse_db(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Pick any file from BLAST DB", "",
            "BLAST DB files (*.nhr *.nin *.nsq);;All files (*)",
        )
        if path:
            base = os.path.splitext(path)[0]
            self.db_edit.setText(base)

    def _on_run(self):
        lines = []
        try:
            tools = discover_blast()
            lines.append(f"blastn:      {tools['blastn']}")
            lines.append(f"makeblastdb: {tools['makeblastdb']}")
            ver = check_version(tools["blastn"])
            lines.append(f"version:     {ver}")
        except BlastError as exc:
            self.output.setPlainText(f"ERROR: {exc}")
            return
        db = self.db_edit.text().strip()
        if db:
            lines.append("")
            lines.append(f"Checking DB prefix: {db}")
            sidecars = {ext: os.path.exists(db + ext) for ext in (".nhr", ".nin", ".nsq")}
            for ext, present in sidecars.items():
                lines.append(f"  {ext}: {'OK' if present else 'MISSING'}")
            if all(sidecars.values()):
                lines.append("DB looks healthy.")
            else:
                lines.append(
                    "Sidecar(s) missing — re-build with Tools → Build BLAST database…"
                )
        self.output.setPlainText("\n".join(lines))


# ---------------------------------------------------------------------------
# Re-export from results.csv
# ---------------------------------------------------------------------------

class ReExportDialog(QDialog):
    """`flexiprobes export` — convert results.csv to csv/tsv/json/bed."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Re-export results.csv")
        self.setModal(True)
        self.resize(540, 220)

        outer = QVBoxLayout(self)
        form = QFormLayout()
        outer.addLayout(form)

        self.src_edit = QLineEdit()
        sh = QHBoxLayout(); sh.setContentsMargins(0, 0, 0, 0)
        sh.addWidget(self.src_edit)
        b1 = QPushButton("Browse…"); b1.clicked.connect(self._browse_src)
        sh.addWidget(b1)
        form.addRow("results.csv:", _wrap(sh))

        self.fmt_combo = QComboBox()
        self.fmt_combo.addItems(["csv", "tsv", "json", "bed"])
        form.addRow("Format:", self.fmt_combo)

        self.out_edit = QLineEdit()
        oh = QHBoxLayout(); oh.setContentsMargins(0, 0, 0, 0)
        oh.addWidget(self.out_edit)
        b2 = QPushButton("Browse…"); b2.clicked.connect(self._browse_out)
        oh.addWidget(b2)
        form.addRow("Output file:", _wrap(oh))

        btns = QDialogButtonBox(QDialogButtonBox.Cancel)
        do_btn = btns.addButton("Export", QDialogButtonBox.AcceptRole)
        do_btn.clicked.connect(self._on_export)
        btns.rejected.connect(self.reject)
        outer.addWidget(btns)

    def _browse_src(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Pick results.csv", "", "CSV (*.csv);;All files (*)",
        )
        if path:
            self.src_edit.setText(path)

    def _browse_out(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Output file", self.out_edit.text() or "",
        )
        if path:
            self.out_edit.setText(path)

    def _on_export(self):
        import csv
        src = self.src_edit.text().strip()
        out = self.out_edit.text().strip()
        fmt = self.fmt_combo.currentText()
        if not src or not os.path.exists(src):
            QMessageBox.warning(self, "Missing results.csv", "Pick a valid file.")
            return
        if not out:
            QMessageBox.warning(self, "Missing output", "Choose an output path.")
            return
        try:
            if fmt in ("csv", "tsv"):
                sep = "," if fmt == "csv" else "\t"
                with open(src) as fin, open(out, "w", newline="") as fout:
                    w = csv.writer(fout, delimiter=sep)
                    for row in csv.reader(fin):
                        w.writerow(row)
            elif fmt == "json":
                with open(src) as fin:
                    rows = list(csv.DictReader(fin))
                with open(out, "w") as fout:
                    json.dump(rows, fout, indent=2)
            elif fmt == "bed":
                with open(src) as fin, open(out, "w") as fout:
                    for row in csv.DictReader(fin):
                        fout.write(
                            f"{row.get('gene_id','')}\t{row.get('start',0)}\t"
                            f"{row.get('end',0)}\t{row.get('probe_id','')}\t"
                            f"{row.get('composite_score',0)}\n"
                        )
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))
            return
        QMessageBox.information(self, "Exported", f"Wrote {out}")
        self.accept()


# ---------------------------------------------------------------------------
# Preferences
# ---------------------------------------------------------------------------

class PreferencesDialog(QDialog):
    """Edit persisted app preferences (currently the default output folder)."""

    def __init__(self, settings, parent=None):
        super().__init__(parent)
        self._settings = settings
        self.setWindowTitle("Preferences")
        self.setModal(True)
        self.resize(560, 180)

        outer = QVBoxLayout(self)
        form = QFormLayout()
        outer.addLayout(form)

        self.out_edit = QLineEdit(settings.default_output_dir())
        oh = QHBoxLayout(); oh.setContentsMargins(0, 0, 0, 0)
        oh.addWidget(self.out_edit)
        b = QPushButton("Browse…"); b.clicked.connect(self._browse_out)
        oh.addWidget(b)
        form.addRow("Default output folder:", _wrap(oh))

        # Read at open, so a folder chosen in New query since the last visit
        # shows here: New query saves its choice through the same setting.
        hint = QLabel(
            "New queries pre-fill this folder. Choosing a different output "
            "folder in New query also makes it the new default. Each run still "
            "gets its own timestamped sub-folder inside it (grouped by project "
            "when set)."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet(f"color: {muted_hex()};")
        outer.addWidget(hint)

        btns = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        btns.accepted.connect(self._on_save)
        btns.rejected.connect(self.reject)
        outer.addWidget(btns)

    def _browse_out(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select default output folder", self.out_edit.text() or "",
        )
        if path:
            self.out_edit.setText(path)

    def _on_save(self):
        path = self.out_edit.text().strip()
        if not path:
            QMessageBox.warning(self, "Preferences", "Choose a default output folder.")
            return
        self._settings.set_default_output_dir(path)
        self.accept()


# ---------------------------------------------------------------------------
# About / version
# ---------------------------------------------------------------------------

def show_about(parent=None):
    try:
        tools = discover_blast()
        blast_info = (
            f"blastn: {tools['blastn']}\n"
            f"makeblastdb: {tools['makeblastdb']}\n"
            f"blastn version: {check_version(tools['blastn'])}"
        )
    except BlastError as exc:
        blast_info = f"BLAST+: unavailable ({exc})"
    # The 10x / IDT disclaimer appeared only in NOTICE, which a GUI-only user
    # never opens. Trademark nominative fair use and the "no patent rights"
    # position both depend on the disclaimer being visible to the people
    # actually using the tool, so it belongs here too. Kept short on purpose;
    # NOTICE has the full text.
    QMessageBox.about(
        parent, "About FlexiProbes",
        f"FlexiProbes {__version__}\n\n"
        f"{blast_info}\n\n"
        "Licensed GPL-3.0-or-later. Bundles NCBI BLAST+ (public domain), "
        "Qt/PySide6 (GPL-3.0 option), primer3-py (GPL-2.0-or-later), MAFFT "
        "(BSD), Biopython (BSD-style) and NumPy (BSD-3-Clause); full "
        "licence texts ship with the application under licenses/.\n\n"
        "INDEPENDENT, UNOFFICIAL TOOL. Not affiliated with, endorsed by, or "
        "supported by 10x Genomics, Inc. or Integrated DNA Technologies, Inc. "
        "\"10x Genomics\", \"Flex\", \"Visium\", \"Chromium\", \"IDT\" and "
        "\"oPools\" are trademarks of their respective owners, used here "
        "descriptively to indicate compatibility. No patent rights are "
        "granted; the Flex probe-hybridisation chemistry may be patented, and "
        "obtaining any licence required for your use is your responsibility.\n\n"
        "Research software — NOT a medical device, and not validated for "
        "diagnostic or clinical use. Verify chemistry constants against the "
        "current 10x technical note, and order files against IDT's current "
        "requirements, before ordering reagents.\n\n"
        "See the bundled NOTICE file for the full text.",
    )


def _wrap(layout):
    w = QWidget()
    w.setLayout(layout)
    return w
