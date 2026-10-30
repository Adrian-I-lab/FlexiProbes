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

"""GUI for the HPC export / import workflow.

The whole-transcriptome method is too memory-heavy for a laptop. These dialogs
let the user run it on a cluster instead:

  • :class:`HpcExportDialog` — turn the run(s) configured in the New-query dialog
    into self-contained PBS/SLURM job bundles (config + script + staged inputs).
  • :class:`ImportHpcDialog` — bring a finished job's output (folder or archive)
    back in as a run, so it opens in the review dialog like a local one.

The heavy lifting lives in ``flexprobe.hpc_export``; these are thin forms.
"""
from __future__ import annotations

import os
import re
import shutil
from typing import Any, Dict, List, Optional, Tuple

from PySide6.QtCore import QEventLoop, Qt, QThread, Signal
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPlainTextEdit, QProgressDialog,
    QPushButton, QRadioButton, QSpinBox, QVBoxLayout, QWidget,
)

from flexprobe.hpc_export import (
    DEFAULT_CHAIN_LINKS, DEFAULT_CORES, DEFAULT_MEM_GB, DEFAULT_WALLTIME,
    HpcResources, import_hpc_output, write_hpc_bundle,
)
from flexprobe.models import RunConfig

from .run_history import RunRecord
from .settings import AppSettings
from .theme import muted_hex


def _wrap(layout) -> QWidget:
    w = QWidget()
    w.setLayout(layout)
    return w


def _slugify(name: str, fallback: str = "flexiprobes") -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", (name or "").strip()).strip("_")
    return slug or fallback


# Resource defaults come from flexprobe.hpc_export so the dialog, the CLI's
# export-hpc cannot drift apart.

# Settings keys for remembering the user's cluster between exports.
_S_SCHEDULER = "hpc_scheduler"
_S_IMAGE = "hpc_image"
_S_PARTITION = "hpc_partition"
_S_ACCOUNT = "hpc_account"
_S_WALLTIME = "hpc_walltime"
_S_MEM = "hpc_mem_gb"
_S_CORES = "hpc_cores"
_S_LINKS = "hpc_links"
_S_MODULE = "hpc_apptainer_module"
_S_STAGE = "hpc_stage_inputs"
_S_DEST = "hpc_dest_dir"


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

class HpcExportDialog(QDialog):
    """Write PBS/SLURM job bundles for one or more configured runs."""

    def __init__(
        self,
        runs: List[Tuple[RunRecord, RunConfig, Optional[str]]],
        settings: Optional[AppSettings] = None,
        parent=None,
    ):
        super().__init__(parent)
        self._runs = runs
        self._settings = settings or AppSettings()
        self.setWindowTitle("Export job bundle for HPC")
        self.setModal(True)
        self.resize(620, 560)

        s = self._settings
        outer = QVBoxLayout(self)

        intro = QLabel(
            f"Writes a self-contained job bundle for <b>{len(runs)} run(s)</b>. "
            "Copy it to your cluster, adjust the script if needed, submit it, then "
            "bring the results back with <b>Tools → Import HPC run…</b>."
        )
        intro.setWordWrap(True)
        outer.addWidget(intro)

        form = QFormLayout()
        outer.addLayout(form)

        # Scheduler
        self.slurm_radio = QRadioButton("SLURM (sbatch)")
        self.pbs_radio = QRadioButton("PBS (qsub)")
        if str(s.get(_S_SCHEDULER, "slurm")) == "pbs":
            self.pbs_radio.setChecked(True)
        else:
            self.slurm_radio.setChecked(True)
        srow = QHBoxLayout(); srow.setContentsMargins(0, 0, 0, 0)
        srow.addWidget(self.slurm_radio); srow.addWidget(self.pbs_radio)
        srow.addStretch(1)
        form.addRow("Scheduler:", _wrap(srow))

        # Container image
        self.image_edit = QLineEdit(str(s.get(_S_IMAGE, "")))
        self.image_edit.setPlaceholderText("/path/on/cluster/flexiprobe.sif")
        ih = QHBoxLayout(); ih.setContentsMargins(0, 0, 0, 0)
        ih.addWidget(self.image_edit)
        ib = QPushButton("Browse…"); ib.clicked.connect(self._browse_image)
        ih.addWidget(ib)
        form.addRow("Container image (.sif):", _wrap(ih))

        # A batch job does not inherit interactively-loaded modules, and on any
        # site where apptainer needs one this is THE way a chain dies: link 1
        # exits in about a second, writes the abort sentinel, and every
        # remaining link no-ops. The APPTAINER_MODULE escape hatch always
        # existed and is exactly what nobody remembers to type.
        self.module_edit = QLineEdit(str(s.get(_S_MODULE, "")))
        self.module_edit.setPlaceholderText(
            "(optional) e.g. apptainer/1.4.1 — leave blank if apptainer is "
            "already on a batch job's PATH"
        )
        self.module_edit.setToolTip(
            "Environment module that puts `apptainer` on PATH inside the job.\n\n"
            "A BATCH JOB DOES NOT INHERIT THE MODULES YOU LOADED "
            "INTERACTIVELY, so running `module load apptainer` in your shell "
            "before submitting does nothing for the job. If your site needs a "
            "module and this is blank, every link of the chain aborts about a "
            "second after it starts.\n\n"
            "Find the name on a login node with:  module avail apptainer\n\n"
            "Leave blank if apptainer is on the default PATH. The job script "
            "also falls back to a bare `module load apptainer`, and you can "
            "still override per submission with "
            "APPTAINER_MODULE=<name> ./submit_chain.sh."
        )
        form.addRow("Container module:", self.module_edit)

        # Resources
        self.walltime_edit = QLineEdit(
            str(s.get(_S_WALLTIME, DEFAULT_WALLTIME)))
        self.walltime_edit.setToolTip(
            "Walltime for ONE CHAIN LINK, not for the whole run. "
            "submit_chain.sh queues several copies of the job, each depending "
            "on the last, and every stage checkpoints — so a link killed on "
            "walltime resumes where it stopped instead of starting over. The "
            "budget you actually get is this number times the link count "
            "below.\n\n"
            f"{DEFAULT_WALLTIME} is deliberately short. A long walltime request "
            "queues behind everything: schedulers backfill short jobs into "
            "gaps a long one cannot fit, and — the case that costs whole weeks "
            "— a job whose walltime would still be running when a scheduled "
            "maintenance window opens is not started at all. It sits queued "
            "until the window closes.\n\n"
            "So if the run needs longer, add links. Do not lengthen this."
        )
        form.addRow("Walltime per link (HH:MM:SS):", self.walltime_edit)

        # The chain link count was hardcoded at 12 inside the generated
        # submit_chain.sh and reachable from neither the dialog nor the CLI —
        # the user had to edit the script or remember to pass an argument at
        # submit time. It is half of the walltime budget, so it belongs here
        # beside the other half.
        self.links_spin = QSpinBox(); self.links_spin.setRange(1, 200)
        self.links_spin.setValue(
            int(s.get(_S_LINKS, DEFAULT_CHAIN_LINKS) or DEFAULT_CHAIN_LINKS))
        self.links_spin.setToolTip(
            "How many chained jobs ./submit_chain.sh queues. Each resumes from "
            "the last checkpoint, so together they outlast any single "
            "walltime limit.\n\n"
            "Over-provisioning is free: once the run finishes, every remaining "
            "link sees the .flexiprobes_DONE sentinel and exits in about a "
            "second. This is the knob to raise for a long run — not the "
            "walltime.\n\n"
            "You can still override it at submit time: ./submit_chain.sh 24"
        )
        form.addRow("Chain links:", self.links_spin)

        # Walltime x links is the number that decides whether the run finishes,
        # and it is the one neither field shows on its own.
        self.budget_label = QLabel()
        self.budget_label.setWordWrap(True)
        self.budget_label.setStyleSheet(f"color: {muted_hex()};")
        form.addRow("", self.budget_label)
        self.walltime_edit.textChanged.connect(self._refresh_budget)
        self.links_spin.valueChanged.connect(self._refresh_budget)
        self._refresh_budget()

        self.mem_spin = QSpinBox(); self.mem_spin.setRange(1, 4096)
        self.mem_spin.setSuffix(" GB")
        self.mem_spin.setValue(int(s.get(_S_MEM, DEFAULT_MEM_GB) or DEFAULT_MEM_GB))
        self.mem_spin.setToolTip(
            f"{DEFAULT_MEM_GB} GB suits an ordinary multi-gene design.\n\n"
            "A WHOLE TRANSCRIPTOME NEEDS FAR MORE. A human cDNA run can be "
            "killed for running out of memory at 64 GB during tile+filter, "
            "and at 128 GB during the screen. This dialog writes a single "
            "job, so for a whole transcriptome request at least 300 GB."
        )
        form.addRow("Memory:", self.mem_spin)

        self.cores_spin = QSpinBox(); self.cores_spin.setRange(1, 512)
        self.cores_spin.setValue(int(s.get(_S_CORES, DEFAULT_CORES) or DEFAULT_CORES))
        self.cores_spin.setToolTip(
            "Passed to the pipeline as BOTH --workers (tiling + filtering) and "
            "--screen-workers (the off-target screen), so the whole run "
            f"parallelises across the allocation. {DEFAULT_CORES} suits an "
            "ordinary design; a whole-transcriptome screen suits about 24. "
            "Screen workers share the reference and index in memory, "
            "so each extra worker adds little by itself, but on a single job "
            "the screen starts from the tiling stage's memory peak. Raise "
            "memory with cores for a transcriptome."
        )
        form.addRow("CPU cores:", self.cores_spin)

        self.partition_edit = QLineEdit(str(s.get(_S_PARTITION, "")))
        self.partition_edit.setPlaceholderText("(optional) SLURM partition / PBS queue")
        form.addRow("Partition / queue:", self.partition_edit)

        self.account_edit = QLineEdit(str(s.get(_S_ACCOUNT, "")))
        self.account_edit.setPlaceholderText("(optional) allocation / account")
        form.addRow("Account:", self.account_edit)

        # Input staging
        self.stage_chk = QCheckBox(
            "Copy the input FASTA into the bundle"
        )
        self.stage_chk.setChecked(bool(s.get(_S_STAGE, True)))
        self.stage_chk.setToolTip(
            "Unchecked, config.json keeps THIS computer's path to the FASTA. "
            "Edit its \"fasta_path\" to the path on the cluster before "
            "submitting, or the job will not find its input."
        )
        form.addRow("Inputs:", self.stage_chk)

        # Optional override for BLAST DB paths. Left blank (the default), the
        # export reuses the staged input FASTA as the DB when they're the same
        # transcriptome (whole-transcriptome case → self-contained bundle), and
        # leaves any genuinely-separate DB path as-is. Deliberately NOT
        # pre-filled with the local paths, which would defeat that reuse.
        self.db_edit = QPlainTextEdit()
        self.db_edit.setPlaceholderText(
            "Optional. Reference path(s) AS THEY EXIST ON THE CLUSTER, one "
            "per line.\nLeave blank to reuse the staged input FASTA "
            "(whole-transcriptome) and keep any other configured paths."
        )
        self.db_edit.setFixedHeight(70)
        configured = self._collect_db_paths()
        if configured:
            self.db_edit.setToolTip(
                "Currently configured DB(s) in the app:\n  "
                + "\n  ".join(configured)
                + "\n\nA DB that is the same file as your input FASTA is reused "
                "from inputs/ automatically. Enter cluster paths above only to "
                "override (e.g. a pre-built DB elsewhere on the cluster)."
            )
        self.db_edit.setToolTip(
            (self.db_edit.toolTip() + "\n\n" if self.db_edit.toolTip() else "")
            + "Each line is a FASTA, or a BLAST DB prefix under the BLAST "
            "engine. Cluster paths outside the bundle also need a matching "
            "--bind for their folder on the exec line of the job script.")
        form.addRow("Override off-target references:", self.db_edit)

        # Destination
        self.dest_edit = QLineEdit(str(s.get(_S_DEST, "")))
        self.dest_edit.setPlaceholderText("Folder to write the bundle(s) into")
        dh = QHBoxLayout(); dh.setContentsMargins(0, 0, 0, 0)
        dh.addWidget(self.dest_edit)
        db = QPushButton("Browse…"); db.clicked.connect(self._browse_dest)
        dh.addWidget(db)
        form.addRow("Destination folder:", _wrap(dh))

        note = QLabel(
            "Rebuild the container from app/packaging/container/ for every "
            "new FlexiProbes version (see its README): the job refuses to run "
            "on an image with a different engine version. With the input staged and used as its own off-target DB "
            "(whole-transcriptome), the bundle is self-contained — copy the whole "
            "folder to the cluster and submit. Separate cluster DBs are referenced, "
            "not copied."
        )
        note.setWordWrap(True)
        note.setStyleSheet(f"color: {muted_hex()};")
        outer.addWidget(note)

        btns = QDialogButtonBox(QDialogButtonBox.Cancel)
        self.export_btn = btns.addButton("Export bundle", QDialogButtonBox.AcceptRole)
        btns.accepted.connect(self._on_export)
        btns.rejected.connect(self.reject)
        outer.addWidget(btns)

        # Open at least as tall as the contents need. resize() above is only
        # a starting size, and word-wrapped help text is allowed to squash
        # below its natural height: this dialog's text overlapped or went
        # undrawn at the fixed size (2026-09-22 screenshot pass).
        self.resize(self.width(), max(self.height(), self.sizeHint().height()))

    # ---- helpers ----

    @staticmethod
    def _walltime_hours(walltime: str):
        """HH:MM:SS -> hours, or None if it is not a walltime we can read.

        Deliberately permissive about the hours field: PBS and Slurm both
        accept values well past 24, and "100:00:00" is a walltime, not an
        error.
        """
        parts = (walltime or "").strip().split(":")
        if not 2 <= len(parts) <= 3:
            return None
        try:
            nums = [int(x) for x in parts]
        except ValueError:
            return None
        if any(x < 0 for x in nums):
            return None
        h, m = nums[0], nums[1]
        sec = nums[2] if len(nums) == 3 else 0
        return h + m / 60.0 + sec / 3600.0

    def _refresh_budget(self) -> None:
        """Show walltime x links, the number that decides if the run finishes."""
        hours = self._walltime_hours(self.walltime_edit.text())
        links = self.links_spin.value()
        if hours is None:
            self.budget_label.setText(
                "Walltime should look like HH:MM:SS (e.g. 24:00:00)."
            )
            return
        total = hours * links
        self.budget_label.setText(
            f"Total resumable budget: {links} x {hours:g} h = "
            f"<b>{total:g} h</b> (~{total / 24:.1f} days). The run resumes "
            f"across links, so this is what has to cover it — not the walltime "
            f"on its own."
        )

    def _collect_db_paths(self) -> List[str]:
        seen: List[str] = []
        for _rec, cfg, _cached in self._runs:
            for p in cfg.screen.db_paths:
                if p and p not in seen:
                    seen.append(p)
        return seen

    def _scheduler(self) -> str:
        return "pbs" if self.pbs_radio.isChecked() else "slurm"

    def _browse_image(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select container image", self.image_edit.text() or "",
            "Container images (*.sif *.img);;All files (*)",
        )
        if path:
            self.image_edit.setText(path)

    def _browse_dest(self):
        path = QFileDialog.getExistingDirectory(
            self, "Choose destination folder", self.dest_edit.text() or "",
        )
        if path:
            self.dest_edit.setText(path)

    def _persist(self) -> None:
        s = self._settings
        s.set(_S_SCHEDULER, self._scheduler())
        s.set(_S_IMAGE, self.image_edit.text().strip())
        s.set(_S_PARTITION, self.partition_edit.text().strip())
        s.set(_S_ACCOUNT, self.account_edit.text().strip())
        s.set(_S_WALLTIME, self.walltime_edit.text().strip())
        s.set(_S_MEM, self.mem_spin.value())
        s.set(_S_CORES, self.cores_spin.value())
        s.set(_S_LINKS, self.links_spin.value())
        s.set(_S_MODULE, self.module_edit.text().strip())
        s.set(_S_STAGE, self.stage_chk.isChecked())
        s.set(_S_DEST, self.dest_edit.text().strip())

    def _on_export(self):
        image = self.image_edit.text().strip()
        dest = self.dest_edit.text().strip()
        if not image:
            QMessageBox.warning(
                self, "Container image",
                "Enter the path to the .sif container image on the cluster.",
            )
            return
        if not dest:
            QMessageBox.warning(
                self, "Destination", "Choose a folder to write the bundle into.",
            )
            return
        if not self.walltime_edit.text().strip():
            QMessageBox.warning(self, "Walltime", "Enter a walltime (HH:MM:SS).")
            return

        remote_dbs = [
            line.strip() for line in self.db_edit.toPlainText().splitlines()
            if line.strip()
        ]
        scheduler = self._scheduler()
        stage = self.stage_chk.isChecked()

        written: List[str] = []
        try:
            for i, (rec, cfg, _cached) in enumerate(self._runs):
                base = _slugify(rec.name or os.path.splitext(
                    os.path.basename(cfg.fasta_path))[0] or f"run{i + 1}")
                # Disambiguate when several runs share a name.
                bundle_dir = os.path.join(dest, base)
                if os.path.exists(bundle_dir) and len(self._runs) > 1:
                    bundle_dir = os.path.join(dest, f"{base}_{i + 1}")
                resources = HpcResources(
                    job_name=base,
                    walltime=self.walltime_edit.text().strip(),
                    mem_gb=self.mem_spin.value(),
                    cores=self.cores_spin.value(),
                    links=self.links_spin.value(),
                    partition=self.partition_edit.text().strip(),
                    account=self.account_edit.text().strip(),
                )
                write_hpc_bundle(
                    cfg,
                    scheduler=scheduler,
                    resources=resources,
                    image=image,
                    dest_dir=bundle_dir,
                    stage_inputs=stage,
                    remote_db_paths=remote_dbs or None,
                    apptainer_module=self.module_edit.text().strip(),
                )
                written.append(bundle_dir)
        except Exception as exc:
            QMessageBox.critical(
                self, "Export failed", f"Could not write the bundle:\n{exc}",
            )
            return

        self._persist()
        listing = "\n".join(f"  • {p}" for p in written)
        QMessageBox.information(
            self, "Bundle exported",
            f"Wrote {len(written)} job bundle(s):\n{listing}\n\n"
            f"Copy each to your cluster and submit the chain from inside "
            f"it:\n\n    ./submit_chain.sh\n\n"
            f"That queues {self.links_spin.value()} linked jobs of "
            f"{self.walltime_edit.text().strip()} each, resuming from "
            f"checkpoints until the run finishes; spare links exit in a "
            f"second. A single "
            f"{'sbatch' if scheduler == 'slurm' else 'qsub'} of "
            f"run_flexiprobes.{'slurm' if scheduler == 'slurm' else 'pbs'} "
            f"also works, but you have to resubmit it yourself if it hits "
            f"walltime.",
        )
        self.accept()


# ---------------------------------------------------------------------------
# Import
# ---------------------------------------------------------------------------

class _ImportWorker(QThread):
    """Runs ``import_hpc_output`` off the UI thread.

    What it does is not cheap: a cluster's ``flexiprobes_output.tar.gz`` holds
    a results.csv that reaches 1.17 GB on this project, and decompressing it
    inline froze the app for the whole extraction with no progress and no way
    out.

    Same shape as ``tools_dialog._BuildDBWorker`` — the other "run one blocking
    library call, report ok/failed" worker in this package.

    Cancellation is honest about its limits: neither ``tarfile`` nor
    ``zipfile`` can be interrupted mid-``extractall``, so Cancel returns the UI
    to the user immediately and the extraction is DISOWNED — it finishes in the
    background and its output directory is then deleted (only ever a directory
    this import itself created under the imported-runs root). Nothing is
    registered in the run history.
    """
    finished_ok = Signal(object)    # the info dict from import_hpc_output
    failed = Signal(str, bool)      # message, is_user_error

    def __init__(self, src: str, parent=None):
        super().__init__(parent)
        self._src = src
        self.disowned = False

    def run(self) -> None:
        try:
            info = import_hpc_output(self._src)
        except (FileNotFoundError, ValueError) as exc:
            self.failed.emit(str(exc), True)
            return
        except Exception as exc:
            self.failed.emit(f"Could not import the run:\n{exc}", False)
            return
        if self.disowned:
            self._discard(info.get("output_dir", ""))
            return
        self.finished_ok.emit(info)

    @staticmethod
    def _discard(output_dir: str) -> None:
        """Remove an extraction this cancelled import created, and nothing else.

        Guarded on the imported-runs root: importing a FOLDER registers it in
        place, and deleting the user's own cluster output because they clicked
        Cancel would be catastrophic.
        """
        from flexprobe.hpc_export import _imported_root
        try:
            root = os.path.realpath(_imported_root())
            target = os.path.realpath(output_dir)
        except OSError:
            return
        if not target or not os.path.isdir(target):
            return
        if os.path.commonpath([root, target]) != root or target == root:
            return
        shutil.rmtree(target, ignore_errors=True)


class ImportHpcDialog(QDialog):
    """Register a finished HPC run (folder or archive) into the run history."""

    def __init__(self, history, on_imported=None, parent=None):
        super().__init__(parent)
        self._history = history
        self._on_imported = on_imported
        # A cancelled import whose extraction is still running. Held so Qt does
        # not destroy a live QThread, and waited on in closeEvent.
        self._pending_import: Optional[_ImportWorker] = None
        self.setWindowTitle("Import HPC run")
        self.setModal(True)
        self.resize(600, 220)

        outer = QVBoxLayout(self)
        intro = QLabel(
            "Select the output a finished cluster job produced — either the "
            "<code>flexiprobes_output.tar.gz</code> archive or the unpacked "
            "<code>output/</code> folder (anything containing results.csv)."
        )
        intro.setWordWrap(True)
        outer.addWidget(intro)

        form = QFormLayout()
        outer.addLayout(form)

        self.src_edit = QLineEdit()
        self.src_edit.setPlaceholderText("Folder or .tar.gz / .zip from the cluster")
        sh = QHBoxLayout(); sh.setContentsMargins(0, 0, 0, 0)
        sh.addWidget(self.src_edit)
        bf = QPushButton("Folder…"); bf.clicked.connect(self._browse_folder)
        ba = QPushButton("Archive…"); ba.clicked.connect(self._browse_archive)
        sh.addWidget(bf); sh.addWidget(ba)
        form.addRow("Source:", _wrap(sh))

        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("(optional) name for this run")
        form.addRow("Run name:", self.name_edit)

        btns = QDialogButtonBox(QDialogButtonBox.Cancel)
        self.import_btn = btns.addButton("Import", QDialogButtonBox.AcceptRole)
        btns.accepted.connect(self._on_import)
        btns.rejected.connect(self.reject)
        outer.addWidget(btns)

    def _browse_folder(self):
        path = QFileDialog.getExistingDirectory(self, "Select output folder", "")
        if path:
            self.src_edit.setText(path)

    def _browse_archive(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select output archive", "",
            "Archives (*.tar.gz *.tgz *.zip);;All files (*)",
        )
        if path:
            self.src_edit.setText(path)

    def _on_import(self):
        src = self.src_edit.text().strip()
        if not src:
            QMessageBox.warning(self, "Source", "Choose a folder or archive.")
            return
        info = self._import_off_thread(src)
        if info is None:
            return

        from .run_import import SOURCE_HPC, find_duplicate, record_from_folder
        # The same run imported twice used to become two rows. Match on the
        # resolved folder, or on the manifest identity for a second archive of
        # the same run (each archive extracts to a fresh folder).
        dup = find_duplicate(
            self._history.list_all(), info["output_dir"],
            info.get("manifest_timestamp", ""), info.get("fasta_sha256", ""))
        if dup is not None:
            # Only ever removes an extraction this import created (guarded on
            # the imported-runs root), never a folder imported in place.
            _ImportWorker._discard(info["output_dir"])
            QMessageBox.information(
                self, "Already imported",
                f"This run is already in the list as “{dup.name}”. Nothing "
                f"was added.",
            )
            return

        name = self.name_edit.text().strip() or info.get("name") or "Imported HPC run"
        # Shared with "Import runs…": fills the FASTA path, the settings Re-run
        # needs (run_config.json), the counts (design_summary.json) and the
        # manifest identity "Locate…" checks against.
        record = record_from_folder(info["output_dir"], source=SOURCE_HPC,
                                    name=name, info=info)
        self._history.append(record)
        if self._on_imported is not None:
            self._on_imported()
        QMessageBox.information(
            self, "Run imported",
            f"Imported “{name}”. Open it from the Completed runs list to review "
            f"and select probes.",
        )
        self.accept()

    def closeEvent(self, event) -> None:
        # Never let Qt destroy a running QThread: "QThread: Destroyed while
        # thread is still running" is a hard crash, and a disowned extraction
        # can easily outlive the user's decision to close this dialog.
        worker = self._pending_import
        if worker is not None and worker.isRunning():
            worker.wait()
        super().closeEvent(event)

    def _import_off_thread(self, src: str) -> Optional[Dict[str, Any]]:
        """Extract/scan the cluster output on a worker thread.

        Returns the info dict, or None when it failed or the user cancelled
        (both already reported). Blocks on a NESTED event loop rather than on
        the extraction, so the progress dialog animates and Cancel is live.
        """
        worker = _ImportWorker(src, self)
        loop = QEventLoop(self)

        # Indeterminate on purpose: tarfile/zipfile expose no progress hook, and
        # a fake percentage on a multi-minute decompression is worse than none.
        progress = QProgressDialog(self)
        progress.setLabelText(
            "Reading the cluster output…\n"
            "(an archive is decompressed first; a 1 GB results.csv takes a "
            "while)")
        progress.setWindowTitle("Import HPC run")
        progress.setRange(0, 0)
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(400)
        progress.setAutoClose(False)
        progress.setAutoReset(False)

        result: Dict[str, Any] = {}

        def _ok(info):
            result["info"] = info
            loop.quit()

        def _failed(msg, is_user_error):
            result["error"] = (msg, is_user_error)
            loop.quit()

        def _cancel():
            # See _ImportWorker: extraction cannot be interrupted, so the thread
            # is disowned and cleans up after itself.
            worker.disowned = True
            result["cancelled"] = True
            loop.quit()

        worker.finished_ok.connect(_ok)
        worker.failed.connect(_failed)
        progress.canceled.connect(_cancel)
        worker.finished.connect(worker.deleteLater)

        worker.start()
        progress.setValue(0)
        loop.exec()
        # QProgressDialog.close() emits canceled() from its own closeEvent, so
        # the handler has to come off FIRST or every successful import reports
        # itself as cancelled on the way out.
        progress.canceled.disconnect(_cancel)
        progress.close()

        if result.get("cancelled"):
            self._pending_import = worker    # keep it alive until it finishes
            return None
        worker.wait()
        if "error" in result:
            QMessageBox.critical(self, "Import failed", result["error"][0])
            return None
        return result.get("info")
