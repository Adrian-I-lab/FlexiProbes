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

"""Sequential run queue. A single background QThread runs at a time."""
from __future__ import annotations

import logging
import os
import threading
import traceback
from collections import deque
from typing import Any, Deque, Dict, Optional, Tuple

from PySide6.QtCore import QObject, QThread, Signal

from flexprobe.auto_export import AutoExportStatus, auto_export
from flexprobe.constants import (
    MIN_PROBE_GAP_BP, MIN_PROBE_GAP_FLOOR_BP, RELAXED_PROBE_GAP_BP,
)
from flexprobe.errors import FlexprobeError
from flexprobe.models import RunConfig
from flexprobe.pipeline import ProgressEvent, run_pipeline
from flexprobe.summary import write_failed_summary, write_ok_summary

from .run_history import RunHistory, RunRecord, now_iso

logger = logging.getLogger("flexprobe_gui.run_manager")

# Sidecar holding the exact RunConfig, written into a run's output dir at submit
# time. It lets an interrupted run (cancelled, crashed, or power-off) be resumed
# byte-for-byte into the SAME folder, where its checkpoints live.
RUN_CONFIG_SIDECAR = "run_config.json"


#: Stage text shown on the home page for each non-exported outcome. Phrased as
#: what the user should do next, because that is the only question a finished
#: run raises. Every non-``EXPORTED`` status leaves the run in
#: ``ready_for_review``: the probes exist in ``results.csv``, so a manual pick
#: in the review dialog still gets the user their order sheets.
_AUTO_EXPORT_STAGE = {
    AutoExportStatus.NO_RESULTS: "no results to auto-export; review manually",
    AutoExportStatus.NOTHING_SELECTED:
        "auto-select picked no probes; review manually",
    AutoExportStatus.EXPORT_FAILED: "auto-export failed; review manually",
}


def _config_sidecar_path(output_dir: str) -> str:
    return os.path.join(output_dir, RUN_CONFIG_SIDECAR)


def load_run_config(output_dir: str) -> Optional[RunConfig]:
    """Load the persisted RunConfig from a run's output dir, or None."""
    path = _config_sidecar_path(output_dir)
    if not os.path.exists(path):
        return None
    try:
        from flexprobe.config_io import load_config
        return load_config(path)
    except Exception:
        logger.warning("Could not load %s from %s", RUN_CONFIG_SIDECAR,
                       output_dir, exc_info=True)
        return None


def is_resumable(record: RunRecord) -> bool:
    """True if ``record`` has saved progress on disk that a resume can reuse.

    Requires the config sidecar (so the exact run can be reconstructed) plus at
    least one checkpoint dir (candidate cache or BLAST chunk cache). A fully
    completed run has its checkpoints cleared, so it is never "resumable".
    """
    od = record.output_dir
    if not od or not os.path.isdir(od):
        return False
    if not os.path.exists(_config_sidecar_path(od)):
        return False
    from flexprobe import screen_cache, candidate_cache
    return (os.path.isdir(candidate_cache.cache_dir(od))
            or os.path.isdir(screen_cache.cache_dir(od)))


class _DesignWorker(QThread):
    progress = Signal(str, int, str)     # run_id, pct, "stage: message"
    finished_ok = Signal(str, dict)      # run_id, pipeline result dict
    failed = Signal(str, str)            # run_id, error message
    cancelled = Signal(str)              # run_id

    def __init__(
        self,
        run_id: str,
        config: RunConfig,
        cached_hits: Optional[Any] = None,
        parent: Optional[QObject] = None,
    ):
        super().__init__(parent)
        self._run_id = run_id
        self._config = config
        self._cached_hits = cached_hits
        self._cancel = threading.Event()

    def request_cancel(self) -> None:
        self._cancel.set()

    def _on_progress(self, ev: ProgressEvent) -> None:
        # pct is None for messages that carry no progress information. The Qt
        # signal is typed int, so -1 is the sentinel for "no change"; the slot
        # leaves the bar alone and updates only the stage text.
        pct = -1 if ev.pct is None else int(ev.pct)
        self.progress.emit(self._run_id, pct, f"{ev.stage}: {ev.message}")

    def run(self) -> None:
        try:
            result = run_pipeline(
                self._config,
                cancel_event=self._cancel,
                progress_callback=self._on_progress,
                screen_hits_by_probe=self._cached_hits,
            )
        except FlexprobeError as exc:
            # An emptied funnel is a RESULT, not a non-event, and without this
            # a GUI run that produced 0 probes left no summary AND no manifest
            # — indistinguishable from a run that never started. The CLI has
            # written this since the 0.3.10 "0 tiles" mislabelling; the GUI
            # path was simply never wired up to the same writer.
            write_failed_summary(self._config.output_dir, exc)
            if self._cancel.is_set():
                self.cancelled.emit(self._run_id)
            else:
                self.failed.emit(self._run_id, str(exc))
            return
        except Exception as exc:  # noqa: BLE001 — reported verbatim below
            write_failed_summary(self._config.output_dir, exc)
            self.failed.emit(self._run_id, traceback.format_exc())
            return

        if self._cancel.is_set():
            self.cancelled.emit(self._run_id)
        else:
            write_ok_summary(self._config.output_dir, result["stats"])
            self.finished_ok.emit(self._run_id, result)


class RunManager(QObject):
    """FIFO queue over a single background worker.

    Submitting a run while another is active enqueues it. When the active
    worker finishes, the next queued run is dispatched automatically. Emits
    ``record_changed`` (run_id) so the UI can refresh its display.
    """

    record_changed = Signal(str)

    def __init__(self, history: RunHistory, parent: Optional[QObject] = None):
        super().__init__(parent)
        self._history = history
        self._queue: Deque[Dict[str, Any]] = deque()
        self._current_worker: Optional[_DesignWorker] = None
        self._current_run_id: Optional[str] = None

    @property
    def history(self) -> RunHistory:
        return self._history

    def submit(
        self,
        record: RunRecord,
        config: RunConfig,
        cached_hits: Optional[Any] = None,
    ) -> None:
        record.status = "queued"
        self._history.append(record)
        self._persist_config(record, config)
        self._queue.append({"record": record, "config": config, "cached_hits": cached_hits})
        self.record_changed.emit(record.run_id)
        self._maybe_start_next()

    def resume(self, run_id: str) -> bool:
        """Re-launch an interrupted run into the SAME output dir.

        Loads the run's persisted config, forces resume + overwrite (so the
        existing checkpoints are reused and any partial final output is
        replaced), and re-queues the existing record. Returns False if the run
        can't be resumed (no record, no output dir, or no saved config).
        """
        rec = self._history.get(run_id)
        if rec is None or not rec.output_dir:
            return False
        cfg = load_run_config(rec.output_dir)
        if cfg is None:
            return False
        cfg.output_dir = rec.output_dir       # never relocate a resume
        cfg.resume_screen = True               # reuse the checkpoints
        cfg.export.overwrite = True           # replace any partial final outputs
        rec.status = "queued"
        rec.error = None
        rec.completed_at = None
        rec.stage = "resuming"
        self._history.update(rec)
        self._queue.append({"record": rec, "config": cfg, "cached_hits": None})
        self.record_changed.emit(rec.run_id)
        self._maybe_start_next()
        return True

    def cancel(self, run_id: str) -> None:
        if self._current_run_id == run_id and self._current_worker is not None:
            self._current_worker.request_cancel()
            return
        for item in list(self._queue):
            if item["record"].run_id == run_id:
                self._queue.remove(item)
                rec = item["record"]
                rec.status = "cancelled"
                rec.completed_at = now_iso()
                self._history.update(rec)
                self.record_changed.emit(run_id)
                return

    # ---- internals ----

    def _persist_config(self, record: RunRecord, config: RunConfig) -> None:
        """Write the run's config sidecar so it can be resumed after a crash."""
        if not config.output_dir:
            return
        try:
            os.makedirs(config.output_dir, exist_ok=True)
            from flexprobe.config_io import save_config
            save_config(config, _config_sidecar_path(config.output_dir))
        except OSError:
            logger.warning("Could not persist %s for run %s",
                           RUN_CONFIG_SIDECAR, record.run_id, exc_info=True)

    def _maybe_start_next(self) -> None:
        if self._current_worker is not None or not self._queue:
            return
        item = self._queue.popleft()
        record: RunRecord = item["record"]
        record.status = "running"
        record.started_at = now_iso()
        self._history.update(record)
        # Persist the now-running record to disk so a crash / power-off mid-run
        # leaves a trace; on next launch it loads back as "interrupted" and can
        # be resumed (its checkpoints are intact in the output dir).
        self._history.flush()

        worker = _DesignWorker(
            record.run_id, item["config"], item["cached_hits"], parent=self,
        )
        worker.progress.connect(self._on_progress)
        worker.finished_ok.connect(self._on_finished_ok)
        worker.failed.connect(self._on_failed)
        worker.cancelled.connect(self._on_cancelled)
        worker.finished.connect(worker.deleteLater)
        self._current_worker = worker
        self._current_run_id = record.run_id
        self.record_changed.emit(record.run_id)
        worker.start()

    def _on_progress(self, run_id: str, pct: int, msg: str) -> None:
        rec = self._history.get(run_id)
        if rec is None:
            return
        if pct >= 0:          # -1 == "no progress information", keep the bar
            rec.progress = pct
        rec.stage = msg
        self._history.update(rec)
        self.record_changed.emit(run_id)

    def _finalize_terminal(self, _run_id: str) -> None:
        # _run_id kept for a uniform call signature; not needed here.
        self._current_worker = None
        self._current_run_id = None
        self._maybe_start_next()

    def _on_finished_ok(self, run_id: str, result: dict) -> None:
        rec = self._history.get(run_id)
        if rec is not None:
            stats = result.get("stats") or {}
            # Discovery-mode runs stop at "ready_for_review" — the user picks
            # their final probes in the review dialog, which then writes the
            # order CSVs and flips status to "complete". Auto-select bypasses
            # the manual review entirely (whole-transcriptome convenience).
            discovery = bool(stats.get("discovery_only"))
            auto_n = int(stats.get("auto_select_per_gene") or 0)
            if discovery and auto_n > 0:
                # The outcome DECIDES the status. It used to be discarded and
                # "complete" assigned unconditionally two lines later, so a run
                # whose exporter raised — or that selected nothing at all — was
                # presented as finished and exported with no order CSVs on disk.
                outcome = self._auto_export(rec, auto_n)
                rec.progress = 100
                if outcome is AutoExportStatus.EXPORTED:
                    rec.status = "complete"
                    rec.stage = f"auto-exported (top {auto_n}/gene)"
                else:
                    rec.status = "ready_for_review"
                    rec.stage = _AUTO_EXPORT_STAGE[outcome]
            else:
                rec.status = "ready_for_review" if discovery else "complete"
                rec.progress = 100
                rec.stage = "ready for review" if discovery else "done"
            rec.completed_at = now_iso()
            rec.stats = stats
            self._history.update(rec)
            self.record_changed.emit(run_id)
        self._finalize_terminal(run_id)

    def _auto_export(self, rec: RunRecord, n_per_gene: int) -> AutoExportStatus:
        """Auto-pick top-N probes per gene and write the order CSVs.

        Returns what happened, because the caller has to decide the run's final
        status from it: only ``EXPORTED`` may be shown as a completed, exported
        run. See :class:`AutoExportStatus`.

        The selection, the export and the wording of every message live in
        ``flexprobe.auto_export``, shared with ``cli._auto_export``. All that is
        left here is the channel — the diagnostics log — and finding the run's
        spacing flags, which the GUI has to read back off disk.
        """
        min_gap, relaxed_gap, floor_gap = self._spacing_gaps(rec)
        result = auto_export(rec.output_dir, n_per_gene, min_gap=min_gap,
                             relaxed_gap=relaxed_gap, floor_gap=floor_gap)
        for line in result.report_lines():
            logger.info("Run %s: %s", rec.run_id, line)
        for line in result.warning_lines():
            logger.warning("Run %s: %s", rec.run_id, line)
        if result.status is AutoExportStatus.EXPORT_FAILED:
            # Non-fatal: the probes are picked and saved, so the review dialog
            # can still write the order sheets by hand. It is NOT a success,
            # though — the caller turns this into "ready_for_review", never
            # "complete". The traceback goes to the log file installed by
            # ``flexprobe_gui.diagnostics`` at startup, which is the only place
            # a packaged .app has to write it.
            logger.error("Auto-export failed for run %s; falling back to "
                         "manual review.", rec.run_id,
                         exc_info=result.error)
        return result.status

    @staticmethod
    def _spacing_gaps(rec: RunRecord) -> Tuple[int, int, int]:
        """The run's own ``--min/--relaxed/--floor-probe-gap``, off its sidecar.

        The CLI is handed these; the GUI has to recover them from the
        ``run_config.json`` written at submit time. When that file is missing or
        unreadable the module defaults are the only thing left — but a run the
        user asked to space differently would then be picked at 150/40/15 and
        the order CSVs would look completely normal, so the substitution is
        logged rather than made silently.
        """
        run_cfg = load_run_config(rec.output_dir)
        if run_cfg is None:
            logger.warning(
                "Run %s: no %s to read the probe spacing from; auto-selecting "
                "at the defaults (%d/%d/%d bp), which may not be what this run "
                "was launched with.", rec.run_id, RUN_CONFIG_SIDECAR,
                MIN_PROBE_GAP_BP, RELAXED_PROBE_GAP_BP, MIN_PROBE_GAP_FLOOR_BP,
            )
            return MIN_PROBE_GAP_BP, RELAXED_PROBE_GAP_BP, MIN_PROBE_GAP_FLOOR_BP
        return (run_cfg.min_probe_gap_bp, run_cfg.relaxed_probe_gap_bp,
                run_cfg.floor_probe_gap_bp)

    def _on_failed(self, run_id: str, message: str) -> None:
        rec = self._history.get(run_id)
        if rec is not None:
            rec.status = "failed"
            rec.error = message
            rec.completed_at = now_iso()
            self._history.update(rec)
            self.record_changed.emit(run_id)
        self._finalize_terminal(run_id)

    def _on_cancelled(self, run_id: str) -> None:
        rec = self._history.get(run_id)
        if rec is not None:
            rec.status = "cancelled"
            rec.completed_at = now_iso()
            self._history.update(rec)
            self.record_changed.emit(run_id)
        self._finalize_terminal(run_id)
