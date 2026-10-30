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

"""Persisted run history for the home page.

Active runs (queued / running) live only in memory; once a run reaches a
terminal status (complete / failed / cancelled) the full list is flushed to
``~/.flexprobe/runs.json`` so the home-page table survives restarts.
"""
from __future__ import annotations

import datetime
import json
import logging
import os
import shutil
import threading
import uuid
from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Optional

from .paths import flexprobe_home

logger = logging.getLogger("flexprobe_gui.run_history")


TERMINAL_STATUSES = {"complete", "failed", "cancelled", "ready_for_review",
                     "interrupted"}
ACTIVE_STATUSES = {"queued", "running"}
# Statuses flushed to disk. "running" is persisted (once, when a run starts) so a
# crash or power-off mid-run leaves a trace; on the next launch :meth:`_load`
# rewrites any such stale "running"/"queued" record as "interrupted" — a
# terminal, resumable state (its checkpoints survive in the output dir).
PERSIST_STATUSES = TERMINAL_STATUSES | {"running"}
# Discovery-mode runs land in ``ready_for_review`` after the pipeline finishes
# and stay there until the user confirms a selection in the review dialog.
# They count as "completed runs" for the home-page table.


@dataclass
class HistoryLoadError:
    """A runs.json that could not be read, and where its remains were put."""

    path: str
    backup_path: Optional[str]
    message: str


#: Called once with a :class:`HistoryLoadError` when a history file cannot be
#: parsed. Left unset by default so importing this module never assumes a UI
#: exists (the CLI, the tests and the HPC tooling all construct RunHistory);
#: ``flexprobe_gui.app`` installs a dialog here at startup.
_load_error_notifier: Optional[Callable[[HistoryLoadError], None]] = None


def set_load_error_notifier(
    fn: Optional[Callable[[HistoryLoadError], None]],
) -> None:
    """Install (or clear, with ``None``) the "your history is damaged" hook."""
    global _load_error_notifier
    _load_error_notifier = fn


def _backup_corrupt_file(path: str) -> Optional[str]:
    """Copy an unreadable ``path`` aside and return where it went.

    ``.bak`` first, then ``.bak.1``, ``.bak.2``… A second corruption must not
    overwrite the first backup: the first one is the copy most likely to still
    hold the user's real history, since each subsequent load starts from a file
    this code has already rewritten.

    COPY, not move. The app has to keep working, which means writing to ``path``
    again, and a move would leave the next reader with no file and no hint that
    anything was ever wrong.
    """
    candidate = path + ".bak"
    n = 0
    while os.path.exists(candidate):
        n += 1
        candidate = f"{path}.bak.{n}"
        if n > 999:                # pathological; overwrite the newest, not the first
            break
    try:
        shutil.copy2(path, candidate)
    except OSError:
        logger.error("Could not preserve the damaged history %s", path,
                     exc_info=True)
        return None
    return candidate


def now_iso() -> str:
    return datetime.datetime.now().isoformat(timespec="seconds")


def history_path() -> str:
    return os.path.join(flexprobe_home(), "runs.json")


@dataclass
class RunRecord:
    run_id: str = field(default_factory=lambda: uuid.uuid4().hex[:12])
    name: str = ""
    created_at: str = field(default_factory=now_iso)
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    modality_label: str = ""
    fasta_path: str = ""
    output_dir: str = ""
    status: str = "queued"          # queued | running | complete | failed | cancelled
    progress: int = 0
    stage: str = ""
    stats: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    config_summary: Dict[str, Any] = field(default_factory=dict)
    # Empty string = not attached to any project.
    project_id: str = ""
    # Full New-query form state (a ``NewQueryDialog.to_draft()`` payload) so
    # the run can be re-opened, pre-filled and re-launched. Empty for runs
    # created before re-run support; those fall back to ``config_summary``.
    query_draft: Dict[str, Any] = field(default_factory=dict)
    # The finished run's identity, copied from its manifest.json (``timestamp``
    # and ``inputs.fasta_sha256``). "Locate…" uses it to prove that a folder
    # found somewhere else is THIS run and not merely one with the same name.
    # Empty for records saved before it was recorded, and for runs that never
    # wrote a manifest. See flexprobe_gui.run_import.
    manifest_timestamp: str = ""
    fasta_sha256: str = ""

    def to_json(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "RunRecord":
        # Forward-compatible: ignore unknown keys.
        known = {f for f in cls.__dataclass_fields__}
        return cls(**{k: v for k, v in data.items() if k in known})


class RunHistory:
    """Thread-safe list of RunRecord with on-disk persistence."""

    def __init__(self, path: Optional[str] = None):
        self._path = path or history_path()
        self._lock = threading.Lock()
        #: Set when the history on disk could not be read. Not cleared by later
        #: writes: it is the record of an event, and the UI has to be able to
        #: tell the user about it whenever it gets round to asking.
        self.load_error: Optional[HistoryLoadError] = None
        self._records: List[RunRecord] = self._load()

    def _load(self) -> List[RunRecord]:
        if not os.path.exists(self._path):
            return []
        try:
            with open(self._path) as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            return self._on_corrupt(f"{type(exc).__name__}: {exc}")
        if not isinstance(data, list):
            # A dict parses fine and then iterates as its KEYS, so the old
            # comprehension below turned it into an empty history without a
            # word — a file that was never a run history at all, silently
            # accepted as an empty one.
            return self._on_corrupt(
                f"expected a list of runs, found {type(data).__name__}")
        records = [RunRecord.from_json(d) for d in data if isinstance(d, dict)]
        # A run still marked active in the file means the app exited before it
        # finished (crash / force-quit / power-off). It can't be live now, so
        # rewrite it as "interrupted" — terminal and resumable.
        for r in records:
            if r.status in ACTIVE_STATUSES:
                r.status = "interrupted"
                if not r.stage:
                    r.stage = "interrupted — resume to continue"
        return records

    def _on_corrupt(self, why: str) -> List[RunRecord]:
        """Preserve a damaged history, say so loudly, and carry on empty.

        Carrying on empty is unavoidable — there is nothing to load — but it
        used to be the WHOLE behaviour: the empty list was written straight back
        over the file by the next terminal run, and three recorded runs became
        zero with nothing kept and nothing said. A half-written runs.json is
        precisely what a power-off during a flush leaves behind, so this is the
        failure that follows the failure the flush exists to survive.
        """
        backup = _backup_corrupt_file(self._path)
        kept = (f"a copy was kept at {backup}" if backup
                else "the file could NOT be copied aside")
        self.load_error = HistoryLoadError(
            path=self._path, backup_path=backup,
            message=(f"The run history at {self._path} could not be read "
                     f"({why}). Previously recorded runs will not be listed; "
                     f"{kept}. Their output folders are untouched."),
        )
        logger.error("Damaged run history %s (%s); %s",
                     self._path, why, kept)
        if _load_error_notifier is not None:
            try:
                _load_error_notifier(self.load_error)
            except Exception:       # noqa: BLE001 — a UI hook must not block startup
                logger.exception("The run-history error notifier itself failed")
        return []

    def _save_locked(self):
        to_save = [r.to_json() for r in self._records if r.status in PERSIST_STATUSES]
        tmp = self._path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(to_save, fh, indent=2)
        os.replace(tmp, self._path)

    def flush(self) -> None:
        """Force a save (used to persist a run the moment it starts running)."""
        with self._lock:
            self._save_locked()

    def append(self, record: RunRecord) -> None:
        with self._lock:
            self._records.append(record)
            if record.status in TERMINAL_STATUSES:
                self._save_locked()

    def update(self, record: RunRecord) -> None:
        with self._lock:
            for i, r in enumerate(self._records):
                if r.run_id == record.run_id:
                    self._records[i] = record
                    break
            else:
                self._records.append(record)
            if record.status in TERMINAL_STATUSES:
                self._save_locked()

    def get(self, run_id: str) -> Optional[RunRecord]:
        with self._lock:
            for r in self._records:
                if r.run_id == run_id:
                    return r
        return None

    def delete(self, run_id: str) -> bool:
        """Drop a run from history. Returns True if a record was removed.

        Only persists when the deleted record was terminal (the on-disk file
        holds terminal runs only); deleting an active run just removes it from
        the in-memory list.
        """
        with self._lock:
            for i, r in enumerate(self._records):
                if r.run_id == run_id:
                    del self._records[i]
                    if r.status in TERMINAL_STATUSES:
                        self._save_locked()
                    return True
            return False

    def delete_many(self, run_ids: Iterable[str]) -> int:
        """Drop several runs in one shot. Returns how many were removed.

        Persists once if *any* removed record was terminal (the on-disk file
        holds terminal runs only), so a bulk delete is a single disk write.
        """
        ids = set(run_ids)
        if not ids:
            return 0
        with self._lock:
            kept: List[RunRecord] = []
            removed = 0
            removed_terminal = False
            for r in self._records:
                if r.run_id in ids:
                    removed += 1
                    if r.status in TERMINAL_STATUSES:
                        removed_terminal = True
                else:
                    kept.append(r)
            if removed:
                self._records = kept
                if removed_terminal:
                    self._save_locked()
            return removed

    def set_project_many(self, run_ids: Iterable[str], project_id: str) -> int:
        """Reassign several runs to ``project_id`` (``""`` = unfiled).

        Returns how many records actually changed. Persists once if any
        changed record is terminal.
        """
        ids = set(run_ids)
        if not ids:
            return 0
        with self._lock:
            changed = 0
            changed_terminal = False
            for r in self._records:
                if r.run_id in ids and r.project_id != project_id:
                    r.project_id = project_id
                    changed += 1
                    if r.status in TERMINAL_STATUSES:
                        changed_terminal = True
            if changed_terminal:
                self._save_locked()
            return changed

    def list_all(self) -> List[RunRecord]:
        """Every record, active and terminal, in insertion order."""
        with self._lock:
            return list(self._records)

    def list_active(self) -> List[RunRecord]:
        """Queued + running, oldest first (FIFO ordering)."""
        with self._lock:
            return sorted(
                (r for r in self._records if r.status in ACTIVE_STATUSES),
                key=lambda r: r.created_at,
            )

    def list_completed(self) -> List[RunRecord]:
        """Terminal-state runs, newest first."""
        with self._lock:
            return sorted(
                (r for r in self._records if r.status in TERMINAL_STATUSES),
                key=lambda r: r.completed_at or r.created_at,
                reverse=True,
            )
