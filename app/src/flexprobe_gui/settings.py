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

"""Persisted application preferences.

A tiny JSON-backed key/value store living next to the run/project history in
``~/.flexprobe/settings.json`` (or under ``$FLEXPROBE_HOME``). Currently holds
the user's default output folder so every new query pre-fills the same place
instead of the hard-coded ``~/flexprobe_out``.
"""
from __future__ import annotations

import json
import os
import tempfile
import threading
from typing import Any, Dict, Optional


def settings_path() -> str:
    base = os.environ.get("FLEXPROBE_HOME") or os.path.join(
        os.path.expanduser("~"), ".flexprobe"
    )
    os.makedirs(base, exist_ok=True)
    return os.path.join(base, "settings.json")


_DEFAULT_OUTPUT_DIR = os.path.join(os.path.expanduser("~"), "flexprobe_out")


class AppSettings:
    """Thread-safe JSON-backed preferences blob.

    The file, not this object's cache, is the source of truth. Two copies of
    the app can be open at once, and each used to hold its own snapshot from
    start-up and rewrite the WHOLE file from it on every save, so the second
    copy silently put back whatever the first had just changed. Every read and
    every write therefore re-reads the file first, and a write only replaces
    the one key it was asked to set.
    """

    def __init__(self, path: Optional[str] = None):
        self._path = path or settings_path()
        self._lock = threading.Lock()
        self._data: Dict[str, Any] = self._load() or {}

    def _load(self) -> Optional[Dict[str, Any]]:
        """The file's contents; {} when absent; None when it cannot be read.

        None is kept distinct from {} so a transient read failure never looks
        like "no settings" and gets written back over the real file.
        """
        if not os.path.exists(self._path):
            return {}
        try:
            with open(self._path) as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError):
            return None
        return data if isinstance(data, dict) else None

    def _refresh_locked(self) -> None:
        data = self._load()
        if data is not None:
            self._data = data

    def _save_locked(self) -> None:
        # A unique temp file per write: a fixed "<path>.tmp" is shared by every
        # process, so two concurrent saves could interleave into one file.
        fd, tmp = tempfile.mkstemp(
            prefix=os.path.basename(self._path) + ".",
            suffix=".tmp", dir=os.path.dirname(self._path) or ".",
        )
        try:
            with os.fdopen(fd, "w") as fh:
                json.dump(self._data, fh, indent=2)
            os.replace(tmp, self._path)
        except BaseException:
            if os.path.exists(tmp):
                os.remove(tmp)
            raise

    def get(self, key: str, default: Any = None) -> Any:
        with self._lock:
            self._refresh_locked()
            return self._data.get(key, default)

    def set(self, key: str, value: Any) -> None:
        with self._lock:
            self._refresh_locked()
            self._data[key] = value
            self._save_locked()

    # ---- Convenience accessors ----

    def default_output_dir(self) -> str:
        """The configured default output base dir, or the built-in fallback."""
        value = self.get("default_output_dir", "")
        return value or _DEFAULT_OUTPUT_DIR

    def set_default_output_dir(self, path: str) -> None:
        self.set("default_output_dir", path.strip())
