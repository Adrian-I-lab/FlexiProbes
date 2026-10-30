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

"""Saved query drafts — resume a half-configured query after restarting.

A draft is just the New-query dialog's form state captured as a flat JSON
dict, stored one-file-per-draft under ``~/.flexprobe/query_drafts/`` (or
``$FLEXPROBE_HOME/query_drafts``). Drafts are keyed by name; saving under an
existing name overwrites it.

The dialog owns the shape of the payload (see ``NewQueryDialog.to_draft`` /
``apply_draft``); this module only handles naming, listing and persistence.
"""
from __future__ import annotations

import json
import os
import threading
from typing import Any, Dict, List, Optional


def drafts_dir() -> str:
    base = os.environ.get("FLEXPROBE_HOME") or os.path.join(
        os.path.expanduser("~"), ".flexprobe"
    )
    d = os.path.join(base, "query_drafts")
    os.makedirs(d, exist_ok=True)
    return d


def _safe_filename(name: str) -> str:
    """Map a user-visible draft name to a filesystem-safe stem."""
    return "".join(c if c.isalnum() or c in "_-." else "_" for c in name).strip("_") or "draft"


class QueryDraftStore:
    """Thread-safe one-file-per-draft store of query form state."""

    def __init__(self, directory: Optional[str] = None):
        self._dir = directory or drafts_dir()
        self._lock = threading.Lock()

    def _path_for(self, name: str) -> str:
        return os.path.join(self._dir, _safe_filename(name) + ".json")

    def list(self) -> List[str]:
        """Return saved draft names (from each file's stored ``name``), sorted."""
        names: List[str] = []
        with self._lock:
            try:
                files = os.listdir(self._dir)
            except OSError:
                return []
            for fn in files:
                if not fn.endswith(".json"):
                    continue
                try:
                    with open(os.path.join(self._dir, fn)) as fh:
                        data = json.load(fh)
                except (OSError, json.JSONDecodeError):
                    continue
                if isinstance(data, dict) and data.get("name"):
                    names.append(str(data["name"]))
        return sorted(names, key=str.lower)

    def save(self, name: str, payload: Dict[str, Any]) -> str:
        name = (name or "").strip()
        if not name:
            raise ValueError("Draft name must be non-empty.")
        record = dict(payload)
        record["name"] = name
        with self._lock:
            path = self._path_for(name)
            tmp = path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(record, fh, indent=2)
            os.replace(tmp, path)
        return name

    def load(self, name: str) -> Optional[Dict[str, Any]]:
        with self._lock:
            path = self._path_for(name)
            if not os.path.exists(path):
                return None
            try:
                with open(path) as fh:
                    data = json.load(fh)
            except (OSError, json.JSONDecodeError):
                return None
        return data if isinstance(data, dict) else None

    def delete(self, name: str) -> bool:
        with self._lock:
            path = self._path_for(name)
            if os.path.exists(path):
                os.remove(path)
                return True
            return False
