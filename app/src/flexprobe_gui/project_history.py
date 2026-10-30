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

"""Project bag — groups multiple runs together.

A Project is just a named bucket: it doesn't own data of its own beyond
``name`` + ``description`` + ``created_at``. Runs are tagged with
``project_id`` in their ``RunRecord``; the project is what the home-page
filter and the bulk-export action key off of.

Persisted to ``~/.flexprobe/projects.json`` so the list survives restarts.
"""
from __future__ import annotations

import datetime
import json
import os
import threading
import uuid
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Iterable, List, Optional


def _now_iso() -> str:
    return datetime.datetime.now().isoformat(timespec="seconds")


def projects_path() -> str:
    base = os.environ.get("FLEXPROBE_HOME") or os.path.join(
        os.path.expanduser("~"), ".flexprobe"
    )
    os.makedirs(base, exist_ok=True)
    return os.path.join(base, "projects.json")


@dataclass
class Project:
    project_id: str = field(default_factory=lambda: uuid.uuid4().hex[:12])
    name: str = ""
    description: str = ""
    created_at: str = field(default_factory=_now_iso)

    def to_json(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Project":
        known = {f for f in cls.__dataclass_fields__}
        return cls(**{k: v for k, v in data.items() if k in known})


class ProjectHistory:
    """Thread-safe JSON-backed project list."""

    def __init__(self, path: Optional[str] = None):
        self._path = path or projects_path()
        self._lock = threading.Lock()
        self._projects: List[Project] = self._load()

    def _load(self) -> List[Project]:
        if not os.path.exists(self._path):
            return []
        try:
            with open(self._path) as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError):
            return []
        return [Project.from_json(d) for d in data if isinstance(d, dict)]

    def _save_locked(self) -> None:
        tmp = self._path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump([p.to_json() for p in self._projects], fh, indent=2)
        os.replace(tmp, self._path)

    def list(self) -> List[Project]:
        with self._lock:
            return sorted(self._projects, key=lambda p: p.name.lower())

    def get(self, project_id: str) -> Optional[Project]:
        with self._lock:
            for p in self._projects:
                if p.project_id == project_id:
                    return p
        return None

    def create(self, name: str, description: str = "") -> Project:
        name = (name or "").strip()
        if not name:
            raise ValueError("Project name must be non-empty.")
        with self._lock:
            # Re-use any existing project with the same name (case-insensitive)
            for p in self._projects:
                if p.name.lower() == name.lower():
                    return p
            project = Project(name=name, description=description.strip())
            self._projects.append(project)
            self._save_locked()
            return project

    def update(
        self,
        project_id: str,
        name: Optional[str] = None,
        description: Optional[str] = None,
    ) -> Optional[Project]:
        """Rename and/or re-describe a project in place.

        ``name``/``description`` left as ``None`` are unchanged. A blank new
        name is rejected; a name colliding (case-insensitively) with a
        *different* project raises ValueError so the caller can warn.
        """
        with self._lock:
            target: Optional[Project] = None
            for p in self._projects:
                if p.project_id == project_id:
                    target = p
                    break
            if target is None:
                return None
            if name is not None:
                name = name.strip()
                if not name:
                    raise ValueError("Project name must be non-empty.")
                for p in self._projects:
                    if p.project_id != project_id and p.name.lower() == name.lower():
                        raise ValueError(
                            f"Another project is already named '{name}'."
                        )
                target.name = name
            if description is not None:
                target.description = description.strip()
            self._save_locked()
            return target

    def delete(self, project_id: str) -> bool:
        with self._lock:
            n_before = len(self._projects)
            self._projects = [p for p in self._projects if p.project_id != project_id]
            if len(self._projects) != n_before:
                self._save_locked()
                return True
            return False

    def delete_many(self, project_ids: Iterable[str]) -> int:
        """Delete several projects at once. Returns how many were removed.

        Like :meth:`delete`, runs attached to the removed projects stay on
        disk and simply lose their project filter. Persists a single time.
        """
        ids = set(project_ids)
        if not ids:
            return 0
        with self._lock:
            n_before = len(self._projects)
            self._projects = [p for p in self._projects if p.project_id not in ids]
            removed = n_before - len(self._projects)
            if removed:
                self._save_locked()
            return removed
