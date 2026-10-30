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

"""Where flexiprobes keeps its per-user state.

One resolution of ``$FLEXPROBE_HOME``, resolved at CALL time. Both halves of
that matter:

* **One.** The same three lines were copy-pasted into ``settings``,
  ``run_history``, ``project_history``, ``query_drafts`` and
  ``flexprobe.hpc_export``. Five copies is five chances for the app to keep
  half its state somewhere the user did not ask for.
* **At call time.** Caching the directory in a module-level constant would bind
  it at import, and the test suite redirects ``FLEXPROBE_HOME`` per test with an
  autouse fixture. A module that read the environment once at import would spend
  the whole suite writing into the developer's real ``~/.flexprobe``.
"""
from __future__ import annotations

import os


def flexprobe_home(*parts: str, create: bool = True) -> str:
    """``$FLEXPROBE_HOME`` (default ``~/.flexprobe``), joined with ``parts``.

    Creates the directory by default, because every caller so far wants to
    write into it moments later. Pass ``create=False`` to resolve the path
    without touching the filesystem.
    """
    base = os.environ.get("FLEXPROBE_HOME") or os.path.join(
        os.path.expanduser("~"), ".flexprobe"
    )
    path = os.path.join(base, *parts) if parts else base
    if create:
        os.makedirs(path, exist_ok=True)
    return path
