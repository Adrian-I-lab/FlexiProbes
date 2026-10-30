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

"""Open a folder in the OS file browser.

Single home for what used to be two byte-identical copies (``_reveal`` in
home_window and ``reveal_in_finder`` in review_dialog)."""

from __future__ import annotations

import os
import subprocess
import sys


def reveal_in_finder(path: str) -> bool:
    """Reveal ``path`` in Finder / Explorer / the default file manager.

    Returns False, doing nothing, for a missing / non-directory path, so a
    caller can tell the user instead of the button silently doing nothing (a
    moved run folder used to look exactly like a broken Open button). True
    means a launcher was asked, not that a window opened. Best-effort: a
    launcher that isn't present just fails silently (``check=False``) rather
    than raising into the UI thread."""
    if not path or not os.path.isdir(path):
        return False
    if sys.platform == "darwin":
        subprocess.run(["open", path], check=False)
    elif os.name == "nt":
        subprocess.run(["explorer", path], check=False)
    else:
        subprocess.run(["xdg-open", path], check=False)
    return True
