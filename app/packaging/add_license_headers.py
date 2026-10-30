#!/usr/bin/env python3
"""Add the GPL notice header to every Python source file under ``app/src/``.

Why
---
GPLv3's "How to Apply These Terms to Your New Programs" recommends attaching
the notice to the start of each source file. Without it, a file lifted out of
the tree carries no indication of its terms, and automated licence scanners
see nothing.

Placement
---------
The block goes in as a **comment**, after any shebang and encoding line and
*before* the module docstring. Comments are not statements, so the docstring
remains the first statement in the module and ``__doc__`` is untouched. Header
insertion therefore cannot disturb documented behaviour, doctests, or
``help()``.

Idempotent: a file that already carries the SPDX identifier is left alone, so
this is safe to re-run and safe to wire into CI.

Usage
-----
    python packaging/add_license_headers.py            # add to app/src
    python packaging/add_license_headers.py --check    # CI: fail if any lack it
    python packaging/add_license_headers.py --exclude src/flexprobe/cli.py
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent.parent
DEFAULT_ROOTS = [APP_DIR / "src"]

SPDX = "SPDX-License-Identifier: GPL-3.0-or-later"

# Copyright holder matches app/NOTICE. Keep the two in sync.
HEADER = f"""\
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
# {SPDX}
"""

SHEBANG_RE = re.compile(r"^#!")
# PEP 263 encoding declaration, valid on line 1 or 2.
CODING_RE = re.compile(r"^[ \t\f]*#.*?coding[:=][ \t]*[-_.a-zA-Z0-9]+")


def insertion_point(lines: list[str]) -> int:
    """Index at which the header should be inserted.

    After a shebang and/or a PEP 263 encoding line, and before everything else
    (notably before the module docstring, which stays the first *statement*).
    """
    idx = 0
    if idx < len(lines) and SHEBANG_RE.match(lines[idx]):
        idx += 1
    # PEP 263 only honours the declaration on the first two physical lines.
    if idx < len(lines) and idx < 2 and CODING_RE.match(lines[idx]):
        idx += 1
    return idx


def has_header(text: str) -> bool:
    return SPDX in text


def add_header(path: Path) -> bool:
    """Insert the header. Returns True if the file was modified."""
    text = path.read_text(encoding="utf-8")
    if has_header(text):
        return False

    lines = text.splitlines(keepends=True)
    idx = insertion_point(lines)

    block = HEADER
    # Keep exactly one blank line between the header and whatever follows.
    if idx < len(lines) and lines[idx].strip():
        block += "\n"

    lines[idx:idx] = [block]
    path.write_text("".join(lines), encoding="utf-8")
    return True


def iter_sources(roots: list[Path], excludes: set[Path]):
    for root in roots:
        for path in sorted(root.rglob("*.py")):
            if "__pycache__" in path.parts:
                continue
            if path.resolve() in excludes:
                continue
            yield path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "roots", nargs="*", type=Path, default=None,
        help="directories to process (default: app/src)",
    )
    ap.add_argument(
        "--check", action="store_true",
        help="do not modify; exit non-zero if any file lacks the header.",
    )
    ap.add_argument(
        "--exclude", action="append", default=[], metavar="PATH",
        help="path to skip, relative to app/ or absolute. Repeatable.",
    )
    args = ap.parse_args()

    roots = args.roots or DEFAULT_ROOTS
    roots = [r if r.is_absolute() else (APP_DIR / r) for r in roots]
    excludes = {
        (Path(e) if Path(e).is_absolute() else APP_DIR / e).resolve()
        for e in args.exclude
    }

    paths = list(iter_sources(roots, excludes))
    if not paths:
        print("no Python sources found", file=sys.stderr)
        return 1

    if args.check:
        missing = [p for p in paths if not has_header(p.read_text(encoding="utf-8"))]
        if missing:
            print(
                f"{len(missing)} of {len(paths)} file(s) lack the GPL header:",
                file=sys.stderr,
            )
            for p in missing:
                print(f"  - {p.relative_to(APP_DIR)}", file=sys.stderr)
            print(
                f"\nRun: {sys.executable} {Path(__file__).relative_to(APP_DIR)}",
                file=sys.stderr,
            )
            return 1
        print(f"OK — all {len(paths)} source file(s) carry the GPL header.")
        return 0

    changed = [p for p in paths if add_header(p)]
    for p in changed:
        print(f"added header: {p.relative_to(APP_DIR)}")
    print(f"\n{len(changed)} added, {len(paths) - len(changed)} already had one.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
