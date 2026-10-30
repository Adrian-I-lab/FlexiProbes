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

"""Which code is the engine, and the stamp a resume checkpoint carries.

A resume checkpoint (``.screen_cache`` or ``.candidate_cache``) holds answers
the engine computed. Reading them back is only honest if the same engine would
compute the same answers now. A cache keyed on inputs alone would let a run
resumed after an engine change read back answers the changed engine would not
give, and an engine change need not touch either the version string or
``CACHE_SCHEMA``.

Each checkpoint therefore carries :func:`current_stamp`. On resume a missing or
different stamp raises :class:`~flexprobe.errors.CacheStampMismatchError`,
naming both stamps. The checkpoint is never deleted or silently recomputed:
the user decides whether to resume with the code that wrote it or to clear it.

Only the ENGINE files count. The GUI, writers, logging and docs keep changing
after the freeze, and a change there cannot change a cached answer, so it must
not block a resume. Every module in the package is classified below, and a test
fails when one is added without being classified, so a new engine module can
never escape the stamp by being forgotten.
"""
from __future__ import annotations

import functools
import hashlib
import os
import platform
from pathlib import Path
from typing import Any, Dict, Optional

#: Modules that can change a probe, a verdict, a score, a count, or anything a
#: checkpoint stores. Paths are relative to the ``flexprobe`` package.
ENGINE_FILES = frozenset({
    "constants.py",        # every tunable threshold
    "models.py",           # RunConfig defaults
    "config_io.py",        # absent keys are filled from defaults on load
    "cli.py",              # argument defaults become RunConfig values
    "io_fasta.py",         # input parsing and ambiguity codes
    "tiling.py",           # windows and reverse complement
    "filters.py",          # candidate filter funnel
    "_utils.py",           # shared Tm and GC helpers
    "reference_qc.py",     # which reference records are screened
    "screen.py",           # off-target search engines
    "kmer_index.py",
    "blast.py",
    "screen_parallel.py",
    "specificity.py",      # verdicts
    "consensus.py",        # multi-strain consensus and scores
    "ranking.py",          # composite score
    "isoform_select.py",   # per-gene selection
    "auto_select.py",
    "auto_export.py",
    "chemistry.py",        # oligo assembly
    "pipeline.py",         # orchestrator
    "null_model.py",       # calibration reported in run stats
    "screen_cache.py",     # what a checkpoint stores and how it is read back
    "candidate_cache.py",
    "engine_identity.py",  # this list
})

#: Modules that cannot change a cached answer. Classified explicitly, not by
#: omission, so that forgetting a new module fails a test instead of passing.
NON_ENGINE_FILES = frozenset({
    "__init__.py",         # version string only; stamped as engine_version
    "errors.py",
    "logging_setup.py",
    "orientation.py",      # warns only; never flips or refuses a record
    "summary.py",
    "prescreen_export.py",
    "hpc_export.py",
    "select_export.py",
    "writers/__init__.py",
    "writers/cellranger_probe_set.py",
    "writers/idt_opools.py",
    "writers/manifest.py",
    "writers/results_table.py",
})

#: Set by the container image to its build ID. Unset outside a container,
#: which the stamp records as ``None``.
CONTAINER_BUILD_ID_ENV = "FLEXIPROBES_CONTAINER_BUILD_ID"

PACKAGE_ROOT = Path(__file__).resolve().parent

#: The committed per-file SHA-256 of :data:`ENGINE_FILES`, plus the
#: combined :func:`engine_sha256`. ``tests/test_engine_freeze.py`` fails when
#: an engine file changes without this file changing in the same commit.
FREEZE_MANIFEST = PACKAGE_ROOT / "ENGINE_FREEZE.json"


@functools.lru_cache(maxsize=1)
def engine_sha256() -> str:
    """SHA-256 over the source of every file in :data:`ENGINE_FILES`.

    Content, not a version string or a git commit: an engine change need not
    change the version string, and a commit is absent in an installed wheel or
    a container. Each file contributes its relative path and its bytes, in sorted
    order. A frozen build that ships no ``.py`` source returns
    ``"unavailable"``, and ``engine_version`` is then the only identity left.
    """
    h = hashlib.sha256()
    found = 0
    for rel in sorted(ENGINE_FILES):
        p = PACKAGE_ROOT / rel
        try:
            data = p.read_bytes()
        except OSError:
            continue
        found += 1
        h.update(rel.encode())
        h.update(b"\0")
        h.update(data)
        h.update(b"\0")
    return h.hexdigest() if found else "unavailable"


def _tool_versions(blastn: bool) -> Dict[str, str]:
    """Versions of the external code a cached answer depends on.

    ``blastn`` only when the BLAST engine produced the checkpoint: a k-mer run
    never invokes it, and stamping a binary that never ran would make a k-mer
    resume depend on whether blastn happens to be on PATH.
    """
    versions: Dict[str, str] = {"python": platform.python_version()}
    try:
        import primer3
        pv = getattr(primer3, "__version__", None)
        if not pv:
            from importlib.metadata import version as _pkg_version
            pv = _pkg_version("primer3-py")
        versions["primer3"] = str(pv)
    except Exception:
        versions["primer3"] = "unknown"
    if blastn:
        from .blast import blastn_version
        versions["blastn"] = blastn_version() or "unknown"
    return versions


def current_stamp(*, blastn: bool = False) -> Dict[str, Any]:
    """The stamp a checkpoint written now would carry."""
    from . import __version__

    build_id: Optional[str] = os.environ.get(CONTAINER_BUILD_ID_ENV) or None
    return {
        "engine_version": __version__,
        "engine_sha256": engine_sha256(),
        "container_build_id": build_id,
        "tool_versions": _tool_versions(blastn),
    }


def write_freeze_manifest(engine_version: Optional[str] = None) -> Path:
    """Write :data:`FREEZE_MANIFEST` from the engine files as they are now.

    Run only to freeze, or to unfreeze on purpose together with an engine
    version bump. ``engine_version`` defaults to the package version.
    """
    import json

    if engine_version is None:
        from . import __version__ as engine_version
    files = {rel: hashlib.sha256((PACKAGE_ROOT / rel).read_bytes()).hexdigest()
             for rel in sorted(ENGINE_FILES)}
    engine_sha256.cache_clear()
    data = {"engine_version": engine_version,
            "engine_sha256": engine_sha256(),
            "files": files}
    with open(FREEZE_MANIFEST, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return FREEZE_MANIFEST
