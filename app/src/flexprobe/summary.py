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

"""``design_summary.json`` — the one file that is written after every run.

A run whose filter funnel empties is a RESULT: the tiles existed and every one
was rejected. Without a summary a consumer that keys off the run's output files
cannot tell that apart from a run that never started, and a real candidate set
that was filtered out would be read as "0 tiles".

This module holds the writers so that both entry points, the CLI and the GUI,
write the same summary and manifest when the funnel empties.
"""
from __future__ import annotations

import json
import logging
import os
from typing import Any, Dict, Optional

logger = logging.getLogger("flexprobe.summary")

SUMMARY_NAME = "design_summary.json"


def _write(outdir: str, payload: Dict[str, Any]) -> Optional[str]:
    """Atomically write the summary. Never raises.

    Called on error paths, where a failure to write the explanation must not
    replace the original error.
    """
    try:
        path = os.path.join(outdir, SUMMARY_NAME)
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(payload, fh, indent=2)
        os.replace(tmp, path)
        return path
    except OSError:
        logger.debug("Could not write %s under %s", SUMMARY_NAME, outdir,
                     exc_info=True)
        return None


def write_ok_summary(outdir: str, stats: Dict[str, Any]) -> Optional[str]:
    """Summary for a run that reached the end of the pipeline without raising.

    ``status`` is ``"ok"`` when the run produced at least one probe pair, or
    ``"no_probes"`` when every stage completed cleanly but the funnel still
    emptied before export — typically every remaining candidate failed
    OFF-TARGET SPECIFICITY (or a downstream pairing/export step) rather than
    the input FILTERS. An emptied *filter* funnel is a different, earlier
    failure mode: ``run_pipeline`` raises before screening even starts, and
    that path is caught and reported by :func:`write_failed_summary` as
    ``"no_candidates"``. Both are real results, not non-events — a run that
    designs zero probes is exactly as reportable as one that designs 40 — so
    neither silently shares ``"ok"`` with a run that actually produced
    something. See ``cli.cmd_design`` for the matching non-zero CLI exit code
    (the GUI path, ``_DesignWorker.run``, has no process exit code to flip,
    but reads this same ``status`` for its own reporting).
    """
    pairs_output = int((stats or {}).get("pairs_output", 0) or 0)
    status = "ok" if pairs_output > 0 else "no_probes"
    return _write(outdir, {"status": status, "stats": stats})


def write_prepared_summary(outdir: str, stats: Dict[str, Any]) -> Optional[str]:
    """Summary for a ``--prepare-only`` run: status ``"prepared"``.

    A deliberately DISTINCT status, not ``"ok"`` with a flag buried in
    ``stats``. Every consumer of this file keys on ``status``, so a prepared
    output directory has to be unreadable as a finished one by the simplest
    check anybody writes. It designed no probes and has no results.csv, and it
    is not a failure either: the run did exactly what it was asked to do, and
    the screen job has yet to run.

    The counterpart on disk is ``constants.PREPARED_MARKER``, written by the
    pipeline beside the candidate checkpoint it certifies. This file says what
    happened; that file is what a shell script tests for.
    """
    return _write(outdir, {"status": "prepared", "stats": stats})


def write_failed_summary(outdir: str, exc: BaseException) -> Optional[str]:
    """Summary for a run that did not reach the end.

    ``status`` is one of:

    * ``no_candidates`` — tiles were generated and the filters rejected all of
      them. ``stats`` carries the real ``tiles_total`` and the per-reason
      rejection counts lifted from ``filter_report.json`` (which the pipeline
      writes before it aborts), so the funnel that produced the empty result is
      legible without reopening the run.
    * ``cancelled`` — the user stopped it.
    * ``failed`` — anything else (bad input, BLAST failure).
    """
    from .errors import CancellationError

    status = "failed"
    stats: Dict[str, Any] = {}
    counts: Dict[str, Any] = {}
    try:
        with open(os.path.join(outdir, "filter_report.json")) as fh:
            counts = json.load(fh).get("counts", {}) or {}
    except (OSError, ValueError):
        counts = {}
    if counts:
        total = int(counts.get("total", 0))
        passed = int(counts.get("passed", 0))
        if total > 0 and passed == 0:
            status = "no_candidates"
        stats = {
            "tiles_total": total,
            # The SAME key names a COMPLETED run writes. If the abort path
            # emitted different keys, a failed run would produce a
            # design_summary.json that is present, looks ordinary and silently
            # has a different schema, so a consumer keyed on the normal names
            # would raise KeyError instead of reading "this run produced
            # nothing". A summary whose shape depends on whether the run
            # succeeded is not a summary. `status` is what says what happened,
            # not the presence or absence of a field.
            "filter_passed": passed,
            "specificity_passed": 0,
            "candidates_passed_filters": passed,   # legacy alias, kept readable
            "pairs_output": 0,
            "filter_stats": counts,
        }
    if isinstance(exc, CancellationError):
        status = "cancelled"
    return _write(outdir, {
        "status": status,
        "error": str(exc),
        "error_type": type(exc).__name__,
        "stats": stats,
    })
