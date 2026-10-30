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

"""Bring run folders into the run list, and find them again after they move.

Two user actions live here. Both are MANUAL: nothing in this module runs at
start-up or on Refresh. The home window calls only :func:`is_folder_missing`
per displayed row, which is a single ``os.path.isdir``.

* **Import runs…** scans a picked folder, and up to
  :data:`flexprobe.hpc_export.RUN_SCAN_DEPTH` levels below it, for finished run
  folders and registers each one.
* **Locate…** points a record whose output folder has gone at the folder's new
  place, after checking that the folder really is that run.

The rules, stated once
----------------------

Badge. Every imported record has ``stats["imported"] = True`` and
``stats["import_source"]`` of ``"hpc"`` or ``"desktop"``. "Import HPC run…"
always records ``"hpc"``. "Import runs…" records ``"hpc"`` only when
:func:`flexprobe.hpc_export.cluster_evidence` finds cluster-only evidence in
the folder, and ``"desktop"`` otherwise. A record with ``imported`` and no
``import_source`` predates the key; only the HPC dialog made those, so it reads
as ``"hpc"``.

Duplicate. A found folder is already in the list when either (a) its resolved
path equals the resolved ``output_dir`` of an existing record, or (b) its
manifest ``timestamp`` is non-empty and it and ``inputs.fasta_sha256`` both
equal an existing record's. (b) catches a copy of a run in a second place. When
(b) matches a record whose own folder is missing, the import relinks that
record to the found folder instead of adding a second row for the same run.

Locate. A folder is accepted for a record only if it holds ``results.csv``, and

* when the record has a stored manifest timestamp: the folder's manifest.json
  has the same timestamp and, where both sides record one, the same FASTA
  SHA-256. Anything else is refused;
* when it has none (saved before the identity was stored): the folder has the
  same name as the old one, and the user is asked to confirm, because a name
  is not proof. A different name is refused.
"""
from __future__ import annotations

import datetime
import json
import os
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

from flexprobe.hpc_export import (
    RESULTS_FILE, RUN_SCAN_DEPTH, describe_run_folder, find_run_folders,
)
from flexprobe.models import RunConfig

from ._nq_build import (
    RUN_MODE_AUTO, RUN_MODE_MANUAL, advanced_summary, specificity_summary,
)
from .modality import MODALITY_PRESETS
from .run_history import RunHistory, RunRecord
from .run_manager import load_run_config

SOURCE_HPC = "hpc"
SOURCE_DESKTOP = "desktop"

MANIFEST_FILE = "manifest.json"

# Locate verdicts.
ACCEPT = "accept"
CONFIRM = "confirm"
REJECT = "reject"

#: config_summary keys the draft map reads that an import deliberately leaves
#: out. ``multi_count`` is the New-query spin box value, which run_config.json
#: does not record; without it apply_draft counts the recorded barcodes, which
#: is the fallback that exists for exactly such records.
SUMMARY_KEYS_NOT_FROM_CONFIG = frozenset({"multi_count"})


# ---------------------------------------------------------------------------
# Small predicates the home table uses
# ---------------------------------------------------------------------------

def import_source(rec: RunRecord) -> str:
    """``"hpc"`` or ``"desktop"`` for an imported run, ``""`` for one made here."""
    stats = rec.stats or {}
    if not stats.get("imported"):
        return ""
    return str(stats.get("import_source") or SOURCE_HPC)


def is_folder_missing(rec: RunRecord) -> bool:
    """True when the record names an output folder that is no longer there.

    A record with no ``output_dir`` at all has nothing to locate, so it is not
    reported as missing.
    """
    return bool(rec.output_dir) and not os.path.isdir(rec.output_dir)


def nearest_existing_parent(path: str) -> str:
    """The closest ancestor of ``path`` that still exists, or ``""``.

    Where a folder picker for a moved run should open: usually the old parent,
    one level up from the folder that went.
    """
    cur = os.path.dirname(os.path.abspath(path)) if path else ""
    while cur:
        if os.path.isdir(cur):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    return ""


# ---------------------------------------------------------------------------
# Identity
# ---------------------------------------------------------------------------

def read_identity(output_dir: str) -> Tuple[str, str]:
    """``(manifest timestamp, input FASTA SHA-256)``; blanks when unreadable."""
    try:
        with open(os.path.join(output_dir, MANIFEST_FILE),
                  encoding="utf-8") as fh:
            man = json.load(fh)
    except (OSError, ValueError):
        return "", ""
    if not isinstance(man, dict):
        return "", ""
    inputs = man.get("inputs") or {}
    sha = inputs.get("fasta_sha256") if isinstance(inputs, dict) else ""
    ts = man.get("timestamp")
    return (ts if isinstance(ts, str) else ""), (sha if isinstance(sha, str) else "")


def stamp_identity(rec: RunRecord) -> bool:
    """Copy the run's manifest identity onto ``rec`` if it has none yet.

    Returns True when the record changed (the caller persists it). Cheap when
    already stamped: no file is opened.
    """
    if rec.manifest_timestamp or not rec.output_dir \
            or not os.path.isdir(rec.output_dir):
        return False
    ts, sha = read_identity(rec.output_dir)
    if not ts:
        return False
    rec.manifest_timestamp = ts
    rec.fasta_sha256 = sha
    return True


def _record_identity(rec: RunRecord) -> Tuple[str, str]:
    if rec.manifest_timestamp:
        return rec.manifest_timestamp, rec.fasta_sha256
    if rec.output_dir and os.path.isdir(rec.output_dir):
        return read_identity(rec.output_dir)
    return "", ""


def _path_key(path: str) -> str:
    return os.path.normcase(os.path.realpath(path)) if path else ""


def find_duplicate(
    records: Iterable[RunRecord],
    output_dir: str,
    timestamp: str,
    fasta_sha256: str,
    identity: Optional[Callable[[RunRecord], Tuple[str, str]]] = None,
) -> Optional[RunRecord]:
    """The existing record this folder duplicates, or None. See module docs."""
    identity = identity or _record_identity
    records = list(records)
    key = _path_key(output_dir)
    for r in records:
        if r.output_dir and _path_key(r.output_dir) == key:
            return r
    if timestamp:
        for r in records:
            if identity(r) == (timestamp, fasta_sha256):
                return r
    return None


# ---------------------------------------------------------------------------
# Building a record from a folder
# ---------------------------------------------------------------------------

def summary_from_config(cfg: RunConfig) -> Dict[str, Any]:
    """``RunRecord.config_summary`` rebuilt from a run's own run_config.json.

    Same keys and value conventions that ``BuildMixin.build_runs`` writes, so
    the existing Re-run path (``HomeWindow._record_to_draft`` ->
    ``summary_to_draft``) restores an imported run the way it restores one made
    here. Without it an imported run re-ran at the form defaults with no FASTA.
    ``tests/test_gui_run_import.py`` fails when the draft map reads a key this
    does not produce (other than :data:`SUMMARY_KEYS_NOT_FROM_CONFIG`).
    """
    f = cfg.filters
    multi = bool(cfg.multi_strain)
    strategy = cfg.consensus_strategy if multi else ""
    return {
        "tile_size": cfg.tile_length,
        "tile_step": cfg.tile_step,
        "gc_min": f.gc_min,
        "gc_max": f.gc_max,
        "enforce_gc_range": f.enforce_gc_range,
        "tm_min": f.tm_min,
        "max_n_per_half": f.max_n_per_half,
        "max_homopolymer_run": f.max_homopolymer_run,
        "require_tn_junction": f.require_tn_junction,
        "tn_junction_mode": f.tn_junction_mode,
        "screen_workers": cfg.screen_workers,
        "workers": cfg.workers,
        "screen_engine": cfg.screen.engine,
        "skip_screen": cfg.screen.skip_screen,
        "primer3": cfg.primer3_scoring,
        "group_by_gene": cfg.group_by_gene,
        "min_isoform_coverage": cfg.min_isoform_coverage,
        "isoform_min_count": cfg.isoform_min_count,
        "screen_top_per_gene": cfg.screen_top_per_gene,
        "extra_on_target_ids": list(cfg.extra_on_target_ids),
        "ranking_weights": (dict(cfg.ranking_weights)
                            if cfg.ranking_weights is not None else None),
        "multiplex": cfg.multiplex,
        "probe_barcode": cfg.probe_barcodes[0] if cfg.probe_barcodes else "",
        "probe_barcodes": list(cfg.probe_barcodes),
        "auto_select_per_gene": cfg.auto_select_per_gene,
        "run_mode": (RUN_MODE_AUTO if cfg.auto_select_per_gene > 0
                     else RUN_MODE_MANUAL),
        "multi_strain": multi,
        "consensus_strategy": strategy,
        "consensus_reference": cfg.consensus_reference,
        "conserved_max_degenerate_per_half": (
            cfg.conserved_max_degenerate_per_half
            if multi and strategy == "conserved" else None
        ),
        "strain_unique_min_mismatches": (
            cfg.strain_unique_min_mismatches if multi else None),
        "strain_unique_per_half": (
            cfg.strain_unique_per_half if multi else None),
        **specificity_summary(cfg.specificity, cfg.screen),
        **advanced_summary(cfg),
    }


def _gui_modality_label(cfg: Optional[RunConfig], fallback: str) -> str:
    """The New-query preset name for the run's modality, else ``fallback``.

    Re-run restores the modality combo by exact text, so the flat manifest
    label (``flex_v2``) restores nothing; the preset name does.
    """
    if cfg is not None:
        m = cfg.modality
        for label, preset in MODALITY_PRESETS.items():
            if (m.name, m.version, m.spatial) == \
                    (preset.name, preset.version, preset.spatial):
                return label
    return fallback


def _local_iso(timestamp: str) -> Optional[datetime.datetime]:
    """A manifest timestamp (UTC, with offset) as naive local time."""
    try:
        dt = datetime.datetime.fromisoformat(timestamp)
    except (TypeError, ValueError):
        return None
    if dt.tzinfo is not None:
        dt = dt.astimezone().replace(tzinfo=None)
    return dt.replace(microsecond=0)


def _manifest_fasta_path(output_dir: str, info: dict) -> str:
    """The manifest's input FASTA as a path the GUI can open.

    Since G16 the manifest names it relative to the run folder, so a relative
    value is joined to ``output_dir``. An absolute one (older runs) is kept.
    """
    fasta = str(info.get("fasta_path") or "")
    if fasta and not os.path.isabs(fasta):
        fasta = os.path.normpath(os.path.join(output_dir, fasta))
    return fasta


def record_from_folder(output_dir: str, *, source: str, name: str = "",
                       info: Optional[dict] = None) -> RunRecord:
    """A completed-run record filled from the folder's own files.

    * name: ``name``, else the folder name (or the FASTA name for an archive's
      random extraction folder);
    * FASTA path and design settings: run_config.json, so Re-run is the run;
    * counts: design_summary.json;
    * completed / started: the manifest timestamp, and that minus the
      manifest's total elapsed seconds. Approximate by the export stage's
      length, since the manifest is written during export;
    * app version: the manifest, kept on ``stats["app_version"]``.
    """
    info = info if info is not None else describe_run_folder(output_dir)
    out = info["output_dir"]
    cfg = load_run_config(out) if info.get("has_run_config") else None
    version = str(info.get("app_version") or "")
    stats = dict(info.get("stats") or {})
    stats["imported"] = True
    stats["import_source"] = source
    if version:
        stats["app_version"] = version
    if source == SOURCE_HPC:
        stage = "imported from HPC"
    else:
        stage = "imported" + (f" (made with FlexiProbes {version})"
                              if version else "")
    rec = RunRecord(
        name=name or info.get("name") or os.path.basename(out),
        modality_label=_gui_modality_label(cfg, info.get("modality_label", "")),
        fasta_path=(cfg.fasta_path if cfg is not None and cfg.fasta_path
                    else _manifest_fasta_path(out, info)),
        output_dir=out,
        status=info.get("status", "ready_for_review"),
        progress=100,
        stage=stage,
        stats=stats,
        config_summary=summary_from_config(cfg) if cfg is not None else {},
        manifest_timestamp=str(info.get("manifest_timestamp") or ""),
        fasta_sha256=str(info.get("fasta_sha256") or ""),
    )
    completed = _local_iso(rec.manifest_timestamp)
    rec.completed_at = (completed.isoformat() if completed is not None
                        else rec.created_at)
    elapsed = info.get("elapsed_seconds")
    if completed is not None and isinstance(elapsed, (int, float)) \
            and elapsed >= 0:
        started = completed - datetime.timedelta(seconds=round(elapsed))
        rec.started_at = started.isoformat()
    return rec


def relink(history: RunHistory, rec: RunRecord, folder: str) -> None:
    """Point ``rec`` at ``folder`` and persist. The caller has matched it."""
    rec.output_dir = os.path.abspath(folder)
    stamp_identity(rec)
    history.update(rec)


# ---------------------------------------------------------------------------
# Import runs…
# ---------------------------------------------------------------------------

@dataclass
class ImportSummary:
    """What one "Import runs…" did, for the summary dialog and the tests."""

    root: str
    imported: List[RunRecord] = field(default_factory=list)
    #: (folder, name of the record already holding it)
    already_present: List[Tuple[str, str]] = field(default_factory=list)
    #: (folder, name): a missing record found at a new place and relinked
    relinked: List[Tuple[str, str]] = field(default_factory=list)
    #: (folder, reason)
    skipped: List[Tuple[str, str]] = field(default_factory=list)

    def message(self) -> str:
        n_new = len(self.imported)
        lines = [f"Imported {n_new} run{'s' if n_new != 1 else ''} "
                 f"from {self.root}."]
        if self.relinked:
            lines.append(f"{len(self.relinked)} run(s) already in the list had "
                         f"a missing folder and now point here.")
        if self.already_present:
            lines.append(f"{len(self.already_present)} already in the list, "
                         f"not added again.")
        if self.skipped:
            lines.append(f"{len(self.skipped)} skipped (see details).")
        if not (n_new or self.relinked or self.already_present
                or self.skipped):
            lines.append(
                f"No run folders were found. A run folder holds results.csv "
                f"and manifest.json, and is looked for up to {RUN_SCAN_DEPTH} "
                f"levels below the chosen folder.")
        return "\n".join(lines)

    def details(self) -> str:
        out: List[str] = []
        for title, rows in (
            ("Imported", [(r.output_dir, r.name) for r in self.imported]),
            ("Relinked", self.relinked),
            ("Already in the list", self.already_present),
            ("Skipped", self.skipped),
        ):
            if rows:
                out.append(f"{title}:")
                out.extend(f"  {a}  ({b})" for a, b in rows)
        return "\n".join(out)


def import_runs(root: str, history: RunHistory,
                max_depth: int = RUN_SCAN_DEPTH) -> ImportSummary:
    """Register every finished run folder at or below ``root``."""
    summary = ImportSummary(root=root)
    if not os.path.isdir(root):
        summary.skipped.append((root, "not a folder"))
        return summary
    runs, bare = find_run_folders(root, max_depth)
    for path in bare:
        summary.skipped.append((
            path, "results.csv but no manifest.json, so it cannot be "
                  "identified; use Import HPC run… for a bare results file"))

    records = history.list_all()
    # Each existing record's identity is read at most once per import, however
    # many folders are compared against it.
    cache: Dict[str, Tuple[str, str]] = {}

    def identity(r: RunRecord) -> Tuple[str, str]:
        if r.run_id not in cache:
            cache[r.run_id] = _record_identity(r)
        return cache[r.run_id]

    for path in runs:
        info = describe_run_folder(path)
        if not info.get("manifest_readable"):
            summary.skipped.append((path, "manifest.json could not be read"))
            continue
        ts = str(info.get("manifest_timestamp") or "")
        sha = str(info.get("fasta_sha256") or "")
        dup = find_duplicate(records, path, ts, sha, identity)
        if dup is not None:
            if is_folder_missing(dup):
                relink(history, dup, path)
                cache.pop(dup.run_id, None)
                summary.relinked.append((path, dup.name))
            else:
                summary.already_present.append((path, dup.name))
            continue
        source = SOURCE_HPC if info.get("cluster_evidence") else SOURCE_DESKTOP
        rec = record_from_folder(path, source=source, info=info)
        history.append(rec)
        records.append(rec)
        summary.imported.append(rec)
    return summary


# ---------------------------------------------------------------------------
# Locate…
# ---------------------------------------------------------------------------

@dataclass
class LocateMatch:
    verdict: str            # ACCEPT | CONFIRM | REJECT
    reason: str = ""


def _same_name(a: str, b: str) -> bool:
    return os.path.normcase(os.path.basename(os.path.normpath(a))) == \
        os.path.normcase(os.path.basename(os.path.normpath(b)))


def match_folder(rec: RunRecord, folder: str) -> LocateMatch:
    """Whether ``folder`` is the moved output folder of ``rec``. See module docs."""
    if not folder or not os.path.isdir(folder):
        return LocateMatch(REJECT, f"{folder or 'That'} is not a folder.")
    if not os.path.isfile(os.path.join(folder, RESULTS_FILE)):
        return LocateMatch(
            REJECT, f"There is no {RESULTS_FILE} in {folder}, so it is not a "
                    f"finished run folder.")
    if rec.manifest_timestamp:
        ts, sha = read_identity(folder)
        if not ts:
            return LocateMatch(
                REJECT, f"{folder} has no readable {MANIFEST_FILE}, so it "
                        f"cannot be checked against this run.")
        if ts != rec.manifest_timestamp:
            return LocateMatch(
                REJECT, f"That folder holds a different run: its "
                        f"{MANIFEST_FILE} was written at {ts}, and this run's "
                        f"at {rec.manifest_timestamp}.")
        if rec.fasta_sha256 and sha and sha != rec.fasta_sha256:
            return LocateMatch(
                REJECT, "That folder's manifest names a different input FASTA "
                        "(its SHA-256 does not match this run's).")
        return LocateMatch(ACCEPT)
    old = os.path.basename(os.path.normpath(rec.output_dir))
    new = os.path.basename(os.path.normpath(folder))
    if not _same_name(rec.output_dir, folder):
        return LocateMatch(
            REJECT, f"This run was recorded before FlexiProbes kept a "
                    f"fingerprint of its manifest, so a folder can only be "
                    f"matched by name, and “{new}” is not “{old}”.")
    return LocateMatch(
        CONFIRM, f"This run was recorded before FlexiProbes kept a fingerprint "
                 f"of its manifest, so “{new}” can only be matched by its name "
                 f"and its {RESULTS_FILE}, which does not prove it is the same "
                 f"run.")


def _moved_roots(old_dir: str, new_dir: str) -> Tuple[str, str]:
    """The two roots left after stripping the trailing components both share.

    ``/a/b/proj/run`` -> ``/x/proj/run`` gives ``(/a/b, /x)``. With nothing in
    common (the folder itself was renamed) the parents are returned.
    """
    old_parts = os.path.normpath(os.path.abspath(old_dir)).split(os.sep)
    new_parts = os.path.normpath(os.path.abspath(new_dir)).split(os.sep)
    n = 0
    while (n < min(len(old_parts), len(new_parts)) - 1
           and os.path.normcase(old_parts[-1 - n])
           == os.path.normcase(new_parts[-1 - n])):
        n += 1
    n = max(n, 1)

    def _join(parts: List[str]) -> str:
        root = os.sep.join(parts) or os.sep
        return root + os.sep if root.endswith(":") else root

    return _join(old_parts[:-n]), _join(new_parts[:-n])


def sibling_candidates(
    history: RunHistory, old_dir: str, new_dir: str, exclude_run_id: str,
) -> List[Tuple[RunRecord, str, LocateMatch]]:
    """Other missing runs that appear to have moved along with ``old_dir``.

    For each missing record two places are tried: ``<new parent>/<its folder
    name>`` (runs that sat side by side and still do), and the same relative
    path under the moved root (``<old root>/other/run2`` becoming
    ``<new root>/other/run2``, which covers a moved ``<base>/<project>/<run>``
    tree). Every candidate must pass :func:`match_folder`; refusals are dropped
    and CONFIRM matches are returned for the caller to show as unverified.
    """
    old_root, new_root = _moved_roots(old_dir, new_dir)
    new_parent = os.path.dirname(os.path.abspath(new_dir))
    found: List[Tuple[RunRecord, str, LocateMatch]] = []
    for r in history.list_completed():
        if r.run_id == exclude_run_id or not is_folder_missing(r):
            continue
        candidates = [os.path.join(
            new_parent, os.path.basename(os.path.normpath(r.output_dir)))]
        try:
            rel = os.path.relpath(r.output_dir, old_root)
        except ValueError:          # a different drive on Windows
            rel = ""
        if rel and not rel.startswith(os.pardir):
            candidates.append(os.path.join(new_root, rel))
        for cand in dict.fromkeys(os.path.normpath(c) for c in candidates):
            if not os.path.isdir(cand):
                continue
            m = match_folder(r, cand)
            if m.verdict != REJECT:
                found.append((r, cand, m))
                break
    return found
