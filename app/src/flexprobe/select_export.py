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

"""Write the IDT oPool + Cell Ranger probe-set CSVs from user selections.

Used after a discovery-mode run: the pipeline has already emitted
``results.csv`` (every candidate, ranked) and ``manifest.json``. The user
reviews candidates in the GUI and picks a final subset per gene; this
module turns that selection into orderable order CSVs.
"""
from __future__ import annotations

import csv
import json
import os
import tempfile
import threading
from collections import OrderedDict
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from .chemistry import BARCODE_PLACEHOLDER, FLEX_V1_PROBE_BARCODES
from .constants import PROBE_ROLE_CONSENSUS, PROBE_ROLE_VARIANT
from .errors import WriterError
from .writers.cellranger_probe_set import (
    cr_header_lines, default_panel_identity, warn_if_transcript_like_gene_ids,
)


SELECTIONS_FILE = "selections.json"
RESULTS_FILE = "results.csv"
MANIFEST_FILE = "manifest.json"


# ---------------------------------------------------------------------------
# Selection persistence
# ---------------------------------------------------------------------------

def selections_path(output_dir: str) -> str:
    return os.path.join(output_dir, SELECTIONS_FILE)


def load_selections(output_dir: str) -> Dict[str, List[str]]:
    """Return ``{gene_id: [probe_id, ...]}`` from disk, or ``{}`` if absent."""
    path = selections_path(output_dir)
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as fh:
            doc = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}
    raw = doc.get("selections", {}) if isinstance(doc, dict) else {}
    return {
        str(g): [str(p) for p in (ids or [])]
        for g, ids in raw.items()
    }


def save_selections(
    output_dir: str,
    selections: Dict[str, List[str]],
    *,
    note: Optional[str] = None,
) -> str:
    """Persist selections + timestamp as ``selections.json``."""
    os.makedirs(output_dir, exist_ok=True)
    path = selections_path(output_dir)
    doc = {
        "selections": {g: list(ids) for g, ids in selections.items()},
        "saved_at": datetime.now().isoformat(timespec="seconds"),
    }
    if note:
        doc["note"] = note
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(doc, fh, indent=2, sort_keys=True)
    os.replace(tmp, path)
    return path


# ---------------------------------------------------------------------------
# Inspect a discovery run
# ---------------------------------------------------------------------------

class ResultsLoadCancelled(Exception):
    """Raised by :func:`load_results_rows` when ``should_cancel()`` goes true.

    A distinct type so a caller that offered the user a Cancel button can tell
    "the user stopped this" apart from "the file is broken" and stay silent
    about the former.
    """


# ---------------------------------------------------------------------------
# results.csv cache
#
# BOUNDED, because a whole-transcriptome results.csv can be over a gigabyte.
# With an unbounded cache every run the user opened in one session would stay
# resident until the process exited.
#
# Eviction policy
# ---------------
#   * Keyed by absolute path; the stored value is (mtime, rows), so a file that
#     changed on disk is a miss even when its path is still resident.
#   * Strict LRU over an OrderedDict: a hit moves its key to the end, and
#     eviction pops from the FRONT (least recently used).
#   * TWO independent caps, both enforced after every insert, evicting until
#     BOTH hold:
#       - _RESULTS_CACHE_MAX_ENTRIES runs, and
#       - _RESULTS_CACHE_MAX_BYTES of *source CSV bytes* summed over residents.
#   * A single file larger than the byte budget is NEVER cached: admitting it
#     would evict every other resident and still blow the budget. Such a run is
#     re-read from disk each time it is opened, which is the deliberate trade —
#     the cache exists to make small and medium runs instant, not to hold a
#     whole-transcriptome results.csv in RAM.
#
# The budget counts bytes ON DISK, not the retained size of the parsed rows: a
# list of dict-of-str inflates the source by roughly 5-10x, so the caps below
# are set well under what the machine can hold. Raise them only with that
# multiplier in mind.
#
# Guarded by a lock: the GUI parses results.csv on a worker QThread while
# the UI thread is live, so two threads really can reach this dict at once.
# ---------------------------------------------------------------------------
_RESULTS_CACHE_MAX_ENTRIES = 4
_RESULTS_CACHE_MAX_BYTES = 64 * 1024 * 1024

# How often the streaming reader stops to look at the cancel flag and push a
# progress update: often enough that Cancel feels immediate on a slow disk,
# rarely enough that the check is noise next to the CSV parse itself.
_CANCEL_POLL_BYTES = 1 << 20

_RESULTS_CACHE: "OrderedDict[str, Tuple[float, int, List[Dict[str, str]]]]" = (
    OrderedDict()
)
_RESULTS_CACHE_LOCK = threading.Lock()


def _cache_get(path: str, mtime: float) -> Optional[List[Dict[str, str]]]:
    with _RESULTS_CACHE_LOCK:
        entry = _RESULTS_CACHE.get(path)
        if entry is None or entry[0] != mtime:
            return None
        _RESULTS_CACHE.move_to_end(path)
        return entry[2]


def _cache_put(path: str, mtime: float, size: int,
               rows: List[Dict[str, str]]) -> None:
    if size > _RESULTS_CACHE_MAX_BYTES:
        # Too big to admit — see the eviction policy above.
        return
    with _RESULTS_CACHE_LOCK:
        _RESULTS_CACHE[path] = (mtime, size, rows)
        _RESULTS_CACHE.move_to_end(path)
        while _RESULTS_CACHE and (
            len(_RESULTS_CACHE) > _RESULTS_CACHE_MAX_ENTRIES
            or sum(e[1] for e in _RESULTS_CACHE.values())
            > _RESULTS_CACHE_MAX_BYTES
        ):
            if len(_RESULTS_CACHE) == 1:
                break        # the sole resident is within budget by admission
            _RESULTS_CACHE.popitem(last=False)


def clear_results_cache() -> None:
    """Drop every cached results.csv. For tests and for freeing memory."""
    with _RESULTS_CACHE_LOCK:
        _RESULTS_CACHE.clear()


def results_cache_paths() -> List[str]:
    """Resident cache keys, least-recently-used first. For tests."""
    with _RESULTS_CACHE_LOCK:
        return list(_RESULTS_CACHE)


def load_results_rows(
    output_dir: str,
    *,
    progress: Optional[Any] = None,
    should_cancel: Optional[Any] = None,
) -> List[Dict[str, str]]:
    """Load every row of ``results.csv`` as a dict, cached by file mtime.

    ``progress`` is called as ``progress(pct)`` with an int 0-100 as the file is
    consumed, and ``should_cancel()`` is polled between chunks; when it returns
    true this raises :class:`ResultsLoadCancelled`. Both are optional and the
    default call is exactly what it always was — the hooks exist so the GUI can
    read a 1.17 GB results.csv on a worker thread with a live progress bar and a
    Cancel that actually stops the read rather than only hiding the dialog.
    """
    path = os.path.join(output_dir, RESULTS_FILE)
    if not os.path.exists(path):
        raise FileNotFoundError(f"results.csv not found in {output_dir}")
    mtime = os.path.getmtime(path)
    size = os.path.getsize(path)
    cached = _cache_get(path, mtime)
    if cached is not None:
        if progress is not None:
            progress(100)
        return cached

    rows: List[Dict[str, str]] = []
    with open(path, newline="") as fh:
        # csv.reader takes any iterable of lines, so wrapping the file lets the
        # read be metered and interrupted WITHOUT buffering the whole file first
        # — fh.tell() is unusable here (Python disables it once __next__ has
        # been called on a text file), so bytes are counted as they go past.
        def _lines():
            seen = 0
            next_check = _CANCEL_POLL_BYTES
            last_pct = -1
            for line in fh:
                seen += len(line)
                if seen >= next_check:
                    next_check = seen + _CANCEL_POLL_BYTES
                    if should_cancel is not None and should_cancel():
                        raise ResultsLoadCancelled(path)
                    if progress is not None and size:
                        pct = min(100, int(100 * seen / size))
                        if pct != last_pct:
                            last_pct = pct
                            progress(pct)
                yield line

        rows.extend(csv.DictReader(_lines()))
    if progress is not None:
        progress(100)
    _cache_put(path, mtime, size, rows)
    return rows


def load_gene_lengths(output_dir: str) -> Dict[str, int]:
    """Pull ``manifest.params.gene_lengths`` from a previous run."""
    path = os.path.join(output_dir, MANIFEST_FILE)
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as fh:
            man = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}
    params = (man.get("params") or {}) if isinstance(man, dict) else {}
    raw = params.get("gene_lengths") or {}
    return {str(g): int(n) for g, n in raw.items() if str(n).isdigit() or isinstance(n, int)}


def load_variant_positions(output_dir: str) -> List[int]:
    """Return the consensus variant positions recorded in the manifest."""
    info = load_modality_info(output_dir)
    raw = info.get("variant_position_list") or []
    out = []
    for p in raw:
        try:
            out.append(int(p))
        except (TypeError, ValueError):
            pass
    return out


def load_modality_info(output_dir: str) -> Dict[str, Any]:
    """Return the manifest's ``params`` dict (modality, multiplex, etc.)."""
    path = os.path.join(output_dir, MANIFEST_FILE)
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as fh:
            man = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}
    return (man.get("params") or {}) if isinstance(man, dict) else {}


def modality_label(params: Dict[str, Any]) -> str:
    """The flat ``"flex_v1"``-style label from a manifest's ``params``.

    Manifests now carry the whole RunConfig, so ``params['modality']`` is the
    structured ModalityDefinition and the flat string moved to
    ``params['modality_label']``. Manifests written before that (including the
    committed validation runs) only have the string under ``params['modality']``
    — read both, newest first.
    """
    label = params.get("modality_label")
    if isinstance(label, str) and label:
        return label
    legacy = params.get("modality")
    if isinstance(legacy, str):
        return legacy
    if isinstance(legacy, dict):        # a manifest with no explicit label
        name = str(legacy.get("name", ""))
        version = str(legacy.get("version", ""))
        return f"{name}_{version}".strip("_")
    return ""


def applied_probe_barcodes(params: Dict[str, Any]) -> List[str]:
    """The probe barcodes the run ACTUALLY applied.

    ``params['probe_barcodes']`` is now the request recorded from the config;
    chemistry can drop it entirely (Flex v2 has no probe barcode), so the
    applied list is recorded separately. Older manifests only carry the applied
    list, under the ``probe_barcodes`` key — read both, newest first.
    """
    for key in ("probe_barcodes_applied", "probe_barcodes"):
        val = params.get(key)
        if isinstance(val, list):
            return [str(v) for v in val]
    return []


# ---------------------------------------------------------------------------
# Export selected probes to the order-format CSVs
# ---------------------------------------------------------------------------

def _atomic_write(path: str, write_fn, overwrite: bool) -> str:
    if os.path.exists(path) and not overwrite:
        raise WriterError(f"Output exists and overwrite=False: {path}")
    dir_name = os.path.dirname(path) or "."
    os.makedirs(dir_name, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=dir_name, suffix=".tmp")
    try:
        with os.fdopen(fd, "w", newline="") as fh:
            write_fn(fh)
        os.replace(tmp, path)
    except Exception:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    return path


def _select_rows(
    rows: List[Dict[str, str]],
    selections: Dict[str, List[str]],
) -> List[Dict[str, str]]:
    """Filter ``rows`` to those matching ``selections``, preserving the
    selection's per-gene order so the user's intent survives into the CSVs.

    A degenerate CONSENSUS row is never orderable — its sequence is a mixed-base
    string that would land on an IDT order sheet as a mixed-base oligo and in the
    Cell Ranger probe set as an invalid sequence. Selecting one is nonetheless
    the natural thing for a reviewer to do, because the parent is what identifies
    the SITE. So a selected parent is replaced by its variants: the concrete ACGT
    oligos that are what you actually order for that site, all of them, since
    they go into the same pool. Nothing that is not pure ACGT leaves here.
    """
    by_key = {(r.get("gene_id", ""), r.get("probe_id", "")): r for r in rows}
    children: Dict[str, List[Dict[str, str]]] = {}
    for r in rows:
        parent = r.get("parent_probe_id") or ""
        if parent and r.get("probe_role") == PROBE_ROLE_VARIANT:
            children.setdefault(parent, []).append(r)
    selected: List[Dict[str, str]] = []
    for gene_id, ids in selections.items():
        for pid in ids:
            row = by_key.get((gene_id, pid))
            if row is None:
                continue
            if row.get("probe_role") == PROBE_ROLE_CONSENSUS:
                selected.extend(children.get(pid, []))
            else:
                selected.append(row)
    return selected


def _resolve_barcode_sequences(labels: List[str]) -> List[Tuple[str, str]]:
    """Return ``[(label, 8-mer), ...]`` for each user-selected BC label.

    Accepts either the canonical ``BC003``-style label or a raw 8-mer (in
    which case the label is shown verbatim).
    """
    pairs: List[Tuple[str, str]] = []
    for raw in labels:
        if not raw:
            continue
        if raw in FLEX_V1_PROBE_BARCODES:
            pairs.append((raw, FLEX_V1_PROBE_BARCODES[raw]))
        else:
            pairs.append((raw, raw))   # assume user passed a raw 8-mer
    return pairs


def export_selected_probes(
    output_dir: str,
    selections: Dict[str, List[str]],
    *,
    overwrite: bool = True,
    write_idt: bool = True,
    write_cellranger: bool = True,
) -> Dict[str, str]:
    """Read ``results.csv`` in ``output_dir`` and write order CSVs.

    Returns a ``{kind: path}`` dict for the files written (``idt`` and
    ``cellranger``). ``write_idt`` / ``write_cellranger`` are the run's
    ``ExportSettings`` switches (``--no-idt`` / ``--no-cellranger``), and the
    auto-select path honours them too.

    Multi-barcode Flex v1 (multiplex): when ``manifest.params.probe_barcodes``
    lists more than one barcode, each selected probe is repeated for every
    barcode/sample. The LHS oligo is shared (no barcode), the RHS placeholder
    ``{BC}`` is substituted with each barcode's 8-mer, and each pool is named
    after its barcode (pool ``BC001`` for sample 1, ``BC008`` for sample 8…).
    """
    rows = load_results_rows(output_dir)
    selected = _select_rows(rows, selections)
    if not selected:
        raise WriterError("No selections supplied (or no matching probes).")
    # Belt and braces on top of _select_rows' substitution. A mixed-base oligo
    # on an order sheet is a real order, of a real (wrong) thing, and money.
    if any(r.get("probe_role") == PROBE_ROLE_CONSENSUS for r in selected):
        raise WriterError(
            "A degenerate consensus probe reached the order export. Consensus "
            "probes are parents, not oligos — order their variants instead."
        )

    paths: Dict[str, str] = {}
    modality = load_modality_info(output_dir)
    cr_version = "v1"
    m_label = modality_label(modality).lower()
    if "v2" in m_label:
        cr_version = "v2"
    # NOTE: chemistry does not change the Cell Ranger CSV's column set. See
    # writers/cellranger_probe_set.py for why `region` is omitted for every
    # chemistry rather than added only for Visium.

    # Multi-barcode expansion (Flex v1). Empty/1-element = single-barcode
    # behaviour (one pool per gene with one LHS + one RHS row per probe).
    barcode_pairs = _resolve_barcode_sequences(applied_probe_barcodes(modality))
    multi_barcode = len(barcode_pairs) > 1

    # IDT oPools — single-barcode = 2 rows per pair; multi-barcode = 1 LHS
    # row + N RHS rows per pair, grouped into per-sample pools.
    def _write_idt(fh):
        w = csv.writer(fh)
        w.writerow(["Pool Name", "Sequence Name", "Sequence"])
        if multi_barcode:
            # One pool per barcode/sample, named after the barcode, containing
            # every selected LHS + the corresponding RHS with that barcode
            # substituted in.
            for bc_label, bc_seq in barcode_pairs:
                pool = bc_label
                # LHS oligos first (shared across barcodes, but each pool
                # needs its own copy for hybridization).
                for r in selected:
                    pid = r.get("probe_id", "")
                    w.writerow([pool, f"{pid}_LHS", r.get("lhs_oligo", "")])
                # Then RHS oligos with the barcode substituted in.
                for r in selected:
                    pid = r.get("probe_id", "")
                    rhs = r.get("rhs_oligo", "").replace(
                        BARCODE_PLACEHOLDER, bc_seq,
                    )
                    w.writerow([pool, f"{pid}_RHS_{bc_label}", rhs])
        else:
            # Single-barcode behaviour: one pool per gene.
            for r in selected:
                gene = r.get("gene_id", "")
                pid = r.get("probe_id", "")
                pool = f"{gene}_pool"
                w.writerow([pool, f"{pid}_LHS", r.get("lhs_oligo", "")])
                w.writerow([pool, f"{pid}_RHS", r.get("rhs_oligo", "")])
    if write_idt:
        paths["idt"] = _atomic_write(
            os.path.join(output_dir, "idt_opools.csv"), _write_idt, overwrite,
        )

    # Cell Ranger — 1 row per pair with the ligated 50-bp probe seq. Header
    # shape (metadata block + included column, no region/gene_name) matches
    # writers/cellranger_probe_set.py — see that module's docstring for why.
    # The FASTA basename recorded in this run's manifest is the closest
    # flexiprobes concept to a 10x panel identity / reference genome.
    fasta_ident = os.path.splitext(os.path.basename(
        modality.get("fasta_path") or "",
    ))[0] or os.path.basename(output_dir.rstrip(os.sep)) or "run"
    panel_name = default_panel_identity(f"flexiprobes custom probe set: {fasta_ident}")

    def _write_cr(fh):
        for line in cr_header_lines(panel_name, fasta_ident):
            fh.write(line + "\n")
        w = csv.writer(fh)
        w.writerow(["gene_id", "probe_seq", "probe_id", "included"])
        gene_ids_written = []
        for idx, r in enumerate(selected):
            gene = r.get("gene_id", "")
            pid = f"{gene}|{cr_version}|{idx:04d}"
            seq = r.get("probe_rc_50bp") or (
                r.get("target_LHS_25bp", "") + r.get("target_RHS_25bp", "")
            )
            gene_ids_written.append(gene)
            w.writerow([gene, seq, pid, "TRUE"])
        warn_if_transcript_like_gene_ids(gene_ids_written)
    if write_cellranger:
        paths["cellranger"] = _atomic_write(
            os.path.join(output_dir, "cellranger_probe_set.csv"), _write_cr,
            overwrite,
        )

    return paths


# ---------------------------------------------------------------------------
# Project-level bulk export
# ---------------------------------------------------------------------------

def export_project_bulk(
    runs: List[Dict[str, Any]],
    output_dir: str,
    *,
    project_name: str = "project",
    overwrite: bool = True,
) -> Dict[str, Any]:
    """Combine selections from many runs into one set of order CSVs.

    ``runs`` is a list of dicts ``{name, output_dir, selections}`` —
    typically each run in a project. Pool names in the combined IDT
    oPools CSV are prefixed with the run name so probes from different
    target organisms / viruses stay distinguishable in the order sheet.

    Returns ``{"paths": {kind: path}, "summary": {...}, "skipped": [...]}``.
    """
    os.makedirs(output_dir, exist_ok=True)

    combined_rows: List[Dict[str, str]] = []
    summary: Dict[str, int] = {}
    skipped: List[str] = []
    cr_version = "v1"

    for run in runs:
        run_name = run.get("name") or "run"
        run_out = run.get("output_dir") or ""
        sel = run.get("selections") or {}
        if not run_out or not os.path.isdir(run_out) or not sel:
            skipped.append(run_name)
            continue
        try:
            rows = load_results_rows(run_out)
        except FileNotFoundError:
            skipped.append(run_name)
            continue
        modality = load_modality_info(run_out)
        m_label = modality_label(modality).lower()
        if "v2" in m_label:
            cr_version = "v2"
        # Multi-barcode expansion happens *per run* — different runs in the
        # same project may have been configured with different barcode sets.
        run_barcodes = _resolve_barcode_sequences(applied_probe_barcodes(modality))
        run_rows = _select_rows(rows, sel)
        emitted = 0
        for r in run_rows:
            if len(run_barcodes) > 1 and BARCODE_PLACEHOLDER in r.get("rhs_oligo", ""):
                for bc_label, bc_seq in run_barcodes:
                    tagged = dict(r)
                    tagged["__run__"] = run_name
                    tagged["__bc_label__"] = bc_label
                    tagged["rhs_oligo"] = r.get("rhs_oligo", "").replace(
                        BARCODE_PLACEHOLDER, bc_seq,
                    )
                    combined_rows.append(tagged)
                    emitted += 1
            else:
                tagged = dict(r)
                tagged["__run__"] = run_name
                tagged["__bc_label__"] = ""
                combined_rows.append(tagged)
                emitted += 1
        summary[run_name] = emitted

    if not combined_rows:
        raise WriterError(
            "Nothing to bulk-export — none of the project's runs had a saved "
            "selection. Open each run, confirm a selection, then re-try."
        )

    def _write_idt(fh):
        w = csv.writer(fh)
        w.writerow(["Pool Name", "Sequence Name", "Sequence"])
        for r in combined_rows:
            run_name = r.get("__run__", "")
            gene = r.get("gene_id", "")
            pid = r.get("probe_id", "")
            bc = r.get("__bc_label__", "")
            # In multi-barcode runs, group probes by barcode (one pool per
            # sample/barcode). In single-barcode runs, one pool per
            # gene. Run name is always prefixed to avoid clashes
            # across viruses in a project.
            if bc:
                pool = f"{run_name}__{bc}"
                rhs_suffix = f"_RHS_{bc}"
            else:
                pool = f"{run_name}__{gene}_pool"
                rhs_suffix = "_RHS"
            w.writerow([pool, f"{run_name}__{pid}_LHS", r.get("lhs_oligo", "")])
            w.writerow([pool, f"{run_name}__{pid}{rhs_suffix}", r.get("rhs_oligo", "")])

    # Header shape matches writers/cellranger_probe_set.py — see that
    # module's docstring. A bulk export spans multiple runs/organisms, so
    # ``project_name`` (not any single run's FASTA) is the panel identity.
    panel_name = default_panel_identity(f"flexiprobes custom probe set: {project_name}")

    def _write_cr(fh):
        for line in cr_header_lines(panel_name, project_name):
            fh.write(line + "\n")
        w = csv.writer(fh)
        w.writerow(["gene_id", "probe_seq", "probe_id", "included"])
        gene_ids_written = []
        for idx, r in enumerate(combined_rows):
            run_name = r.get("__run__", "")
            gene = r.get("gene_id", "")
            pid = f"{run_name}__{gene}|{cr_version}|{idx:04d}"
            seq = r.get("probe_rc_50bp") or (
                r.get("target_LHS_25bp", "") + r.get("target_RHS_25bp", "")
            )
            gene_ids_written.append(gene)
            w.writerow([f"{run_name}__{gene}", seq, pid, "TRUE"])
        warn_if_transcript_like_gene_ids(gene_ids_written)

    paths = {
        "idt": _atomic_write(
            os.path.join(output_dir, f"{project_name}_idt_opools.csv"),
            _write_idt, overwrite,
        ),
        "cellranger": _atomic_write(
            os.path.join(output_dir, f"{project_name}_cellranger_probe_set.csv"),
            _write_cr, overwrite,
        ),
    }

    # Bulk-export manifest summarizing what went in.
    manifest_path = os.path.join(output_dir, f"{project_name}_bulk_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump({
            "project": project_name,
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "runs_included": summary,
            "runs_skipped": skipped,
            "total_pairs": len(combined_rows),
        }, fh, indent=2)

    return {"paths": paths, "summary": summary, "skipped": skipped}


# ---------------------------------------------------------------------------
# Overlap detection (used by the review-dialog visualization)
# ---------------------------------------------------------------------------

def find_overlaps(
    intervals: List[Dict[str, Any]],
    *,
    min_gap: int = 0,
) -> List[tuple]:
    """Return pairs of (original) indices that are closer than ``min_gap`` bp.

    With ``min_gap=0`` (default) this returns pairs whose ``[start, end)``
    intervals truly overlap. With ``min_gap>0`` it also returns pairs that
    are separated by fewer than ``min_gap`` bp — used to flag chosen probes
    that sit too close together (see ``constants.MIN_PROBE_GAP_BP``).

    Sweep-line over intervals sorted by start. The per-step rebuild of the
    active list costs O(|active|), so the total is O(N log N + K) where K is the
    work proportional to overlap density (the number of flagged pairs scanned).
    For sparse, well-spaced probe layouts |active| stays small so this is close
    to O(N log N); for pathologically dense overlap it degrades toward O(N²) —
    which is unavoidable, since K itself can be O(N²). Still far better than an
    unconditional O(N²) double loop. ``intervals`` is a list of dicts with at
    least ``start`` and ``end``.
    """
    order = sorted(
        range(len(intervals)),
        key=lambda i: (int(intervals[i]["start"]), int(intervals[i]["end"])),
    )
    overlaps: List[tuple] = []
    active: List[tuple] = []   # (end, original_index) of still-open intervals
    for idx in order:
        s = int(intervals[idx]["start"])
        e = int(intervals[idx]["end"])
        # Drop intervals that close more than min_gap before this one starts.
        active = [(ae, ai) for (ae, ai) in active if ae > s - min_gap]
        # Everything still active overlaps or is within min_gap of this one.
        for _ae, ai in active:
            overlaps.append((min(ai, idx), max(ai, idx)))
        active.append((e, idx))
    return overlaps
