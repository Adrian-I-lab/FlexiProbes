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

"""On-disk checkpoint for the post-cap candidate set (resumable runs).

Tiling + filtering a whole transcriptome generates hundreds of millions of tiles
before the off-target BLAST even starts. When an HPC job is killed by its
walltime and resubmitted, [`screen_cache`](screen_cache.py) already lets BLAST pick
up where it stopped — but only *after* the pipeline has re-tiled and re-filtered
everything from scratch, which can cost a large slice of each restart's wall
time.

This module checkpoints the candidate list that survives Stage 3b (the pre-screen
``screen_top_per_gene`` cap) — the exact set of probes that go into BLAST. On a
resubmission with an unchanged input + config, the pipeline reloads it and skips
tiling/filtering/capping entirely, jumping straight to BLAST (which then resumes
from its own chunk cache).

A ``_meta.json`` fingerprint (input FASTA identity + the whole ``RunConfig``
minus a documented list of provably irrelevant fields) guards correctness
exactly like ``screen_cache``: change ``tile_step``, any filter, the tile size,
the cap, or the input and the stale checkpoint is ignored and the candidates are
regenerated. The cache dir is dot-prefixed and removed on a fully successful run,
so it never lands in the results tarball.

The two caches are deliberately kept in *separate* directories so that
``screen_cache.open_cache`` wiping a stale BLAST cache can never delete a freshly
written candidate checkpoint (and vice versa).
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
from typing import Any, Dict, Iterable, List, Optional, Tuple

from .models import ProbeCandidate
from .errors import CacheStampMismatchError

CACHE_DIRNAME = ".candidate_cache"
_META = "_meta.json"
_DATA = "candidates.jsonl"

# Bump whenever the ON-DISK format OR the fingerprint derivation changes, so a
# checkpoint written by older code can never be matched by newer code (it is
# folded into the digest, exactly as ``screen_cache.CACHE_SCHEMA`` is).
#   1 = allowlist fingerprint over a hand-picked subset of the config. An
#       allowlist can miss a field that changes the candidate set (a pre-tiling
#       option, or a filter such as filters.reject_junction_base_n), and then
#       a stale checkpoint returns probes the flag should have excluded.
#   2 = denylist fingerprint derived from the whole RunConfig
#   3 = ``_meta.json`` carries an engine stamp
#       (:func:`flexprobe.engine_identity.current_stamp`). A checkpoint with a
#       missing or different stamp raises CacheStampMismatchError on load and
#       is left in place, so candidates tiled and filtered by other engine
#       code are never reloaded and never silently discarded
CACHE_SCHEMA = 3

# The ONLY config fields excluded from the fingerprint, each because it cannot
# change which candidates the tile+filter+cap stage produces. Everything else —
# including any field added to RunConfig in future — is folded in by DEFAULT, so
# the failure mode of forgetting to update this module is a redundant cache miss
# (slow) rather than a stale checkpoint silently answering for a config it was
# not built under (wrong).
#
# Entries are dotted paths into ``config_to_dict``; naming a nested block drops
# the whole block. Justification per entry:
#   output_dir          - where results are written; the cache lives inside it
#   export              - output-writer toggles; run after the checkpoint
#   screen              - DB paths / skip flags; the screen runs after the checkpoint
#   specificity         - off-target thresholds; applied to hits, post-checkpoint
#   extra_on_target_ids - specificity input only
#   min_isoform_coverage- gate applied after BLAST
#   discovery_only      - which files get written
#   auto_select_per_gene- post-pipeline selection step
#   primer3_scoring     - chemistry/ranking stage; prescreen_score never uses it
#   multiplex           - chemistry (which RHS handle), post-checkpoint
#   probe_barcodes      - chemistry/export, post-checkpoint
#   seed                - provenance only; the pipeline has no RNG
#   workers             - parallelism only; the parallel path is required to
#                         produce a candidate set identical to the serial one
#   screen_workers      - parallelism for the OFF-TARGET SCREEN, a stage that
#                         runs entirely AFTER this checkpoint is written. It
#                         cannot reach the tile+filter+cap stage at all, so
#                         folding it in would re-tile a whole transcriptome
#                         (hours on human) to reproduce a byte-identical
#                         candidate set the moment anyone raised it.
#   resume_screen        - whether caching happens at all, not what it contains
#   screen_cache_mode    - the BLAST cache's on-disk format, not the candidates
#   emit_prescreen_candidates / emit_prescreen_precap_rows
#                       - write out the candidate set as an analysis artefact;
#                         reporting only, and the checkpoint is what they report
#                         ON. Folding them in would make turning the artefact on
#                         invalidate a perfectly good checkpoint mid-chain and
#                         re-tile a whole transcriptome to produce the same set.
#   keep_candidate_cache- whether this cache survives the end of the run
#   prepare_only        - WHERE the run stops, never which candidates it built.
#                         It is read once, after this checkpoint has already
#                         been saved (pipeline.run_pipeline), and the tile,
#                         filter and cap stages cannot see it at all.
#
#                         THIS ONE IS LOAD-BEARING. The split
#                         chain is one bundle and one config.json: job A adds
#                         --prepare-only on the command line and job B does
#                         not, which is the only difference between them. In
#                         the fingerprint, that difference alone gives the two
#                         jobs different digests, so job B misses job A's
#                         checkpoint and re-tiles the whole transcriptome
#                         inside the smaller allocation the split just gave it
#                         — i.e. the split would silently cost more than it
#                         saved. tests/test_prepare_only.py
#                         (test_prepare_only_does_not_change_the_checkpoint_
#                         fingerprint) is the guard that keeps them equal.
FINGERPRINT_IGNORED = frozenset({
    "output_dir",
    "export",
    "screen",
    "specificity",
    "extra_on_target_ids",
    "min_isoform_coverage",
    "discovery_only",
    "auto_select_per_gene",
    "primer3_scoring",
    "multiplex",
    "probe_barcodes",
    "seed",
    "workers",
    "screen_workers",
    "resume_screen",
    "screen_cache_mode",
    "emit_prescreen_candidates",
    "emit_prescreen_precap_rows",
    "keep_candidate_cache",
    "prepare_only",
})

# Only the fields needed to reconstruct a pre-screen candidate. Everything else
# (specificity verdict, isoform coverage, chemistry-derived oligos, ranking) is
# filled in by the downstream stages, so persisting it would just bloat the file.
_PERSIST_FIELDS = (
    "gene_id", "probe_id", "sequence", "start", "end",
    "gc_content", "gc_left", "gc_right", "homopolymer_max",
    # Set by the filter funnel, which a resumed run skips — so it has to be
    # persisted or a resumed results.csv would differ from a clean one. Ranking
    # also reads it (the TN-junction preference), so losing it would silently
    # change probe order on resume, not just one output column.
    "junction_base",
)


def cache_dir(output_dir: str) -> str:
    return os.path.join(output_dir, CACHE_DIRNAME)


def fingerprint_payload(config) -> Dict[str, Any]:
    """The config subset the fingerprint is derived from (whole config minus
    :data:`FINGERPRINT_IGNORED`). Exposed so a test can enumerate exactly what
    is and is not folded in."""
    from .config_io import config_to_dict

    return {
        k: v for k, v in config_to_dict(config).items()
        if k not in FINGERPRINT_IGNORED
    }


def _stable(o: Any) -> Any:
    """JSON fallback: order-independent for sets, string-coerced otherwise."""
    if isinstance(o, (set, frozenset)):
        return sorted(str(x) for x in o)
    return str(o)


def fingerprint(config) -> str:
    """Stable hash of the input FASTA + every config field affecting the candidate set.

    Derived by DENYLIST, not allowlist: the whole ``RunConfig`` is serialized and
    only the fields in :data:`FINGERPRINT_IGNORED` are dropped. That inversion is
    the point. An allowlist of hand-picked values silently ignores every other
    field, including pre-tiling options and filters such as
    ``filters.reject_junction_base_n``. The documented HPC workflow is "job hits
    walltime, resubmit with a flag added". Under an allowlist, adding such a flag
    gives an IDENTICAL digest, so the stale checkpoint is reloaded, tiling and
    filtering are skipped, and the results come back as if the flag had been
    honoured. A silent wrong answer, which the denylist prevents.

    With the denylist, a config field added in future is folded in by default:
    the worst case is an unnecessary cache miss, never a stale hit.

    The FASTA is identified cheaply by size + mtime — it does not change between
    a job and its resubmission in the same run dir, and re-hashing a multi-GB
    transcriptome on every restart would be wasteful. (The manifest records the
    real SHA-256 of the input, once, for the provenance question.)
    """
    h = hashlib.sha256()
    h.update(f"schema={CACHE_SCHEMA}\0".encode())
    fasta_path = getattr(config, "fasta_path", "") or ""
    try:
        st = os.stat(fasta_path)
        h.update(f"fasta={os.path.basename(fasta_path)}\0".encode())
        h.update(f"size={st.st_size}\0mtime={st.st_mtime_ns}\0".encode())
    except OSError:
        # No stat (e.g. path not yet present) — the bare path is still in the
        # config payload below, so a later mismatch is at worst a conservative
        # cache miss, never a false hit.
        h.update(b"fasta=<unstatable>\0")
    h.update(json.dumps(
        fingerprint_payload(config), sort_keys=True, default=_stable,
    ).encode())
    return h.hexdigest()


def load(
    output_dir: str, fp: str, *, stamp: Optional[Dict[str, Any]] = None,
) -> Optional[Tuple[List[ProbeCandidate], Dict[str, Any]]]:
    """Return ``(candidates, meta)`` if a matching checkpoint exists, else None.

    Raises :class:`~flexprobe.errors.CacheStampMismatchError`, leaving the
    checkpoint in place, when its engine stamp is missing or differs from
    ``stamp`` (default :func:`flexprobe.engine_identity.current_stamp`).

    Returns None when the cache is absent, the fingerprint differs (stale), or
    the data file is missing/torn (a write killed mid-flush) — in every case the
    caller regenerates the candidates from scratch. ``meta`` carries the cached
    headline counts (``tiles_total``, ``filter_passed``) so a resumed run can
    still report coherent stats without re-tiling.
    """
    cdir = cache_dir(output_dir)
    meta_path = os.path.join(cdir, _META)
    data_path = os.path.join(cdir, _DATA)
    try:
        with open(meta_path, encoding="utf-8") as fh:
            meta = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None
    if stamp is None:
        from .engine_identity import current_stamp
        stamp = current_stamp()
    stamp = dict(stamp)
    cached_stamp = meta.get("stamp") if isinstance(meta, dict) else None
    if cached_stamp != stamp:
        # Checked before the fingerprint: other engine code must stop the run,
        # not quietly fall through to a re-tile that overwrites the checkpoint.
        raise CacheStampMismatchError(cdir, cached_stamp, stamp)
    if meta.get("fingerprint") != fp:
        return None
    if not meta.get("complete"):
        # The save was interrupted before it committed — don't trust a partial.
        return None
    out: List[ProbeCandidate] = []
    try:
        with open(data_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                d = json.loads(line)
                out.append(ProbeCandidate(**d))
    except (OSError, json.JSONDecodeError, TypeError):
        return None
    return out, meta


def save(
    output_dir: str, fp: str, candidates: Iterable[ProbeCandidate],
    *, tiles_total: int = 0, filter_passed: int = 0,
    stamp: Optional[Dict[str, Any]] = None,
) -> None:
    """Atomically checkpoint ``candidates`` for fingerprint ``fp``.

    ``stamp`` (default :func:`flexprobe.engine_identity.current_stamp`) is
    written into ``_meta.json`` so :func:`load` can refuse a checkpoint from
    other engine code.

    The data file is streamed to a temp path and renamed into place, then a
    ``_meta.json`` carrying the fingerprint, the headline counts, and a
    ``complete`` flag is committed last. ``load`` only trusts a checkpoint whose
    meta says ``complete``, so a job killed mid-write is treated as no checkpoint
    (the candidates regenerate) rather than as a truncated one.
    """
    if stamp is None:
        from .engine_identity import current_stamp
        stamp = current_stamp()
    stamp = dict(stamp)
    cdir = cache_dir(output_dir)
    # Start clean so a stale (different-fingerprint) checkpoint never lingers.
    # Only reached after load() has accepted or found no checkpoint, so this
    # never removes one with a different engine stamp.
    shutil.rmtree(cdir, ignore_errors=True)
    os.makedirs(cdir, exist_ok=True)

    data_path = os.path.join(cdir, _DATA)
    data_tmp = data_path + ".tmp"
    valid = set(_PERSIST_FIELDS)
    with open(data_tmp, "w", encoding="utf-8") as fh:
        for c in candidates:
            row = {k: getattr(c, k) for k in _PERSIST_FIELDS if hasattr(c, k)}
            # Defensive: a future field rename would otherwise silently drop data.
            assert set(row).issubset(valid)
            fh.write(json.dumps(row))
            fh.write("\n")
    os.replace(data_tmp, data_path)

    meta_path = os.path.join(cdir, _META)
    meta_tmp = meta_path + ".tmp"
    with open(meta_tmp, "w", encoding="utf-8") as fh:
        json.dump({
            "fingerprint": fp, "complete": True, "stamp": stamp,
            "tiles_total": int(tiles_total), "filter_passed": int(filter_passed),
        }, fh)
    os.replace(meta_tmp, meta_path)


def clear(output_dir: str) -> None:
    """Remove the cache dir (called after a fully successful run)."""
    shutil.rmtree(cache_dir(output_dir), ignore_errors=True)
