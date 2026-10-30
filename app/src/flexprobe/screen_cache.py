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

"""On-disk checkpoint cache for the chunked BLAST stage (resumable runs).

A whole-transcriptome off-target screen is split into thousands of BLAST chunks
and can outlast an HPC walltime. This cache lets an interrupted run resume: each
chunk's parsed hits are written to ``<output>/.screen_cache/chunk_NNNNN.json`` as
it completes, so a resubmitted job reloads finished chunks and only BLASTs the
ones that are missing.

A ``_meta.json`` fingerprint (candidate count + a hash of every probe id + the
DB + thresholds + the reference's own content hash) guards correctness: if the
current run's candidates differ from the cached ones — e.g. ``tile_step``
changed, or the reference FASTA at the same path was edited — the stale cache
is wiped and BLAST starts fresh, so results from two different configs (or two
different reference files at one path) can never be mixed.

The cache dir is dot-prefixed and removed on a fully successful run, so it exists
only while a run is mid-flight and never lands in the results tarball.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

from .constants import REFERENCE_STRAND_BOTH
from .errors import CacheStampMismatchError

CACHE_DIRNAME = ".screen_cache"
_META = "_meta.json"

# Bump whenever the ON-DISK chunk format changes. This is folded into the
# fingerprint, so a cache written in an older format can never be read back
# by newer code: the digest won't match and the stale dir is wiped. Every
# change listed below is a silent-wrong-answer hazard if an old chunk were
# reused, which is what the schema exists to prevent.
#   1 = raw hits  {probe_id: [hit, ...]}
#   2 = verdicts  {probe_id: {verdict field: v}}  (see VERDICT_FIELDS)
#   3 = full-half mismatch recount. blast.BLASTN_OUTFMT_COLS carries
#       ``qseq``/``sseq`` (a `hits`-mode chunk without them would be re-scored
#       by the conservative path while claiming to be a normal run), and
#       VERDICT_FIELDS carries ``offtarget_recount_basis`` (without it a
#       verdict chunk restores every probe to the "not_screened" default).
#   4 = junction-duplex gate. A `hits`-mode chunk carries one junction-duplex
#       record per probe (screen.KmerScreen._junction_scan), and VERDICT_FIELDS
#       carries offtarget_clean_left_bp / _right_bp. A chunk without the
#       junction record would report the gate as "not measured" for every
#       probe, so the gate would never fire.
#   5 = degenerate probe halves. A probe half carrying an 'N' or an IUPAC code
#       is screened fail-closed and recorded as `unscreenable_half`. A schema-4
#       chunk has no such record and would restore a hit list (or a verdict)
#       that looks complete, so the pair would keep a pass it never earned.
#   6 = ``_meta.json`` carries an engine stamp
#       (:func:`flexprobe.engine_identity.current_stamp`). A cache keyed on
#       inputs only would let a run resumed after an engine fix read back
#       answers the fixed engine would not give (for example, BLAST hits
#       capped by an HSP limit the current engine no longer applies). A missing
#       or different stamp raises CacheStampMismatchError; see open_cache.
CACHE_SCHEMA = 6


# What a checkpoint stores. See RunConfig.screen_cache_mode for the trade-off.
MODE_VERDICTS = "verdicts"   # scored fields; size independent of hit volume
MODE_HITS = "hits"           # raw hits; re-scorable for free, but scales badly
CACHE_MODES = (MODE_VERDICTS, MODE_HITS)

# Every field `specificity.batch_evaluate()` writes onto a ProbeCandidate.
# This IS the resume contract: anything the scorer sets but we fail to persist
# would silently revert to its model default on resume. `isoform_hit_ids` is a
# list (not a scalar) — downstream isoform selection reads it without a FASTA.
VERDICT_FIELDS: Sequence[str] = (
    "specificity_pass",
    "specificity_status",
    "specificity_margin",
    # Not cosmetic: without it a resumed chunk would come back claiming
    # "not_screened" provenance for a verdict that WAS screened, or — worse —
    # silently lose a "conservative" marking and read as exactly recounted.
    "offtarget_recount_basis",
    "offtarget_min_qcov",
    "offtarget_worst_mismatches",
    "offtarget_worst_junction_mm",
    # The junction-duplex evidence. Without these a resumed chunk would report
    # -1 ("never measured") for probes that WERE measured, which is the same
    # class of silent provenance loss as offtarget_recount_basis above.
    "offtarget_clean_left_bp",
    "offtarget_clean_right_bp",
    "offtarget_nearmiss_count",
    "best_offtarget_pident",
    "best_offtarget_bitscore",
    "best_offtarget_evalue",
    "isoform_total",
    "isoform_hits",
    "isoform_coverage",
    "isoform_hit_ids",
)


def cache_dir(output_dir: str) -> str:
    return os.path.join(output_dir, CACHE_DIRNAME)


def fingerprint(
    probe_ids: Iterable[str],
    *,
    db_arg: str,
    chunk_size: int,
    perc_identity: Optional[float] = None,
    qcov_hsp_perc: Optional[float] = None,
    reference_strand: str = REFERENCE_STRAND_BOTH,
    spec_params: Optional[Mapping[str, Any]] = None,
    cache_mode: str = MODE_VERDICTS,
    engine: str = "blast",
    ref_content_hash: str = "",
) -> str:
    """Stable hash of the candidate set + engine parameters + cache format.

    O(n) over the probe ids; identical inputs (same FASTA + config) reproduce the
    same digest across runs, so chunk K always maps to the same probes.

    ``cache_mode`` participates because the two modes write different on-disk
    formats under the same schema: flipping the knob must invalidate an existing
    cache, never silently misread it.

    ``engine`` participates for the same reason, and it is the more dangerous
    case: two engines searching the same DB with the same candidates produce
    DIFFERENT hits (one local and trimmed, one end-to-end), so a digest blind to
    the engine would let a BLAST cache be read back as a k-mer run's answers and
    vice versa. Cached hits are the *input* to scoring, so that substitution
    would be silent and would change verdicts.

    ``ref_content_hash`` closes the gap the two paragraphs above leave open:
    ``db_arg`` is a PATH, not the reference's content. Re-run against the same
    path after the reference FILE changed underneath it (edited in place, or a
    different FASTA copied over the old one) and ``db_arg`` is byte-identical,
    so without this the stale cache would be read back and silently reused —
    the exact hazard ``engine`` and ``cache_mode`` above exist to prevent, just
    for the database's identity instead of its format. Callers pass whatever
    content signature they already have cheaply on hand (``Screen.
    fingerprint_parts``: the k-mer engine reuses the SHA-256 it computed to
    load its index; the BLAST engine hashes the source FASTA it resolves,
    falling back to the on-disk DB files when no FASTA is present). Defaults to
    ``""`` rather than being required, so a caller that genuinely cannot
    determine a content signature degrades to path-only keying instead of
    raising, the same "unknown, not fabricated" convention
    ``_blast_db_sha256s`` follows.

    ``perc_identity`` / ``qcov_hsp_perc`` default to None so an engine with no
    server-side pushdown need not invent values for them; engines pass whatever
    identifies their own search via ``Screen.fingerprint_parts``.

    ``spec_params`` MUST be supplied in ``verdicts`` mode and MUST be omitted in
    ``hits`` mode. A verdict is the *output* of specificity scoring, so it is only
    valid for the exact thresholds that produced it. Raw hits are scoring-
    agnostic, so keeping them out of the digest is what preserves the
    re-score-without-re-BLAST property that ``hits`` mode exists to offer.

    ``reference_strand`` participates for the same reason ``engine`` does: it
    changes which hits an engine reports (a "plus" run never even searches the
    minus strand), so a resumed run under a different setting must not read
    back another setting's answers.
    """
    h = hashlib.sha256()
    h.update(f"schema={CACHE_SCHEMA}\0mode={cache_mode}\0".encode())
    h.update(f"engine={engine}\0".encode())
    h.update(f"refhash={ref_content_hash}\0".encode())
    h.update(f"db={db_arg}\0chunk={chunk_size}\0".encode())
    h.update(f"pid={perc_identity}\0qcov={qcov_hsp_perc}\0".encode())
    h.update(f"strand={reference_strand}\0".encode())
    if spec_params is not None:
        # sort_keys => stable across dict ordering; default=str handles sets.
        h.update(json.dumps(spec_params, sort_keys=True, default=_stable).encode())
        h.update(b"\0")
    n = 0
    for pid in probe_ids:
        h.update(pid.encode())
        h.update(b"\0")
        n += 1
    h.update(f"n={n}".encode())
    return h.hexdigest()


def _stable(o: Any) -> Any:
    """JSON fallback: make sets order-independent so the digest is stable."""
    if isinstance(o, (set, frozenset)):
        return sorted(str(x) for x in o)
    return str(o)


def open_cache(
    output_dir: str, fp: str, *, stamp: Optional[Mapping[str, Any]] = None,
) -> str:
    """Prepare the cache dir for a run with fingerprint ``fp``.

    ``stamp`` identifies the engine writing the cache (default
    :func:`flexprobe.engine_identity.current_stamp` with no blastn). An existing
    cache is checked against it FIRST. A cache with no readable stamp, or a
    different one, raises :class:`~flexprobe.errors.CacheStampMismatchError`
    and is left exactly as it is: answers from other engine code are never
    reused, and never silently thrown away either.

    With the stamp matching, a matching fingerprint keeps the cache (its chunks
    will be reused). A different fingerprint means the same engine was given
    different inputs (the documented "resubmit with a flag changed" workflow),
    so the stale chunks are removed and a fresh cache created.
    Returns the cache directory path.
    """
    if stamp is None:
        from .engine_identity import current_stamp
        stamp = current_stamp()
    stamp = dict(stamp)
    cdir = cache_dir(output_dir)
    meta_path = os.path.join(cdir, _META)
    if os.path.isdir(cdir):
        try:
            with open(meta_path, encoding="utf-8") as fh:
                meta = json.load(fh)
        except (OSError, json.JSONDecodeError):
            meta = {}
        cached_stamp = meta.get("stamp") if isinstance(meta, dict) else None
        if cached_stamp != stamp:
            raise CacheStampMismatchError(cdir, cached_stamp, stamp)
        if meta.get("fingerprint") == fp:
            return cdir            # resume: keep the cached chunks
        shutil.rmtree(cdir, ignore_errors=True)   # same engine, new inputs
    os.makedirs(cdir, exist_ok=True)
    tmp = meta_path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump({"fingerprint": fp, "stamp": stamp}, fh)
    os.replace(tmp, meta_path)
    return cdir


def _chunk_path(cdir: str, ci: int) -> str:
    return os.path.join(cdir, f"chunk_{ci:05d}.json")


def load_chunk(cdir: str, ci: int) -> Optional[Dict[str, List[Dict[str, Any]]]]:
    """Return the cached hits for chunk ``ci`` (probe_id -> hits), or None."""
    path = _chunk_path(cdir, ci)
    if not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        # A torn write (killed mid-flush) — treat as not-done so it re-runs.
        return None


def save_chunk(cdir: str, ci: int, hits: Dict[str, List[Dict[str, Any]]]) -> None:
    """Atomically persist a chunk's hits (probes with >=1 hit only).

    Schema-1 path, kept for the NON-streaming pipeline, which merges raw hits
    into one dict and genuinely needs them. Do not use for whole-transcriptome
    runs: at 3.18M candidates this averaged 1.89 GB/chunk => ~1.2 TB over 635
    chunks. Streaming runs use save_verdicts() instead.
    """
    path = _chunk_path(cdir, ci)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(hits, fh)
    os.replace(tmp, path)


def save_verdicts(cdir: str, ci: int, candidates: Sequence[Any]) -> None:
    """Atomically persist a chunk's SCORED verdicts (schema 2).

    Stores the ~11 fields batch_evaluate() writes per candidate instead of the
    raw hits that produced them: ~2500x smaller (1.89 GB -> ~1-2 MB/chunk), and
    it is the only thing resume actually needs.

    Every candidate in the chunk is written, including probes with zero hits —
    they still carry a real verdict (specificity_pass=True via the clean path),
    and load_verdicts() uses full coverage to detect a torn write.
    """
    payload = {
        c.probe_id: {f: getattr(c, f) for f in VERDICT_FIELDS}
        for c in candidates
    }
    path = _chunk_path(cdir, ci)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh)
    os.replace(tmp, path)


def load_verdicts(cdir: str, ci: int) -> Optional[Dict[str, Dict[str, Any]]]:
    """Return cached verdicts for chunk ``ci`` (probe_id -> fields), or None."""
    path = _chunk_path(cdir, ci)
    if not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        # A torn write (killed mid-flush) — treat as not-done so it re-runs.
        return None


def apply_verdicts(
    candidates: Sequence[Any], verdicts: Mapping[str, Mapping[str, Any]]
) -> Optional[int]:
    """Restore cached verdicts onto ``candidates``; return the count passing.

    Returns None — meaning "cache miss, re-BLAST this chunk" — if ANY candidate
    is absent or any field is missing. A partially-applied chunk would silently
    leave some probes at their model defaults (specificity_pass=False), quietly
    dropping real probes from the final set; re-running the chunk is cheap by
    comparison. Nothing is mutated unless the whole chunk validates.
    """
    staged = []
    for c in candidates:
        v = verdicts.get(c.probe_id)
        if v is None or any(f not in v for f in VERDICT_FIELDS):
            return None
        staged.append((c, v))
    passed = 0
    for c, v in staged:
        for f in VERDICT_FIELDS:
            setattr(c, f, v[f])
        if c.specificity_pass:
            passed += 1
    return passed


def clear(output_dir: str) -> None:
    """Remove the cache dir (called after a fully successful run)."""
    shutil.rmtree(cache_dir(output_dir), ignore_errors=True)
