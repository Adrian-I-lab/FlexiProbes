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

"""Robust FASTA parser with validation and ambiguity-code support."""
from __future__ import annotations
import logging
import os
import re
from typing import Dict, Iterator, List, Optional, Tuple
from .errors import FastaParseError

logger = logging.getLogger("flexprobe.io_fasta")

_IUPAC = set("ACGTURYSWKMBDHVNacgturyswkmbdhvn")
# Alignment gap characters ('-' and '.') are tolerated at ingest so aligned /
# consensus FASTAs (e.g. multi-strain segments, MSA-derived references) parse
# without aborting the whole run. Gaps are KEPT in the sequence so tile
# coordinates stay aligned to the input; the filter funnel then rejects any
# tile that contains a gap (see filters._filter_body), so a gap can never leak
# into a designed probe.
_GAP = set("-.")
_ALLOWED = _IUPAC | _GAP
_HEADER_RE = re.compile(r'^>\s*(\S+)(.*)')

# parse_fasta reads the whole file into memory. Whole-transcriptome FASTAs are
# legitimately large, so a big file is only a warning, not an error — it flags
# the common "I dragged the wrong file" mistake. Callers that need a hard
# ceiling (e.g. the GUI on a constrained machine) pass an explicit max_bytes.
_WARN_BYTES = 500 * 1024 * 1024  # 500 MB


def iter_fasta(
    path: str,
    *,
    max_bytes: Optional[int] = None,
    split_header_on_pipe: bool = False,
    check_duplicates: bool = True,
) -> Iterator[Tuple[str, str, str]]:
    """Stream a FASTA file one record at a time.

    Yields (gene_id, description, sequence) tuples, validating each record as it
    is produced, so a consumer that processes records independently never needs
    the whole file resident. ``parse_fasta`` is the list-materialising wrapper
    for callers (gene-aware / multi-strain) that genuinely need random access.

    ``max_bytes`` refuses a file larger than the limit up front; otherwise files
    over ~500 MB only warn.

    The record id is the FULL first whitespace-delimited token of the header, as
    the FASTA convention specifies. It is deliberately NOT truncated at ``|``:
    doing so collapses every legacy NCBI-style header (``>gi|12345|ref|NM_...``)
    onto the id ``gi``, and because the pipeline keys sequences by id in a dict,
    the records then silently overwrite one another — N transcripts in, one gene
    designed against, no warning. ``split_header_on_pipe=True`` restores the old
    behaviour for inputs where the pipe-split really is wanted; it is opt-in
    precisely because it can merge records.

    ``check_duplicates`` makes the stream fail closed on the second occurrence of
    an id. :func:`parse_fasta` turns it off and does its own full-file check so
    it can report *every* collision at once rather than only the first.
    """
    try:
        size = os.path.getsize(path)
    except OSError as exc:
        raise FastaParseError(f"Cannot open FASTA file: {path}") from exc
    if max_bytes is not None and size > max_bytes:
        raise FastaParseError(
            f"FASTA file is {size / 1e6:.1f} MB, exceeding the "
            f"{max_bytes / 1e6:.1f} MB limit ({path}). Split the input into "
            f"smaller files or raise the limit."
        )
    if size > _WARN_BYTES:
        logger.warning(
            "Large FASTA (%.1f MB): %s", size / 1e6, path,
        )

    def _finish(gid, desc, parts):
        seq = ''.join(parts)
        if not seq:
            raise FastaParseError(f"Empty sequence for {gid}")
        bad = set(seq) - _ALLOWED
        if bad:
            raise FastaParseError(f"Invalid characters {bad} in sequence for {gid}")
        if _GAP & set(seq):
            logger.warning(
                "Sequence for %s contains alignment gap(s); gap-containing "
                "tiles will be skipped during filtering.", gid,
            )
        return (gid, desc, seq)

    current_id: str | None = None
    current_desc: str = ""
    seq_parts: List[str] = []
    seen: Dict[str, int] = {}
    try:
        with open(path, 'r') as fh:
            for lineno, raw_line in enumerate(fh, 1):
                line = raw_line.strip()
                if not line or line.startswith('#'):
                    continue
                if line.startswith('>'):
                    if current_id is not None:
                        yield _finish(current_id, current_desc, seq_parts)
                    m = _HEADER_RE.match(line)
                    if not m:
                        raise FastaParseError(
                            f"Malformed header at line {lineno}: {line!r}")
                    token = m.group(1)
                    current_id = token.split('|')[0] if split_header_on_pipe else token
                    if not current_id:
                        raise FastaParseError(
                            f"Empty record id at line {lineno}: {line!r}")
                    if check_duplicates:
                        if current_id in seen:
                            raise FastaParseError(
                                f"Duplicate FASTA record id {current_id!r} at line "
                                f"{lineno} (already seen at line {seen[current_id]}). "
                                f"Ids must be unique — the pipeline keys sequences by "
                                f"id, so duplicates would silently overwrite each "
                                f"other."
                            )
                        seen[current_id] = lineno
                    current_desc = m.group(2).strip()
                    seq_parts = []
                else:
                    if current_id is None:
                        raise FastaParseError(
                            f"Sequence data before first header at line {lineno}")
                    seq_parts.append(line)
    except OSError as exc:
        raise FastaParseError(f"Cannot open FASTA file: {path}") from exc
    if current_id is not None:
        yield _finish(current_id, current_desc, seq_parts)


def parse_fasta(
    path: str,
    *,
    max_bytes: Optional[int] = None,
    split_header_on_pipe: bool = False,
) -> List[Tuple[str, str, str]]:
    """Parse a FASTA file into a list of (gene_id, description, sequence).

    The id is the full first whitespace-delimited token of the header (see
    :func:`iter_fasta` for why it is not split on ``|`` unless asked).

    Raises :class:`FastaParseError` if two records share an id, listing every
    collision. The pipeline stores sequences in a ``{id: seq}`` dict, so without
    this a duplicate would silently overwrite its predecessor and the run would
    design against fewer targets than the user supplied — while still reporting
    the full input count.

    Thin materialising wrapper over :func:`iter_fasta`; the gene-aware and
    multi-strain modes need every record resident, so this list is the memory
    floor for those runs (single-target/streaming consumers can use iter_fasta).
    """
    records = list(iter_fasta(
        path, max_bytes=max_bytes,
        split_header_on_pipe=split_header_on_pipe,
        # Collect everything first so the error can name EVERY collision, not
        # just the first one a user would then have to fix and re-run to find.
        check_duplicates=False,
    ))
    if not records:
        raise FastaParseError(f"No records found in {path}")
    _reject_duplicate_ids(records, path, split_header_on_pipe)
    return records


def _reject_duplicate_ids(
    records: List[Tuple[str, str, str]],
    path: str,
    split_header_on_pipe: bool,
) -> None:
    """Raise ``FastaParseError`` naming every id that appears more than once."""
    counts: Dict[str, int] = {}
    for gid, _desc, _seq in records:
        counts[gid] = counts.get(gid, 0) + 1
    dups = sorted((gid, n) for gid, n in counts.items() if n > 1)
    if not dups:
        return
    shown = ", ".join(f"{gid!r} x{n}" for gid, n in dups[:20])
    if len(dups) > 20:
        shown += f", ... ({len(dups) - 20} more)"
    hint = (
        " Header ids were truncated at '|' (split_header_on_pipe=True); turn "
        "that off to keep the full header token."
        if split_header_on_pipe else ""
    )
    raise FastaParseError(
        f"Duplicate record id(s) in {path}: {shown}. "
        f"{len(records)} records collapse to {len(counts)} unique ids, so "
        f"sequences would silently overwrite each other.{hint}"
    )


def fasta_to_dict(path: str, *, split_header_on_pipe: bool = False) -> Dict[str, str]:
    """Convenience: parse FASTA and return {gene_id: sequence}."""
    return {
        gid: seq for gid, _desc, seq
        in parse_fasta(path, split_header_on_pipe=split_header_on_pipe)
    }
