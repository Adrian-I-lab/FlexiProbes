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

"""Small shared helpers used by more than one module."""
from __future__ import annotations

import math

from .constants import DEGENERATE_BASES


def max_homopolymer(seq: str) -> int:
    """Return the length of the longest single-base run in ``seq``.

    Case-insensitive; returns 0 for an empty sequence.
    """
    if not seq:
        return 0
    best = run = 1
    for i in range(1, len(seq)):
        if seq[i].upper() == seq[i - 1].upper():
            run += 1
            if run > best:
                best = run
        else:
            run = 1
    return best


def gc_fraction(seq: str) -> float:
    """Return the G+C fraction of ``seq`` (0.0 for empty). Case-insensitive.

    Counts literal G and C only, so an IUPAC code contributes ZERO however
    G/C-rich its resolutions are: an 'S' (G or C) is G+C with certainty and
    still counts 0, and an 'N' dilutes a 25-mer's GC by 4% per base. That is
    the honest answer for a mixed-base string — there is no single GC for a
    sequence that is several sequences — but it is the wrong number to judge a
    degenerate probe by. A consensus parent therefore reports a GC RANGE over
    its concrete variants (``gc_lhs_min/max``, ``gc_rhs_min/max``) rather than
    a GC computed from its own string; see ``models.ProbeCandidate`` and
    ``filters._filter_body``.
    """
    if not seq:
        return 0.0
    return sum(1 for c in seq.upper() if c in "GC") / len(seq)


def degenerate_count(seq: str) -> int:
    """Number of non-ACGT (IUPAC ambiguity or 'N') bases in ``seq``.

    Counts every code in ``constants.DEGENERATE_BASES``, not just literal 'N'.
    A consensus built from disagreeing records emits two-way codes (M/R/W/S/Y/K)
    far more often than 'N', so an 'N'-only count measures almost nothing on a
    real strain alignment.

    This runs twice per candidate in the filter funnel — hundreds of millions of
    times on a whole-transcriptome run — where it replaced a single C-level
    ``str.count('N')``. Hence the fast path: four C-level counts settle the
    overwhelmingly common "this half is pure ACGT" case without touching the
    per-character loop, which then only runs on the rare ambiguous half. The
    two-step form is not an approximation; both branches count exactly the
    members of ``DEGENERATE_BASES``.
    """
    s = (seq or "").upper()
    if len(s) == s.count("A") + s.count("C") + s.count("G") + s.count("T"):
        return 0
    return sum(1 for c in s if c in DEGENERATE_BASES)


def safe_tm(seq: str) -> float:
    """Return ``primer3.calc_tm(seq)`` when possible, otherwise NaN.

    Returns NaN (rather than 0.0) in these cases:
      • primer3 isn't installed,
      • the sequence is empty,
      • the sequence contains any base outside ACGT after uppercasing
        (primer3 rejects ambiguity codes like N).

    NaN is the honest sentinel: 0.0 looked like a real (impossibly low) Tm.
    Inputs are uppercased defensively because input FASTA may be lowercase
    and primer3 wants uppercase.
    """
    try:
        import primer3
    except ImportError:
        return float("nan")
    s = (seq or "").upper()
    if not s:
        return float("nan")
    if any(c not in "ACGT" for c in s):
        return float("nan")
    from .constants import TM_DNA_CONC, TM_DNTP_CONC, TM_DV_CONC, TM_MV_CONC
    try:
        # Explicit salt conditions (documented as adjustable; primer3 PCR
        # defaults, not the unpublished Flex buffer). See constants.TM_MV_CONC.
        return float(primer3.calc_tm(
            s, mv_conc=TM_MV_CONC, dv_conc=TM_DV_CONC,
            dntp_conc=TM_DNTP_CONC, dna_conc=TM_DNA_CONC,
        ))
    except Exception:
        return float("nan")


def is_nan(value) -> bool:
    """True if ``value`` is a float NaN. Safe for non-float inputs."""
    return isinstance(value, float) and math.isnan(value)


def fmt_tm(value) -> str:
    """Render a Tm for CSV/UI: ``n/a`` for NaN/None, ``%.2f`` otherwise."""
    if value is None or is_nan(value):
        return "n/a"
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return "n/a"
