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

"""Exhaustive ground truth for the off-target screen.

Every other test in this suite checks the screen against its own conventions:
that a rule is applied where it should be, that a fallback lands in the safe
direction, that a column is written. None of them can say whether the answer is
*right*, because they all reason in the same terms the implementation does — and
four separate defects in the specificity path survived exactly that blind spot,
each one shipping a plausible number for an unknown one.

This module answers the different question. It finds every position in a
reference where a 25-mer sits within four mismatches, by construction rather
than by heuristic:

    a 25-mer carrying <= 4 mismatches must contain at least one of its five
    5 nt blocks intact (pigeonhole), so indexing every 5-mer of the reference
    and verifying each candidate offset by a FULL-LENGTH comparison cannot miss
    a hit

No local alignment, no trimming, no score threshold, nothing reconstructed. It
is slow and it is the point of reference the BLAST path is judged against — the
same oracle that adjudicated the removal of the bounded pushdown floor, where
BLAST-plus-recount agreed with truth on 119 of 120 disputed probes.

Marked ``slow`` and deselected by default. Run it with ``pytest -m slow``
whenever the specificity path changes.
"""
from __future__ import annotations

import random

import pytest

from flexprobe.specificity import full_half_mismatch_detail

HALF = 25
MAX_MM = 4
BLOCK = HALF // (MAX_MM + 1)          # 5 nt, the pigeonhole block size


def _exhaustive(reference: str, query: str, max_mm: int = MAX_MM):
    """Every (offset, mismatches) with mismatches <= max_mm. Ground truth."""
    n = len(query)
    out = []
    for i in range(len(reference) - n + 1):
        mm = sum(1 for a, b in zip(reference[i:i + n], query) if a != b)
        if mm <= max_mm:
            out.append((i, mm))
    return out


def _pigeonhole(reference: str, query: str, max_mm: int = MAX_MM):
    """Same answer, via 5-mer seeding — the method, not the definition."""
    n = len(query)
    index: dict = {}
    for i in range(len(reference) - BLOCK + 1):
        index.setdefault(reference[i:i + BLOCK], []).append(i)
    cands = set()
    for b in range(max_mm + 1):
        off = b * BLOCK
        for pos in index.get(query[off:off + BLOCK], ()):
            start = pos - off
            if 0 <= start <= len(reference) - n:
                cands.add(start)
    out = []
    for start in sorted(cands):
        mm = sum(1 for a, b in zip(reference[start:start + n], query) if a != b)
        if mm <= max_mm:
            out.append((start, mm))
    return out


@pytest.mark.slow
def test_pigeonhole_seeding_finds_exactly_what_brute_force_finds():
    """The oracle's own method must equal the oracle's definition.

    If seeding and brute force ever disagree, every verdict this module
    adjudicates is worthless — so it is checked before it is trusted.
    """
    rng = random.Random(20260805)
    for _ in range(25):
        ref = "".join(rng.choice("ACGT") for _ in range(4000))
        # Plant a query derived from the reference so hits actually exist.
        start = rng.randrange(0, len(ref) - HALF)
        q = list(ref[start:start + HALF])
        for pos in rng.sample(range(HALF), rng.randint(0, MAX_MM)):
            q[pos] = rng.choice([b for b in "ACGT" if b != q[pos]])
        query = "".join(q)
        assert _pigeonhole(ref, query) == _exhaustive(ref, query)


@pytest.mark.slow
def test_recount_matches_ground_truth_on_planted_off_targets():
    """The screen's mismatch count must equal the true full-length count.

    Builds a subject with a known number of planted mismatches, hands the
    recount the HSP a local aligner would report for it — deliberately TRIMMED
    to the longest clean core, which is what blastn-short actually emits — and
    requires the recounted answer to be the planted truth, not the geometry of
    the fragment it was given.
    """
    rng = random.Random(4578)

    class _Subject:
        def __init__(self, seq): self.seq = seq
        def window_clamped(self, sseqid, start, end):
            lo, hi = max(1, start), min(len(self.seq), end)
            win = self.seq[lo - 1:hi] if hi >= lo else ""
            return win, (end - start + 1) - len(win)

    for planted in range(0, MAX_MM + 1):
        for _ in range(40):
            query = "".join(rng.choice("ACGT") for _ in range(HALF))
            subj = list(query)
            for pos in rng.sample(range(HALF), planted):
                subj[pos] = rng.choice([b for b in "ACGT" if b != subj[pos]])
            subject = "".join(subj)
            # The longest clean run — what BLAST would report as the core.
            best_lo = best_len = 0
            run_lo = 0
            for i in range(HALF + 1):
                if i == HALF or query[i] != subject[i]:
                    if i - run_lo > best_len:
                        best_lo, best_len = run_lo, i - run_lo
                    run_lo = i + 1
            if best_len == 0:
                continue
            hit = {
                "qseqid": "p__L", "sseqid": "S",
                "pident": 100.0, "length": best_len,
                "qstart": best_lo + 1, "qend": best_lo + best_len,
                "sstart": best_lo + 1, "send": best_lo + best_len,
                "qseq": query[best_lo:best_lo + best_len],
                "sseq": subject[best_lo:best_lo + best_len],
            }
            n, basis = full_half_mismatch_detail(
                hit, HALF, query_seq=query, subjects=_Subject(subject))
            assert basis == "exact"
            assert n == planted, (
                f"recount said {n}, truth is {planted}; the HSP it was handed "
                f"covered only {best_len}/{HALF} nt"
            )
