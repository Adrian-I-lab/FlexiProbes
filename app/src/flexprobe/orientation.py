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

"""Input-orientation sanity check: is this FASTA the sense strand?

WHY THIS EXISTS
---------------
Every probe this tool designs is built from the input sequence read as the
SENSE mRNA — see the chemistry module's header diagram. Both halves of a pair
are reverse complements of the input window, so they hybridise to the input
strand and to nothing else.

That makes the input's orientation a silent correctness precondition. Nothing
downstream can notice it is wrong: GC content, Tm, homopolymer runs, secondary
structure and the off-target screen are all strand-symmetric or are computed
against the same (wrong) strand, so a reverse-complemented input passes every
filter in the pipeline and yields a full, well-scored, orderable probe set that
detects nothing in the assay.

The flagship application makes this a live trap rather than a theoretical one.
Oropouche virus is a negative-sense RNA virus: its genome (vRNA) is the
COMPLEMENT of the mRNA the probes must bind, and GenBank segment records are
deposited in either convention depending on the submission. A user who
downloads the genome strand and designs against it gets a panel that looks
perfect and reports nothing.

WHAT THIS DOES
--------------
Finds the longest open reading frame on each strand of each record and warns —
never refuses — when the reverse strand's is substantially longer, which is the
signature of a sequence supplied in genome/vRNA orientation.

It is a heuristic and it is deliberately quiet. A UTR-only fragment, a short
amplicon or a genuinely non-coding target has no meaningful ORF on either
strand, and this module says nothing about those rather than crying wolf. See
:data:`ORF_RATIO_WARN` for the threshold and the measurements behind it.
"""
from __future__ import annotations

import logging
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .tiling import reverse_complement

logger = logging.getLogger(__name__)

_STOPS = frozenset(("TAA", "TAG", "TGA"))
_CODON_RE = re.compile("...")

# ---------------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------------

#: Warn when the reverse strand's longest ORF is at least this many times the
#: forward strand's.
#:
#: Calibrated on orthobunyavirus segment records (complete and near-complete
#: coding segments) read in both orientations. Correctly oriented records sit
#: near or below 1, and the same records reverse-complemented sit at roughly 5
#: to 40, so the two populations are well separated and 2.0 sits in the gap
#: with room on both sides. It is also the
#: right order of magnitude a priori: under a null of non-coding sequence the
#: two strands' longest ORFs are draws from the same distribution and the ratio
#: concentrates near 1, so a factor of two is already well outside the noise.
#:
#: Left at 2.0 rather than tightened toward the lowest ratio seen on a flipped
#: record: the ratio is the SENSITIVE half of the test and :data:`MIN_ORF_NT`
#: is the specific half. Raising this would buy almost nothing on false
#: positives (which MIN_ORF_NT controls) while narrowing the margin on a
#: genuinely flipped input whose two strands happen to be closer together.
#:
#: Some records do not separate at all. A short, ambiguity-rich consensus of an
#: orthobunyavirus S segment, whose N and NSs reading frames overlap, scores
#: close to 1 in both orientations and is silent either way. That is the
#: intended failure mode: such a record does not carry the signal, and this
#: check declines to guess rather than firing at random on it.
ORF_RATIO_WARN = 2.0

#: A reverse-strand ORF shorter than this never triggers a warning, whatever the
#: ratio. This floor, not the ratio, is what does the false-positive work: a
#: 2x ratio between two ORFs that are both noise is itself noise, and short
#: non-coding sequence produces exactly that.
#:
#: 600 nt = 200 codons. Calibrated on 4,500 simulated non-coding records
#: (lengths 1-30 kb x GC 0.30/0.40/0.50, uniform random), counting how often
#: the pair (floor, ratio >= 2.0) fires on sequence that has no coding strand
#: at all:
#:
#:     floor  300 nt -> 1.09%      floor  600 nt -> 0.09%
#:     floor  450 nt -> 0.38%      floor  750 nt -> 0.00%
#:
#: 600 is chosen over 750 so that a short viral segment supplied backwards is
#: still caught: a ~900 nt segment can leave a reverse-strand ORF of only
#: about 700 nt, which a floor of 750 would miss. 600 costs 0.09% on noise.
#:
#: The practical consequence to know about: a genuinely short CDS — under ~600
#: nt — supplied backwards will NOT be flagged. That is the deliberate side of
#: "warn, never refuse": this check earns its keep on whole viral segments and
#: transcripts, and stays quiet where it cannot tell.
MIN_ORF_NT = 600

#: Records shorter than this are skipped outright — they cannot carry an ORF
#: over MIN_ORF_NT on either strand, so scanning them only costs time. Tied to
#: MIN_ORF_NT so the two can never drift apart.
MIN_RECORD_NT = MIN_ORF_NT

# ---------------------------------------------------------------------------
# Cost control
# ---------------------------------------------------------------------------
# A whole-transcriptome input is hundreds of thousands of records and hundreds
# of megabases; a six-frame scan of all of it would be a visible tax on every
# run for a check that is aimed at hand-assembled viral inputs. So the scan is
# SAMPLED, and says so in the log when it samples.
#
# Sampling is sound for what this detects: an input in the wrong orientation is
# wrong as a whole (someone downloaded the genome strand, or reverse-
# complemented a file), so a spread sample of a few dozen records finds it just
# as reliably as scanning every one. It is NOT designed to find a single flipped
# record inside an otherwise-correct transcriptome, and does not claim to.

#: Most records examined. Small inputs — every viral panel this tool is aimed
#: at — are well under this and are scanned in full.
MAX_RECORDS = 64

#: Most bases read from any one record. Long enough that a full viral segment or
#: a long transcript is scanned end to end; short enough that a scaffold-sized
#: record cannot dominate the cost.
MAX_BASES = 30_000

#: Most warnings emitted before they are collapsed into a single count. A wholly
#: reverse-complemented input would otherwise produce one line per record.
MAX_WARNINGS = 5


def longest_orf(seq: str) -> int:
    """Length in nt of the longest ORF in any of *seq*'s three forward frames.

    An ORF is an ATG followed by in-frame codons up to (and including) the first
    stop, or up to the end of the sequence if no in-frame stop occurs before it.
    Running off the 3' end counts: a deposited genome segment is frequently a
    near-complete CDS whose stop codon is outside the record, and requiring a
    terminated ORF would score those at zero, which would read a multi-kb
    segment with a completely open reading frame as having no ORF at all.

    Ambiguity codes (N and the IUPAC degenerate set) are treated as neither a
    start nor a stop, so a codon containing one neither opens nor closes an ORF.
    That is the conservative direction: ambiguity then lengthens the apparent
    ORFs on BOTH strands roughly equally, which pushes the strand ratio toward
    1 — toward silence — and never toward a false alarm. The alternative
    (breaking an ORF at every ambiguous codon) destroys the signal on an
    N-containing consensus, cutting a multi-kb ORF to a few tens of nt.
    """
    seq = seq.upper()
    n = len(seq)
    best = 0
    for frame in range(3):
        codons = _CODON_RE.findall(seq[frame:])
        start: Optional[int] = None      # codon index of the open ATG
        for i, cod in enumerate(codons):
            if cod in _STOPS:
                if start is not None:
                    best = max(best, (i + 1 - start) * 3)
                    start = None
            elif start is None and cod == "ATG":
                start = i
        if start is not None:
            # Open at the 3' end: count to the end of the sequence.
            best = max(best, n - (frame + start * 3))
    return best


def orf_strand_ratio(seq: str) -> Tuple[int, int]:
    """``(forward_orf_nt, reverse_orf_nt)`` for one record."""
    return longest_orf(seq), longest_orf(reverse_complement(seq))


def _sample(
    records: Sequence[Tuple[str, str]], limit: int,
) -> Tuple[List[Tuple[str, str]], bool]:
    """Up to *limit* records, spread evenly across the input order.

    Even stride rather than the first *limit*: a transcriptome FASTA is usually
    sorted (by chromosome, or by gene id), so the leading records are not a
    representative sample of anything, and taking a stride costs nothing.
    """
    if len(records) <= limit:
        return list(records), False
    stride = len(records) / float(limit)
    return [records[int(i * stride)] for i in range(limit)], True


def check_orientation(
    records: Iterable[Tuple[str, str]],
    *,
    ratio_warn: float = ORF_RATIO_WARN,
    min_orf_nt: int = MIN_ORF_NT,
) -> List[Dict[str, object]]:
    """Warn about records whose dominant ORF is on the reverse strand.

    Parameters
    ----------
    records
        ``(record_id, sequence)`` pairs, as supplied by the user.

    Returns
    -------
    list of dict
        One entry per SUSPECT record — ``record_id``, ``forward_orf_nt``,
        ``reverse_orf_nt``, ``ratio``. Empty (the overwhelmingly common case)
        means nothing looked wrong, which is NOT the same as "the orientation
        was verified": a record with no meaningful ORF on either strand is
        reported as fine because this check has nothing to say about it.

    Emits a ``logger.warning`` per suspect record as it goes, capped at
    :data:`MAX_WARNINGS` with a summary line if there are more. Never raises and
    never refuses a run: it is a heuristic, some legitimate targets are
    non-coding, and the user is the one who knows what they supplied.
    """
    candidates = [
        (rid, seq) for rid, seq in records if len(seq) >= MIN_RECORD_NT
    ]
    if not candidates:
        return []
    sampled, was_sampled = _sample(candidates, MAX_RECORDS)
    if was_sampled:
        logger.info(
            "Input-orientation check: sampling %d of %d records (evenly "
            "spread). A whole input supplied in the wrong orientation is "
            "detected from a sample; a single reversed record among correct "
            "ones may not be.",
            len(sampled), len(candidates),
        )

    suspects: List[Dict[str, object]] = []
    for rid, seq in sampled:
        fwd, rev = orf_strand_ratio(seq[:MAX_BASES])
        if rev < min_orf_nt:
            continue                       # nothing worth calling an ORF
        if rev < ratio_warn * max(fwd, 1):
            continue                       # forward strand holds its own
        suspects.append({
            "record_id": rid,
            "forward_orf_nt": fwd,
            "reverse_orf_nt": rev,
            "ratio": rev / max(fwd, 1),
        })
        if len(suspects) <= MAX_WARNINGS:
            logger.warning(
                "Record %r looks like it may be in the WRONG ORIENTATION: its "
                "longest open reading frame is %d nt on the REVERSE strand but "
                "only %d nt on the strand you supplied (%.1fx). Probes are "
                "designed against the input strand read as sense mRNA, so "
                "probes designed against this record will NOT hybridise to the "
                "target mRNA — the sequence may be in genome/vRNA orientation "
                "(as a negative-sense RNA virus genome is). Supply the "
                "mRNA/cDNA sense strand, i.e. the reverse complement of this "
                "record. If this target really is non-coding or antisense, "
                "ignore this — the run is not affected and nothing was changed.",
                rid, rev, fwd, rev / max(fwd, 1),
            )
    if len(suspects) > MAX_WARNINGS:
        logger.warning(
            "... and %d further record(s) with the same reverse-strand ORF "
            "signature. The whole input may be in genome/vRNA orientation.",
            len(suspects) - MAX_WARNINGS,
        )
    return suspects
