"""The input-orientation check.

The trap this guards: probes are built as reverse complements of the input
window, so they bind the strand the user supplied and nothing else. A FASTA in
genome/vRNA orientation — the natural download for a negative-sense RNA virus
such as Oropouche, this tool's flagship target — passes every filter, the
chemistry and the off-target screen, and produces a complete, well-scored,
orderable panel that detects nothing in the assay. Nothing else in the pipeline
can notice.

These tests pin the three behaviours that matter: it fires on a reversed coding
sequence, it is silent on the correct orientation, and it is silent on input
with no coding strand at all. The last is not a nicety — a check that cried
wolf on UTR fragments and non-coding targets would be turned off, and then it
would not be there for the case it exists for.
"""
from __future__ import annotations

import logging
import random

from flexprobe.orientation import (
    MIN_ORF_NT,
    ORF_RATIO_WARN,
    check_orientation,
    longest_orf,
    orf_strand_ratio,
)
from flexprobe.tiling import reverse_complement

# All 61 sense codons. Drawing uniformly from the full set is what makes the
# fixture honest: the reverse strand then accumulates chance stops at the rate
# real coding sequence does, so "the reverse strand is quiet" is a property of
# the construction rather than something the fixture was rigged for. A
# hand-picked codon subset can produce a sequence with NO stops on either
# strand (a 1.0 ratio), which is not what a CDS looks like.
_SENSE_CODONS = [
    a + b + c
    for a in "ACGT" for b in "ACGT" for c in "ACGT"
    if a + b + c not in ("TAA", "TAG", "TGA")
]


def _synthetic_cds(n_codons: int = 500, seed: int = 11) -> str:
    """ATG + *n_codons* uniform sense codons + TAA. Deterministic."""
    rng = random.Random(seed)
    body = "".join(rng.choice(_SENSE_CODONS) for _ in range(n_codons))
    return "ATG" + body + "TAA"


def _synthetic_noncoding(length: int = 3000, seed: int = 23) -> str:
    """Sequence with no coding strand — the "say nothing" case.

    Built by inserting an in-frame stop codon every few codons in all three
    frames' worth of positions, which is what a real UTR or intron looks like
    to this check: chance ORFs of a few dozen codons on both strands.
    """
    rng = random.Random(seed)
    out = []
    while len("".join(out)) < length:
        out.append("".join(rng.choice("ACGT") for _ in range(rng.randint(9, 30))))
        out.append(rng.choice(["TAA", "TAG", "TGA"]))
        out.append(rng.choice(["TTA", "CTA", "TCA"]))   # stops on the reverse
    return "".join(out)[:length]


# ---------------------------------------------------------------------------
# longest_orf itself
# ---------------------------------------------------------------------------

def test_longest_orf_counts_the_stop_codon():
    assert longest_orf("ATG" + "AAA" * 4 + "TAA") == 18


def test_longest_orf_counts_an_orf_that_runs_off_the_three_prime_end():
    """A deposited genome segment is often a CDS with no in-frame stop inside
    the record. Scoring those at ZERO would read a 6.7 kb L segment with a
    completely open frame as having no ORF at all.
    """
    seq = "ATG" + "AAA" * 200                    # no stop anywhere
    assert longest_orf(seq) == len(seq)


def test_longest_orf_ignores_ambiguity_codes_rather_than_breaking_on_them():
    """N must neither open nor close an ORF.

    Breaking at every ambiguous codon destroys the signal on CORRECTLY oriented
    N-containing consensus inputs (a multi-kb ORF collapses to a few codons).
    Treating N as neutral
    lengthens both strands roughly equally, which biases toward silence.
    """
    clean = "ATG" + "AAA" * 30 + "TAA"
    with_n = "ATG" + "AAA" * 10 + "ANA" + "AAA" * 19 + "TAA"
    assert longest_orf(with_n) == longest_orf(clean)


# ---------------------------------------------------------------------------
# The three behaviours that matter
# ---------------------------------------------------------------------------

def test_fires_on_a_reverse_complemented_coding_sequence(caplog):
    cds = _synthetic_cds()
    flipped = reverse_complement(cds)

    with caplog.at_level(logging.WARNING, logger="flexprobe.orientation"):
        suspects = check_orientation([("SEG_L_vRNA", flipped)])

    assert [s["record_id"] for s in suspects] == ["SEG_L_vRNA"]
    assert float(suspects[0]["ratio"]) >= ORF_RATIO_WARN
    text = caplog.text
    assert "SEG_L_vRNA" in text
    # The message has to say what it COSTS, not just that something looks odd.
    assert "will NOT hybridise" in text
    assert "genome/vRNA orientation" in text


def test_silent_on_the_correct_orientation(caplog):
    cds = _synthetic_cds()
    with caplog.at_level(logging.WARNING, logger="flexprobe.orientation"):
        assert check_orientation([("SEG_L_mRNA", cds)]) == []
    assert caplog.text == ""


def test_silent_on_a_non_coding_input(caplog):
    """No meaningful ORF on either strand: say nothing rather than guess."""
    seq = _synthetic_noncoding()
    fwd, rev = orf_strand_ratio(seq)
    assert max(fwd, rev) < MIN_ORF_NT, (
        f"fixture is not actually non-coding (fwd={fwd}, rev={rev})"
    )
    with caplog.at_level(logging.WARNING, logger="flexprobe.orientation"):
        assert check_orientation([("UTR_FRAGMENT", seq)]) == []
        assert check_orientation([("UTR_FRAGMENT_RC", reverse_complement(seq))]) == []
    assert caplog.text == ""


def test_silent_on_a_short_amplicon_however_lopsided_the_ratio():
    """The absolute ORF floor, not the ratio, is what keeps this quiet.

    A 400 bp amplicon can easily show a 3x strand ratio between two ORFs that
    are both noise. MIN_ORF_NT is what stops that being a warning.
    """
    amplicon = _synthetic_cds(n_codons=60)       # 186 nt ORF, well under the floor
    assert len(amplicon) < MIN_ORF_NT
    assert check_orientation([("AMPLICON", reverse_complement(amplicon))]) == []


# ---------------------------------------------------------------------------
# Cost control
# ---------------------------------------------------------------------------

def test_samples_rather_than_scanning_a_whole_transcriptome(caplog):
    """A whole-transcriptome run must not pay for a six-frame scan of everything.

    The check is documented as sampling; this pins that it SAYS so, because a
    silent sample would let a reader over-read an empty result as "every record
    was checked".
    """
    rng = random.Random(5)
    records = [
        (f"ENST{i:011d}", "".join(rng.choice("ACGT") for _ in range(1200)))
        for i in range(500)
    ]
    with caplog.at_level(logging.INFO, logger="flexprobe.orientation"):
        check_orientation(records)
    assert "sampling" in caplog.text


def test_records_below_the_length_floor_are_skipped():
    short = [("TINY", "ATG" + "AAA" * 20)]
    assert check_orientation(short) == []
