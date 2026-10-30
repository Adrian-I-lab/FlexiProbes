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

"""Minus-strand phantom off-targets.

At ``reference_strand="both"``, both engines count a reference hit on the
MINUS strand (a reverse-complement match) as an off-target. A probe half is
single-stranded DNA that hybridises to sense-strand reference sequence, so
against a sense cDNA reference a minus-only hit is sequence the probe cannot
physically bind: a "phantom" off-target that rejects a good probe for no
biological reason.

These tests build a reverse-complement decoy and run it against the real
engines: a reference carrying ONLY the reverse complement of a probe's target
window must stop being reported as an off-target hit once the caller has
opted into ``reference_strand="plus"`` — which is only correct for a
sense-oriented cDNA/transcriptome reference, so it is NOT the default (see
``constants.REFERENCE_STRAND_DEFAULT``). A control decoy carrying the SAME
strand must still be caught, in every mode, on every engine — the fix must
never hide a genuine off-target.
"""
from __future__ import annotations

import random
import shutil

import pytest

from flexprobe.blast import build_db, discover_blast
from flexprobe.constants import REFERENCE_STRAND_PLUS
from flexprobe.models import ProbeCandidate
from flexprobe.screen import BlastScreen, KmerScreen
from flexprobe.specificity import evaluate_specificity, split_junction_records

TARGET_ID = "TARGETGENE"

needs_blast = pytest.mark.skipif(
    shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
    reason="ncbi-blast+ (blastn/makeblastdb) not installed",
)


def _rand(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _cand(seq, pid="p1"):
    return ProbeCandidate(gene_id=TARGET_ID, probe_id=pid, sequence=seq,
                          start=0, end=len(seq))


def _kmer_screen(tmp_path, fasta_text, name, **kw):
    fa = tmp_path / name
    fa.write_text(fasta_text)
    s = KmerScreen(cache_dir=str(tmp_path / f"idx_{name}"), mmap=False, **kw)
    s.prepare([str(fa)])
    return s


def _blast_screen(tmp_path, fasta_text, name, **kw):
    fa = tmp_path / name
    fa.write_text(fasta_text)
    db = str(tmp_path / f"{name}.db")
    build_db(discover_blast()["makeblastdb"], str(fa), db, title="t")
    s = BlastScreen(**kw)
    s.prepare([db])
    return s


def _decoys(rng):
    """A 50 nt probe (an exact tile, 0 mismatches to itself) and two decoy
    references, built exactly as the verify doc reproduced the bug: 200 bp of
    unrelated random filler either side of the probe or its reverse
    complement, under a record id distinct from the probe's own gene.
    """
    probe = _rand(50, rng)
    filler = _rand(200, rng)
    rc_only = f">DECOY_RC\n{filler}{_revcomp(probe)}{filler}\n"
    same_strand = f">DECOY_SENSE\n{filler}{probe}{filler}\n"
    return probe, rc_only, same_strand


# ---------------------------------------------------------------------------
# k-mer engine — half search (screen.py:457, the "quiet" channel)
# ---------------------------------------------------------------------------

def test_kmer_half_search_ignores_an_rc_only_decoy_under_plus(tmp_path):
    rng = random.Random(1001)
    probe, rc_only, _ = _decoys(rng)

    both = _kmer_screen(tmp_path, rc_only, "both.fa")
    plus = _kmer_screen(tmp_path, rc_only, "plus.fa",
                        reference_strand=REFERENCE_STRAND_PLUS)

    hits_both, _ = split_junction_records(
        both.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"])
    hits_plus, _ = split_junction_records(
        plus.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"])

    assert hits_both, "both (default) must still find the minus-strand hit"
    assert not hits_plus, (
        "plus must not report a hit on sequence the probe cannot bind")


def test_kmer_half_search_still_catches_a_same_strand_decoy_under_plus(tmp_path):
    """Control: the fix must never hide a genuine (same-strand) off-target."""
    rng = random.Random(1002)
    probe, _, same_strand = _decoys(rng)

    plus = _kmer_screen(tmp_path, same_strand, "plus.fa",
                        reference_strand=REFERENCE_STRAND_PLUS)
    hits, _ = split_junction_records(
        plus.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"])
    assert hits, "a genuine same-strand off-target must still be caught"


# ---------------------------------------------------------------------------
# k-mer engine — junction scan (screen.py:547, the DOMINANT channel: it has
# no per-site strand field, so no post-filter of the hit list could fix this)
# ---------------------------------------------------------------------------

def test_junction_scan_ignores_an_rc_only_decoy_under_plus(tmp_path):
    rng = random.Random(2001)
    probe, rc_only, _ = _decoys(rng)

    both = _kmer_screen(tmp_path, rc_only, "both.fa")
    plus = _kmer_screen(tmp_path, rc_only, "plus.fa",
                        reference_strand=REFERENCE_STRAND_PLUS)

    _, jd_both = split_junction_records(
        both.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"])
    _, jd_plus = split_junction_records(
        plus.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"])

    assert jd_both and jd_both[0]["sites"], (
        "both (default) must still find the minus-strand junction site")
    assert jd_plus and not jd_plus[0]["sites"], (
        "plus must not report a junction the probe cannot hybridise to, let "
        "alone ligate across")


def test_junction_scan_still_catches_a_same_strand_decoy_under_plus(tmp_path):
    rng = random.Random(2002)
    probe, _, same_strand = _decoys(rng)

    plus = _kmer_screen(tmp_path, same_strand, "plus.fa",
                        reference_strand=REFERENCE_STRAND_PLUS)
    _, jd = split_junction_records(
        plus.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"])
    assert jd and jd[0]["sites"], (
        "a genuine same-strand junction duplex must still be caught")


# ---------------------------------------------------------------------------
# BLAST engine (blast.py:396-415 / screen.py's BlastScreen.search_chunk).
# No junction scan exists on this engine (verify doc §5), so only the
# alignment-hit channel is exercised here.
# ---------------------------------------------------------------------------

@needs_blast
def test_blast_engine_ignores_an_rc_only_decoy_under_plus(tmp_path):
    rng = random.Random(3001)
    probe, rc_only, _ = _decoys(rng)

    # The pushdown floors are passed because the PIPELINE always passes them,
    # and without them this test exercises a configuration that never ships.
    # It matters under the complete-coverage settings: word_size 5 with dust
    # and soft masking off finds incidental 6 nt and 12 nt matches on the plus
    # strand, so an unfiltered BlastScreen does not return an empty hit list
    # even though
    # the reverse-complement decoy is genuinely suppressed. Both of those junk
    # hits sit at 24% and 48% query coverage against a 56% floor, so the real
    # screen drops them server-side and never sees them.
    # NOTE the argument names: BlastScreen takes min_pident / min_aln /
    # derive_qcov and has a **_ignored catch-all, so passing perc_identity or
    # qcov_hsp_perc here is silently swallowed and changes nothing.
    from flexprobe.models import SpecificityThresholds
    from flexprobe.pipeline import _derive_qcov_hsp_perc, _off_target_pushdown
    _pid, _aln = _off_target_pushdown(SpecificityThresholds(), 25)
    floors = dict(min_pident=_pid, min_aln=_aln, derive_qcov=_derive_qcov_hsp_perc)

    both = _blast_screen(tmp_path, rc_only, "both", **floors)
    plus = _blast_screen(tmp_path, rc_only, "plus",
                         reference_strand=REFERENCE_STRAND_PLUS, **floors)

    # A probe with zero hits is simply absent from the dict (search_chunk only
    # ever does hits.setdefault(...).append(...)), not present with an empty
    # list — which is itself part of what "plus fully suppresses it" means.
    hits_both = both.search_chunk([_cand(probe)], 0, str(tmp_path)).get("p1", [])
    hits_plus = plus.search_chunk([_cand(probe)], 0, str(tmp_path)).get("p1", [])

    assert hits_both, "both (default, blastn's own default) must still find it"
    assert not hits_plus, (
        "-strand plus must stop blastn from even searching the minus strand")


@needs_blast
def test_blast_engine_still_catches_a_same_strand_decoy_under_plus(tmp_path):
    rng = random.Random(3002)
    probe, _, same_strand = _decoys(rng)

    plus = _blast_screen(tmp_path, same_strand, "plus",
                         reference_strand=REFERENCE_STRAND_PLUS)
    hits = plus.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"]
    assert hits, "a genuine same-strand off-target must still be caught"


# ---------------------------------------------------------------------------
# Verdict level — the thing a user actually experiences: with
# reference_strand="plus" a probe hit only by an RC decoy passes, and a probe
# that should be rejected still is.
# ---------------------------------------------------------------------------

def test_evaluate_specificity_stops_rejecting_for_an_rc_only_decoy(tmp_path):
    """Today's default ("both") keeps rejecting — that is the deliberately
    conservative, unchanged behaviour a genome or negative-sense-viral
    reference still needs. Opting into "plus" for a reference asserted to be
    sense cDNA is what fixes the phantom rejection."""
    rng = random.Random(4001)
    probe, rc_only, _ = _decoys(rng)

    both = _kmer_screen(tmp_path, rc_only, "both.fa")
    c_both = _cand(probe)
    evaluate_specificity(
        c_both, both.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"],
        {TARGET_ID},
    )
    assert c_both.specificity_pass is False
    assert c_both.specificity_status == "failed_off_target"

    plus = _kmer_screen(tmp_path, rc_only, "plus.fa",
                        reference_strand=REFERENCE_STRAND_PLUS)
    c_plus = _cand(probe)
    evaluate_specificity(
        c_plus, plus.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"],
        {TARGET_ID},
    )
    assert c_plus.specificity_pass is True


def test_evaluate_specificity_still_rejects_a_same_strand_decoy_under_plus(tmp_path):
    rng = random.Random(4002)
    probe, _, same_strand = _decoys(rng)

    plus = _kmer_screen(tmp_path, same_strand, "plus.fa",
                        reference_strand=REFERENCE_STRAND_PLUS)
    c = _cand(probe)
    evaluate_specificity(
        c, plus.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"],
        {TARGET_ID},
    )
    assert c.specificity_pass is False
    assert c.specificity_status == "failed_off_target"
