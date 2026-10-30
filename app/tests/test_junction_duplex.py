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

"""The junction-duplex gate: OUR rule, on top of 10x's count rule.

10x's off-target sentence is a pure COUNT and says nothing about WHERE the
mismatches sit. RTL needs both halves to hybridise AND the ligase to seal the
nick between them, so a pair whose mismatches against another gene are all
terminal still presents a long perfectly matched duplex ACROSS the junction: it
hybridises, it ligates, and it puts reads on the wrong gene — while satisfying
the count rule comfortably.

These tests pin three separate claims:

  * WHAT IS MEASURED — the unbroken matched run on each side of the nick, per
    side and never summed (``junction_clean_duplex``).
  * THAT IT IS REACHABLE AT ALL — the dangerous probe carries 6+ mismatches per
    half and therefore produces NO HIT in a list capped at 5, so the scan has to
    seed on the junction rather than post-process the hits.
  * WHAT IT DOES TO A VERDICT — rejection under its own status, only ever on
    top of the count rule, and inherited by degenerate consensus parents.
"""
from __future__ import annotations

import random

import pytest

from flexprobe.constants import (
    OFF_TARGET_JUNCTION_DUPLEX_BP, OFF_TARGET_JUNCTION_DUPLEX_MIN_BP,
)
from flexprobe.constants import PROBE_ROLE_CONSENSUS, PROBE_ROLE_VARIANT
from flexprobe.models import ProbeCandidate
from flexprobe.screen import KmerScreen
from flexprobe.specificity import (
    JUNCTION_DUPLEX_KEY, STATUS_FAILED_JUNCTION_DUPLEX, evaluate_specificity,
    junction_clean_duplex, split_junction_records,
)

HALF = 25


def _rand(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _flip(base):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


def _mutate(seq, idxs):
    out = list(seq)
    for i in idxs:
        out[i] = _flip(out[i])
    return "".join(out)


def _cand(seq, pid="p1", gene="TARGETGENE"):
    return ProbeCandidate(gene_id=gene, probe_id=pid, sequence=seq,
                          start=0, end=len(seq))


def _jd_record(*sites, probe="p1"):
    """A hand-built junction-duplex record, as the k-mer screen emits it."""
    return {"qseqid": f"{probe}__J", JUNCTION_DUPLEX_KEY: True,
            "sites": [list(s) for s in sites]}


# ---------------------------------------------------------------------------
# What is measured
# ---------------------------------------------------------------------------

def test_clean_run_is_counted_outward_from_the_nick_on_each_side():
    q = "A" * 50
    # Mismatches at query offsets 15 and 34 leave 25-16 = 9 clean bases to the
    # left of the nick (offsets 16..24) and 9 to the right (25..33).
    subject = "C" * 16 + "A" * 18 + "C" * 16
    assert junction_clean_duplex(q, subject) == (9, 9)


def test_a_probe_can_be_perfect_on_one_side_and_still_be_measured_per_side():
    """15+3 is NOT 9+9, and the rule must not be able to confuse them.

    The ligase needs a stable duplex on EACH side of the nick; three matched
    bases on one side do not hold it there. A total-length form would have
    scored 15+3 (=18) as worse than 9+9 (=18 too, and dangerous), which is
    measuring homology rather than ligatability.
    """
    q = "A" * 50
    subject = "C" * 10 + "A" * 18 + "C" * 22       # 15 left, 3 right
    assert junction_clean_duplex(q, subject) == (15, 3)


def test_an_ambiguous_base_ends_the_run_on_either_side():
    """An N is not a matched base pair, so no duplex extends through one."""
    q = "A" * 50
    assert junction_clean_duplex(q, "A" * 20 + "N" + "A" * 29) == (4, 25)
    # Identical characters are not enough: N against N still pairs with nothing.
    assert junction_clean_duplex("A" * 24 + "NN" + "A" * 24,
                                 "A" * 24 + "NN" + "A" * 24) == (0, 0)


def test_a_subject_that_runs_out_ends_the_run_rather_than_extending_it():
    """Absent sequence is KNOWN, not unknown — the same reading as
    ``specificity._flank_mismatches`` gives a flank off a transcript end."""
    q = "A" * 50
    assert junction_clean_duplex(q, "A" * 30) == (25, 5)


def test_the_nick_position_is_explicit_and_defaults_to_the_midpoint():
    q = "A" * 50
    assert junction_clean_duplex(q, q) == (25, 25)
    assert junction_clean_duplex(q, q, 10) == (10, 40)


# ---------------------------------------------------------------------------
# What it does to a verdict
# ---------------------------------------------------------------------------

def _verdict(left, right, *, sseqid="OTHERGENE", **kw):
    c = _cand("A" * 50)
    evaluate_specificity(c, [_jd_record((sseqid, left, right))],
                         {"TARGETGENE"}, **kw)
    return c


def test_nine_and_nine_is_rejected_under_its_own_status():
    c = _verdict(9, 9)
    assert c.specificity_pass is False
    assert c.specificity_status == STATUS_FAILED_JUNCTION_DUPLEX
    # NOT the count rule's status. The whole claim of this gate is that it
    # rejects probes 10x's own sentence would have shipped, and a reader has to
    # be able to separate the two populations in results.csv.
    assert c.specificity_status != "failed_off_target"
    assert (c.offtarget_clean_left_bp, c.offtarget_clean_right_bp) == (9, 9)


def test_eight_and_eight_passes():
    """One base below the threshold on both sides, and the gate stays silent —
    so the rejection above is the threshold acting, not the gate firing on
    anything it is shown."""
    c = _verdict(8, 8)
    assert c.specificity_pass is True
    assert c.specificity_status == "passed"
    assert (c.offtarget_clean_left_bp, c.offtarget_clean_right_bp) == (8, 8)


def test_fifteen_and_three_passes_because_the_rule_is_per_side():
    """Longer TOTAL run than 9+9, and it must still pass.

    18 matched bases either way. The gate is per side because the ligase is:
    three bases on the RHS cannot hold a duplex there, so no product forms and
    the probe generates no false signal.
    """
    c = _verdict(15, 3)
    assert c.specificity_pass is True
    assert (c.offtarget_clean_left_bp, c.offtarget_clean_right_bp) == (15, 3)


def test_the_threshold_is_configurable_and_zero_disables_the_gate():
    assert _verdict(9, 9, off_target_junction_duplex_bp=0).specificity_pass
    assert _verdict(8, 8, off_target_junction_duplex_bp=8).specificity_pass is False
    assert _verdict(9, 9, off_target_junction_duplex_bp=11).specificity_pass


def test_the_evidence_is_reported_even_when_the_gate_is_off():
    """Same contract as offtarget_worst_junction_mm: reported, not only gated.

    A run that did not act on the number should still let a reviewer see it.
    """
    c = _verdict(20, 20, off_target_junction_duplex_bp=0)
    assert c.specificity_pass is True
    assert (c.offtarget_clean_left_bp, c.offtarget_clean_right_bp) == (20, 20)


def test_a_site_on_the_probes_own_target_never_counts():
    """A probe is tiled from its target, so its own gene is a perfect 25+25
    site by construction. Counting it would reject every probe ever designed."""
    c = _verdict(25, 25, sseqid="TARGETGENE")
    assert c.specificity_pass is True
    assert (c.offtarget_clean_left_bp, c.offtarget_clean_right_bp) == (0, 0)


def test_an_excused_isoform_never_counts_either():
    c = _cand("A" * 50)
    evaluate_specificity(
        c, [_jd_record(("NM_999999.1", 25, 25))], {"TARGETGENE"},
        extra_on_target_ids={"NM_999999.1"},
    )
    assert c.specificity_pass is True


def test_the_worst_non_target_site_is_the_one_reported():
    c = _cand("A" * 50)
    evaluate_specificity(
        c,
        [_jd_record(("TARGETGENE", 25, 25), ("GENE_A", 6, 6),
                    ("GENE_B", 12, 7), ("GENE_C", 8, 8))],
        {"TARGETGENE"},
    )
    # Largest min(left, right) wins: GENE_C's 8 beats GENE_B's 7 even though
    # GENE_B has the longer total run. The minimum is what the gate reads.
    assert (c.offtarget_clean_left_bp, c.offtarget_clean_right_bp) == (8, 8)


def test_no_junction_record_means_unmeasured_not_zero():
    """-1, and it must never collapse into a measured 0.

    The BLAST engine emits no junction scan at all, and neither does a hit cache
    written before this existed. Reporting 0 there would claim a clean junction
    for a probe nothing looked at.
    """
    c = _cand("A" * 50)
    evaluate_specificity(c, [], {"TARGETGENE"})
    assert c.offtarget_clean_left_bp == -1
    assert c.offtarget_clean_right_bp == -1
    assert c.specificity_pass is True

    # And a scan that ran and found nothing is a DIFFERENT answer.
    c2 = _cand("A" * 50)
    evaluate_specificity(c2, [_jd_record()], {"TARGETGENE"})
    assert (c2.offtarget_clean_left_bp, c2.offtarget_clean_right_bp) == (0, 0)


def test_the_count_rule_still_decides_first():
    """This gate only ever rejects MORE than 10x's rule, never instead of it.

    A pair the count rule already rejects keeps ``failed_off_target``, so
    ``failed_junction_duplex`` names exactly one population: probes the vendor's
    own sentence would have shipped.
    """
    from conftest import aligned_hit
    q = "A" * HALF
    close = "A" * HALF                      # 0 mismatches on both halves
    hits = [aligned_hit("OTHERGENE", q, close, 1, HALF, half=h) for h in "LR"]
    c = _cand("A" * 50)
    evaluate_specificity(c, hits + [_jd_record(("OTHERGENE", 25, 25))],
                         {"TARGETGENE"})
    assert c.specificity_pass is False
    assert c.specificity_status == "failed_off_target"
    # The evidence is still recorded on the way out — a rejected probe reports
    # what its junction looked like, like every other set-before-return field.
    assert (c.offtarget_clean_left_bp, c.offtarget_clean_right_bp) == (25, 25)


def test_a_junction_record_is_never_scored_as_an_alignment():
    """It has no qseq/sseq and no half. Scored as a hit it would read as a
    0-mismatch off-target on BOTH halves and reject every probe."""
    aln, jd = split_junction_records(
        [{"qseqid": "p1__L", "sseqid": "X"}, _jd_record(("Y", 9, 9))]
    )
    assert len(aln) == 1 and len(jd) == 1
    c = _cand("A" * 50)
    # Only the junction record, at a harmless length: nothing may reject here.
    evaluate_specificity(c, [_jd_record(("OTHERGENE", 6, 6))], {"TARGETGENE"})
    assert c.specificity_pass is True


# ---------------------------------------------------------------------------
# That it is reachable at all — the scan, end to end
# ---------------------------------------------------------------------------

def _screen(tmp_path, fasta_text, name="ref.fa"):
    fa = tmp_path / name
    fa.write_text(fasta_text)
    s = KmerScreen(cache_dir=str(tmp_path / f"idx_{name}"), mmap=False)
    s.prepare([str(fa)])
    return s


def _dangerous_pair(rng):
    """A probe and an off-target that the COUNT RULE CANNOT SEE.

    Six mismatches in each half, all pushed out to the oligo ends, plus one more
    on each side to trim the junction-spanning clean run to exactly 9+9. Eight
    mismatches per half is well past ``KmerScreen.MAX_MISMATCH_REPORTED`` (5),
    so this off-target produces no hit whatsoever — which is why the scan cannot
    be a post-processing pass over the hit list.
    """
    probe = _rand(50, rng)
    off = _mutate(probe, list(range(0, 6)) + list(range(44, 50)))
    off = _mutate(off, [15, 34])            # trim the clean run to 9 + 9
    return probe, off


def test_the_dangerous_off_target_produces_no_hit_and_is_still_caught(tmp_path):
    """The case the whole feature exists for.

    Both halves are 8 mismatches from the off-target, so 10x's ">= 5 mismatches
    in at least one half" is satisfied twice over and the hit list is empty of
    it. The junction is nonetheless perfectly matched for 9 bases each side, so
    the pair hybridises and ligates on a gene it was not designed for.
    """
    rng = random.Random(101)
    probe, off = _dangerous_pair(rng)
    filler = _rand(400, rng)
    s = _screen(tmp_path,
                f">TARGETGENE\n{filler}{probe}{filler}\n"
                f">OTHERGENE\n{filler}{off}{filler}\n")
    hits = s.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"]
    aln, jd = split_junction_records(hits)

    # 1) The off-target is invisible to the count rule: every alignment is on
    #    the probe's own gene.
    assert aln, "the probe's own locus must still be found"
    assert {h["sseqid"] for h in aln} == {"TARGETGENE"}

    # 2) The junction scan sees it, per side.
    assert len(jd) == 1
    assert ["OTHERGENE", 9, 9] in jd[0]["sites"]

    # 3) And the verdict follows.
    c = _cand(probe)
    evaluate_specificity(c, hits, {"TARGETGENE"})
    assert c.specificity_status == STATUS_FAILED_JUNCTION_DUPLEX
    assert (c.offtarget_clean_left_bp, c.offtarget_clean_right_bp) == (9, 9)

    # Control: with the gate off, the same probe and the same hits pass — so
    # the rejection above is this gate and not something else in the screen.
    c2 = _cand(probe)
    evaluate_specificity(c2, hits, {"TARGETGENE"},
                         off_target_junction_duplex_bp=0)
    assert c2.specificity_pass is True


def test_the_scan_searches_the_minus_strand_too(tmp_path):
    """This is the mechanism, exercised at the DEFAULT (``reference_strand=
    "both"``) setting the pipeline still ships: pin that the scan is even
    CAPABLE of finding a minus-strand site, and that when it does, the counts
    come back in the PROBE's frame, not the reference's.

    Whether that site is a REAL off-target is a separate, conditional
    question. For a sense-oriented cDNA reference, OTHERGENE's minus-strand match means
    OTHERGENE's own sense sequence is ~identical to the probe (RC of RC),
    which is NOT complementary to it, so no duplex forms there and the site
    is a phantom. It IS a genuine off-target when the reference is a genome
    (OTHERGENE really transcribes off the opposite strand) or a negative-
    sense-viral record. Both are still-supported use cases this codebase
    cannot currently tell apart from a sense cDNA file by inspection, which
    is why "both" stays the default and ``reference_strand="plus"`` (see
    constants.REFERENCE_STRAND_DEFAULT) is an explicit opt-in, not this."""
    rng = random.Random(202)
    probe = _rand(50, rng)
    # Asymmetric on purpose — 9 clean bases 5' of the nick and 19 3' of it — so
    # a frame that came back reversed would be caught rather than look right.
    off = _mutate(probe, list(range(0, 6)) + [15] + list(range(44, 50)))
    filler = _rand(400, rng)
    s = _screen(tmp_path,
                f">TARGETGENE\n{filler}{probe}{filler}\n"
                f">OTHERGENE\n{filler}{_revcomp(off)}{filler}\n")
    hits = s.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"]
    sites = dict((n, (a, b)) for n, a, b in
                 split_junction_records(hits)[1][0]["sites"])
    # Not (19, 9): reversing the frame must not swap the probe's own 5' and 3'.
    assert sites["OTHERGENE"] == (9, 19)


def test_every_site_the_scan_returns_clears_the_seeding_floor(tmp_path):
    """The seed is an exact 12-mer centred on the nick, so 6+6 is the floor by
    construction — which is why a threshold below it is refused rather than
    silently run at 6."""
    rng = random.Random(303)
    ref = _rand(4000, rng)
    s = _screen(tmp_path, f">REC1\n{ref}\n")
    probe = ref[1000:1050]
    hits = s.search_chunk([_cand(probe)], 0, str(tmp_path))["p1"]
    sites = split_junction_records(hits)[1][0]["sites"]
    assert sites
    for _name, left, right in sites:
        assert left >= OFF_TARGET_JUNCTION_DUPLEX_MIN_BP
        assert right >= OFF_TARGET_JUNCTION_DUPLEX_MIN_BP


def test_a_run_stops_at_a_record_boundary(tmp_path):
    """Two transcripts are not one molecule. A duplex cannot walk from the end
    of one record into the start of the next, whatever the concatenated
    reference looks like in memory."""
    rng = random.Random(404)
    probe = _rand(50, rng)
    s = _screen(tmp_path,
                f">A\n{_rand(300, rng)}{probe[:30]}\n"
                f">B\n{probe[30:]}{_rand(300, rng)}\n")
    hits = s.search_chunk([_cand(probe, gene="A")], 0, str(tmp_path))["p1"]
    for _name, left, right in split_junction_records(hits)[1][0]["sites"]:
        # Record A holds only the first 30 bases of the window, so at most 5
        # can sit 3' of the nick.
        assert right <= 5


def test_an_n_across_the_nick_is_unmeasured_rather_than_clean(tmp_path):
    """The seed cannot be built, so the scan never searched. Saying "found
    nothing" for a window nobody looked at is the failure this avoids."""
    rng = random.Random(505)
    ref = _rand(2000, rng)
    probe = ref[500:524] + "N" + ref[525:550]
    s = _screen(tmp_path, f">REC1\n{ref}\n")
    hits = s.search_chunk([_cand(probe)], 0, str(tmp_path)).get("p1", [])
    assert split_junction_records(hits)[1] == []
    c = _cand(probe)
    evaluate_specificity(c, hits, {"REC1"})
    assert c.offtarget_clean_left_bp == -1


# ---------------------------------------------------------------------------
# Degenerate consensus parents inherit the verdict
# ---------------------------------------------------------------------------

def _variant(pid, parent, **kw):
    return ProbeCandidate(
        gene_id="g", probe_id=pid, sequence="A" * 50, start=0, end=50,
        probe_role=PROBE_ROLE_VARIANT, parent_probe_id=parent, **kw)


def test_a_consensus_parent_inherits_a_junction_duplex_rejection():
    """The parent's string is never screened — the concrete variants are — and
    every variant goes into the same pool, so one ligatable variant
    contaminates it. The parent must name which variant and WHY."""
    from flexprobe.pipeline import _rollup_degenerate_specificity

    parent = ProbeCandidate(gene_id="g", probe_id="g_t000000",
                            sequence="A" * 50, start=0, end=50,
                            probe_role=PROBE_ROLE_CONSENSUS, n_variants=2)
    good = _variant("g_t000000_v00", "g_t000000", specificity_pass=True,
                    specificity_status="passed", specificity_margin=9.0,
                    offtarget_clean_left_bp=6, offtarget_clean_right_bp=6)
    bad = _variant("g_t000000_v01", "g_t000000", specificity_pass=False,
                   specificity_status=STATUS_FAILED_JUNCTION_DUPLEX,
                   specificity_margin=0.0,
                   offtarget_clean_left_bp=9, offtarget_clean_right_bp=14)

    assert _rollup_degenerate_specificity([parent, good, bad]) == 0
    assert parent.specificity_pass is False
    assert "g_t000000_v01" in parent.specificity_status
    # The reason survives the rollup: a parent that just said "off_target"
    # would send the reviewer looking for a mismatch count that is not the
    # thing that rejected it.
    assert STATUS_FAILED_JUNCTION_DUPLEX in parent.specificity_status
    # The pool is only as specific as its least specific member, on this axis
    # too: the parent reports the WORST child's junction, not the first one's.
    assert (parent.offtarget_clean_left_bp,
            parent.offtarget_clean_right_bp) == (9, 14)


def test_a_parent_whose_variants_all_pass_still_reports_the_worst_junction():
    from flexprobe.pipeline import _rollup_degenerate_specificity

    parent = ProbeCandidate(gene_id="g", probe_id="g_t000000",
                            sequence="A" * 50, start=0, end=50,
                            probe_role=PROBE_ROLE_CONSENSUS, n_variants=2)
    kids = [
        _variant("g_t000000_v00", "g_t000000", specificity_pass=True,
                 specificity_status="passed", specificity_margin=4.0,
                 offtarget_clean_left_bp=6, offtarget_clean_right_bp=6),
        _variant("g_t000000_v01", "g_t000000", specificity_pass=True,
                 specificity_status="passed", specificity_margin=9.0,
                 offtarget_clean_left_bp=8, offtarget_clean_right_bp=7),
    ]
    assert _rollup_degenerate_specificity([parent] + kids) == 1
    assert parent.specificity_pass is True
    assert (parent.offtarget_clean_left_bp,
            parent.offtarget_clean_right_bp) == (8, 7)


def test_an_unmeasured_variant_does_not_read_as_a_clean_one():
    """-1 is not 0. A child the scan could not measure must not pull the
    parent's evidence down to "measured, and nothing found"."""
    from flexprobe.pipeline import _rollup_degenerate_specificity

    parent = ProbeCandidate(gene_id="g", probe_id="g_t000000",
                            sequence="A" * 50, start=0, end=50,
                            probe_role=PROBE_ROLE_CONSENSUS, n_variants=2)
    kids = [
        _variant("g_t000000_v00", "g_t000000", specificity_pass=True,
                 specificity_status="passed", specificity_margin=4.0),
        _variant("g_t000000_v01", "g_t000000", specificity_pass=True,
                 specificity_status="passed", specificity_margin=9.0,
                 offtarget_clean_left_bp=7, offtarget_clean_right_bp=7),
    ]
    assert _rollup_degenerate_specificity([parent] + kids) == 1
    assert (parent.offtarget_clean_left_bp,
            parent.offtarget_clean_right_bp) == (7, 7)

    # And a parent NONE of whose children were measured stays unmeasured.
    parent2 = ProbeCandidate(gene_id="g", probe_id="g_t000001",
                             sequence="A" * 50, start=0, end=50,
                             probe_role=PROBE_ROLE_CONSENSUS, n_variants=1)
    _rollup_degenerate_specificity([parent2, _variant(
        "g_t000001_v00", "g_t000001", specificity_pass=True,
        specificity_status="passed")])
    assert parent2.offtarget_clean_left_bp == -1


# ---------------------------------------------------------------------------
# Configuration surface
# ---------------------------------------------------------------------------

def test_the_default_is_on_at_the_measured_optimum():
    from flexprobe.models import SpecificityThresholds
    assert OFF_TARGET_JUNCTION_DUPLEX_BP == 9
    assert SpecificityThresholds().off_target_junction_duplex_bp == 9


def test_a_threshold_below_the_scans_floor_is_refused_not_clamped(tmp_path):
    """Running the gate at 6 while the config says 3 would make the knob and
    its effect disagree silently — which is what the recount basis exists to
    prevent elsewhere."""
    from flexprobe.errors import ConfigError
    from flexprobe.models import RunConfig, SpecificityThresholds

    cfg = RunConfig(fasta_path="x.fa", output_dir=str(tmp_path),
                    specificity=SpecificityThresholds(
                        off_target_junction_duplex_bp=3))
    with pytest.raises(ConfigError, match="junction"):
        cfg.validate()

    # 0 is fine — that is how you ask for 10x's published rule and nothing else.
    cfg.specificity.off_target_junction_duplex_bp = 0
    cfg.validate()


def test_cli_flags_follow_the_off_target_same_subject_pattern():
    from flexprobe.cli import _build_run_config, build_parser

    p = build_parser()
    base = ["design", "--fasta", "x.fa"]
    assert _build_run_config(p.parse_args(base)).specificity \
        .off_target_junction_duplex_bp == OFF_TARGET_JUNCTION_DUPLEX_BP
    assert _build_run_config(
        p.parse_args(base + ["--off-target-junction-duplex", "11"])
    ).specificity.off_target_junction_duplex_bp == 11
    for off in (["--no-off-target-junction-duplex"],
                ["--off-target-junction-duplex", "0"]):
        assert _build_run_config(
            p.parse_args(base + off)
        ).specificity.off_target_junction_duplex_bp == 0


def test_results_csv_carries_both_columns_for_every_probe(tmp_path):
    """Emitted for every row, not only for rejected ones — and -1 for a probe
    whose junction was never measured, which is not the same as 0."""
    import csv

    from flexprobe.writers.results_table import write_results_table

    measured = _cand("A" * 50, pid="p_measured")
    measured.offtarget_clean_left_bp = 9
    measured.offtarget_clean_right_bp = 12
    out = str(tmp_path / "results.csv")
    write_results_table([measured, _cand("A" * 50, pid="p_unmeasured")], out)
    rows = {r["probe_id"]: r for r in csv.DictReader(open(out))}
    assert rows["p_measured"]["offtarget_clean_left_bp"] == "9"
    assert rows["p_measured"]["offtarget_clean_right_bp"] == "12"
    assert rows["p_unmeasured"]["offtarget_clean_left_bp"] == "-1"
    assert rows["p_unmeasured"]["offtarget_clean_right_bp"] == "-1"


def test_the_verdict_survives_the_resume_cache():
    """A resumed chunk must not report -1 for a probe that WAS measured — the
    same class of silent provenance loss as a lost recount basis."""
    from flexprobe import screen_cache

    for f in ("offtarget_clean_left_bp", "offtarget_clean_right_bp"):
        assert f in screen_cache.VERDICT_FIELDS
    # And the on-disk format changed, so a chunk written before the scan
    # existed can never be read back as one that carries it.
    assert screen_cache.CACHE_SCHEMA >= 4


# ---------------------------------------------------------------------------
# B5: the manifest must say whether the gate actually ran, and BLAST users
# must be told out loud, not left to infer it from -1 in results.csv.
# ---------------------------------------------------------------------------

import shutil as _shutil

_HAS_BLAST = (_shutil.which("blastn") is not None
             and _shutil.which("makeblastdb") is not None)

# A 120-bp target, same shape as the BLAST e2e fixtures elsewhere in the suite.
_ENGINE_TARGET = (
    "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGTTACCGATGCATGC"
    "ATGCATGCATGCATGCATGCAACGTACGTACGTACGTACGTTGCATGCATG"
    "CATGCATGCATG"
)
_UNRELATED_OFFTARGET = "TTTTTTTTAAAAAAAACCCCCCCCGGGGGGGG" * 4


def _engine_cfg(fa, out, db, engine):
    from flexprobe.models import (
        ExportSettings, ModalityDefinition, RunConfig, ScreenResources,
    )
    return RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=False, db_paths=[db], engine=engine),
        export=ExportSettings(overwrite=True), discovery_only=True,
    )


def test_manifest_junction_duplex_measured_true_under_kmer(tmp_path):
    """The default engine actually runs KmerScreen._junction_scan, so the
    manifest must say the gate was measured."""
    from flexprobe.pipeline import run_pipeline
    from flexprobe.screen import ENGINE_KMER

    fa = tmp_path / "target.fasta"
    fa.write_text(f">GENE1\n{_ENGINE_TARGET}\n")
    dbfa = tmp_path / "off.fasta"
    dbfa.write_text(f">UNRELATED\n{_UNRELATED_OFFTARGET}\n")
    out = tmp_path / "out"

    result = run_pipeline(_engine_cfg(fa, out, str(dbfa), ENGINE_KMER))
    assert result["manifest"].params["junction_duplex_measured"] is True
    assert result["manifest"].params["screen_engine"] == ENGINE_KMER


def test_blast_engine_with_the_gate_on_is_refused(tmp_path):
    """BLAST + a live junction gate is now a CONFIG ERROR, not a warning.

    Warn-and-continue is the wrong shape: the manifest would say
    ``junction_duplex_measured: false`` and a warning would go to the log, but
    the option is documented with a default of 9 and no engine caveat, so the
    run would still SAY it applied a 9 bp gate. A log line the reader never
    sees does not make that true, and every downstream artefact (results.csv,
    any methods text written from it) would carry the unearned claim. So the
    combination is refused with an
    actionable message; ``blast.py``'s module docstring says why BLAST cannot
    supply the records instead. No BLAST binary is needed to test a refusal
    that happens before anything runs.
    """
    from flexprobe.errors import ConfigError
    from flexprobe.screen import ENGINE_BLAST, ENGINE_KMER

    cfg = _engine_cfg(tmp_path / "t.fa", tmp_path / "out", "db", ENGINE_BLAST)
    with pytest.raises(ConfigError) as ei:
        cfg.validate()
    msg = str(ei.value)
    assert "off_target_junction_duplex_bp" in msg
    # Both ways out, named: switch engine, or turn the gate off on purpose.
    assert ENGINE_KMER in msg and "0" in msg


@pytest.mark.skipif(not _HAS_BLAST,
                    reason="ncbi-blast+ (makeblastdb) not installed")
def test_manifest_junction_duplex_measured_false_under_blast(tmp_path):
    """With the gate explicitly off, a BLAST run still records that it was not
    measured — the manifest must never imply a gate the engine cannot run."""
    from flexprobe.blast import build_db, discover_blast
    from flexprobe.models import SpecificityThresholds
    from flexprobe.pipeline import run_pipeline
    from flexprobe.screen import ENGINE_BLAST

    fa = tmp_path / "target.fasta"
    fa.write_text(f">GENE1\n{_ENGINE_TARGET}\n")
    dbfa = tmp_path / "off.fasta"
    dbfa.write_text(f">UNRELATED\n{_UNRELATED_OFFTARGET}\n")
    db = str(tmp_path / "offdb")
    build_db(discover_blast()["makeblastdb"], str(dbfa), db, title="off")
    out = tmp_path / "out"

    cfg = _engine_cfg(fa, out, db, ENGINE_BLAST)
    cfg.specificity = SpecificityThresholds(off_target_junction_duplex_bp=0)
    result = run_pipeline(cfg)

    assert result["manifest"].params["junction_duplex_measured"] is False
    assert result["manifest"].params["screen_engine"] == ENGINE_BLAST


@pytest.mark.skipif(not _HAS_BLAST,
                    reason="ncbi-blast+ (blastn/makeblastdb) not installed")
def test_no_warning_when_gate_already_disabled_under_blast(tmp_path, caplog):
    """--off-target-junction-duplex 0 already means no gate under EITHER
    engine, so warning about BLAST specifically would be noise, not signal."""
    import logging

    from flexprobe.blast import build_db, discover_blast
    from flexprobe.models import SpecificityThresholds
    from flexprobe.pipeline import run_pipeline
    from flexprobe.screen import ENGINE_BLAST

    fa = tmp_path / "target.fasta"
    fa.write_text(f">GENE1\n{_ENGINE_TARGET}\n")
    dbfa = tmp_path / "off.fasta"
    dbfa.write_text(f">UNRELATED\n{_UNRELATED_OFFTARGET}\n")
    db = str(tmp_path / "offdb")
    build_db(discover_blast()["makeblastdb"], str(dbfa), db, title="off")
    out = tmp_path / "out"

    cfg = _engine_cfg(fa, out, db, ENGINE_BLAST)
    cfg.specificity = SpecificityThresholds(off_target_junction_duplex_bp=0)
    with caplog.at_level(logging.WARNING, logger="flexprobe.pipeline"):
        run_pipeline(cfg)
    assert not any("junction-duplex" in r.getMessage().lower()
                  for r in caplog.records)
