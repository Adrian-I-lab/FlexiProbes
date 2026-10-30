"""Direct tests for the specificity engine."""
from conftest import FakeSubjects, aligned_hit
from flexprobe.models import ProbeCandidate
from flexprobe.specificity import (
    StrainRule, batch_evaluate, classify_hit, evaluate_specificity, id_matches,
)


def _hit(sseqid, pident, length):
    return {"sseqid": sseqid, "pident": pident, "length": length}


# --- id_matches: the boundary-aware target matcher --------------------------

def test_id_matches_exact_and_delimiter():
    assert id_matches("OROV_BeAn19991", "OROV_BeAn19991")
    assert id_matches("OROV_BeAn19991|Large|ACC123", "OROV_BeAn19991")
    assert id_matches("KR013237.1", "KR013237")


def test_id_matches_rejects_underscore_prefix():
    # The old startswith() would wrongly match these; the boundary rule must not.
    assert not id_matches("OROV_BeAn19991", "OROV")
    assert not id_matches("TRUBV2", "TRUB")
    assert not id_matches("", "TRUB")
    assert not id_matches("TRUB", "")


# --- classify_hit -----------------------------------------------------------

def test_classify_on_off_ignore():
    tids = {"TRUB"}
    assert classify_hit(_hit("TRUB|Small|ACC", 99.0, 25), tids) == "on_target"
    # on-target but below identity/length floor -> ignore
    assert classify_hit(_hit("TRUB|Small|ACC", 85.0, 25), tids) == "ignore"
    # different gene, above off-target floor -> off_target
    assert classify_hit(_hit("HUMAN_GENE|x", 95.0, 24), tids) == "off_target"
    # different gene, below off-target floor -> ignore
    assert classify_hit(_hit("HUMAN_GENE|x", 70.0, 18), tids) == "ignore"


# --- evaluate_specificity ---------------------------------------------------

def _cand(seq="A" * 50):
    # ``sequence`` is load-bearing since the full-half recount slices the probe's
    # two 25 nt halves out of it to compare the positions BLAST never aligned.
    return ProbeCandidate(probe_id="p1", sequence=seq)


def test_passes_with_on_target_no_off_target():
    c = _cand()
    evaluate_specificity(c, [_hit("TRUB|S", 99.0, 25)], {"TRUB"})
    assert c.specificity_pass and c.specificity_status == "passed"


def test_rejected_on_off_target():
    c = _cand()
    evaluate_specificity(
        c, [_hit("TRUB|S", 99.0, 25), _hit("HUMAN|x", 92.0, 24)], {"TRUB"}
    )
    assert not c.specificity_pass and c.specificity_status == "failed_off_target"


def test_off_target_screen_passes_clean_probe():
    # The fix: BLAST is an off-target screen. A probe with no off-target hit
    # PASSES by default even with no on-target hit (its on-target binding is
    # given — it was tiled from the target). The old behaviour wrongly rejected
    # this as "no_on_target_hit", producing zero probes.
    c = _cand()
    evaluate_specificity(c, [_hit("OTHER|x", 70.0, 18)], {"TRUB"})  # sub-threshold
    assert c.specificity_pass and c.specificity_status == "passed"
    c2 = _cand()
    evaluate_specificity(c2, [], {"TRUB"})       # no hits at all
    assert c2.specificity_pass


def test_specificity_margin_discriminates_passing_probes():
    # Audit #4: a passing probe's margin should reflect how close its nearest
    # near-miss off-target came to the disqualification threshold, so it
    # actually discriminates (the old best_on - best_off was ~0 for all).
    from flexprobe.specificity import off_target_reject_pident
    # The margin is expressed against whatever off-target rule is in force. In
    # the default full-half mode the disqualification boundary is 4 mismatches
    # over 25 nt == 84% identity, not the legacy 80% local-alignment floor.
    threshold = off_target_reject_pident()
    assert threshold == 84.0
    # No near-miss off-target at all -> full margin (== threshold).
    clean = _cand()
    evaluate_specificity(clean, [], {"TRUB"})
    assert clean.specificity_margin == threshold
    # A full-length non-target hit just below threshold (79% over 25bp) is a
    # near-miss: small margin, and strictly smaller than the clean probe's.
    close = _cand()
    evaluate_specificity(close, [_hit("OTHER|x", 79.0, 25)], {"TRUB"})
    assert close.specificity_pass
    assert 0.0 <= close.specificity_margin < clean.specificity_margin
    # A short hit (below the length floor) whose UNALIGNED positions are shown
    # to be mismatched is not a real specificity concern and must NOT shrink the
    # margin. Note what supplies that showing: the subject window. A geometry
    # formula that ASSUMED the unaligned positions were mismatches would also
    # get this case right, but the same assumption, applied to a trimmed core,
    # certifies non-specific probes as clean. Here the assumption happens to be
    # true, and it is checked rather than presumed.
    query = "ACGTACGTACGTACGTACGTACGTA"
    subject = query[:15] + "TGCATGCATG"          # 10 unaligned, all mismatched
    short = _cand(seq=query + "A" * 25)
    hit = aligned_hit("OTHER|x", query, subject, 1, 15, half="L", probe="p1")
    assert hit["length"] == 15 and hit["pident"] == 100.0
    evaluate_specificity(short, [hit], {"TRUB"},
                         subjects=FakeSubjects({"OTHER|x": subject}))
    assert short.specificity_margin == threshold
    assert short.offtarget_recount_basis == "exact"

    # ...and with no way to check, it must be treated as disqualifying. The
    # conservative fallback counts the 10 unaligned positions as MATCHES, so the
    # hit reads as 0 mismatches. Both halves are given the same hit here because
    # the 10x rule is a disjunction — one dirty half never rejects a pair.
    # Over-rejecting is the safe direction; the old behaviour's failure mode was
    # the opposite.
    unknown = _cand(seq=query * 2)
    r_hit = aligned_hit("OTHER|x", query, subject, 1, 15, half="R", probe="p1")
    evaluate_specificity(unknown, [hit, r_hit], {"TRUB"})
    assert unknown.specificity_pass is False
    assert unknown.offtarget_recount_basis == "conservative"


def test_strict_mode_requires_on_target():
    c = _cand()
    evaluate_specificity(c, [_hit("OTHER|x", 70.0, 18)], {"TRUB"},
                         require_on_target=True)
    assert not c.specificity_pass and c.specificity_status == "no_on_target_hit"


def test_exact_set_requires_all_strains():
    tids = {"OROV_BeAn19991", "OROV_PEIAM4578"}
    # only one strain hit -> EXACT_SET fails (strict mode only)
    c = _cand()
    evaluate_specificity(c, [_hit("OROV_BeAn19991|L", 99.0, 25)], tids,
                         strain_rule=StrainRule.EXACT_SET, require_on_target=True)
    assert not c.specificity_pass
    # both strains hit -> passes
    c2 = _cand()
    evaluate_specificity(
        c2,
        [_hit("OROV_BeAn19991|L", 99.0, 25), _hit("OROV_PEIAM4578|L", 99.0, 25)],
        tids, strain_rule=StrainRule.EXACT_SET, require_on_target=True,
    )
    assert c2.specificity_pass


def test_batch_evaluate_counts():
    cands = [_cand(), _cand()]
    cands[1].probe_id = "p2"
    hits = {
        "p1": [_hit("TRUB|S", 99.0, 25)],
        "p2": [_hit("OTHER|x", 99.0, 25)],   # off-target only
    }
    n = batch_evaluate(cands, hits, {"TRUB"})
    assert n == 1


# ---------------------------------------------------------------------------
# The full-half mismatch RECOUNT
#
# Every sequence and coordinate below was produced by running blastn-short
# against a DB of public OROV GenBank records (KP052850.1, PQ073184.1), NOT
# constructed to make a point:
#
#   $ blastn -task blastn-short -query <LHS half> -db <OROV_PEIAM4578|Large> \
#           -evalue 1000 -word_size 7 -max_hsps 1
#   ...__L  OROV_PEIAM4578|Large|PQ073184.1  100.000  11  0 0  15 25  3719 3729
#
#   $ blastdbcmd -db ... -entry 'OROV_PEIAM4578|Large|PQ073184.1' \
#           -range 3705-3729 -outfmt %s
#   GGAGGTTTTTGCTTACAGCAGTAGG
#
# The probe half is the LHS 25 nt of the 50 bp window starting at position 3761
# of OROV_BeAn19991|Large|KP052850.1.
# ---------------------------------------------------------------------------

OROV_PROBE_ID = "OROV_BeAn19991|Large|KP052850.1_t003761"
OROV_SUBJECT = "OROV_PEIAM4578|Large|PQ073184.1"
# The 50 bp target window; its two halves are what BLAST is queried with.
OROV_TILE = "GGAGATTTCTACTCACAGCAGTAGGAGACTGTGCCTATCTAGGGCCTTAT"
OROV_LHS = OROV_TILE[:25]           # GGAGATTTCTACTCACAGCAGTAGG
OROV_RHS = OROV_TILE[25:]           # AGACTGTGCCTATCTAGGGCCTTAT
# Subject bases opposite each half, retrieved with blastdbcmd. LHS differs at
# query positions 5, 9, 11, 14 -> 4 mismatches. RHS differs at 1, 16, 19 -> 3.
OROV_SUBJ_LHS = "GGAGGTTTTTGCTTACAGCAGTAGG"
OROV_SUBJ_RHS = "GGACTGTGCCTATCTGGGCCCTTAT"


def _orov_lhs_hit():
    """The REAL trimmed HSP: 11 nt at 100%, covering query 15-25 only."""
    return aligned_hit(OROV_SUBJECT, OROV_LHS, OROV_SUBJ_LHS, 15, 25,
                       subject_offset=3705, half="L", probe=OROV_PROBE_ID)


def _orov_rhs_hit():
    """The REAL RHS HSP: 24 nt at 91.667%, covering query 2-25."""
    return aligned_hit(OROV_SUBJECT, OROV_RHS, OROV_SUBJ_RHS, 2, 25,
                       subject_offset=3730, half="R", probe=OROV_PROBE_ID)


def _orov_subjects():
    # Padded so the two windows sit at their REAL plus-strand coordinates
    # (3705-3729 and 3730-3754). Placing them at offset 1 instead would make the
    # flank lookups run off the end and quietly take the conservative path — the
    # test would still "pass" while testing nothing.
    return FakeSubjects(
        {OROV_SUBJECT: "N" * 3704 + OROV_SUBJ_LHS + OROV_SUBJ_RHS}
    )


def test_off_target_full_half_underestimates_mismatches():
    """A hit 4 mismatches away must be read as 4, whatever BLAST reported.

    blastn-short (-penalty -3) reports the highest-scoring LOCAL alignment. For
    a 25-mer sitting 4 spread mismatches from a subject that is an 11 nt exact
    core (score 11), not the full-length alignment (21 - 12 = 9), so BLAST
    emits the core. Inferring the count from that geometry gives
    ``25 - round(11 * 100/100) = 14``, which is 10 mismatches past the truth and
    certifies the probe specific when it can ligate on the other strain.

    Nothing here is synthetic: the HSP is the one blastn actually emitted for
    OROV probe ...KP052850.1_t003761 against OROV_PEIAM4578|Large, and the
    subject window is what blastdbcmd returned for it.
    """
    from flexprobe.specificity import full_half_mismatches

    hit = _orov_lhs_hit()
    # Guard the guard: if this drifts, the assertion below stops testing the
    # trimmed-HSP case and silently becomes a full-coverage one.
    assert (hit["length"], hit["pident"]) == (11, 100.0)
    assert (hit["qstart"], hit["qend"]) == (15, 25), "HSP must be end-trimmed"

    n = full_half_mismatches(hit, 25, query_seq=OROV_LHS,
                             subjects=_orov_subjects())
    assert n == 4
    assert n != 14, "14 is the old geometry inference — the defect itself"


def test_recount_is_exact_and_fallback_is_conservative():
    """The same trimmed HSP, counted both ways, must differ — and differ DOWN.

    Exact: 4 mismatches (the truth). Conservative: 0, because the 14 query
    positions BLAST never aligned are counted as MATCHES. Both are below the
    5-mismatch floor, so both reject — which is the point. The fallback errs
    toward calling a probe non-specific, the exact inversion of the old bug,
    which could only ever certify a non-specific probe as clean.
    """
    from flexprobe.specificity import (
        RECOUNT_CONSERVATIVE, RECOUNT_EXACT, full_half_mismatch_detail,
    )
    hit = _orov_lhs_hit()

    assert full_half_mismatch_detail(
        hit, 25, query_seq=OROV_LHS, subjects=_orov_subjects(),
    ) == (4, RECOUNT_EXACT)

    # No retrieval handle at all.
    assert full_half_mismatch_detail(hit, 25, query_seq=OROV_LHS) == (
        0, RECOUNT_CONSERVATIVE)
    # Retrieval available but the subject isn't in the DB — same answer, and it
    # must NOT silently fall back to the geometry formula's 14.
    assert full_half_mismatch_detail(
        hit, 25, query_seq=OROV_LHS, subjects=FakeSubjects({}),
    ) == (0, RECOUNT_CONSERVATIVE)

    # The conservative count is a LOWER BOUND: never above the exact one.
    for hit, query, subj in ((_orov_lhs_hit(), OROV_LHS, OROV_SUBJ_LHS),
                             (_orov_rhs_hit(), OROV_RHS, OROV_SUBJ_RHS)):
        exact, _ = full_half_mismatch_detail(
            hit, 25, query_seq=query, subjects=_orov_subjects())
        cons, _ = full_half_mismatch_detail(hit, 25, query_seq=query)
        assert cons <= exact
        assert exact == sum(1 for a, b in zip(query, subj) if a != b)


def test_recount_handles_a_minus_strand_hsp():
    """A minus-strand HSP walks the subject BACKWARDS from the query.

    Real case: the reverse complement of the same OROV LHS half, blasted against
    the same subject, returns qstart=1 qend=11 sstart=3729 send=3719 (sstart >
    send => minus). The truth is still 4 mismatches, so getting the orientation
    wrong — reading the flank forwards, or forgetting to complement it — shows up
    immediately as a different number.
    """
    from flexprobe.specificity import RECOUNT_EXACT, full_half_mismatch_detail

    def _rc(s):
        return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    query = _rc(OROV_LHS)               # CCTACTGCTGTGAGTAGAAATCTCC
    hit = {
        "qseqid": "rc__L", "sseqid": OROV_SUBJECT,
        "pident": 100.0, "length": 11, "mismatch": 0, "gapopen": 0,
        "qstart": 1, "qend": 11, "sstart": 3729, "send": 3719,
        "evalue": 0.022, "bitscore": 22.3,
        "qseq": "CCTACTGCTGT", "sseq": "CCTACTGCTGT",
    }
    # Subject placed so that plus-strand 3705..3729 holds the LHS window.
    subjects = FakeSubjects({OROV_SUBJECT: "N" * 3704 + OROV_SUBJ_LHS})
    assert hit["sstart"] > hit["send"], "fixture must be a minus-strand HSP"
    assert full_half_mismatch_detail(
        hit, 25, query_seq=query, subjects=subjects) == (4, RECOUNT_EXACT)


def test_trimmed_hsp_now_rejects_the_probe_it_used_to_certify():
    """End to end: the OROV probe the old rule passed must now fail.

    Both halves of ...KP052850.1_t003761 sit within 4 mismatches of
    OROV_PEIAM4578|Large (LHS 4, RHS 3), so under the both-halves rule the pair
    can hybridise AND ligate there and must be rejected. The old geometry read
    the LHS as 14 mismatches, leaving only one dirty half, and passed it. This is
    one of the 143/182.
    """
    hits = [_orov_lhs_hit(), _orov_rhs_hit()]

    c = ProbeCandidate(probe_id=OROV_PROBE_ID, sequence=OROV_TILE)
    evaluate_specificity(c, hits, {"OROV_BeAn19991"}, subjects=_orov_subjects())
    assert c.specificity_pass is False
    assert c.specificity_status == "failed_off_target"
    assert c.offtarget_recount_basis == "exact"

    # Watch it fail the other way: with the LHS hit alone only one half is
    # dirty, so the both-halves rule passes the pair. That is what the old
    # geometry effectively did to EVERY probe like this one.
    one_half = ProbeCandidate(probe_id=OROV_PROBE_ID, sequence=OROV_TILE)
    evaluate_specificity(one_half, [_orov_lhs_hit()], {"OROV_BeAn19991"},
                         subjects=_orov_subjects())
    assert one_half.specificity_pass is True


def test_recount_basis_is_recorded_per_probe():
    """A run must never be readable as exactly-recounted when it was bounded."""
    hits = [_orov_lhs_hit(), _orov_rhs_hit()]

    exact = ProbeCandidate(probe_id=OROV_PROBE_ID, sequence=OROV_TILE)
    subs = _orov_subjects()
    evaluate_specificity(exact, hits, {"OROV_BeAn19991"}, subjects=subs)
    assert exact.offtarget_recount_basis == "exact"
    assert (subs.n_exact, subs.n_fallback) == (2, 0)

    bounded = ProbeCandidate(probe_id=OROV_PROBE_ID, sequence=OROV_TILE)
    empty = FakeSubjects({})
    evaluate_specificity(bounded, hits, {"OROV_BeAn19991"}, subjects=empty)
    assert bounded.offtarget_recount_basis == "conservative"
    assert empty.n_fallback == 2

    # Legacy local-alignment mode makes no full-half count at all, and must say
    # so rather than claiming either.
    legacy = ProbeCandidate(probe_id=OROV_PROBE_ID, sequence=OROV_TILE)
    evaluate_specificity(legacy, hits, {"OROV_BeAn19991"},
                         off_target_full_half=False, subjects=_orov_subjects())
    assert legacy.offtarget_recount_basis == "legacy_local"

    # Never screened stays never screened.
    assert ProbeCandidate().offtarget_recount_basis == "not_screened"


# ---------------------------------------------------------------------------
# offtarget_min_qcov — how much of the half BLAST actually looked at
# ---------------------------------------------------------------------------

def test_hit_query_coverage_matches_blast_qcovhsp_definition():
    """Derived, not requested from BLAST — so it must match BLAST's own value.

    Checked against real `qcovhsp` output on 20,817 hits with zero
    disagreements; this pins the formula so a future edit cannot drift from it.
    """
    from flexprobe.specificity import hit_query_coverage
    assert hit_query_coverage({"qstart": 1, "qend": 25}, 25) == 1.0
    assert hit_query_coverage({"qstart": 10, "qend": 24}, 25) == 0.6   # 15/25
    assert hit_query_coverage({"qstart": 1, "qend": 11}, 25) == 0.44   # 11/25
    # No usable range is "unknown", never a fabricated 0.0 that would poison a
    # minimum.
    assert hit_query_coverage({}, 25) is None
    assert hit_query_coverage({"qstart": 1, "qend": 25}, 0) is None


def test_min_qcov_reports_the_trimmed_core_not_the_full_length_hit():
    """The whole point: a verdict resting on a short core must say so.

    A probe with one full-length off-target near-miss and one 11 nt trimmed core
    must report 0.44, not 1.0 — the minimum is what tells you how much of the
    half went unexamined.
    """
    from flexprobe.models import ProbeCandidate
    from flexprobe.specificity import evaluate_specificity
    c = ProbeCandidate(gene_id="g", probe_id="p", sequence="A" * 50,
                       start=0, end=50)
    hits = [
        {"qseqid": "p__L", "sseqid": "OTHER1", "pident": 84.0, "length": 25,
         "qstart": 1, "qend": 25, "sstart": 1, "send": 25,
         "evalue": 1e-3, "bitscore": 30.0},
        {"qseqid": "p__R", "sseqid": "OTHER2", "pident": 100.0, "length": 11,
         "qstart": 1, "qend": 11, "sstart": 1, "send": 11,
         "evalue": 1e-2, "bitscore": 22.0},
    ]
    evaluate_specificity(c, hits, {"TARGET"})
    assert c.offtarget_min_qcov == 0.44


def test_min_qcov_is_nan_when_nothing_off_target_was_seen():
    """No non-target hits at all is the STRONGEST pass, not coverage of zero.

    If this ever reads 0.0 the two become indistinguishable, which is exactly
    the confusion the field exists to remove.
    """
    import math
    from flexprobe.models import ProbeCandidate
    from flexprobe.specificity import evaluate_specificity
    c = ProbeCandidate(gene_id="g", probe_id="p", sequence="A" * 50,
                       start=0, end=50)
    evaluate_specificity(c, [], {"TARGET"})
    assert c.specificity_pass is True
    assert math.isnan(c.offtarget_min_qcov)


# ---------------------------------------------------------------------------
# Junction-proximal mismatches on an OFF-target read backwards: protective
# ---------------------------------------------------------------------------

def _off_hit(half, qstart, qseq, sseq):
    return {"qseqid": f"p__{half}", "sseqid": "OTHER", "pident": 84.0,
            "length": len(qseq), "qstart": qstart, "qend": qstart + len(qseq) - 1,
            "sstart": 1, "send": len(qseq), "evalue": 1e-3, "bitscore": 30.0,
            "qseq": qseq, "sseq": sseq}


def test_junction_zone_is_beside_the_nick_not_the_middle_of_the_half():
    """LHS: the LAST 4 bases. RHS: the FIRST 4. The nick sits between them.

    Getting this wrong would measure the middle of each half — a region where a
    mismatch is thermodynamically ordinary and enzymatically irrelevant.
    """
    from flexprobe.specificity import hit_junction_mismatches
    # LHS half, full 25 nt aligned, single mismatch at query position 25 (the
    # base abutting the nick).
    q = "A" * 25
    s = "A" * 24 + "C"
    assert hit_junction_mismatches(_off_hit("L", 1, q, s), 25, 4) == 1
    # Same mismatch moved to the OUTER end of the LHS half: not protective.
    s2 = "C" + "A" * 24
    assert hit_junction_mismatches(_off_hit("L", 1, q, s2), 25, 4) == 0
    # RHS half: the first base abuts the nick.
    assert hit_junction_mismatches(_off_hit("R", 1, q, s2), 25, 4) == 1
    assert hit_junction_mismatches(_off_hit("R", 1, q, s), 25, 4) == 0


def test_junction_mismatches_are_unknown_not_zero_when_unmeasurable():
    """Unknown must never render as 0 — 0 means "nothing protective"."""
    from flexprobe.specificity import hit_junction_mismatches
    # HSP that stops well short of the nick: the zone was never examined.
    assert hit_junction_mismatches(_off_hit("L", 1, "A" * 10, "A" * 10), 25, 4) is None
    # No half attribution at all.
    h = _off_hit("L", 1, "A" * 25, "A" * 25); h["qseqid"] = "p"
    assert hit_junction_mismatches(h, 25, 4) is None


def test_closest_off_target_profile_is_recorded_on_the_candidate():
    """Two probes, same flat mismatch count, very different real risk."""
    from flexprobe.models import ProbeCandidate
    from flexprobe.specificity import evaluate_specificity
    q = "A" * 25
    junction_mm = "A" * 21 + "CCCC"      # 4 mismatches beside the nick
    end_mm = "CCCC" + "A" * 21           # 4 mismatches at the outer end
    for subject, expect in ((junction_mm, 4), (end_mm, 0)):
        c = ProbeCandidate(gene_id="g", probe_id="p", sequence="A" * 50,
                           start=0, end=50)
        evaluate_specificity(c, [_off_hit("L", 1, q, subject)], {"TARGET"})
        assert c.offtarget_worst_mismatches == 4
        assert c.offtarget_worst_junction_mm == expect


def test_flank_off_the_end_of_the_subject_counts_as_mismatches_not_unknown():
    """Running off a transcript end is KNOWN, not unknown.

    There is no subject sequence past the end of the record, so the probe cannot
    base-pair there — those positions are mismatches. Treating them as unknown
    would send the hit down the conservative path, which counts them as
    MATCHES, makes an off-target look closer than it is, and rejects a probe
    that is fine. It is not an edge case: any hit within a half-length of a
    transcript end has one, and on a human transcriptome screen that is enough
    to reject most of an otherwise passing panel.
    """
    from flexprobe.specificity import full_half_mismatch_detail

    class _Subj:
        """A retrievable record of 30 nt — short enough to overhang."""
        seq = "ACGT" * 7 + "AC"          # 30 nt

        def window_clamped(self, sseqid, start, end):
            lo, hi = max(1, start), min(len(self.seq), end)
            win = self.seq[lo - 1:hi] if hi >= lo else ""
            return win, (end - start + 1) - len(win)

    # Query positions 17-25 aligned to subject 7-15: reconstructing the full
    # 25 nt half needs subject -9..15, so 9 positions lie off the record start.
    hit = {"qseqid": "p__L", "sseqid": "S", "pident": 100.0, "length": 9,
           "qstart": 17, "qend": 25, "sstart": 7, "send": 15,
           "qseq": "ACGTACGTA", "sseq": "ACGTACGTA"}
    n, basis = full_half_mismatch_detail(
        hit, 25, query_seq="A" * 25, subjects=_Subj())
    assert basis == "exact", "an overhang is known, not unknown"
    assert n >= 9, "the 9 unpairable positions must count as mismatches"
