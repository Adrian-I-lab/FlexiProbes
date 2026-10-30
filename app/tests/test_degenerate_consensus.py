"""Degenerate consensus probes: the tolerance rule, the fan-out, the guards.

The parent/variant model lets a pan-strain window carry IUPAC ambiguity where
the input records disagree. Everything here exists because some part of that is
easy to get subtly, silently wrong — the gap-vs-ambiguity conflation above all.
"""
import csv
import os

import pytest

from flexprobe.auto_select import auto_pick_per_gene
from flexprobe.consensus import (
    StrainAlignment, conserved_start_positions, degenerate_conserved_windows,
    resolve_window_per_record,
)
from flexprobe.constants import (
    PROBE_ROLE_CONSENSUS, PROBE_ROLE_PLAIN, PROBE_ROLE_VARIANT,
)
from flexprobe.errors import WriterError
from flexprobe.filters import evaluate_candidate, new_filter_counters
from flexprobe.models import (
    ScreenResources, ExportSettings, FilterThresholds, ModalityDefinition,
    ProbeCandidate, ProbePair, RunConfig,
)
from flexprobe.pipeline import _rollup_degenerate_specificity, run_pipeline
from flexprobe.select_export import _select_rows
from flexprobe.writers.cellranger_probe_set import write_cellranger_probe_set
from flexprobe.writers.idt_opools import write_idt_opools

from conftest import write_fasta


# A 50-mer with no long homopolymer and a mixed GC, used as the "both records
# agree" baseline. Variants are introduced into a copy at chosen offsets.
BASE_50 = "ACGTGACTGACTGGATCCAGTACGTTCAGGATCCAGTACGATCGTTAGCA"
assert len(BASE_50) == 50

_SWAP = {"A": "G", "C": "T", "G": "A", "T": "C"}


def _aln(*seqs, ids=None):
    """A minimal StrainAlignment over pre-aligned, equal-length strings."""
    ids = ids or [f"rec{i}" for i in range(len(seqs))]
    return StrainAlignment(
        fragment_name="frag", strain_ids=list(ids),
        aligned_seqs=list(seqs), consensus=seqs[0],
    )


def _mutate(seq, *offsets):
    """Return ``seq`` with a transition substituted at each offset."""
    chars = list(seq)
    for off in offsets:
        chars[off] = _SWAP[chars[off]]
    return "".join(chars)


# ---------------------------------------------------------------------------
# The tolerance rule
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("offset", [21, 24, 25, 28])
def test_junction_clearance_rejects_variants_inside_the_band(offset):
    """No variant may sit within 4 bp either side of the ligation junction.

    For a 50-mer that is offsets 21..28 inclusive. 24 and 25 straddle the nick
    itself; 21 and 28 are the band edges, which the pre-0.5.0 two-base check
    (``seq[half-1:half+1]``) did not look at.
    """
    aln = _aln(BASE_50, _mutate(BASE_50, offset))
    assert degenerate_conserved_windows(aln, max_degenerate_per_half=1) == []


@pytest.mark.parametrize("offset", [20, 29])
def test_junction_clearance_accepts_variants_just_outside_the_band(offset):
    """20 and 29 are one base beyond the band and must be admitted.

    Paired with the test above this pins both edges: shifting the band by one in
    either direction fails one of the two.
    """
    aln = _aln(BASE_50, _mutate(BASE_50, offset))
    windows = degenerate_conserved_windows(aln, max_degenerate_per_half=1)
    assert [w.start for w in windows] == [0]
    assert windows[0].degenerate_positions == [offset]


def test_junction_band_is_derived_from_tile_length_not_a_fixed_25():
    """A 40-mer's junction is at 20, so its band is 16..23 — not 21..28.

    The codebase has already been bitten once by assuming a 25 nt half (audit
    H1). With a hard-coded band, offset 18 would pass here and offset 25 fail;
    both assertions below invert.
    """
    base = BASE_50[:40]
    inside = _aln(base, _mutate(base, 18))
    outside = _aln(base, _mutate(base, 25))
    assert degenerate_conserved_windows(
        inside, tile_length=40, max_degenerate_per_half=1) == []
    assert len(degenerate_conserved_windows(
        outside, tile_length=40, max_degenerate_per_half=1)) == 1


def test_budget_is_per_half_not_per_window():
    """Two variants in ONE half breach a budget of 1; one in each does not.

    Both windows carry two variant positions, so a per-WINDOW budget of 1 would
    reject both and a per-window budget of 2 would accept both. Only a per-half
    rule separates them — which is the right unit, because the halves are
    separate oligos that hybridise independently.
    """
    same_half = _aln(BASE_50, _mutate(BASE_50, 4, 10))
    split = _aln(BASE_50, _mutate(BASE_50, 4, 40))
    assert degenerate_conserved_windows(same_half, max_degenerate_per_half=1) == []
    assert len(degenerate_conserved_windows(split, max_degenerate_per_half=1)) == 1
    # ... and a budget of 2 admits the two-in-one-half window.
    assert len(degenerate_conserved_windows(
        same_half, max_degenerate_per_half=2)) == 1


def test_indel_window_is_rejected_even_though_its_degenerate_count_is_one():
    """*** THE TRAP. *** A gap column is NOT a one-position variant.

    ``_consensus_from_aligned`` emits 'N' for a gap, and 'N' is in
    DEGENERATE_BASES — so a budget implemented by counting non-ACGT bases in the
    consensus sees "1 degenerate position in the LHS, within budget" and admits
    the window. It must not. An indel means one record genuinely has a base the
    other does not: a probe there does not mispair, it BULGES, and everything
    after it is frame-shifted. The window is not designable as a pan-record probe
    at any tolerance.

    Note what this window looks like to a naive implementation: exactly ONE
    non-ACGT consensus base, at offset 4, nowhere near the junction, zero true
    IUPAC codes. It is indistinguishable from a legitimate one-variant window
    unless you look at the records rather than at the consensus.
    """
    gapped = BASE_50[:4] + "-" + BASE_50[5:]
    aln = _aln(BASE_50, gapped)
    assert degenerate_conserved_windows(aln, max_degenerate_per_half=1) == []
    # And the same window with a real substitution at the same offset IS taken,
    # so the rejection is about the gap and not about offset 4.
    assert len(degenerate_conserved_windows(
        _aln(BASE_50, _mutate(BASE_50, 4)), max_degenerate_per_half=1)) == 1


def test_literal_n_in_an_input_record_is_rejected():
    """An 'N' the user supplied is unknown sequence, not a two-way ambiguity.

    Keeping it out is what lets the fan-out claim that every degenerate position
    has a real, concrete base in every record — and it is what the pre-0.5.0
    ``is_window_conserved`` did too.
    """
    aln = _aln(BASE_50, BASE_50[:4] + "N" + BASE_50[5:])
    assert degenerate_conserved_windows(aln, max_degenerate_per_half=1) == []


def test_tolerance_zero_reproduces_the_legacy_conserved_predicate():
    """Budget 0 must be exact identity — the default, and the compatibility bar.

    Checked against ``conserved_start_positions``, which is a separate, older
    implementation, so this is a genuine cross-check rather than a restatement.
    """
    aln = _aln(
        BASE_50 + BASE_50 + BASE_50,
        BASE_50 + _mutate(BASE_50, 10) + BASE_50,
    )
    starts = [w.start for w in degenerate_conserved_windows(
        aln, max_degenerate_per_half=0)]
    assert starts == conserved_start_positions(aln, tile_length=50)
    assert all(w.degenerate_positions == []
               for w in degenerate_conserved_windows(aln))


def test_zero_variant_window_carries_no_degenerate_annotation():
    """An identical window stays a plain candidate: no parent, no children."""
    aln = _aln(BASE_50, BASE_50)
    (w,) = degenerate_conserved_windows(aln, max_degenerate_per_half=1)
    assert w.degenerate_positions == [] and w.iupac_codes == []
    assert w.sequence == BASE_50


# ---------------------------------------------------------------------------
# Per-record fan-out
# ---------------------------------------------------------------------------

def test_fan_out_is_per_record_not_the_cartesian_product():
    """Two codes in one window give TWO variants, not four.

    The product of the IUPAC codes would invent sequences no input organism
    carries — here the two cross terms (record A's base at one position with
    record B's at the other), which would be synthesised, screened and ordered
    for nothing.
    """
    a = BASE_50
    b = _mutate(BASE_50, 4, 40)
    (w,) = degenerate_conserved_windows(_aln(a, b), max_degenerate_per_half=1)
    resolutions = resolve_window_per_record(
        w.sequence, [a, b], ["A", "B"], 0)
    assert len(resolutions) == 2
    assert {seq for seq, _ids in resolutions} == {a, b}


def test_records_resolving_identically_collapse_to_one_variant():
    """Three records, two of them identical here, give two variants — and the
    surviving variant lists BOTH the records it covers."""
    a, c = BASE_50, _mutate(BASE_50, 4)
    (w,) = degenerate_conserved_windows(
        _aln(a, a, c, ids=["A", "B", "C"]), max_degenerate_per_half=1)
    resolutions = resolve_window_per_record(w.sequence, [a, a, c],
                                            ["A", "B", "C"], 0)
    assert [ids for _seq, ids in resolutions] == [["A", "B"], ["C"]]


# ---------------------------------------------------------------------------
# Filters: the degenerate budget, the junction band, and the Tm hole
# ---------------------------------------------------------------------------

def _thresholds(**kw):
    base = dict(gc_min=0.2, gc_max=0.9, max_n_per_half=1, max_homopolymer_run=6,
                require_tn_junction=False, tm_min=None, tm_max=None)
    base.update(kw)
    return FilterThresholds(**base)


def _candidate(seq, **kw):
    return ProbeCandidate(gene_id="g", probe_id="g_t000000", sequence=seq,
                          start=0, end=len(seq), **kw)


def test_degenerate_budget_counts_iupac_codes_not_just_n():
    """M/R/W/S/Y/K count against ``max_n_per_half``, not only 'N'.

    Two two-way codes in one half is over a budget of 1, although
    ``left.count('N')`` would score that half as 0.
    """
    seq = "R" + "Y" + BASE_50[2:]
    pc = _candidate(seq)
    counters = new_filter_counters()
    assert evaluate_candidate(pc, _thresholds(), counters) is False
    assert "n_count_left" in pc.filter_fail_reasons


def test_junction_band_filter_catches_an_iupac_code_off_the_nick():
    """An 'M' at offset 22 fails: inside the clearance band, off the nick.

    The old check inspected exactly ``seq[24:26]`` and only for a literal 'N',
    so it saw neither the position nor the code.
    """
    seq = BASE_50[:22] + "M" + BASE_50[23:]
    pc = _candidate(seq)
    counters = new_filter_counters()
    assert evaluate_candidate(pc, _thresholds(reject_junction_base_n=True),
                              counters) is False
    assert "junction_base" in pc.filter_fail_reasons


def _parent(seq, **kw):
    """A consensus parent with one degenerate position, outside the band."""
    base = dict(probe_role=PROBE_ROLE_CONSENSUS, n_variants=2,
                degenerate_positions=[4], iupac_codes=["R"])
    base.update(kw)
    return _candidate(seq, **base)


def test_a_degenerate_parent_cannot_bypass_the_tm_gate():
    """*** The most dangerous hole in the feature. ***

    ``safe_tm`` returns NaN for any non-ACGT input (primer3 raises on N and on
    every IUPAC code), and the filter treats NaN as "cannot judge -> pass". Since
    GC was demoted to a ranking term, Tm is THE duplex-stability gate — so a
    parent scored on its own mixed-base string is an UNGATED probe.

    The parent below carries an LHS Tm range whose floor is well under the gate.
    It must fail. Score the parent's own string instead and both halves come back
    NaN and it sails through.
    """
    seq = BASE_50[:4] + "R" + BASE_50[5:]
    cold = _parent(seq, tm_lhs_min=40.0, tm_lhs_max=42.0,
                   tm_rhs_min=65.0, tm_rhs_max=66.0)
    counters = new_filter_counters()
    assert evaluate_candidate(cold, _thresholds(tm_min=55.0), counters) is False
    assert "tm_range_left" in cold.filter_fail_reasons

    # Control: the identical parent with an in-band range passes, so the
    # rejection above is the Tm gate and not the mixed base itself.
    warm = _parent(seq, tm_lhs_min=60.0, tm_lhs_max=61.0,
                   tm_rhs_min=65.0, tm_rhs_max=66.0)
    assert evaluate_candidate(warm, _thresholds(tm_min=55.0),
                              new_filter_counters()) is True


def test_a_parent_without_a_variant_aggregate_fails_closed():
    """No variants, or no Tm aggregate at all, is a fan-out bug — reject.

    This is what makes the gate impossible to skip rather than merely wired up:
    a parent that arrives unaggregated cannot ride through on the NaN
    "cannot judge -> pass" convention the way an ordinary candidate does.
    """
    seq = BASE_50[:4] + "R" + BASE_50[5:]
    counters = new_filter_counters()
    orphan = _parent(seq, n_variants=0)
    assert evaluate_candidate(orphan, _thresholds(tm_min=55.0), counters) is False
    assert "degenerate_unresolved" in orphan.filter_fail_reasons

    unaggregated = _parent(seq)          # n_variants=2 but every Tm still NaN
    assert evaluate_candidate(unaggregated, _thresholds(tm_min=55.0),
                              new_filter_counters()) is False
    assert "degenerate_unresolved" in unaggregated.filter_fail_reasons


def test_parent_gc_is_gated_on_the_variant_range_not_its_own_string():
    """'S' is G+C with certainty yet ``gc_fraction`` counts it 0.

    A parent gated on its own string would use that deflated number. Gated on
    the range, a parent whose variants straddle the GC band fails.
    """
    seq = BASE_50[:4] + "S" + BASE_50[5:]
    pc = _parent(seq, gc_lhs_min=0.30, gc_lhs_max=0.95,
                 gc_rhs_min=0.50, gc_rhs_max=0.55,
                 tm_lhs_min=60.0, tm_lhs_max=61.0,
                 tm_rhs_min=60.0, tm_rhs_max=61.0)
    counters = new_filter_counters()
    assert evaluate_candidate(
        pc, _thresholds(enforce_gc_range=True, gc_min=0.4, gc_max=0.7),
        counters) is False
    assert "gc_range_left" in pc.filter_fail_reasons


# ---------------------------------------------------------------------------
# BLAST rollup: all children must pass
# ---------------------------------------------------------------------------

def test_rollup_fails_the_parent_when_exactly_one_variant_fails():
    """One off-target variant contaminates the pool, so the parent fails —
    and the parent says WHICH variant failed and why, rather than a bare
    "off_target" that sends the reviewer hunting."""
    parent = _candidate(BASE_50, probe_role=PROBE_ROLE_CONSENSUS, n_variants=2)
    good = ProbeCandidate(
        gene_id="g", probe_id="g_t000000_v00", sequence=BASE_50, start=0, end=50,
        probe_role=PROBE_ROLE_VARIANT, parent_probe_id="g_t000000",
        specificity_pass=True, specificity_status="pass", specificity_margin=9.0)
    bad = ProbeCandidate(
        gene_id="g", probe_id="g_t000000_v01", sequence=BASE_50, start=0, end=50,
        probe_role=PROBE_ROLE_VARIANT, parent_probe_id="g_t000000",
        specificity_pass=False, specificity_status="off_target",
        specificity_margin=1.0)

    assert _rollup_degenerate_specificity([parent, good, bad]) == 0
    assert parent.specificity_pass is False
    assert "g_t000000_v01" in parent.specificity_status
    assert "off_target" in parent.specificity_status
    # The pool is only as specific as its least specific member.
    assert parent.specificity_margin == 1.0

    # Control: flip the one failure and the same parent passes, so the test
    # above is measuring the rollup and not a parent that can never pass.
    bad.specificity_pass = True
    bad.specificity_status = "pass"
    assert _rollup_degenerate_specificity([parent, good, bad]) == 1
    assert parent.specificity_pass is True


def test_rollup_fails_a_parent_whose_variants_all_vanished():
    parent = _candidate(BASE_50, probe_role=PROBE_ROLE_CONSENSUS, n_variants=2,
                        specificity_pass=True)
    assert _rollup_degenerate_specificity([parent]) == 0
    assert parent.specificity_pass is False


# ---------------------------------------------------------------------------
# Parents are never orderable
# ---------------------------------------------------------------------------

def _pair(role):
    lhs = _candidate(BASE_50, probe_role=role)
    lhs.probe_rc = BASE_50
    rhs = _candidate(BASE_50, probe_role=role)
    return ProbePair(lhs=lhs, rhs=rhs)


def test_a_parent_never_reaches_the_idt_writer(tmp_path):
    """An IDT oPool row is an order: real money, for an oligo that binds
    nothing in particular. The writer refuses rather than trusting upstream."""
    out = str(tmp_path / "idt.csv")
    with pytest.raises(WriterError, match="never orderable"):
        write_idt_opools([_pair(PROBE_ROLE_CONSENSUS)], out, overwrite=True)
    assert not os.path.exists(out)
    # A plain pair writes fine, so the guard is not simply refusing everything.
    write_idt_opools([_pair(PROBE_ROLE_PLAIN)], out, overwrite=True)
    assert os.path.exists(out)


def test_a_parent_never_reaches_the_cellranger_writer(tmp_path):
    """Cell Ranger reads probe_seq as literal nucleotides; a mixed-base string
    there is an invalid record, not a permissive probe."""
    out = str(tmp_path / "cr.csv")
    with pytest.raises(WriterError, match="never orderable"):
        write_cellranger_probe_set([_pair(PROBE_ROLE_CONSENSUS)], out,
                                   overwrite=True)
    write_cellranger_probe_set([_pair(PROBE_ROLE_PLAIN)], out, overwrite=True)
    assert os.path.exists(out)


def test_cellranger_writer_also_rejects_a_mixed_base_sequence(tmp_path):
    """Belt and braces: caught on the characters as well as on the role, since
    the Cell Ranger probe_seq carries no handle that could excuse one."""
    pair = _pair(PROBE_ROLE_PLAIN)
    pair.lhs.probe_rc = BASE_50[:4] + "R" + BASE_50[5:]
    with pytest.raises(WriterError):
        write_cellranger_probe_set([pair], str(tmp_path / "cr.csv"),
                                   overwrite=True)


def test_selection_of_a_parent_exports_its_variants_instead():
    """Picking the parent is the natural review gesture — it is what names the
    SITE. It must turn into the concrete oligos, all of them, and never into a
    mixed-base row."""
    rows = [
        {"gene_id": "g", "probe_id": "p", "probe_role": PROBE_ROLE_CONSENSUS,
         "parent_probe_id": ""},
        {"gene_id": "g", "probe_id": "p_v00", "probe_role": PROBE_ROLE_VARIANT,
         "parent_probe_id": "p"},
        {"gene_id": "g", "probe_id": "p_v01", "probe_role": PROBE_ROLE_VARIANT,
         "parent_probe_id": "p"},
    ]
    selected = _select_rows(rows, {"g": ["p"]})
    assert [r["probe_id"] for r in selected] == ["p_v00", "p_v01"]
    assert not any(r["probe_role"] == PROBE_ROLE_CONSENSUS for r in selected)


def test_auto_select_picks_a_family_as_one_orderable_probe():
    """One family = one slot, and that slot yields the ORDERABLE ids.

    Worth being precise about what the pre-fix defect was, because the obvious
    guess is wrong. The old code did NOT double-count the family against the
    per-gene cap: a parent and its variants share the same [start, end), so the
    ``min_gap`` spacing check already rejected the children as overlapping the
    parent, and the family consumed exactly one slot either way. Measured on
    these rows the old code returns ``{'g': ['p']}`` at a cap of 1 and
    ``{'g': ['p', 'q']}`` at 2.

    The real defect is that the one slot was spent on the PARENT — which is not
    orderable — so the site contributed nothing at all to the order sheet while
    still consuming the budget of a site that would have. The fix makes the
    parent the cap unit (it carries the score and the spacing) but emits its
    variants' ids, which are what you actually order.
    """
    rows = [
        {"gene_id": "g", "probe_id": "p", "probe_role": PROBE_ROLE_CONSENSUS,
         "parent_probe_id": "", "start": "0", "end": "50",
         "composite_score": "0.9"},
        {"gene_id": "g", "probe_id": "p_v00", "probe_role": PROBE_ROLE_VARIANT,
         "parent_probe_id": "p", "start": "0", "end": "50",
         "composite_score": "0.9"},
        {"gene_id": "g", "probe_id": "p_v01", "probe_role": PROBE_ROLE_VARIANT,
         "parent_probe_id": "p", "start": "0", "end": "50",
         "composite_score": "0.9"},
        {"gene_id": "g", "probe_id": "q", "probe_role": PROBE_ROLE_PLAIN,
         "parent_probe_id": "", "start": "500", "end": "550",
         "composite_score": "0.8"},
    ]
    # Cap of 1: the family is picked, and it yields both variants (which are
    # ordered together) while consuming exactly one slot — so `q` is not picked.
    picked = auto_pick_per_gene(rows, n_per_gene=1, min_gap=10)
    assert picked == {"g": ["p_v00", "p_v01"]}
    # Cap of 2: the family plus the plain probe.
    picked2 = auto_pick_per_gene(rows, n_per_gene=2, min_gap=10)
    assert picked2 == {"g": ["p_v00", "p_v01", "q"]}


# ---------------------------------------------------------------------------
# End-to-end
# ---------------------------------------------------------------------------

# Two 300 bp "strains" differing at three positions, each clear of a junction
# for at least one tiling window.
_STRAIN_A = (BASE_50 + "TTGACCAGTGCATCAGGTACCATGGCATCGATCAGTTGCAAGCTTGACCA"
             + BASE_50 + "GGATCCAGTTACGCATCAGGATCCAGTACGATCGGTACCATGCATGCATA"
             + BASE_50 + "CATGGCATCGATCAGTTGCAAGCTTGACCAGTGCATCAGGTACCATGGCA")
_STRAIN_B = _mutate(_STRAIN_A, 4, 140, 260)


def _e2e_cfg(fa, out, **kw):
    base = dict(
        fasta_path=fa, output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        filters=FilterThresholds(gc_min=0.20, gc_max=0.90, max_n_per_half=1,
                                 max_homopolymer_run=6,
                                 require_tn_junction=False),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        multi_strain=True, consensus_strategy="conserved",
        discovery_only=True, tile_length=50, tile_step=1,
    )
    base.update(kw)
    return RunConfig(**base)


def _rows(out):
    with open(os.path.join(out, "results.csv")) as fh:
        return list(csv.DictReader(fh))


def test_tolerance_zero_is_unchanged_end_to_end(tmp_path):
    """*** THE COMPATIBILITY BAR. *** Explicit 0 must reproduce a prior run.

    The DEFAULT is now 1, so "unchanged" is a promise about the explicit
    setting, not about leaving it alone. A run pinned at 0 must still produce
    exact-identity conserved windows only: every new column present, every one
    blank, no parents and no variants.
    """
    fa = write_fasta(tmp_path / "s.fa", [("A", _STRAIN_A), ("B", _STRAIN_B)])
    out = str(tmp_path / "zero")
    run_pipeline(_e2e_cfg(fa, out, conserved_max_degenerate_per_half=0))
    rows = _rows(out)
    assert rows and all(r["probe_role"] == "" for r in rows)
    assert all(r["parent_probe_id"] == "" for r in rows)
    assert all(not (set(r["mrna_target_50bp"]) - set("ACGT")) for r in rows)


def test_the_default_tolerance_is_one_not_zero(tmp_path):
    """0 is not a neutral default — it reads as biology.

    On real strains, exact identity can leave a whole segment with no
    designable window, which then reads as "this segment has no pan-strain
    probe". One variant position per half recovers usable windows.
    """
    from flexprobe.models import RunConfig
    assert RunConfig().conserved_max_degenerate_per_half == 1

def test_degenerate_run_emits_parents_and_variants(tmp_path):
    """The fan-out end to end: a parent plus one variant per record, ids derived
    from the parent's, and no parent among the emitted pairs."""
    fa = write_fasta(tmp_path / "s.fa", [("A", _STRAIN_A), ("B", _STRAIN_B)])
    out = str(tmp_path / "out")
    result = run_pipeline(
        _e2e_cfg(fa, out, conserved_max_degenerate_per_half=1))
    rows = _rows(out)
    parents = [r for r in rows if r["probe_role"] == PROBE_ROLE_CONSENSUS]
    variants = [r for r in rows if r["probe_role"] == PROBE_ROLE_VARIANT]
    assert parents, "tolerance 1 should admit some degenerate windows"
    assert len(variants) == 2 * len(parents)

    for p in parents:
        kids = [r for r in variants if r["parent_probe_id"] == p["probe_id"]]
        assert len(kids) == int(p["n_variants"]) == 2
        # Ids derive from the parent's, so they are stable across runs — the
        # pipeline is resumable and its caches are keyed on them.
        assert {k["probe_id"] for k in kids} == {
            f"{p['probe_id']}_v00", f"{p['probe_id']}_v01"}
        # Every variant is concrete ACGT; the parent is not.
        for k in kids:
            assert set(k["mrna_target_50bp"]) <= set("ACGT")
            assert k["resolves_to"] in ("A", "B")
        assert not set(p["mrna_target_50bp"]) <= set("ACGT")
        # The parent carries its variants' aggregate, not blanks.
        assert float(p["tm_lhs_min"]) <= float(p["tm_lhs_max"])
        assert p["reduced_confidence"] == "True"
        assert p["degenerate_positions"] and p["iupac_codes"]

    # No parent became a ProbePair, so none can reach an order writer.
    assert all(pair.lhs.probe_role != PROBE_ROLE_CONSENSUS
               for pair in result["pairs"])
    parent_ids = {p["probe_id"] for p in parents}
    assert not (parent_ids & {pair.lhs.probe_id for pair in result["pairs"]})


def test_degenerate_run_never_places_ambiguity_near_the_junction(tmp_path):
    """Every degenerate position emitted by a real run clears the junction."""
    fa = write_fasta(tmp_path / "s.fa", [("A", _STRAIN_A), ("B", _STRAIN_B)])
    out = str(tmp_path / "out")
    run_pipeline(_e2e_cfg(fa, out, conserved_max_degenerate_per_half=1))
    for r in _rows(out):
        if r["probe_role"] != PROBE_ROLE_CONSENSUS:
            continue
        offsets = [int(x) for x in r["degenerate_positions"].split(",") if x]
        assert offsets
        assert not any(21 <= o <= 28 for o in offsets)
        assert len([o for o in offsets if o < 25]) <= 1
        assert len([o for o in offsets if o >= 25]) <= 1
