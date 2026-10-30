"""Junction clearance for CONCRETE mismatches, and the mode-scoping that makes
it safe.

``filters`` has always rejected an IUPAC ambiguity code sitting in the
ligation-junction clearance band. That catches only what is written into the
probe STRING. A probe that is pure ACGT can still mispair a strain at offset 22
— the ligase will not seal across that — and until now nothing gated it; the
mismatch earned a dTm penalty in a report that gates nothing.

CG000621 Rev D asks for both in one sentence:

  "Avoid common single nucleotide polymorphisms (SNPs) and potential mismatches
   at the ligation junction. Refer to the UCSC Genome Browser and the Single
   Nucleotide Polymorphism Database (dbSNP). If avoiding SNPs is not possible,
   SNPs and mismatches should be at least four bp away from the ligation
   junction."

The dangerous half of this feature is not the rule, it is the SCOPE. Against a
record the probe must bind, a junction mismatch is fatal. Against a record the
probe is designed to DISCRIMINATE from, the same mismatch is the deliverable.
So every mode gets its own test, and the strain-unique one is the one that
matters: if it ever starts rejecting, discriminatory design is silently dead.
"""
import csv
import json
import os

import pytest

from flexprobe.consensus import StrainAlignment, score_probe_against_all_strains
from flexprobe.filters import (
    apply_filters, candidate_on_target_records, evaluate_candidate,
    new_filter_counters, on_target_record_ids,
)
from flexprobe.models import (
    ExportSettings, FilterThresholds, ModalityDefinition, ProbeCandidate,
    RunConfig, ScreenResources,
)
from flexprobe.pipeline import _stream_filter, run_pipeline

from conftest import write_fasta


# Mixed GC, no run longer than two, and deliberately NOT an 'A' at offset 25 —
# the TN gate is switched off below so that only the junction rules can speak.
BASE_50 = "ACGTGACTGACTGGATCCAGTACGTTCAGGATCCAGTACGATCGTTAGCA"
assert len(BASE_50) == 50

# Transitions, so a mutated record stays plain ACGT.
_SWAP = {"A": "G", "C": "T", "G": "A", "T": "C"}

# The band under test: half=25, clearance=4 -> offsets 21..28 inclusive.
IN_BAND = 22
OUT_OF_BAND = 10


def _mutate(seq, *offsets):
    chars = list(seq)
    for off in offsets:
        chars[off] = _SWAP[chars[off]]
    return "".join(chars)


def _permissive():
    """Every gate except the junction rules turned off or opened wide.

    A junction test that fails because of GC or a homopolymer proves nothing,
    so the only reasons this funnel can produce are 'junction_base' and
    'junction_mismatch'.
    """
    return FilterThresholds(
        gc_min=0.0, gc_max=1.0, enforce_gc_range=False,
        tm_min=None, tm_max=None,
        max_n_per_half=50, max_homopolymer_run=99,
        require_tn_junction=False,
    )


def _strain_binding(probe_seq, records):
    """The per-record mismatch annotation, built exactly as the pipeline builds
    it (``pipeline._annotate_strain_binding``) from the same scorer, so these
    tests exercise the real dict shape and not a hand-written stand-in.

    ``records`` is ``[(record_id, aligned_sequence)]``.
    """
    aln = StrainAlignment(
        fragment_name="frag",
        strain_ids=[sid for sid, _s in records],
        aligned_seqs=[s for _sid, s in records],
        consensus=records[0][1],
    )
    scores = score_probe_against_all_strains(
        aln, 0, tile_length=len(probe_seq), probe_seq=probe_seq)
    return [
        {
            "strain_id": s.strain_id,
            "total_mismatches": s.total_mismatches,
            "lhs_mismatches": s.lhs_mismatches,
            "rhs_mismatches": s.rhs_mismatches,
            "junction_mismatches": s.junction_mismatches,
            "mismatch_positions": s.mismatch_positions,
            "estimated_dTm": s.estimated_dTm,
            "likely_binds": s.likely_binds,
        }
        for s in scores
    ]


def _candidate(probe_seq, records=None, covered=None):
    pc = ProbeCandidate(
        probe_id="t", sequence=probe_seq, start=0, end=len(probe_seq))
    if records:
        pc.strain_binding = _strain_binding(probe_seq, records)
    if covered:
        pc.covered_record_ids = list(covered)
    return pc


def _run(pc, scope, thresholds=None):
    counters = new_filter_counters()
    thresholds = thresholds or _permissive()
    evaluate_candidate(pc, thresholds, counters, on_target_records=scope)
    return counters


def _scope(strategy="conserved", ids=("A", "B"), **kw):
    return on_target_record_ids(
        RunConfig(multi_strain=True, consensus_strategy=strategy, **kw),
        list(ids))


# ---------------------------------------------------------------------------
# Mode 1 — conserved: a mismatch to ANY record is disqualifying
# ---------------------------------------------------------------------------

def test_conserved_rejects_a_concrete_junction_mismatch():
    """The gap this closes. The probe is plain ACGT — no ambiguity code for the
    existing rule to see — and it mispairs the second record at offset 22."""
    records = [("A", BASE_50), ("B", _mutate(BASE_50, IN_BAND))]
    pc = _candidate(BASE_50, records)

    counters = _run(pc, _scope("conserved"))

    assert pc.passed_filters is False
    assert "junction_mismatch" in pc.filter_fail_reasons
    assert counters["junction_mismatch"] == 1
    # The existing rule cannot see this: the probe carries no degenerate base.
    assert "junction_base" not in pc.filter_fail_reasons
    assert counters["junction_base"] == 0


def test_conserved_allows_a_mismatch_outside_the_band():
    """The rule is a band, not a blanket ban on any strain disagreement."""
    records = [("A", BASE_50), ("B", _mutate(BASE_50, OUT_OF_BAND))]
    pc = _candidate(BASE_50, records)

    counters = _run(pc, _scope("conserved"))

    assert pc.passed_filters is True, pc.filter_fail_reasons
    assert counters["junction_mismatch"] == 0


@pytest.mark.parametrize("offset,rejected", [
    (20, False),   # just outside the band
    (21, True),    # first excluded base on the LHS
    (24, True),    # LHS base adjacent to the nick
    (25, True),    # RHS base adjacent to the nick
    (28, True),    # last excluded base on the RHS
    (29, False),   # just outside the band
])
def test_the_band_is_offsets_21_to_28_inclusive(offset, rejected):
    """Pins the exact band, INCLUDING the off-by-one that CG000621 does not
    resolve (see ``filters._junction_band``). 10x says "at least four bp away";
    whether the base adjacent to the nick is 0 or 1 bp away decides whether
    this band is the rule or one base stricter than it. flexiprobes takes the
    stricter reading deliberately, so this test exists to make any future
    narrowing of the band a visible, argued change rather than a quiet one."""
    records = [("A", BASE_50), ("B", _mutate(BASE_50, offset))]
    pc = _candidate(BASE_50, records)

    _run(pc, _scope("conserved"))

    assert ("junction_mismatch" in pc.filter_fail_reasons) is rejected


# ---------------------------------------------------------------------------
# Mode 2 — strain-unique: the discrimination must survive
# ---------------------------------------------------------------------------

def test_strain_unique_keeps_a_discriminating_junction_mismatch():
    """*** THE TEST THAT MATTERS. ***

    In strain-unique mode the probe IS the target strain's literal sequence and
    its mismatches against the other strains are the entire product. A rule
    scoped to "any record" would reject exactly the windows the run exists to
    find — and would do it silently, since the probe looks fine and the counter
    would just read as a normal rejection.
    """
    records = [("A", BASE_50), ("B", _mutate(BASE_50, IN_BAND))]
    pc = _candidate(BASE_50, records)          # probe == record A's sequence

    counters = _run(pc, _scope("strain_unique", strain_unique_target=0))

    # Guard against a vacuous pass: the mismatch really is there, in the band,
    # against record B — it is simply not B's job to be bound.
    by_id = {d["strain_id"]: d for d in pc.strain_binding}
    assert IN_BAND in by_id["B"]["mismatch_positions"]
    assert by_id["A"]["mismatch_positions"] == []

    assert pc.passed_filters is True, pc.filter_fail_reasons
    assert counters["junction_mismatch"] == 0


def test_strain_unique_still_rejects_a_mismatch_against_its_own_target():
    """Scoping down is not the same as switching off. A probe that mispairs the
    strain it is supposed to be diagnostic FOR is still broken."""
    records = [("A", BASE_50), ("B", _mutate(BASE_50, IN_BAND))]
    pc = _candidate(records[1][1], records)    # probe == B's sequence...

    counters = _run(pc, _scope("strain_unique", strain_unique_target=0))
    #                                            ...but A is the target

    assert "junction_mismatch" in pc.filter_fail_reasons
    assert counters["junction_mismatch"] == 1


# ---------------------------------------------------------------------------
# Mode 3 — single sequence: inert
# ---------------------------------------------------------------------------

def test_single_sequence_run_is_untouched():
    """No multi-strain context, no scope, no rule. A single-sequence candidate
    also carries no ``strain_binding`` at all, so there are two independent
    reasons this can never fire — both asserted, because the pipeline's
    single-sequence path relies on each of them separately."""
    cfg = RunConfig(multi_strain=False, consensus_strategy="conserved")
    assert on_target_record_ids(cfg, ["A"]) is None

    pc = _candidate(BASE_50)
    assert pc.strain_binding == []
    counters = _run(pc, on_target_record_ids(cfg, ["A"]))

    assert pc.passed_filters is True, pc.filter_fail_reasons
    assert counters["junction_mismatch"] == 0


def test_strain_binding_without_a_scope_stays_inert():
    """Fail-open, deliberately. A caller that hands over candidates but no scope
    (a downstream tool, a bare ``apply_filters``) must not have "every record"
    assumed for it — that assumption is precisely what breaks strain-unique."""
    records = [("A", BASE_50), ("B", _mutate(BASE_50, IN_BAND))]
    pc = _candidate(BASE_50, records)

    counters = apply_filters([pc], _permissive())

    assert pc.passed_filters is True, pc.filter_fail_reasons
    assert counters["junction_mismatch"] == 0


# ---------------------------------------------------------------------------
# The existing degenerate-base rule is unchanged
# ---------------------------------------------------------------------------

def test_degenerate_base_rule_is_unchanged():
    """Same band, same reason string, same counter — and it still fires with no
    strain context whatsoever, which is what every pre-existing run relied on."""
    seq = BASE_50[:IN_BAND] + "R" + BASE_50[IN_BAND + 1:]
    pc = _candidate(seq)

    counters = _run(pc, None)

    assert pc.filter_fail_reasons == ["junction_base"]
    assert counters["junction_base"] == 1
    assert counters["junction_mismatch"] == 0


def test_a_degenerate_parent_is_not_double_counted():
    """An IUPAC code reads as "any of these bases" when scoring a record
    (consensus._binding_score_for_window), so a consensus parent's own ambiguity
    is NOT a mismatch. The parent is rejected once, by the rule that owns it."""
    assert BASE_50[IN_BAND] == "C"
    records = [("A", BASE_50), ("B", _mutate(BASE_50, IN_BAND))]
    assert records[1][1][IN_BAND] == "T"
    parent = BASE_50[:IN_BAND] + "Y" + BASE_50[IN_BAND + 1:]   # Y = C or T
    pc = _candidate(parent, records)

    counters = _run(pc, _scope("conserved"))

    assert pc.filter_fail_reasons == ["junction_base"]
    assert counters["junction_mismatch"] == 0


def test_the_two_rules_count_separately_when_both_apply():
    """A code that does NOT cover the record's base is both things at once: an
    ambiguity in the band and a real mispair against a record we must bind. The
    report has to be able to say so, which is why 'junction_mismatch' is not
    folded into 'junction_base'."""
    probe = BASE_50[:IN_BAND] + "M" + BASE_50[IN_BAND + 1:]    # M = A or C
    other = BASE_50[:IN_BAND] + "G" + BASE_50[IN_BAND + 1:]    # G is not in M
    pc = _candidate(probe, [("A", BASE_50), ("B", other)])

    counters = _run(pc, _scope("conserved"))

    assert pc.filter_fail_reasons == ["junction_base", "junction_mismatch"]
    assert counters["junction_base"] == 1
    assert counters["junction_mismatch"] == 1


def test_the_switch_turns_only_its_own_rule_off():
    records = [("A", BASE_50), ("B", _mutate(BASE_50, IN_BAND))]
    thr = _permissive()
    thr.reject_junction_mismatch = False
    pc = _candidate(BASE_50, records)

    counters = _run(pc, _scope("conserved"), thr)

    assert pc.passed_filters is True, pc.filter_fail_reasons
    assert counters["junction_mismatch"] == 0

    # ...and the degenerate rule is still live under the same thresholds.
    degen = _candidate(BASE_50[:IN_BAND] + "R" + BASE_50[IN_BAND + 1:])
    assert _run(degen, _scope("conserved"), thr)["junction_base"] == 1
    assert degen.filter_fail_reasons == ["junction_base"]


def test_it_defaults_on():
    assert FilterThresholds().reject_junction_mismatch is True


# ---------------------------------------------------------------------------
# The mode -> scope mapping itself
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cfg_kw,expected", [
    # Single sequence: no scope at all.
    (dict(multi_strain=False), None),
    # Pan-strain: the probe claims to bind every record, so every record counts.
    (dict(multi_strain=True, consensus_strategy="conserved"), {"A", "B", "C"}),
    # Discriminatory: ONLY the target. This is the inversion to guard.
    (dict(multi_strain=True, consensus_strategy="strain_unique",
          strain_unique_target=1), {"B"}),
    (dict(multi_strain=True, consensus_strategy="strain_unique",
          strain_unique_target=2), {"C"}),
    # Out of range is a config error the pipeline raises on separately; it must
    # not silently degrade to "every record".
    (dict(multi_strain=True, consensus_strategy="strain_unique",
          strain_unique_target=9), None),
    # Constitutive at the default quorum == "shared by all", so all records.
    (dict(multi_strain=True, consensus_strategy="constitutive",
          isoform_min_count=0), {"A", "B", "C"}),
    (dict(multi_strain=True, consensus_strategy="constitutive",
          isoform_min_count=3), {"A", "B", "C"}),
    # Below that quorum the window is designed for ONE isoform group and the
    # records outside it are not targets. WHICH group differs per window, so
    # there is no run-level answer and this still returns None — but that no
    # longer means the rule is skipped: the scope now arrives per candidate on
    # ProbeCandidate.covered_record_ids instead. See the
    # candidate_on_target_records tests below for the half this one cannot see.
    (dict(multi_strain=True, consensus_strategy="constitutive",
          isoform_min_count=2), None),
])
def test_on_target_record_ids_per_mode(cfg_kw, expected):
    got = on_target_record_ids(RunConfig(**cfg_kw), ["A", "B", "C"])
    assert (None if got is None else set(got)) == expected


def test_on_target_record_ids_needs_records():
    cfg = RunConfig(multi_strain=True, consensus_strategy="conserved")
    assert on_target_record_ids(cfg, []) is None


# ---------------------------------------------------------------------------
# Mode 3 — constitutive below its quorum: the scope is the WINDOW's isoform
# group, and it can only arrive per candidate
# ---------------------------------------------------------------------------
#
# Each constitutive window is the shared sequence of its own largest isoform
# group, so "which records must this bind" is a property of the WINDOW, not of
# the run. on_target_record_ids has no answer to give and returns None; before
# the group was carried onto the candidate that None meant the rule was skipped
# outright on this mode. These pin the half that closes it.

def test_a_candidate_group_overrides_the_absent_run_scope():
    """The gap this closes. Run scope is None (partial quorum), so without the
    candidate's own group the rule could not run at all."""
    assert candidate_on_target_records(
        _candidate(BASE_50, covered=["A", "B"]), None) == frozenset({"A", "B"})


def test_a_candidate_without_a_group_keeps_the_run_scope():
    """Every other mode. The field is empty, so the run-level scope decides and
    behaviour is exactly what it was before per-candidate scoping existed."""
    run = frozenset({"A", "B", "C"})
    assert candidate_on_target_records(_candidate(BASE_50), run) == run
    assert candidate_on_target_records(_candidate(BASE_50), None) is None


def test_a_candidate_group_narrows_a_wider_run_scope():
    """Candidate-first precedence, stated explicitly. Nothing populates the
    field alongside a wider run scope today (only constitutive sets it, and at
    the default quorum the two sets are equal), but the precedence is the whole
    contract of the field, so it is pinned rather than left to inference."""
    assert candidate_on_target_records(
        _candidate(BASE_50, covered=["A"]),
        frozenset({"A", "B", "C"})) == frozenset({"A"})


def test_partial_quorum_keeps_a_window_its_group_does_not_cover():
    """*** THE INVERSION GUARD FOR THIS MODE. *** C is not in the isoform group,
    so the probe is never claimed to bind it and its junction mismatch is not a
    defect — it is the isoform this window was allowed to miss. Rejecting here
    would delete exactly the windows a partial quorum exists to find."""
    records = [("A", BASE_50), ("B", BASE_50),
               ("C", _mutate(BASE_50, IN_BAND))]
    counters = _run(_candidate(BASE_50, records, covered=["A", "B"]), None)
    assert counters['junction_mismatch'] == 0
    assert counters['passed'] == 1


def test_partial_quorum_rejects_a_mismatch_against_a_COVERED_isoform():
    """The same candidate and the same mismatch, with C inside the group. Now
    the probe does claim to bind it, so the junction-proximal mispair is fatal
    and the gate is armed rather than merely present."""
    records = [("A", BASE_50), ("B", BASE_50),
               ("C", _mutate(BASE_50, IN_BAND))]
    counters = _run(_candidate(BASE_50, records, covered=["A", "C"]), None)
    assert counters['junction_mismatch'] == 1
    assert counters['passed'] == 0


def test_a_covered_group_still_respects_the_band():
    """Scoping changed WHICH records are consulted, not where the band is: a
    mismatch against a covered isoform outside the clearance band still passes.
    """
    records = [("A", BASE_50), ("B", _mutate(BASE_50, OUT_OF_BAND))]
    counters = _run(_candidate(BASE_50, records, covered=["A", "B"]), None)
    assert counters['junction_mismatch'] == 0
    assert counters['passed'] == 1


def test_conserved_variants_are_untouched_by_per_candidate_scoping():
    """Why this rides on its own field instead of ``resolves_to``.

    A conserved-mode degenerate variant already carries ``resolves_to`` — the
    records whose sequence it matches — and it deliberately MISMATCHES its
    siblings' records. Had the per-candidate scope read that column, every such
    variant would have silently narrowed from "must bind every record" to "must
    bind its own", turning a mismatch the conserved gate is supposed to catch
    into one it ignores. ``covered_record_ids`` is empty on those rows, so the
    run-level scope still decides and the rejection below still happens.
    """
    records = [("A", BASE_50), ("B", _mutate(BASE_50, IN_BAND))]
    pc = _candidate(BASE_50, records)
    pc.probe_role = "variant"
    pc.resolves_to = ["A"]            # populated, and must NOT scope anything
    assert candidate_on_target_records(pc, _scope("conserved")) == frozenset(
        {"A", "B"})
    counters = _run(pc, _scope("conserved"))
    assert counters['junction_mismatch'] == 1


# ---------------------------------------------------------------------------
# Plumbing: the scope has to survive the streaming filter and the pipeline
# ---------------------------------------------------------------------------

def test_stream_filter_forwards_the_scope():
    records = [("A", BASE_50), ("B", _mutate(BASE_50, IN_BAND))]

    passing, counters, _rej = _stream_filter(
        [_candidate(BASE_50, records)], _permissive(),
        on_target_records=_scope("conserved"))
    assert passing == [] and counters["junction_mismatch"] == 1

    passing, counters, _rej = _stream_filter(
        [_candidate(BASE_50, records)], _permissive())
    assert len(passing) == 1 and counters["junction_mismatch"] == 0


_STRAIN_A = (BASE_50 + "TTGACCATGGACTTCAGGATCCAGTACGATCGTTAGCAACGTGACTGACT"
             + BASE_50)
_STRAIN_B = _mutate(_STRAIN_A, 60, 130)


def test_a_conserved_run_still_designs_probes_and_reports_the_counter(tmp_path):
    """Non-regression, end to end. On every path reachable today the new rule is
    a backstop that does not fire — conserved windows are gap-free and keep
    their variant positions out of this band already
    (``consensus.degenerate_conserved_windows`` rules 1 and 3, applied with the
    SAME clearance). So a real multi-strain run must be unchanged, and the
    counter must still be present and zero rather than absent."""
    fa = write_fasta(tmp_path / "s.fa", [("A", _STRAIN_A), ("B", _STRAIN_B)])
    out = str(tmp_path / "out")
    run_pipeline(RunConfig(
        fasta_path=fa, output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        filters=FilterThresholds(gc_min=0.20, gc_max=0.90, max_n_per_half=1,
                                 max_homopolymer_run=6,
                                 require_tn_junction=False),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        multi_strain=True, consensus_strategy="conserved",
        discovery_only=True, tile_length=50, tile_step=1,
    ))

    with open(os.path.join(out, "results.csv")) as fh:
        rows = list(csv.DictReader(fh))
    assert rows

    with open(os.path.join(out, "filter_report.json")) as fh:
        report = json.load(fh)
    assert report["counts"]["junction_mismatch"] == 0
    # The threshold travels with the report, so a run records the funnel that
    # produced its own counts.
    assert report["thresholds"]["reject_junction_mismatch"] is True


# A third isoform that differs from A/B at two offsets, so the largest group at
# those windows is {A, B} — a STRICT subset, which is the only shape that makes
# a partial quorum mean anything.
_ISO_C = _mutate(_STRAIN_A, 60, 130)


def test_a_partial_quorum_run_carries_its_isoform_group_end_to_end(tmp_path):
    """The plumbing, end to end: consensus.constitutive_start_positions knows
    each window's covered isoform indices, and until now run_pipeline dropped
    them on the floor (``[s for s, _w, _idx in triples]``). They now reach the
    candidate, the filter scope and results.csv.

    The counter must still be ZERO. Windows covered by only {A, B} mismatch C in
    the junction band, and C is exactly the isoform a quorum of 2 permits them
    to miss — so a non-zero count here means the scope was widened to every
    record and partial-quorum design is broken.
    """
    fa = write_fasta(tmp_path / "iso.fa", [
        ("A", _STRAIN_A), ("B", _STRAIN_A), ("C", _ISO_C)])
    out = str(tmp_path / "out")
    run_pipeline(RunConfig(
        fasta_path=fa, output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        filters=FilterThresholds(gc_min=0.20, gc_max=0.90, max_n_per_half=1,
                                 max_homopolymer_run=6,
                                 require_tn_junction=False),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        multi_strain=True, consensus_strategy="constitutive",
        isoform_min_count=2,          # PARTIAL: 2 of 3
        discovery_only=True, tile_length=50, tile_step=1,
    ))

    with open(os.path.join(out, "results.csv")) as fh:
        rows = list(csv.DictReader(fh))
    assert rows
    groups = {r["covered_record_ids"] for r in rows}
    # Both shapes must be present, or the fixture stopped testing what it says:
    # windows all three isoforms share, and windows only {A, B} do.
    assert "A;B;C" in groups
    assert "A;B" in groups
    # Never empty on this mode — an empty group would silently disable the gate.
    assert "" not in groups

    with open(os.path.join(out, "filter_report.json")) as fh:
        report = json.load(fh)
    assert report["counts"]["junction_mismatch"] == 0


def test_resolves_to_stays_blank_on_a_constitutive_run(tmp_path):
    """The column this feature deliberately did NOT reuse. ``resolves_to`` is
    documented as the record ids a VARIANT row matches, and these are plain
    rows (probe_role ""); folding the isoform group into it would have put a
    non-blank value on a row no downstream reader expects one on, including the
    GUI's lane logic. Pinned so a later "simplification" cannot quietly merge
    the two columns."""
    fa = write_fasta(tmp_path / "iso.fa", [
        ("A", _STRAIN_A), ("B", _STRAIN_A), ("C", _ISO_C)])
    out = str(tmp_path / "out")
    run_pipeline(RunConfig(
        fasta_path=fa, output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        filters=FilterThresholds(gc_min=0.20, gc_max=0.90, max_n_per_half=1,
                                 max_homopolymer_run=6,
                                 require_tn_junction=False),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        multi_strain=True, consensus_strategy="constitutive",
        isoform_min_count=2, discovery_only=True, tile_length=50, tile_step=1,
    ))
    with open(os.path.join(out, "results.csv")) as fh:
        rows = list(csv.DictReader(fh))
    assert rows
    assert {r["resolves_to"] for r in rows} == {""}
    assert {r["probe_role"] for r in rows} == {""}


def test_a_conserved_run_leaves_the_group_column_blank(tmp_path):
    """Non-regression for every other mode: the new column exists but is empty,
    so the run-level scope is still the only one in play."""
    fa = write_fasta(tmp_path / "s.fa", [("A", _STRAIN_A), ("B", _STRAIN_B)])
    out = str(tmp_path / "out")
    run_pipeline(RunConfig(
        fasta_path=fa, output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        filters=FilterThresholds(gc_min=0.20, gc_max=0.90, max_n_per_half=1,
                                 max_homopolymer_run=6,
                                 require_tn_junction=False),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        multi_strain=True, consensus_strategy="conserved",
        discovery_only=True, tile_length=50, tile_step=1,
    ))
    with open(os.path.join(out, "results.csv")) as fh:
        rows = list(csv.DictReader(fh))
    assert rows
    assert {r["covered_record_ids"] for r in rows} == {""}
