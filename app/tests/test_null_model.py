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
# flexiprobes. If not, see <https://www.gnu.org/licenses/>.
"""Null-model calibration of the off-target gate.

The point of the calibration is that the number it reports is trustworthy, so
these tests pin the two things that could quietly make it wrong: a shuffle that
does not actually preserve dinucleotide composition (which would understate how
often noise finds a BLAST match), and a sample that is not really bounded
(which would make the stage unaffordable at whole-transcriptome scale).
"""
import json
import random

import pytest

from flexprobe import null_model as nm
from flexprobe.models import ProbeCandidate, RunConfig
from flexprobe.specificity import RECOUNT_CONSERVATIVE, RECOUNT_EXACT


SEQS = [
    "ACGTACGTTTGGCCATATATCGCGATTACAGGCATTAGCCATGACTGACT",
    "AAAAACCCCCGGGGGTTTTTACACACGTGTGTACGTACGTAAGGCCTTAA",
    "GCGCGCGCATATATATCGCGCGCGTATATATAGCGCGCGCATATATATCG",
]


# ---------------------------------------------------------------------------
# The shuffle
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("seq", SEQS)
def test_shuffle_preserves_length_composition_and_dinucleotides(seq):
    """The whole reason for an Eulerian shuffle rather than random.shuffle.

    BLAST word seeding and duplex stability both track dinucleotide
    composition, so a null that only preserved base composition would find
    fewer matches than real sequence for a reason that has nothing to do with
    specificity — and would understate the gate's false-positive rate.
    """
    rng = random.Random(0)
    d0, _ = nm._count_and_lists(seq)
    for _ in range(50):
        sh = nm.dinucleotide_shuffle(seq, rng)
        assert len(sh) == len(seq)
        assert sorted(sh) == sorted(seq)
        assert nm._count_and_lists(sh)[0] == d0
        assert nm.preserves_dinucleotides(seq, sh)


def test_shuffle_preserves_both_ends():
    """First and last base are fixed points of the Eulerian construction.

    Load-bearing here: the 50-mer is split into two 25 nt halves and each half
    is BLASTed, so a shuffle that moved the ends would change what the junction
    looks like relative to the real window.
    """
    rng = random.Random(1)
    for seq in SEQS:
        for _ in range(20):
            sh = nm.dinucleotide_shuffle(seq, rng)
            assert sh[0] == seq[0] and sh[-1] == seq[-1]


def test_shuffle_actually_shuffles():
    """A "null" that returns its input would report a false-positive rate of 0."""
    rng = random.Random(2)
    seq = SEQS[0]
    assert any(nm.dinucleotide_shuffle(seq, rng) != seq for _ in range(20))


def test_shuffle_degrades_gracefully_on_non_acgt():
    """Ambiguity codes must not raise — consensus windows carry N and IUPAC.

    The fallback keeps length and composition but not dinucleotide counts, so
    preserves_dinucleotides is False and the caller counts it as degraded
    rather than silently treating it as a proper null.
    """
    rng = random.Random(3)
    seq = "ACGTNNNNACGTACGTRYKMACGTACGTACGTACGTACGTACGTACGTAC"
    sh = nm.dinucleotide_shuffle(seq, rng)
    assert len(sh) == len(seq)
    assert sorted(sh) == sorted(seq)


# ---------------------------------------------------------------------------
# Bounding the cost
# ---------------------------------------------------------------------------
def _cands(n):
    return [
        ProbeCandidate(gene_id="g", probe_id=f"p{i:04d}",
                       sequence=SEQS[i % len(SEQS)], start=i, end=i + 50)
        for i in range(n)
    ]


def test_sample_is_bounded_regardless_of_run_size():
    """This cap is the only reason the stage is affordable at 35M candidates."""
    rng = random.Random(0)
    assert len(nm.sample_windows(_cands(100_000), max_windows=200, rng=rng)) == 200


def test_sample_returns_everything_when_under_the_cap():
    rng = random.Random(0)
    assert len(nm.sample_windows(_cands(7), max_windows=200, rng=rng)) == 7


def test_sample_is_random_not_truncated():
    """Truncation would draw every window from the first gene / the 5' end.

    Their composition is not the run's composition, so the calibration would be
    of a region rather than of the run.
    """
    rng = random.Random(0)
    picked = nm.sample_windows(_cands(1000), max_windows=50, rng=rng)
    assert [c.probe_id for c in picked] != [c.probe_id for c in _cands(50)]


def test_sample_is_reproducible_for_a_seed():
    a = nm.sample_windows(_cands(1000), max_windows=50, rng=random.Random(42))
    b = nm.sample_windows(_cands(1000), max_windows=50, rng=random.Random(42))
    assert [c.probe_id for c in a] == [c.probe_id for c in b]


def test_null_candidates_are_namespaced_and_keep_their_gene():
    """gene_id must survive: on-target classification depends on it, so a null
    window has to be judged by exactly the rule the real window was judged by.
    The probe_id must NOT survive unchanged, or a null row could be mistaken
    for a design."""
    rng = random.Random(0)
    real = _cands(3)
    nulls, degraded = nm.build_null_candidates(real, shuffles=4, rng=rng)
    assert len(nulls) == 12
    assert degraded == 0
    assert all(c.gene_id == "g" for c in nulls)
    assert all("__null" in c.probe_id for c in nulls)
    assert len({c.probe_id for c in nulls}) == 12


# ---------------------------------------------------------------------------
# The reported statistics
# ---------------------------------------------------------------------------
def test_summarise_reports_the_rejection_fraction():
    nulls = _cands(10)
    for i, c in enumerate(nulls):
        c.specificity_pass = i >= 3            # 3 of 10 rejected
        c.specificity_status = "passed" if i >= 3 else "failed_off_target"
        # Exact recount is a PRECONDITION of reporting a rate, not an incidental
        # detail of the fixture -- see the measurability tests below.
        c.offtarget_recount_basis = RECOUNT_EXACT
    out = nm.summarise(nulls, n_real_windows=40, n_sampled=10, shuffles=1,
                       degraded=0, seed=42)
    assert out["measurable"] is True
    assert out["null_rejected"] == 3
    assert out["false_positive_rate"] == pytest.approx(0.3)
    assert out["null_status_counts"]["failed_off_target"] == 3
    assert out["windows_in_run"] == 40


def test_wilson_interval_matches_the_standard_formula():
    """Pin `_wilson_interval` against precomputed values (the standard
    Wilson-score formula, phat=0.1 / n=10 / z=1.96), so a future refactor
    can't silently drift onto the wrong statistic (e.g. the naive Wald
    interval, which is narrower and can leave [0, 1])."""
    lo, hi = nm._wilson_interval(0.1, 10)
    assert lo == pytest.approx(0.017876213095072868, abs=1e-9)
    assert hi == pytest.approx(0.4041500267952385, abs=1e-9)
    # Degenerate n -> the widest honest answer, not a crash.
    assert nm._wilson_interval(0.1, 0) == (0.0, 1.0)


def test_ci_uses_windows_sampled_and_is_wider_than_a_naive_flat_binomial_ci():
    """The sampling is windows_sampled windows x shuffles_per_window shuffles,
    and shuffles of the SAME window are correlated (pseudo-replication) --
    treating null_windows_screened as a flat i.i.d. binomial n understates the
    true uncertainty. This is the regression guard: the reported interval must
    be computed on windows_sampled (here 10), not null_windows_screened (here
    200), and must therefore be visibly wider than the naive interval a flat
    binomial n=200 would give at the SAME observed rate."""
    windows_sampled, shuffles = 10, 20
    nulls = _cands(windows_sampled * shuffles)          # 200 null candidates
    for i, c in enumerate(nulls):
        c.specificity_pass = i >= 20                    # 20 of 200 rejected
        c.specificity_status = "passed" if i >= 20 else "failed_off_target"
        c.offtarget_recount_basis = RECOUNT_EXACT

    out = nm.summarise(nulls, n_real_windows=500, n_sampled=windows_sampled,
                       shuffles=shuffles, degraded=0, seed=42)

    assert out["measurable"] is True
    assert out["false_positive_rate"] == pytest.approx(0.1)
    ci = out["false_positive_rate_ci95"]
    assert ci["effective_n"] == windows_sampled          # NOT 200
    assert ci["lower"] == pytest.approx(0.017876213095072868, abs=1e-6)
    assert ci["upper"] == pytest.approx(0.4041500267952385, abs=1e-6)

    # The clustering correction must actually do something: a flat binomial
    # interval on the same total N (200) is a much tighter (wrong) claim.
    naive_lo, naive_hi = nm._wilson_interval(
        out["false_positive_rate"], out["null_windows_screened"])
    reported_width = ci["upper"] - ci["lower"]
    naive_width = naive_hi - naive_lo
    assert reported_width > naive_width, (
        "the cluster-conservative CI must be wider than the naive flat-n CI, "
        "or the clustering correction isn't doing anything"
    )
    assert ci["lower"] <= naive_lo and ci["upper"] >= naive_hi


def test_summarise_handles_an_empty_population():
    out = nm.summarise([], n_real_windows=0, n_sampled=0, shuffles=20,
                       degraded=0, seed=42)
    assert out["false_positive_rate"] == 0.0
    assert out["null_windows_screened"] == 0


def test_write_report_is_valid_json(tmp_path):
    payload = nm.summarise(_cands(2), n_real_windows=2, n_sampled=2,
                           shuffles=1, degraded=0, seed=42)
    path = nm.write_report(str(tmp_path), payload)
    assert path.endswith(nm.REPORT_FILE)
    with open(path) as fh:
        assert json.load(fh)["seed"] == 42
    # Written atomically — no .tmp left behind for a consumer to trip over.
    assert not (tmp_path / (nm.REPORT_FILE + ".tmp")).exists()


# ---------------------------------------------------------------------------
# Config wiring
# ---------------------------------------------------------------------------
def test_calibration_is_off_by_default():
    """Turning it on changes every run's runtime and output surface, so it is
    an explicit opt-in rather than a new default."""
    assert RunConfig(fasta_path="x", output_dir="y").null_model_shuffles == 0


def test_config_records_the_calibration_settings():
    cfg = RunConfig(fasta_path="x", output_dir="y", null_model_shuffles=20,
                    null_model_max_windows=50, null_model_seed=7)
    assert cfg.null_model_shuffles == 20
    assert cfg.null_model_max_windows == 50
    assert cfg.null_model_seed == 7


# ---------------------------------------------------------------------------
# Measurability: the calibration refuses to report a rate it cannot stand behind
#
# These pin the fix for the defect where the null screen ran BLAST while the run
# used the k-mer engine. BLAST recounts off-target mismatches by reading subject
# bases back through `blastdbcmd -entry`; against a k-mer run's plain-FASTA
# references there is no accession index, so every chance core scored as a
# 0-mismatch off-target and the calibration reported a FABRICATED 100%
# false-positive rate. The CLI refuses to start a BLAST run in that state. This
# path had no equivalent guard, so it published the number instead of refusing.
# ---------------------------------------------------------------------------
def test_summarise_refuses_to_report_when_a_window_was_not_recounted_exactly():
    """One inexact window is enough. A conservative recount over-rejects and
    never over-passes, so it can only push the rate UP -- the exact direction
    that manufactures the "this gate rejects on database size" conclusion."""
    nulls = _cands(10)
    for c in nulls:
        c.specificity_pass = False
        c.specificity_status = "failed_off_target"
        c.offtarget_recount_basis = RECOUNT_EXACT
    nulls[0].offtarget_recount_basis = RECOUNT_CONSERVATIVE

    out = nm.summarise(nulls, n_real_windows=40, n_sampled=10, shuffles=1,
                       degraded=0, seed=42)

    assert out["measurable"] is False
    # Absent, not zero and not a sentinel: a reader doing .get() must receive
    # None rather than something that can be plotted as a measurement.
    assert "false_positive_rate" not in out
    # No interval either -- there is nothing to bound around a suppressed
    # point estimate, and a present CI would wrongly imply one was measured.
    assert "false_positive_rate_ci95" not in out
    assert "1 of 10" in out["unmeasurable_reason"]
    # The evidence is still reported, so the run can be diagnosed from the file.
    assert out["null_recount_basis_counts"][RECOUNT_CONSERVATIVE] == 1
    assert out["null_rejected"] == 10


def test_the_default_recount_basis_alone_makes_a_run_unmeasurable():
    """A candidate that never reached the specificity evaluator keeps the model
    default. That is 'no evidence', which must not be read as 'exact'."""
    nulls = _cands(4)                       # untouched: basis == "not_screened"
    out = nm.summarise(nulls, n_real_windows=4, n_sampled=4, shuffles=1,
                       degraded=0, seed=42)
    assert out["measurable"] is False
    assert "false_positive_rate" not in out
    assert "false_positive_rate_ci95" not in out


def test_summarise_reports_the_basis_breakdown_when_measurable():
    nulls = _cands(5)
    for c in nulls:
        c.specificity_pass = True
        c.specificity_status = "passed"
        c.offtarget_recount_basis = RECOUNT_EXACT
    out = nm.summarise(nulls, n_real_windows=5, n_sampled=5, shuffles=1,
                       degraded=0, seed=42)
    assert out["measurable"] is True
    assert out["null_recount_basis_counts"] == {RECOUNT_EXACT: 5}
    assert out["false_positive_rate"] == 0.0
    # measurable => the interval is present too, bounding a rate of 0.0 with
    # windows_sampled=5 as the effective n (not the naive n=5 here either,
    # since shuffles=1 makes them equal -- see the dedicated CI test above).
    assert out["false_positive_rate_ci95"]["effective_n"] == 5
    assert out["false_positive_rate_ci95"]["lower"] == 0.0


def test_an_unmeasurable_payload_still_serialises(tmp_path):
    """The report is written either way -- refusing to report a RATE is not the
    same as writing no file, and the reason has to reach the user somehow."""
    out = nm.summarise(_cands(3), n_real_windows=3, n_sampled=3, shuffles=1,
                       degraded=0, seed=42)
    path = nm.write_report(str(tmp_path), out)
    with open(path) as fh:
        payload = json.load(fh)
    assert payload["measurable"] is False
    assert "false_positive_rate" not in payload


# ---------------------------------------------------------------------------
# The null screen must use the RUN's engine
# ---------------------------------------------------------------------------
def test_null_screen_inherits_the_runs_engine(tmp_path, monkeypatch):
    """The whole point of the calibration is comparability, so screening the
    nulls with a different engine than the run measures nothing. This omission
    is what produced the fabricated 100% rate; pin it so it cannot come back."""
    from flexprobe import pipeline as pl
    from flexprobe import specificity as sp

    seen = {}

    def _fake_screen(passing, db_paths, output_dir, **kw):
        seen.update(kw)
        return {}

    monkeypatch.setattr(pl, "_run_screen", _fake_screen)
    monkeypatch.setattr(sp, "batch_evaluate",
                        lambda *a, **k: None)

    cfg = RunConfig(fasta_path="x", output_dir=str(tmp_path),
                    null_model_shuffles=2, null_model_max_windows=3)
    cfg.screen.engine = "kmer"
    cfg.screen.exclude_alt_scaffolds = True

    nulls = _cands(3)
    for c in nulls:
        c.specificity_pass = True
        c.specificity_status = "passed"
        c.offtarget_recount_basis = RECOUNT_EXACT

    pl._run_null_model(
        cfg, nulls, str(tmp_path),
        spec_kwargs={}, target_ids=set(),
        half_len=25, emit=lambda *a, **k: None, cancel=None,
    )

    assert seen["engine"] == "kmer"
    # The reference-record scoping has to ride along too: a null screened
    # against a different record set is not comparable either.
    assert seen["exclude_non_primary"] is True
