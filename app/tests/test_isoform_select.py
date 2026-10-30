"""Tests for coverage-maximising per-gene probe selection."""
from __future__ import annotations

import textwrap

from flexprobe.isoform_select import (
    GeneIndex, build_gene_index, coverage_pick, gene_index_from_manifest,
    mismatches_for_identity, results_transcripts,
)


def test_mismatches_for_identity_matches_blast_rule():
    # Audit #26: 90% identity over a 50-mer == 5 mismatches, so the FASTA
    # membership path can be aligned with the engine's BLAST >=90% rule.
    assert mismatches_for_identity(50, 0.90) == 5
    assert mismatches_for_identity(25, 0.90) == 2   # rounds 2.5 -> 2
    assert mismatches_for_identity(50, 1.0) == 0


def _row(probe_id, transcript, start, target, score, **extra):
    r = {
        "probe_id": probe_id, "gene_id": transcript, "start": str(start),
        "end": str(start + len(target)), "mrna_target_50bp": target,
        "composite_score": str(score), "specificity_pass": "True",
        "reduced_confidence": "False",
    }
    r.update({k: str(v) for k, v in extra.items()})
    return r


def _index(seqs_by_transcript, gene_of, symbols=None):
    group_members = {}
    for t, g in gene_of.items():
        group_members.setdefault(g, set()).add(t)
    return GeneIndex(dict(gene_of), group_members,
                     {t: s.upper() for t, s in seqs_by_transcript.items()},
                     dict(symbols or {}))


def test_set_cover_picks_complementary_probes():
    # Gene G has 3 isoforms. No single probe binds all three, but two probes
    # together do. Coverage selection should find the covering pair.
    a = "A" * 30
    b = "C" * 30
    common = "G" * 30
    t1 = a + common          # isoform 1: has A-window + common
    t2 = b + common          # isoform 2: has B-window + common
    t3 = common              # isoform 3: only common
    idx = _index({"t1": t1, "t2": t2, "t3": t3},
                 {"t1": "G", "t2": "G", "t3": "G"})
    # probe_A binds t1 only; probe_B binds t2 only; probe_C (common) binds all 3.
    rows = [
        _row("pA", "t1", 0, a, 0.9),
        _row("pB", "t2", 0, b, 0.8),
        _row("pC", "t1", 30, common, 0.5),
    ]
    res = coverage_pick(rows, idx, n_per_gene=2, min_gap=0)
    chosen = {p for ids in res.selections.values() for p in ids}
    # Two probes, full coverage: the only covering pair within 2 is {pC, anything}
    # — pC alone covers all 3, so it must be picked first; the gene is satisfied.
    assert "pC" in chosen
    assert res.summary["flagged"] == 0


def test_flags_incomplete_coverage():
    # Isoform t3 has a unique window no probe targets -> can't be covered.
    common = "G" * 30
    t1 = "A" * 30 + common
    t2 = "A" * 30 + common
    t3 = "T" * 30            # totally distinct, unreachable
    idx = _index({"t1": t1, "t2": t2, "t3": t3},
                 {"t1": "G", "t2": "G", "t3": "G"})
    rows = [_row("pA", "t1", 0, "A" * 30, 0.9)]
    res = coverage_pick(rows, idx, n_per_gene=3, min_gap=0)
    assert len(res.flags) == 1
    flag = res.flags[0]
    assert flag.gene == "G"
    assert "t3" in flag.missing_isoforms
    assert "incomplete_isoform_coverage" in flag.reasons


def test_spacing_enforced_across_shared_isoform():
    # Two probes on the same isoform closer than min_gap must not both be chosen.
    seq = "ACGT" * 50  # 200 bp
    idx = _index({"t1": seq}, {"t1": "G"})
    w1 = seq[0:30]
    w2 = seq[20:50]   # starts 20 bp after w1 -> overlap/too close at min_gap=25
    rows = [_row("p1", "t1", 0, w1, 0.9), _row("p2", "t1", 20, w2, 0.8)]
    res = coverage_pick(rows, idx, n_per_gene=2, min_gap=25)
    chosen = {p for ids in res.selections.values() for p in ids}
    assert chosen == {"p1"}              # only one fits the spacing rule
    assert any("only_1_spaced_probes" in f.reasons for f in res.flags)


def test_fills_to_n_when_coverage_saturated():
    # One probe already covers the single isoform; with spaced alternatives the
    # gene should still receive up to n_per_gene probes.
    seq = "ACGTACGTAC" * 30  # 300 bp, varied enough for distinct windows
    idx = _index({"t1": seq}, {"t1": "G"})
    rows = [
        _row("p1", "t1", 0, seq[0:30], 0.9),
        _row("p2", "t1", 100, seq[100:130], 0.8),
        _row("p3", "t1", 200, seq[200:230], 0.7),
    ]
    res = coverage_pick(rows, idx, n_per_gene=3, min_gap=25)
    chosen = {p for ids in res.selections.values() for p in ids}
    assert chosen == {"p1", "p2", "p3"}
    assert res.summary["flagged"] == 0


def test_max_mismatches_recovers_divergent_isoform():
    common = "G" * 30
    t1 = "ACGTACGTACGTACGTACGTACGTACGTAC" + common
    # t2 differs by 1 base in the targeted window.
    t2 = "ACGTACGTACGTAAGTACGTACGTACGTAC" + common
    idx = _index({"t1": t1, "t2": t2}, {"t1": "G", "t2": "G"})
    target = t1[:30]
    rows = [_row("p1", "t1", 0, target, 0.9)]
    exact = coverage_pick(rows, idx, n_per_gene=1, min_gap=0, max_mismatches=0)
    assert exact.summary["flagged"] == 1           # misses t2 under exact match
    fuzzy = coverage_pick(rows, idx, n_per_gene=1, min_gap=0, max_mismatches=1)
    assert fuzzy.summary["flagged"] == 0           # 1-mismatch recovers t2


def test_build_gene_index_from_fasta(tmp_path):
    fa = tmp_path / "mini.fa"
    fa.write_text(textwrap.dedent("""\
        >ENST1.1 cdna gene:ENSG1.1 gene_symbol:FOO
        ACGTACGTAC
        >ENST2.1 cdna gene:ENSG1.1 gene_symbol:FOO
        ACGTACGTAC
        >ENST3.1 cdna gene:ENSG2.1 gene_symbol:BAR
        TTTTTTTTTT
    """))
    idx = build_gene_index(str(fa), needed_transcripts={"ENST1.1"})
    assert idx.gene_of["ENST2.1"] == "ENSG1.1"
    assert idx.group_members["ENSG1.1"] == {"ENST1.1", "ENST2.1"}
    assert idx.symbol_of["ENSG1.1"] == "FOO"
    # Sequences loaded for the needed gene's siblings, not the other gene.
    assert "ENST1.1" in idx.seqs and "ENST2.1" in idx.seqs
    assert "ENST3.1" not in idx.seqs


def test_results_transcripts():
    rows = [{"gene_id": "t1"}, {"gene_id": "t2"}, {"gene_id": "t1"}, {"gene_id": ""}]
    assert results_transcripts(rows) == {"t1", "t2"}


# --------------------------------------------------------------------------- #
# FASTA-free path: membership from the persisted isoform_hit_ids column
# --------------------------------------------------------------------------- #
def _col_row(probe_id, transcript, start, hit_ids, score):
    return {
        "probe_id": probe_id, "gene_id": transcript, "start": str(start),
        "end": str(start + 30), "mrna_target_50bp": "",   # no target needed
        "composite_score": str(score), "specificity_pass": "True",
        "reduced_confidence": "False", "isoform_hit_ids": ";".join(hit_ids),
    }


def test_column_mode_no_fasta_needed():
    # Sequence-free index (as built from a manifest's gene_isoforms map).
    idx = gene_index_from_manifest({"G": ["t1", "t2", "t3"]})
    # Membership comes straight from the column. pA covers t1+t2; pB covers t3.
    rows = [
        _col_row("pA", "t1", 0, ["t1", "t2"], 0.9),
        _col_row("pB", "t3", 200, ["t3"], 0.8),
    ]
    res = coverage_pick(rows, idx, n_per_gene=2, min_gap=25)
    chosen = {p for ids in res.selections.values() for p in ids}
    assert chosen == {"pA", "pB"}            # together cover t1,t2,t3
    assert res.summary["flagged"] == 0       # full coverage, FASTA never touched


def test_column_mode_flags_missing_isoform_from_manifest_universe():
    # Manifest says the gene has 3 isoforms, but probes only reach 2 of them.
    idx = gene_index_from_manifest({"G": ["t1", "t2", "t3"]})
    rows = [_col_row("pA", "t1", 0, ["t1", "t2"], 0.9)]
    res = coverage_pick(rows, idx, n_per_gene=3, min_gap=25)
    assert len(res.flags) == 1
    assert res.flags[0].missing_isoforms == ["t3"]
    assert "incomplete_isoform_coverage" in res.flags[0].reasons


def test_preferred_gap_prefers_roomier_spacing():
    # One isoform; three same-transcript probes equally cover it. With a probe
    # already at start 0, a candidate 60 bp away (>= preferred) should be chosen
    # over one only 30 bp away (>= min but < preferred), all else equal.
    idx = gene_index_from_manifest({"G": ["t1"]})
    # 30-bp probes: anchor [0,30). near [60,90) -> 30 bp gap (orange, >=25 <50);
    # far [90,120) -> 60 bp gap (green, >=50).
    rows = [
        _col_row("anchor", "t1", 0, ["t1"], 0.99),   # picked first (top score)
        _col_row("near", "t1", 60, ["t1"], 0.90),    # 30 bp gap -> orange tier
        _col_row("far", "t1", 90, ["t1"], 0.90),     # 60 bp gap -> green tier
    ]
    res = coverage_pick(rows, idx, n_per_gene=2, min_gap=25, preferred_gap=50)
    chosen = [p for ids in res.selections.values() for p in ids]
    assert "anchor" in chosen and "far" in chosen and "near" not in chosen
