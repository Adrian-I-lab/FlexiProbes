"""Gene-aware (group-by-gene) isoform-coverage scoring + gating.

Covers the whole-transcriptome path where input transcripts are grouped by
gene so sibling isoforms are on-target (not off-target) and each probe is
scored on how many of its gene's isoforms it binds.

The first group are pure-Python unit tests (no BLAST). The e2e test at the
bottom runs the real pipeline + bundled blastn and is skipped when blast+
isn't available.
"""
import shutil

import pytest

from flexprobe.models import ProbeCandidate
from flexprobe.screen import ENGINE_BLAST
from flexprobe.ranking import (
    DEFAULT_ISOFORM_COVERAGE_WEIGHT, RankingWeights, compute_composite_score,
)
from flexprobe.specificity import batch_evaluate
from flexprobe.pipeline import _build_gene_groups


# Three isoforms of gene A, one of gene B.
GENE_OF = {
    "ENST_A1": "G_A", "ENST_A2": "G_A", "ENST_A3": "G_A",
    "ENST_B1": "G_B",
}
GROUP_MEMBERS = {
    "G_A": {"ENST_A1", "ENST_A2", "ENST_A3"},
    "G_B": {"ENST_B1"},
}


def _hit(sseqid, pident=100.0, length=25):
    return {"sseqid": sseqid, "pident": pident, "length": length, "qseqid": "q"}


def _cand(gene_id):
    return ProbeCandidate(gene_id=gene_id, probe_id=f"{gene_id}_t0", sequence="A" * 50)


# --- gene grouping ---------------------------------------------------------

def test_build_gene_groups_parses_ensembl_gene_field():
    records = [
        ("ENST_A1", "cdna chromosome:GRCh38:1 gene:G_A gene_symbol:Aaa", "ACGT"),
        ("ENST_A2", "cdna gene:G_A gene_symbol:Aaa", "ACGT"),
        ("ENST_B1", "cdna gene:G_B", "ACGT"),
        ("LOOSE",   "no gene field here", "ACGT"),   # falls back to its own id
    ]
    gene_of, members = _build_gene_groups(records)
    assert gene_of == {"ENST_A1": "G_A", "ENST_A2": "G_A",
                       "ENST_B1": "G_B", "LOOSE": "LOOSE"}
    assert members["G_A"] == {"ENST_A1", "ENST_A2"}
    assert members["LOOSE"] == {"LOOSE"}


# --- coverage + off-target classification ----------------------------------

def test_full_isoform_coverage_all_siblings_hit():
    c = _cand("ENST_A1")
    hits = {c.probe_id: [_hit("ENST_A1"), _hit("ENST_A2"), _hit("ENST_A3")]}
    passed = batch_evaluate([c], hits, set(GENE_OF),
                            gene_of=GENE_OF, group_members=GROUP_MEMBERS)
    assert passed == 1
    assert c.specificity_pass is True
    assert c.isoform_total == 3
    assert c.isoform_hits == 3
    assert c.isoform_coverage == 1.0


def test_partial_isoform_coverage():
    c = _cand("ENST_A1")
    # Binds only one *other* sibling; its own transcript always counts too.
    hits = {c.probe_id: [_hit("ENST_A2")]}
    batch_evaluate([c], hits, set(GENE_OF),
                   gene_of=GENE_OF, group_members=GROUP_MEMBERS)
    assert c.isoform_total == 3
    assert c.isoform_hits == 2                 # A2 (hit) + A1 (own)
    assert c.isoform_coverage == pytest.approx(2 / 3)
    assert c.specificity_pass is True          # siblings are never off-target


def test_cross_gene_hit_is_off_target():
    c = _cand("ENST_A1")
    # A sibling hit (fine) plus a different gene's transcript (off-target).
    hits = {c.probe_id: [_hit("ENST_A2"), _hit("ENST_B1")]}
    batch_evaluate([c], hits, set(GENE_OF),
                   gene_of=GENE_OF, group_members=GROUP_MEMBERS)
    assert c.specificity_pass is False
    assert c.specificity_status == "failed_off_target"
    # Coverage is still annotated even on a failed probe.
    assert c.isoform_total == 3
    assert c.isoform_hits == 2                 # A1 (own) + A2


def test_single_isoform_gene_is_full_coverage():
    c = _cand("ENST_B1")
    passed = batch_evaluate([c], {c.probe_id: []}, set(GENE_OF),
                            gene_of=GENE_OF, group_members=GROUP_MEMBERS)
    assert passed == 1
    assert c.isoform_total == 1
    assert c.isoform_coverage == 1.0


def test_non_gene_aware_leaves_coverage_default():
    # Without gene_of/group_members the candidate keeps the neutral defaults.
    c = _cand("ENST_A1")
    batch_evaluate([c], {c.probe_id: [_hit("ENST_A1")]}, {"ENST_A1"})
    assert c.isoform_total == 1
    assert c.isoform_coverage == 1.0


# --- ranking ----------------------------------------------------------------

def test_isoform_coverage_raises_score():
    weights = RankingWeights(isoform_coverage=DEFAULT_ISOFORM_COVERAGE_WEIGHT)
    high = ProbeCandidate(gene_id="g", probe_id="hi", gc_content=0.5,
                          gc_left=0.5, gc_right=0.5, isoform_coverage=1.0)
    low = ProbeCandidate(gene_id="g", probe_id="lo", gc_content=0.5,
                         gc_left=0.5, gc_right=0.5, isoform_coverage=0.2)
    s_high = compute_composite_score(high, seq_length=100, weights=weights)
    s_low = compute_composite_score(low, seq_length=100, weights=weights)
    assert s_high > s_low
    # Scores are normalised by the weight total, so the gap is the coverage
    # weight times the coverage difference, divided by the total. Neither
    # candidate has strain_binding, so mismatch_placement is inapplicable and
    # its weight is out of that total (see ranking.compute_composite_score).
    denom = weights.total() - weights.mismatch_placement
    assert s_high - s_low == pytest.approx(
        DEFAULT_ISOFORM_COVERAGE_WEIGHT * (1.0 - 0.2) / denom
    )
    # Normalisation keeps the score within [0, 1] even with the extra weight.
    assert 0.0 <= s_high <= 1.0


def test_isoform_weight_zero_is_no_op():
    # Default weights (isoform_coverage=0) ignore coverage entirely.
    a = ProbeCandidate(gc_content=0.5, gc_left=0.5, gc_right=0.5, isoform_coverage=1.0)
    b = ProbeCandidate(gc_content=0.5, gc_left=0.5, gc_right=0.5, isoform_coverage=0.0)
    assert compute_composite_score(a, 100) == compute_composite_score(b, 100)


# --- end-to-end with real blastn -------------------------------------------

@pytest.mark.skipif(
    shutil.which("blastn") is None and shutil.which("makeblastdb") is None,
    reason="ncbi-blast+ not installed",
)
def test_pipeline_group_by_gene_e2e(tmp_path):
    from flexprobe.blast import build_db, discover_blast
    from flexprobe.models import (
        ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
        SpecificityThresholds,
    )
    from flexprobe.pipeline import run_pipeline

    # Two identical isoforms of gene A (so every A probe binds both) + a
    # distinct gene B. The transcriptome is its own BLAST DB.
    iso = (
        "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGTTACCGATGCATGC"
        "ATGCATGCATGCATGCATGCAACGTACGTACGTACGTACGTTGCATGCATG"
        "CATGCATGCATG"
    )
    geneB = (
        "TTGGAATTCCGGAATTCCAATTGGCCAATTGGCCTTAAGGCCTTAAGGCCT"
        "TAAGGTTCCAAGGTTCCAAGGAATTCCGGTTAACCGGTTAACCAATTGGCCA"
    )
    fa = tmp_path / "transcriptome.fasta"
    fa.write_text(
        f">ENST_A1 cdna gene:GENEA\n{iso}\n"
        f">ENST_A2 cdna gene:GENEA\n{iso}\n"
        f">ENST_B1 cdna gene:GENEB\n{geneB}\n"
    )
    db = str(tmp_path / "txdb")
    build_db(discover_blast()["makeblastdb"], str(fa), db, title="tx")
    out = tmp_path / "out"
    cfg = RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=False, db_paths=[db],
                             engine=ENGINE_BLAST),
        # The junction-duplex gate is k-mer-only, so a BLAST config now has
        # to disable it explicitly (models.RunConfig.validate refuses the
        # combination rather than reporting a gate that never ran).
        specificity=SpecificityThresholds(off_target_junction_duplex_bp=0),
        export=ExportSettings(overwrite=True), discovery_only=True,
        group_by_gene=True,
    )
    result = run_pipeline(cfg)
    stats = result["stats"]
    assert stats["group_by_gene"] is True

    by_gene = {}
    for pair in result["pairs"]:
        by_gene.setdefault(pair.lhs.gene_id, []).append(pair.lhs)
    # Gene A's two identical isoforms: every probe binds both -> total 2, full
    # coverage, and crucially NOT rejected as off-target (siblings excused).
    a_probes = by_gene.get("ENST_A1", []) + by_gene.get("ENST_A2", [])
    assert a_probes, "expected probes from gene A isoforms"
    for c in a_probes:
        assert c.isoform_total == 2
        assert c.isoform_coverage == 1.0
        assert c.specificity_pass is True
    # Gene B is a single isoform -> total 1, coverage 1.0.
    for c in by_gene.get("ENST_B1", []):
        assert c.isoform_total == 1
        assert c.isoform_coverage == 1.0

    # results.csv (the user-facing artifact) carries the new columns.
    import csv
    import json
    with open(out / "results.csv") as fh:
        rows = list(csv.DictReader(fh))
    assert {"isoform_coverage", "isoform_hits", "isoform_total",
            "isoform_hit_ids"} <= set(rows[0])
    a_rows = [r for r in rows if r["gene_id"] in ("ENST_A1", "ENST_A2")]
    assert a_rows and all(r["isoform_total"] == "2" for r in a_rows)
    assert all(float(r["isoform_coverage"]) == 1.0 for r in a_rows)
    # The persisted isoform_hit_ids name the actual bound transcripts, and their
    # count matches isoform_hits — so selection needs no FASTA.
    for r in a_rows:
        ids = [t for t in r["isoform_hit_ids"].split(";") if t]
        assert len(ids) == int(r["isoform_hits"])
        assert set(ids) <= {"ENST_A1", "ENST_A2"}

    # The manifest records each gene's full isoform set (the FASTA-free universe).
    with open(out / "manifest.json") as fh:
        man = json.load(fh)
    gene_isoforms = man["params"]["gene_isoforms"]
    assert set(gene_isoforms.get("GENEA", [])) == {"ENST_A1", "ENST_A2"}


def test_pipeline_group_by_gene_e2e_kmer_engine(tmp_path):
    """The same gene-aware assertions on the DEFAULT engine, and with no BLAST.

    Phase 2 of the k-mer migration claimed on-target classification "comes free"
    — >=90% identity over >=21 bp is <=2 mismatches over 25, already inside the
    <=4 search, so require_on_target, the strain rules and isoform coverage need
    no second search and no new code. That claim now sits on the DEFAULT path,
    so it is asserted end to end rather than argued.

    Deliberately NOT skipped when blastn is absent: needing no BLAST at all is
    half the point, and a test that quietly skips would hide a regression on the
    engine most users now run.
    """
    from flexprobe.models import (
        ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
    )
    from flexprobe.pipeline import run_pipeline
    from flexprobe.screen import ENGINE_KMER

    iso = (
        "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGTTACCGATGCATGC"
        "ATGCATGCATGCATGCATGCAACGTACGTACGTACGTACGTTGCATGCATG"
        "CATGCATGCATG"
    )
    geneB = (
        "TTGGAATTCCGGAATTCCAATTGGCCAATTGGCCTTAAGGCCTTAAGGCCT"
        "TAAGGTTCCAAGGTTCCAAGGAATTCCGGTTAACCGGTTAACCAATTGGCCA"
    )
    fa = tmp_path / "transcriptome.fasta"
    fa.write_text(
        f">ENST_A1 cdna gene:GENEA\n{iso}\n"
        f">ENST_A2 cdna gene:GENEA\n{iso}\n"
        f">ENST_B1 cdna gene:GENEB\n{geneB}\n"
    )
    out = tmp_path / "out"
    cfg = RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        # The k-mer engine indexes the FASTA directly — no database to build.
        screen=ScreenResources(skip_screen=False, db_paths=[str(fa)],
                             engine=ENGINE_KMER),
        export=ExportSettings(overwrite=True), discovery_only=True,
        group_by_gene=True,
    )
    result = run_pipeline(cfg)
    assert result["stats"]["group_by_gene"] is True

    by_gene = {}
    for pair in result["pairs"]:
        by_gene.setdefault(pair.lhs.gene_id, []).append(pair.lhs)

    a_probes = by_gene.get("ENST_A1", []) + by_gene.get("ENST_A2", [])
    assert a_probes, "expected probes from gene A isoforms"
    for c in a_probes:
        # Sibling isoforms must be excused, not counted as off-target.
        assert c.isoform_total == 2
        assert c.isoform_coverage == 1.0
        assert c.specificity_pass is True
    for c in by_gene.get("ENST_B1", []):
        assert c.isoform_total == 1

    # And the property the engine exists for, on results.csv — the surface that
    # actually reaches a user. `pair.lhs.offtarget_recount_basis` now agrees
    # with it: the halves are copies of the scored candidate (pipeline._pair_half)
    # rather than fresh ProbeCandidates built from a curated field list, so a
    # fully screened run no longer reads as "not_screened" in memory. That
    # agreement is pinned column-for-column in test_provenance.py::
    # test_pair_halves_carry_the_full_provenance_of_their_results_row.
    import csv
    with open(out / "results.csv") as fh:
        rows = list(csv.DictReader(fh))
    assert rows
    assert {r["offtarget_recount_basis"] for r in rows} == {"exact"}
    assert {p.lhs.offtarget_recount_basis for p in result["pairs"]} == {"exact"}
    assert {r["screen_engine"] for r in rows} == {ENGINE_KMER}
    # Full coverage by construction — no reconstructed positions anywhere.
    assert {r["offtarget_min_qcov"] for r in rows} <= {"1.0000", "n/a"}
