"""Isoform-collapse: constitutive (quorum) window selection + end-to-end run."""
import csv
import importlib.util

import pytest

from flexprobe.consensus import StrainAlignment, constitutive_start_positions

_HAS_BIO = importlib.util.find_spec("Bio") is not None

X = "ACGTACGTAC" * 6                       # 60 bp
GAP10 = "-" * 10 + X[10:]                   # same as X but first 10 bp gapped


def _align(seqs):
    return StrainAlignment(
        fragment_name="gene", strain_ids=[f"iso{i}" for i in range(len(seqs))],
        aligned_seqs=list(seqs), consensus=X,
    )


def test_constitutive_all_vs_quorum():
    aln = _align([X, X, GAP10])             # iso2 lacks the first 10 bp

    # Strict (all 3): only the window fully inside the shared region qualifies.
    strict = constitutive_start_positions(aln, tile_length=50, tile_step=1,
                                          min_isoforms=3)
    assert [s for s, _w, _idx in strict] == [10]
    start, shared, covered = strict[0]
    assert shared == X[10:60] and covered == (0, 1, 2)

    # Quorum of 2: every window is shared by at least iso0 & iso1.
    quorum = constitutive_start_positions(aln, tile_length=50, tile_step=1,
                                          min_isoforms=2)
    assert [s for s, _w, _idx in quorum] == list(range(0, 11))
    # the windows before col 10 are covered by exactly iso0 & iso1
    assert quorum[0][2] == (0, 1)


def test_constitutive_default_is_all():
    aln = _align([X, X, GAP10])
    default = constitutive_start_positions(aln, tile_length=50, tile_step=1)
    assert [s for s, _w, _idx in default] == [10]   # None -> all records


@pytest.mark.skipif(not _HAS_BIO, reason="Biopython needed to align isoforms")
def test_isoform_collapse_end_to_end(tmp_path):
    from flexprobe.models import (
        ScreenResources, ModalityDefinition, RunConfig,
    )
    from flexprobe.pipeline import run_pipeline

    # A 60-bp constitutive block shared by all isoforms, with different
    # flanking "exons" (one isoform skips the 5' exon entirely).
    shared = "CAGTGCGGCTCAATCCAAGTAATCAATGGATCAACTATCGGAACTGGCAGT"  # 50 bp
    shared += "ACGTACGTAC"                                            # 60 bp total
    iso1 = "TTTTGGGGCCCCAAAA" + shared + "GACTGACTGA"
    iso2 = shared + "CCTTGGAACC"                       # 5' exon skipped
    iso3 = "AGAGAGCTCT" + shared
    fa = tmp_path / "isoforms.fasta"
    fa.write_text(f">GENEX-201\n{iso1}\n>GENEX-202\n{iso2}\n>GENEX-203\n{iso3}\n")

    out = tmp_path / "out"
    cfg = RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=True),
        multi_strain=True, consensus_strategy="constitutive",
        isoform_min_count=0,            # all isoforms
    )
    cfg.filters.require_tn_junction = False
    res = run_pipeline(cfg)
    assert res["stats"]["pairs_output"] > 0

    rows = list(csv.DictReader(open(out / "results.csv")))
    assert rows
    # Every constitutive probe's window must be present in ALL three isoforms.
    for r in rows:
        w = r["mrna_target_50bp"].upper()
        assert w in iso1 and w in iso2 and w in iso3


@pytest.mark.skipif(not _HAS_BIO, reason="Biopython needed to align records")
def test_consensus_reference_picks_which_record_supplies_variant_bases(tmp_path):
    """The template decides the base at every position the records disagree on.

    A consensus carries an ambiguity code wherever the inputs differ, but an
    orderable probe must be concrete ACGT — so those positions resolve to ONE
    record. That means a single run only ever emits that record's flavour of a
    variant window and the other records' equally valid probes are never
    generated at all; recovering them requires re-running with a different
    template. Measured on the two OROV L segments, template 0 and template 1
    give 65 and 115 probes respectively and neither is a superset.
    """
    import csv
    from flexprobe.models import ScreenResources, ModalityDefinition, RunConfig
    from flexprobe.pipeline import run_pipeline
    base = ("CAGTGCGGCTCAATCCAAGTAATCAATGGATCAACTATCGGAACTGGCAGT"
            "ACGTACGTACGGATCAACTATCGGAACTGGCAGTACGTACGTAC")
    # One SNP, far enough inside that many windows span it.
    a = base
    b = base[:60] + ("T" if base[60] != "T" else "G") + base[61:]
    fa = tmp_path / "two.fasta"
    fa.write_text(f">REC_A\n{a}\n>REC_B\n{b}\n")

    def _run(ref_idx):
        out = tmp_path / f"out{ref_idx}"
        cfg = RunConfig(
            fasta_path=str(fa), output_dir=str(out),
            modality=ModalityDefinition(name="flex", version="v2"),
            screen=ScreenResources(skip_screen=True),
            multi_strain=True, consensus_strategy="constitutive",
            isoform_min_count=1,           # every window, no conservation filter
            consensus_reference=ref_idx,
        )
        cfg.filters.require_tn_junction = False
        cfg.filters.enforce_gc_range = False
        run_pipeline(cfg)
        return {r["mrna_target_50bp"].upper()
                for r in csv.DictReader(open(out / "results.csv"))}

    win_a, win_b = _run(0), _run(1)
    assert win_a and win_b
    # No N or IUPAC code may survive into an orderable probe either way.
    for w in win_a | win_b:
        assert set(w) <= set("ACGT")
    # The two templates genuinely differ — this is the whole point of the field.
    assert win_a != win_b, "consensus_reference had no effect on the output"
    # Each template's probes carry ITS base at the variant position, so neither
    # run is a superset of the other and only the union is complete.
    assert win_a - win_b and win_b - win_a


@pytest.mark.skipif(not _HAS_BIO, reason="Biopython needed to align records")
def test_consensus_reference_out_of_range_falls_back_to_first_record(tmp_path):
    """An index past the last record must not crash or silently rotate."""
    import csv
    from flexprobe.models import ScreenResources, ModalityDefinition, RunConfig
    from flexprobe.pipeline import run_pipeline
    seq = ("CAGTGCGGCTCAATCCAAGTAATCAATGGATCAACTATCGGAACTGGCAGT"
           "ACGTACGTACGGATCAACTATCGGAACTGGCAGTACGTACGTAC")
    fa = tmp_path / "two.fasta"
    fa.write_text(f">REC_A\n{seq}\n>REC_B\n{seq[:60]}T{seq[61:]}\n")

    def _run(ref_idx, tag):
        out = tmp_path / tag
        cfg = RunConfig(
            fasta_path=str(fa), output_dir=str(out),
            modality=ModalityDefinition(name="flex", version="v2"),
            screen=ScreenResources(skip_screen=True),
            multi_strain=True, consensus_strategy="constitutive",
            isoform_min_count=1, consensus_reference=ref_idx,
        )
        cfg.filters.require_tn_junction = False
        cfg.filters.enforce_gc_range = False
        run_pipeline(cfg)
        return {r["mrna_target_50bp"].upper()
                for r in csv.DictReader(open(out / "results.csv"))}

    assert _run(99, "oor") == _run(0, "first")
