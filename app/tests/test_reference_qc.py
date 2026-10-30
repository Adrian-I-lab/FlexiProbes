"""Reference composition QC: classification, exclusion, tally, and silence.

What these pin down, in the order the failures would hurt:

1. An ALT/PATCH record must be recognised and, by default, kept OUT of the
   searchable reference. On Ensembl's public GRCh38 cdna.all those records are
   7.6% of bases, and they are where spurious perfect-50bp off-target calls
   come from: good probes rejected for matching a second copy of their own gene.
2. A reference WITHOUT Ensembl headers must degrade to "unknown", change
   nothing, and say nothing. Most runs are small custom references (viral
   panels), and a run that warns about scaffolds it does not have is worse than
   no check at all.
3. Pseudogenes must be counted and NEVER dropped: a transcribed pseudogene
   missing from the reference is an off-target the screen cannot see, i.e. a
   bad probe that passes.
"""
import json
import logging
import os
from pathlib import Path

import numpy as np
import pytest

from flexprobe import kmer_index as K
from flexprobe import reference_qc as RQ


# Real header shapes, copied from Ensembl Homo_sapiens.GRCh38.cdna.all.fa.
PRIMARY_HDR = (
    "ENST00000415118.1 cdna chromosome:GRCh38:14:22438547:22438554:1 "
    "gene:ENSG00000223997.1 gene_biotype:TR_D_gene "
    "transcript_biotype:TR_D_gene gene_symbol:TRDD1 "
    "description:T cell receptor delta diversity 1 "
    "[Source:HGNC Symbol;Acc:HGNC:12254]"
)
PATCH_HDR = (
    "ENST00000632585.1 cdna scaffold:GRCh38:HG2290_PATCH:6058:6610:1 "
    "gene:ENSG00000282172.1 gene_biotype:IG_V_gene "
    "transcript_biotype:IG_V_gene gene_symbol:IGKV5-2"
)
ALT_HDR = (
    "ENST00000611780.1 cdna scaffold:GRCh38:HSCHR6_MHC_APD_CTG1:"
    "29942470:29945883:1 gene:ENSG00000277099.1 gene_biotype:protein_coding"
)

# Real header shapes, copied from Ensembl
# Rattus_norvegicus.mRatBN7.2.cdna.all.fa. Rat's location field uses the region
# type ``primary_assembly``, not ``chromosome`` — Ensembl spells it that way for
# every one of that file's 46,856 records. Before this was recognised the whole
# rat reference classified NON_PRIMARY, which left ``should_exclude`` correctly
# dropping nothing (it requires PRIMARY > 0) but made ``warn_if_searching_alt``
# tell the user that 46,856 records "are being SEARCHED as off-target subjects"
# and to remove the flag that would exclude them. Every one of those records is
# on the primary assembly, so the advice was the opposite of the truth.
RAT_PRIMARY_HDR = (
    "ENSRNOT00000094517.1 cdna primary_assembly:mRatBN7.2:15:26770499:"
    "26770801:1 gene:ENSRNOG00000064161.1 gene_biotype:TR_V_gene "
    "transcript_biotype:TR_V_gene"
)
# 69 of that file's records sit on an unplaced scaffold that Ensembl still
# labels ``primary_assembly``. The NAME check is what keeps them non-primary,
# and it is the reason the region type alone must not decide.
RAT_UNPLACED_HDR = (
    "ENSRNOT00000048256.5 cdna primary_assembly:mRatBN7.2:JACYVU010000238.1:"
    "73278:127359:-1 gene:ENSRNOG00000063914.1 gene_biotype:protein_coding "
    "transcript_biotype:protein_coding gene_symbol:LOC120099628 "
    "description:uncharacterized LOC120099628 "
    "[Source:RGD Symbol;Acc:41017459]"
)


def _hdr(tid, region, biotype, region_type="chromosome"):
    return (f"{tid} cdna {region_type}:GRCh38:{region}:100:200:1 "
            f"gene:ENSG{tid} gene_biotype:{biotype} "
            f"transcript_biotype:{biotype} description:made up")


def _seq(n, seed=0):
    """A deterministic non-repeating sequence.

    Non-repeating matters: a tiled motif matches itself at every period, so a
    "one hit vs two" assertion about record exclusion would be swamped by the
    sequence's own internal repeats rather than measuring anything.
    """
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), size=n))


@pytest.fixture(autouse=True)
def _clean_registry():
    """The measured-composition registry is process-global; keep tests apart."""
    RQ.clear_registry()
    yield
    RQ.clear_registry()


# ---------------------------------------------------------------------------
# (1) classification
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("region,expected", [
    ("1", RQ.PRIMARY),
    ("9", RQ.PRIMARY),
    ("22", RQ.PRIMARY),
    ("X", RQ.PRIMARY),
    ("Y", RQ.PRIMARY),
    ("MT", RQ.PRIMARY),
    # Mouse goes to 19 and this must keep working unchanged for it.
    ("19", RQ.PRIMARY),
    # Rat goes to 20.
    ("20", RQ.PRIMARY),
    # Rat's unplaced scaffolds, which Ensembl also labels primary_assembly.
    ("JACYVU010000238.1", RQ.UNPLACED),
    ("MU150189.1", RQ.UNPLACED),
    # Unplaced and unlocalised: unique sequence, RETAINED.
    ("KI270728.1", RQ.UNPLACED),
    ("GL000009.2", RQ.UNPLACED),
    # Copies of primary sequence: still excluded.
    ("HG1524_PATCH", RQ.ALT_PATCH),
    ("HSCHR6_MHC_APD_CTG1", RQ.ALT_PATCH),
])
def test_chromosome_name_decides_primary(region, expected):
    assembly, _bio = RQ.classify_header(_hdr("T1", region, "protein_coding"))
    assert assembly == expected


def test_a_scaffold_region_type_is_never_primary():
    """Even a scaffold NAMED like a chromosome is not the primary assembly."""
    assembly, _ = RQ.classify_header(
        _hdr("T1", "1", "protein_coding", region_type="scaffold"))
    assert assembly != RQ.PRIMARY
    # UNPLACED, not ALT_PATCH: nothing in the name identifies it as an
    # alternate representation, and an unrecognised record is retained.
    assert assembly == RQ.UNPLACED


def test_real_ensembl_headers():
    assert RQ.classify_header(PRIMARY_HDR) == (RQ.PRIMARY, "TR_D_gene")
    assert RQ.classify_header(PATCH_HDR) == (RQ.NON_PRIMARY, "IG_V_gene")
    assert RQ.classify_header(ALT_HDR) == (RQ.NON_PRIMARY, "protein_coding")
    assert RQ.seq_region_of(PATCH_HDR) == "HG2290_PATCH"
    assert RQ.seq_region_of(ALT_HDR) == "HSCHR6_MHC_APD_CTG1"


def test_real_rat_headers_use_the_primary_assembly_region_type():
    """Rat's whole cdna.all is ``primary_assembly:``, and it is primary.

    Measured on Ensembl ``Rattus_norvegicus.mRatBN7.2.cdna.all.fa``: 46,856
    records, region type ``primary_assembly`` on every single one, zero
    ``chromosome`` and zero ``scaffold``. Treating that as non-primary made the
    ALT/PATCH warning fire against a file with no ALT or PATCH record in it.
    """
    assert RQ.classify_header(RAT_PRIMARY_HDR) == (RQ.PRIMARY, "TR_V_gene")


def test_a_rat_unplaced_scaffold_is_still_not_primary():
    """The NAME check survives the region-type widening.

    Ensembl labels rat's unplaced scaffolds ``primary_assembly`` too, so region
    type alone cannot decide. 69 of the file's 46,856 records are these.
    """
    assembly, biotype = RQ.classify_header(RAT_UNPLACED_HDR)
    assert (assembly, biotype) == (RQ.UNPLACED, "protein_coding")
    assert RQ.seq_region_of(RAT_UNPLACED_HDR) == "JACYVU010000238.1"


# ---------------------------------------------------------------------------
# (1b) ALT/PATCH vs unplaced: kept separate
# ---------------------------------------------------------------------------
#
# They must not share one NON_PRIMARY class that is dropped wholesale. They
# are opposite kinds of record and the errors they cause are opposite:
#
#   ALT/PATCH  a near-verbatim COPY of primary sequence. Leaving it in
#              manufactures a perfect off-target for a probe against the real
#              gene, so a GOOD PROBE IS REJECTED (false negative).
#   unplaced   unique sequence that could not be placed on a chromosome. It
#              duplicates nothing. Removing it hides a real off-target, so a
#              BAD PROBE PASSES (false positive).
#
# Measured on Ensembl's public Homo_sapiens.GRCh38.cdna.all.fa, and this is
# the measurement the split rests on: of the 19,713 ALT/PATCH records, the 27
# HG*_NOVEL_TEST ones alone put 23 second copies of a chromosomal gene into the
# reference; of the 35 unplaced records, 0 duplicate a primary gene symbol.

@pytest.mark.parametrize("name", [
    "HSCHR6_MHC_APD_CTG1",       # ALT haplotype
    "HSCHR13_1_CTG5",
    "HG1524_PATCH",              # GRC FIX patch
    "HG2290_PATCH",
    "HG142_HG150_NOVEL_TEST",    # GRC NOVEL patch, no "PATCH" in the name
    "HG151_NOVEL_TEST",
    "CHR_HSCHR1_1_CTG3_alt",     # Ensembl's other ALT spelling
])
def test_alt_and_patch_names_are_recognised(name):
    assert RQ._is_alt_or_patch(name), name
    assert RQ.classify_header(_hdr("T1", name, "protein_coding",
                                   region_type="scaffold"))[0] == RQ.ALT_PATCH


@pytest.mark.parametrize("name", [
    "KI270728.1",            # human unplaced
    "GL000009.2",            # human unlocalised
    "JH584295.1",            # mouse
    "GL456210.1",
    "JACYVU010000238.1",     # rat
    "MU150189.1",
    "HG",                    # not HG<digits>_
    "HGSOMETHING",
    "SOMETHING_WEIRD",       # unrecognised: retained, not dropped
])
def test_unplaced_and_unrecognised_names_are_not_alt_or_patch(name):
    """The conservative default: only a RECOGNISED alternate is excluded."""
    assert not RQ._is_alt_or_patch(name), name


def test_a_pseudogene_on_a_primary_chromosome_is_retained(tmp_path):
    """List item 1."""
    fa = tmp_path / "p1.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(300)}\n"
        f">{_hdr('ENST02', '2', 'processed_pseudogene')}\n{_seq(200, seed=4)}\n"
        f">{_hdr('ENST03', 'HG1524_PATCH', 'protein_coding')}\n{_seq(150, seed=5)}\n"
    )
    ref = K.load_reference(str(fa), exclude_non_primary=True)
    assert any("processed_pseudogene" in n for n in ref.names)


def test_a_pseudogene_on_an_unplaced_scaffold_is_retained(tmp_path):
    """List item 2, applied to unplaced records.

    The module's stated position is that pseudogenes are reported and NEVER
    excluded, because a transcribed pseudogene absent from the reference is an
    off-target the screen cannot see. An exclusion that ignores biotype would
    drop a pseudogene on a non-primary record anyway. Ensembl's rat cDNA has 5
    such records.
    """
    fa = tmp_path / "p2.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(300)}\n"
        f">{_hdr('ENST02', 'HG1524_PATCH', 'protein_coding')}\n{_seq(150, seed=5)}\n"
        f">{_hdr('ENST03', 'JACYVU010000238.1', 'processed_pseudogene')}\n"
        f"{_seq(200, seed=6)}\n"
    )
    ref = K.load_reference(str(fa), exclude_non_primary=True)
    assert ref.composition.excluded_non_primary is True
    kept = " ".join(ref.names)
    assert "JACYVU010000238.1" in kept
    assert "HG1524_PATCH" not in kept


def test_a_protein_coding_transcript_on_an_unplaced_scaffold_is_retained(tmp_path):
    """List item 3. 64 of rat's 69 unplaced records are protein_coding."""
    fa = tmp_path / "p3.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(300)}\n"
        f">{_hdr('ENST02', 'HG1524_PATCH', 'protein_coding')}\n{_seq(150, seed=5)}\n"
        f">{_hdr('ENST03', 'MU150189.1', 'protein_coding')}\n{_seq(200, seed=7)}\n"
    )
    ref = K.load_reference(str(fa), exclude_non_primary=True)
    assert len(ref.names) == 2
    assert any("MU150189.1" in n for n in ref.names)
    assert ref.composition.records[RQ.UNPLACED] == 1
    assert ref.composition.records[RQ.ALT_PATCH] == 1


def test_an_unplaced_scaffold_is_actually_searchable_not_merely_counted(tmp_path):
    """Counting a record is not keeping it. Search for sequence only it has."""
    unique = _seq(300, seed=21)
    fa = tmp_path / "p4.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(300, seed=22)}\n"
        f">{_hdr('ENST02', 'HG1524_PATCH', 'protein_coding')}\n{_seq(300, seed=23)}\n"
        f">{_hdr('ENST03', 'GL000009.2', 'protein_coding')}\n{unique}\n"
    )
    ref = K.load_reference(str(fa), exclude_non_primary=True)
    idx = K.KmerIndex(ref.seq, 12)
    pos, mm = K.search(ref.seq, idx, K.encode_seq(unique[100:150]), 0,
                       (0, 12, 25, 37), 2)
    assert len(pos) == 1 and int(mm[0]) == 0
    assert "GL000009.2" in ref.names[int(ref.record_of(pos)[0])]


def test_a_reference_with_no_alt_or_patch_excludes_nothing(tmp_path):
    """Mouse and rat: many primary records, a few unplaced, zero ALT/PATCH.

    Nothing to exclude means the exclusion must report itself as not having
    run, so the run reuses the full index rather than building a second,
    identical one under the primary-only cache name.
    """
    fa = tmp_path / "mouse_like.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(300)}\n"
        f">{_hdr('ENST02', 'GL456210.1', 'protein_coding')}\n{_seq(120, seed=8)}\n"
    )
    ref = K.load_reference(str(fa), exclude_non_primary=True)
    assert len(ref.names) == 2
    assert ref.composition.excluded_non_primary is False
    assert ref.composition.records[RQ.ALT_PATCH] == 0
    assert ref.composition.records[RQ.UNPLACED] == 1


# ---------------------------------------------------------------------------
# (1c) a reference-build change must not silently reclassify the whole file
# ---------------------------------------------------------------------------

def test_a_reference_with_ensembl_headers_but_no_primary_shouts(tmp_path, caplog):
    """List item 9, in the shape the rat bug actually took.

    46,856 rat records carried parseable Ensembl location fields and NONE
    classified primary, because only the region type "chromosome" was
    recognised. Nothing failed and nothing was dropped. A whole reference
    changing class is exactly what a reference-build change causes, and it has
    to be loud.
    """
    fa = tmp_path / "unrecognised.fa"
    h1 = _hdr("ENST01", "JACYVU010000238.1", "protein_coding",
              region_type="supercontig")
    h2 = _hdr("ENST02", "MU150189.1", "protein_coding",
              region_type="supercontig")
    fa.write_text(
        f">{h1}\n{_seq(300)}\n"
        f">{h2}\n{_seq(300, seed=9)}\n"
    )
    comp = RQ.scan_fasta(str(fa))
    assert comp.records[RQ.PRIMARY] == 0
    assert comp.ensembl_style is True
    with caplog.at_level(logging.WARNING):
        msg = RQ.warn_if_no_primary(comp)
    assert msg is not None
    assert "NONE classified as primary" in msg
    assert "supercontig" in msg          # says WHY, not just THAT


def test_a_healthy_reference_does_not_shout_about_primary(tmp_path):
    """The control: the same guard must stay silent on a normal reference."""
    fa = tmp_path / "healthy.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(300)}\n"
        f">{_hdr('ENST02', 'HG1524_PATCH', 'protein_coding')}\n{_seq(200, seed=10)}\n"
    )
    assert RQ.warn_if_no_primary(RQ.scan_fasta(str(fa))) is None


def test_a_viral_reference_never_shouts_about_primary(tmp_path):
    """A viral panel has no location fields at all, so the guard must not fire."""
    fa = tmp_path / "orov.fa"
    fa.write_text(">OROV_L_BeAn19991\n" + _seq(300) + "\n")
    assert RQ.warn_if_no_primary(RQ.scan_fasta(str(fa))) is None


# ---------------------------------------------------------------------------
# (1d) the report says how much went, and why
# ---------------------------------------------------------------------------

def test_the_manifest_reports_each_exclusion_reason_separately(tmp_path):
    """List item 10.

    A single "non_primary" number cannot distinguish "duplicate copies removed"
    from "unique scaffolds removed", which is the whole reason those two were
    separated.
    """
    fa = tmp_path / "reasons.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(300)}\n"
        f">{_hdr('ENST02', 'HG1524_PATCH', 'protein_coding')}\n{_seq(200, seed=11)}\n"
        f">{_hdr('ENST03', 'HSCHR6_MHC_APD_CTG1', 'lncRNA')}\n{_seq(100, seed=12)}\n"
        f">{_hdr('ENST04', 'GL000009.2', 'protein_coding')}\n{_seq(50, seed=13)}\n"
        f">unclassifiable_viral_record\n{_seq(25, seed=14)}\n"
    )
    man = K.load_reference(str(fa), exclude_non_primary=True) \
        .composition.to_manifest()
    assert man["alt_patch_records"] == 2
    assert man["alt_patch_bases"] == 300
    assert man["alt_patch_seq_regions"] == 2
    assert man["unplaced_records"] == 1
    assert man["unplaced_bases"] == 50
    assert man["unplaced_seq_regions"] == 1
    assert man["primary_records"] == 1
    assert man["unclassified_records"] == 1
    # The old union stays reportable so runs either side of the split compare.
    assert man["non_primary_records"] == 3
    assert man["non_primary_bases"] == 350
    # And the classes are named, so the record does not depend on knowing the
    # release's rule.
    assert man["excluded_classes"] == ["alt_patch"]
    assert set(man["retained_classes"]) == {"primary", "unplaced", "unknown"}
    # Kept + excluded accounts for every record and every base.
    assert man["records_searched"] + man["alt_patch_records"] == man["records"]
    assert man["bases_searched"] + man["alt_patch_bases"] == man["bases"]
    json.dumps(man)


def test_the_two_cache_variants_cannot_collide(tmp_path):
    """A ".primary" index built before the split held FEWER records than a
    primary-only index does now, so the name had to change with the meaning."""
    prefix = K.KmerIndex.cache_prefix(str(tmp_path / "ref.fa"), 12,
                                      cache_dir=str(tmp_path),
                                      primary_only=True)
    assert prefix.endswith("ref.noaltpatch.k12")
    full = K.KmerIndex.cache_prefix(str(tmp_path / "ref.fa"), 12,
                                    cache_dir=str(tmp_path),
                                    primary_only=False)
    assert full.endswith("ref.k12")
    assert prefix != full


def test_a_rat_reference_does_not_claim_it_is_searching_alt_scaffolds(tmp_path):
    """The user-visible symptom, end to end.

    A reference in rat's proportions (many primary_assembly chromosome records,
    a few unplaced scaffolds) must not be told that 100% of its bases are
    ALT/PATCH copies being searched as off-target subjects.
    """
    fa = tmp_path / "rat_like.fa"
    chr20 = _hdr("ENSRNOT02", "20", "protein_coding",
                 region_type="primary_assembly")
    fa.write_text(
        f">{RAT_PRIMARY_HDR}\n{_seq(300, seed=11)}\n"
        f">{chr20}\n{_seq(300, seed=12)}\n"
        f">{RAT_UNPLACED_HDR}\n{_seq(20, seed=13)}\n"
    )
    comp = RQ.scan_fasta(str(fa))
    assert comp.records[RQ.PRIMARY] == 2
    assert comp.records[RQ.UNPLACED] == 1
    assert comp.records[RQ.ALT_PATCH] == 0
    assert comp.records[RQ.UNKNOWN] == 0
    comp.finalise(excluded_non_primary=RQ.should_exclude(comp, True))
    # Rat's cdna.all carries no ALT or PATCH record at all, so there is nothing
    # to exclude and the unplaced scaffolds stay searchable.
    assert comp.excluded_non_primary is False
    assert RQ.warn_if_searching_alt(comp) is None
    assert RQ.warn_if_no_primary(comp) is None


@pytest.mark.parametrize("header", [
    "OROV_L_BeAn19991",
    "NC_045512.2 Severe acute respiratory syndrome coronavirus 2",
    "gene1|transcript1|something:else",
    "MN_123 description:has a colon: in it",
    "",
])
def test_non_ensembl_headers_degrade_to_unknown(header):
    """A viral or user FASTA has none of these fields, and that is normal."""
    assert RQ.classify_header(header) == (RQ.UNKNOWN, RQ.BIOTYPE_UNKNOWN)
    assert RQ.seq_region_of(header) == ""


def test_biotype_is_read_even_without_a_location_field():
    assembly, biotype = RQ.classify_header("T1 gene_biotype:lncRNA")
    assert (assembly, biotype) == (RQ.UNKNOWN, "lncRNA")


def test_pseudogene_matching_covers_the_spellings_ensembl_uses():
    for bio in ("processed_pseudogene", "unprocessed_pseudogene",
                "transcribed_unitary_pseudogene", "IG_V_pseudogene",
                "rRNA_pseudogene", "polymorphic_pseudogene"):
        assert RQ.is_pseudogene(bio), bio
    for bio in ("protein_coding", "lncRNA", "miRNA", RQ.BIOTYPE_UNKNOWN):
        assert not RQ.is_pseudogene(bio), bio


# ---------------------------------------------------------------------------
# (2) exclusion actually removes records from the searchable reference
# ---------------------------------------------------------------------------

def _mixed_fasta(tmp_path):
    """Two primary records and one PATCH record that COPIES the first.

    The copy is the whole point: an ALT/PATCH record is a near-verbatim
    duplicate of a primary region, so leaving it in manufactures a perfect
    off-target for a probe against the real gene.
    """
    gene = _seq(300, seed=1)
    other = _seq(300, seed=2)
    fa = tmp_path / "mixed.fa"
    fa.write_text(
        f">{_hdr('ENST01', '13', 'protein_coding')}\n{gene}\n"
        f">{_hdr('ENST02', '7', 'lncRNA')}\n{other}\n"
        f">{_hdr('ENST03', 'HG1524_PATCH', 'protein_coding')}\n{gene}\n"
    )
    return str(fa), gene


def test_exclusion_drops_the_patch_record_from_the_reference(tmp_path):
    fa, _gene = _mixed_fasta(tmp_path)
    ref = K.load_reference(fa, exclude_non_primary=True)
    assert len(ref.names) == 2
    assert not any("PATCH" in n for n in ref.names)
    assert ref.composition is not None
    assert ref.composition.excluded_non_primary is True
    # The tally still describes the FILE, not the surviving subset — otherwise
    # the manifest could not say what was left out.
    assert ref.composition.total_records == 3
    assert ref.composition.records[RQ.NON_PRIMARY] == 1
    assert ref.composition.records_searched == 2


def test_the_dropped_record_is_no_longer_findable(tmp_path):
    """The point of the exclusion: one copy of the gene, not two."""
    fa, gene = _mixed_fasta(tmp_path)
    query = K.encode_seq(gene[100:150])

    kept = K.load_reference(fa, exclude_non_primary=False)
    idx = K.KmerIndex(kept.seq, 12)
    pos, mm = K.search(kept.seq, idx, query, 0, (0, 12, 25, 37), 2)
    assert len(pos) == 2 and set(mm.tolist()) == {0}

    pruned = K.load_reference(fa, exclude_non_primary=True)
    idx2 = K.KmerIndex(pruned.seq, 12)
    pos2, _ = K.search(pruned.seq, idx2, query, 0, (0, 12, 25, 37), 2)
    assert len(pos2) == 1
    assert "PATCH" not in pruned.names[int(pruned.record_of(pos2)[0])]


def test_record_boundaries_stay_correct_after_a_drop(tmp_path):
    """Positions must map to the record they are really in, not to an offset
    computed from records that are no longer there."""
    fa, _ = _mixed_fasta(tmp_path)
    ref = K.load_reference(fa, exclude_non_primary=True)
    for i, name in enumerate(ref.names):
        assert ref.name_of(int(ref.starts[i])) == name
    assert len(ref.starts) == len(ref.names)


def test_exclusion_can_be_disabled(tmp_path):
    fa, _ = _mixed_fasta(tmp_path)
    ref = K.load_reference(fa, exclude_non_primary=False)
    assert len(ref.names) == 3
    assert ref.composition is not None
    assert ref.composition.excluded_non_primary is False


# ---------------------------------------------------------------------------
# (3) graceful degradation — the case most runs are in
# ---------------------------------------------------------------------------

def test_a_viral_reference_is_untouched_and_silent(tmp_path, caplog):
    fa = tmp_path / "orov.fa"
    fa.write_text(">OROV_L_BeAn19991\n" + _seq(300) + "\n"
                  ">OROV_M_BeAn19991\n" + _seq(300, seed=3) + "\n")
    with caplog.at_level(logging.WARNING):
        ref = K.load_reference(str(fa), exclude_non_primary=True)
    assert len(ref.names) == 2
    assert ref.composition is not None
    assert ref.composition.ensembl_style is False
    assert ref.composition.excluded_non_primary is False
    assert ref.composition.records[RQ.UNKNOWN] == 2
    assert [r for r in caplog.records if r.levelno >= logging.WARNING] == []


def test_a_reference_of_only_alt_records_is_never_emptied(tmp_path):
    """"Has Ensembl headers" is not enough to justify dropping records — the
    file must also have a primary assembly to fall back on."""
    fa = tmp_path / "alt_only.fa"
    fa.write_text(
        f">{_hdr('ENST01', 'HG1524_PATCH', 'protein_coding')}\n{_seq(200)}\n"
        f">{_hdr('ENST02', 'HSCHR6_MHC_APD_CTG1', 'lncRNA')}\n{_seq(200)}\n"
    )
    ref = K.load_reference(str(fa), exclude_non_primary=True)
    assert len(ref.names) == 2
    assert ref.composition.excluded_non_primary is False


# ---------------------------------------------------------------------------
# (4) the tally itself
# ---------------------------------------------------------------------------

def test_biotype_tally_counts_records_and_bases(tmp_path):
    fa = tmp_path / "biotypes.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(300)}\n"
        f">{_hdr('ENST02', '2', 'protein_coding')}\n{_seq(150)}\n"
        f">{_hdr('ENST03', '3', 'processed_pseudogene')}\n{_seq(100)}\n"
        f">{_hdr('ENST04', '4', 'lncRNA')}\n{_seq(50)}\n"
        f">{_hdr('ENST05', 'HG1_PATCH', 'protein_coding')}\n{_seq(400)}\n"
    )
    ref = K.load_reference(str(fa), exclude_non_primary=True)
    comp = ref.composition

    assert comp.total_records == 5
    assert comp.total_bases == 300 + 150 + 100 + 50 + 400
    assert comp.records[RQ.PRIMARY] == 4
    assert comp.bases[RQ.PRIMARY] == 600
    assert comp.records[RQ.NON_PRIMARY] == 1
    assert comp.bases[RQ.NON_PRIMARY] == 400
    assert comp.non_primary_base_fraction == pytest.approx(400 / 1000)

    # The biotype view describes what the screen SEARCHED, so the excluded
    # PATCH record's 400 protein_coding bases are not in it.
    assert comp.biotype_bases() == {
        "protein_coding": 450, "processed_pseudogene": 100, "lncRNA": 50,
    }
    assert comp.biotype_records() == {
        "protein_coding": 2, "processed_pseudogene": 1, "lncRNA": 1,
    }
    assert comp.pseudogene_bases() == 100
    assert comp.bases_searched == 600


def test_biotype_view_widens_when_nothing_is_excluded(tmp_path):
    fa = tmp_path / "biotypes.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(300)}\n"
        f">{_hdr('ENST05', 'HG1_PATCH', 'protein_coding')}\n{_seq(400)}\n"
    )
    comp = K.load_reference(str(fa), exclude_non_primary=False).composition
    assert comp.biotype_bases() == {"protein_coding": 700}
    assert comp.bases_searched == 700


def test_pseudogenes_are_reported_but_never_excluded(tmp_path):
    """The design position here is the opposite of the ALT/PATCH one: a
    transcribed pseudogene missing from the reference is an off-target the
    screen cannot see, i.e. a bad probe that passes."""
    fa = tmp_path / "pseudo.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(300)}\n"
        f">{_hdr('ENST02', '2', 'processed_pseudogene')}\n{_seq(200)}\n"
    )
    ref = K.load_reference(str(fa), exclude_non_primary=True)
    assert len(ref.names) == 2
    man = ref.composition.to_manifest()
    assert man["pseudogene_bases"] == 200
    assert man["pseudogene_base_fraction_of_searched"] == pytest.approx(0.4)


def test_manifest_payload_is_json_serialisable_and_complete(tmp_path):
    fa, _ = _mixed_fasta(tmp_path)
    man = K.load_reference(fa).composition.to_manifest()
    round_tripped = json.loads(json.dumps(man))
    for key in ("records", "bases", "primary_records", "primary_bases",
                "non_primary_records", "non_primary_bases",
                "non_primary_base_fraction", "non_primary_seq_regions",
                "unclassified_records", "unclassified_bases",
                "ensembl_style_headers", "non_primary_excluded",
                "records_searched", "bases_searched", "pseudogene_bases",
                "gene_biotype_bases", "gene_biotype_records"):
        assert key in round_tripped, key
    assert round_tripped["non_primary_seq_regions"] == 1


def test_scan_fasta_matches_what_the_loader_measured(tmp_path):
    """The manifest's fallback for a run that never built a k-mer reference
    must not disagree with the loader about what is in the file."""
    fa, _ = _mixed_fasta(tmp_path)
    loaded = K.load_reference(fa, exclude_non_primary=False).composition
    scanned = RQ.scan_fasta(fa)
    assert scanned.total_records == loaded.total_records
    assert scanned.total_bases == loaded.total_bases
    assert scanned.records == loaded.records
    assert scanned.bases == loaded.bases
    assert scanned.non_primary_seq_regions == loaded.non_primary_seq_regions


def test_scan_of_a_missing_file_is_empty_not_fatal(tmp_path):
    comp = RQ.scan_fasta(str(tmp_path / "nope.fa"))
    assert comp.total_records == 0
    assert comp.ensembl_style is False


def test_the_loader_publishes_what_it_measured(tmp_path):
    fa, _ = _mixed_fasta(tmp_path)
    assert RQ.measured(fa) is None
    K.load_reference(fa)
    assert RQ.measured(fa) is not None
    assert RQ.measured(fa).total_records == 3


# ---------------------------------------------------------------------------
# (5) the warning: once per run, only when it means something
# ---------------------------------------------------------------------------

def test_one_warning_when_alt_records_are_being_searched(tmp_path, caplog):
    fa, _ = _mixed_fasta(tmp_path)
    with caplog.at_level(logging.WARNING):
        K.load_reference(fa, exclude_non_primary=False)
    warnings = [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert len(warnings) == 1, "one line per run, not per record"
    msg = warnings[0].getMessage()
    assert "ALT/PATCH" in msg and "good probes rejected" in msg


def test_no_warning_when_the_records_were_excluded(tmp_path, caplog):
    fa, _ = _mixed_fasta(tmp_path)
    with caplog.at_level(logging.WARNING):
        K.load_reference(fa, exclude_non_primary=True)
    assert [r for r in caplog.records if r.levelno >= logging.WARNING] == []


def test_no_warning_below_the_meaningful_fraction(tmp_path, caplog):
    """One unplaced scaffold in a big reference is not actionable news."""
    fa = tmp_path / "trace.fa"
    fa.write_text(
        f">{_hdr('ENST01', '1', 'protein_coding')}\n{_seq(20000)}\n"
        f">{_hdr('ENST02', 'KI270728.1', 'protein_coding')}\n{_seq(50)}\n"
    )
    with caplog.at_level(logging.WARNING):
        comp = K.load_reference(str(fa), exclude_non_primary=False).composition
    assert comp.non_primary_base_fraction < RQ.ALT_PATCH_WARN_BASE_FRACTION
    assert [r for r in caplog.records if r.levelno >= logging.WARNING] == []


# ---------------------------------------------------------------------------
# (6) the cached index must not be shared between the two record sets
# ---------------------------------------------------------------------------

def test_primary_only_index_gets_its_own_cache_name(tmp_path):
    """Dropping records does not change the FASTA's hash, so the NAME is what
    keeps a full-file index from being read back as a primary-only one — which
    would report hits at positions belonging to other records."""
    fa, _ = _mixed_fasta(tmp_path)
    full = K.KmerIndex.cache_prefix(fa, 12, primary_only=False)
    pruned = K.KmerIndex.cache_prefix(fa, 12, primary_only=True)
    assert full != pruned
    assert pruned.endswith(".noaltpatch.k12")


def test_build_or_load_keeps_the_two_variants_apart(tmp_path):
    fa, _ = _mixed_fasta(tmp_path)
    idx_pruned, ref_pruned = K.KmerIndex.build_or_load(
        fa, 12, exclude_non_primary=True)
    idx_full, ref_full = K.KmerIndex.build_or_load(
        fa, 12, exclude_non_primary=False)
    assert len(ref_pruned.names) == 2 and len(ref_full.names) == 3
    assert len(idx_pruned.pos) < len(idx_full.pos)

    # And each is reused rather than rebuilt: a cache hit returns an index whose
    # positions match the reference it was built for.
    again, ref_again = K.KmerIndex.build_or_load(
        fa, 12, exclude_non_primary=True)
    assert len(again.pos) == len(idx_pruned.pos)
    assert ref_again.names == ref_pruned.names
    assert np.array_equal(np.asarray(again.pos), np.asarray(idx_pruned.pos))


def test_staging_the_wrong_index_variant_warns_and_names_both(tmp_path, caplog):
    """A miss caused by the WRONG VARIANT must say so, at WARNING.

    Shipping a prebuilt index beside a reference is a supported workflow, and
    the easy mistake is to stage the full-file index for a run that -- because
    ``--exclude-alt-scaffolds`` is on by default -- wants the primary-only one.
    That mistake turns a seconds-long run into a full index rebuild, and the
    operator must be able to read the cause off the log without diffing
    filenames by eye.
    """
    fa, _ = _mixed_fasta(tmp_path)
    # Stage only the full-file variant, then ask for the default (primary-only).
    K.KmerIndex.build_or_load(fa, 12, exclude_non_primary=False)
    wanted = K.KmerIndex.cache_prefix(fa, 12, primary_only=True)
    assert not os.path.exists(wanted + ".meta.json")

    with caplog.at_level(logging.WARNING, logger="flexprobe.kmer_index"):
        K.KmerIndex.build_or_load(fa, 12, exclude_non_primary=True)

    # Filter by logger: reference_qc also warns here (the staged full-file
    # index legitimately searches ALT/PATCH records), and that is not this test.
    warnings = [r.getMessage() for r in caplog.records
                if r.levelno >= logging.WARNING
                and r.name == "flexprobe.kmer_index"]
    assert warnings, "a wrong-variant miss must warn, not rebuild silently"
    msg = warnings[0]
    assert "MISS" in msg
    # Both prefixes named, so the fix is "rename/restage this file", not a hunt.
    assert os.path.basename(wanted) in msg
    assert os.path.basename(
        K.KmerIndex.cache_prefix(fa, 12, primary_only=False)) in msg
    # And it must say what it costs, or the reader will not act on it.
    assert "7 GB" in msg


def test_a_changed_reference_misses_quietly_and_says_which_field(tmp_path,
                                                                 caplog):
    """The ordinary miss stays at INFO. Only the wrong-variant case is loud.

    A reference that genuinely changed SHOULD rebuild, and warning about it
    every time would train the operator to ignore the warning that matters.
    """
    fa, _ = _mixed_fasta(tmp_path)
    K.KmerIndex.build_or_load(fa, 12, exclude_non_primary=True)
    Path(fa).write_text(Path(fa).read_text() + ">extra\nACGTACGTACGTACGT\n")

    with caplog.at_level(logging.INFO, logger="flexprobe.kmer_index"):
        K.KmerIndex.build_or_load(fa, 12, exclude_non_primary=True)

    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
    infos = " ".join(r.getMessage() for r in caplog.records)
    assert "reference has changed" in infos


# ---------------------------------------------------------------------------
# (7) the manifest — the reference is characterised, not merely identified
# ---------------------------------------------------------------------------

def _key(off, tmp_path):
    """Manifest blocks are keyed by the path relative to the output dir."""
    import os
    return os.path.relpath(str(off), str(tmp_path / "out"))


def _design_run(tmp_path, off_fasta, **screen_over):
    """A real k-mer pipeline run against ``off_fasta``. No blastn needed."""
    from flexprobe.models import (
        ExportSettings, ModalityDefinition, RunConfig, ScreenResources,
    )
    from flexprobe.pipeline import run_pipeline
    from flexprobe.screen import ENGINE_KMER

    fa = tmp_path / "target.fa"
    fa.write_text(f">GENE1\n{_seq(600, seed=11)}\n")
    out = tmp_path / "out"
    run_pipeline(RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(db_paths=[str(off_fasta)], engine=ENGINE_KMER,
                               **screen_over),
        export=ExportSettings(overwrite=True), discovery_only=True,
    ))
    return json.loads((out / "manifest.json").read_text())


def test_the_manifest_records_what_the_reference_contained(tmp_path):
    """``db_sha256s`` says WHICH reference; this says what was IN it. Two runs
    can cite the same hash and still have incomparable rejection rates."""
    off = tmp_path / "off.fa"
    off.write_text(
        f">{_hdr('ENST01', '13', 'protein_coding')}\n{_seq(3000, seed=21)}\n"
        f">{_hdr('ENST02', '2', 'processed_pseudogene')}\n{_seq(500, seed=22)}\n"
        f">{_hdr('ENST03', 'HG1524_PATCH', 'protein_coding')}\n"
        f"{_seq(3000, seed=21)}\n"
    )
    manifest = _design_run(tmp_path, off)

    comp = manifest["reference_composition"][_key(off, tmp_path)]
    assert comp["records"] == 3
    assert comp["non_primary_records"] == 1
    assert comp["non_primary_excluded"] is True
    assert comp["records_searched"] == 2
    assert comp["ensembl_style_headers"] is True
    # Pseudogenes stay in the searched set and are visible in the breakdown.
    assert comp["gene_biotype_records"]["processed_pseudogene"] == 1
    assert comp["pseudogene_bases"] == 500
    # Composition sits ALONGSIDE identity. A k-mer run reads the FASTA
    # directly and has no BLAST index files, so `db_sha256s` digests the FASTA
    # itself. It must never be {} for a k-mer run, or the manifest carries no
    # reference identity.
    import hashlib
    assert len(comp["sha256"]) == 64
    assert manifest["db_sha256s"] == {
        _key(off, tmp_path): hashlib.sha256(off.read_bytes()).hexdigest()}
    assert manifest["db_sha256s"][_key(off, tmp_path)] == comp["sha256"]


def test_a_viral_reference_is_recorded_as_unclassified(tmp_path):
    """The common case: recorded honestly, with nothing excluded."""
    off = tmp_path / "orov.fa"
    off.write_text(f">OROV_L\n{_seq(3000, seed=31)}\n")
    manifest = _design_run(tmp_path, off)

    comp = manifest["reference_composition"][_key(off, tmp_path)]
    assert comp["ensembl_style_headers"] is False
    assert comp["unclassified_records"] == 1
    assert comp["non_primary_excluded"] is False
    assert comp["records_searched"] == 1


def test_a_skipped_screen_claims_nothing_about_a_reference(tmp_path):
    off = tmp_path / "off.fa"
    off.write_text(f">OROV_L\n{_seq(3000, seed=31)}\n")
    manifest = _design_run(tmp_path, off, skip_screen=True)
    assert manifest["reference_composition"] == {}
