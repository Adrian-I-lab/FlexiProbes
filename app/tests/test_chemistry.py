"""Chemistry layer: the strongest correctness signal in the project.

Verifies LHS/RHS oligo assembly against the EGFP example published in
10x Genomics CG000839 Tables 4 & 5, plus the other modalities' handles.
"""
from flexprobe.chemistry import (
    BARCODE_PLACEHOLDER, CHEMISTRIES, assemble_pair, get_chemistry,
)
from flexprobe.models import ModalityDefinition


# EGFP target window from CG000839: 5' half + 3' half (sense).
EGFP_50 = "GATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACC"


def test_flex_v2_egfp_matches_cg000839():
    ap = assemble_pair(EGFP_50, CHEMISTRIES["flex_v2"], multiplex=True)
    # CG000839 Table 4 (LHS) and Table 5 (RHS), EGFP-LHS-1 / EGFP-RHS-1.
    assert ap.target_lhs == "GGTAGTGGTCGGCGAGCTGCACGCT"
    assert ap.target_rhs == "GCCGTCCTCGATGTTGTGGCGGATC"
    assert ap.lhs_oligo == "CCTTGGCACCCGAGAATTCCA" + ap.target_lhs
    assert ap.rhs_oligo == "/5Phos/" + ap.target_rhs + "CCCATATAAGAAA"


def test_flex_v2_singleplex_uses_pcs1():
    ap = assemble_pair(EGFP_50, CHEMISTRIES["flex_v2"], multiplex=False)
    assert ap.rhs_oligo.endswith("CGGTCCTAGCAA")  # pCS1, 4-Samples Kit


def test_visium_polyA_tail():
    ap = assemble_pair(EGFP_50, CHEMISTRIES["visium"], multiplex=True)
    assert ap.rhs_oligo.endswith("A" * 30)


def test_flex_v1_multiplex_barcode_embedded():
    ap = assemble_pair(EGFP_50, CHEMISTRIES["flex_v1"], multiplex=True,
                       probe_barcode="BC003")
    # Constant Seq + NN + BC003(AGTAGGCT) + pCS1
    assert "ACGCGGTTAGCACGTANN" + "AGTAGGCT" + "CGGTCCTAGCAA" in ap.rhs_oligo


def test_flex_v1_placeholder_kept_for_multibarcode():
    ap = assemble_pair(EGFP_50, CHEMISTRIES["flex_v1"], multiplex=True,
                       leave_barcode_placeholder=True)
    assert BARCODE_PLACEHOLDER in ap.rhs_oligo


def test_tn_junction_flag():
    ok = EGFP_50[:25] + "A" + EGFP_50[26:]
    bad = EGFP_50[:25] + "C" + EGFP_50[26:]
    assert assemble_pair(ok, CHEMISTRIES["flex_v2"], multiplex=True).tn_junction_ok
    assert not assemble_pair(bad, CHEMISTRIES["flex_v2"], multiplex=True).tn_junction_ok


def test_get_chemistry_routing():
    assert get_chemistry(ModalityDefinition(name="flex", version="v2")).key == "flex_v2"
    assert get_chemistry(ModalityDefinition(name="flex", version="v1")).key == "flex_v1"
    assert get_chemistry(ModalityDefinition(spatial=True)).key == "visium"
