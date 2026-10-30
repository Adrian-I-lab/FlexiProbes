"""Greedy auto-pick + selection export (incl. multi-barcode expansion)."""
import csv
import os

import json

import pytest

from flexprobe.auto_select import auto_pick_per_gene
from flexprobe.chemistry import BARCODE_PLACEHOLDER
from flexprobe.errors import WriterError
from flexprobe.select_export import (
    export_project_bulk, export_selected_probes, find_overlaps, save_selections,
)


def test_xlsx_export_write_only_round_trips(tmp_path):
    # Audit #43: results are exported to xlsx via openpyxl write_only (streamed).
    pytest.importorskip("openpyxl")
    import openpyxl
    from flexprobe.models import ProbeCandidate
    from flexprobe.writers.results_table import write_results_table
    c = ProbeCandidate(probe_id="p1", sequence="ACGT" * 12 + "AC", gc_content=0.5)
    out = str(tmp_path / "results.csv")
    write_results_table([c], out, write_xlsx=True, overwrite=True)
    xlsx = out.rsplit(".", 1)[0] + ".xlsx"
    assert os.path.exists(xlsx)
    ws = openpyxl.load_workbook(xlsx).active
    rows = list(ws.iter_rows(values_only=True))
    assert ws.title == "Probes"
    assert rows[0][0] == "rank" and rows[1][1] == "p1"


def test_auto_pick_spaced_top_n():
    # Spacing walks the ladder: 150 bp, then 40, then 15. This gene offers only
    # A2 and A5 at 150 or at 40 (A3 starts 20 bp after A2 ends), which is under
    # the 3-probe quota, so it descends to the 15 bp rung and A3 comes back in.
    rows = [
        {"gene_id": "A", "probe_id": "A1", "start": "0", "end": "50", "composite_score": "0.9"},
        {"gene_id": "A", "probe_id": "A2", "start": "30", "end": "80", "composite_score": "0.95"},
        {"gene_id": "A", "probe_id": "A3", "start": "100", "end": "150", "composite_score": "0.85"},
        {"gene_id": "A", "probe_id": "A4", "start": "60", "end": "110", "composite_score": "0.80"},
        {"gene_id": "A", "probe_id": "A5", "start": "200", "end": "250", "composite_score": "0.75"},
    ]
    sel = auto_pick_per_gene(rows, n_per_gene=3)
    assert sel == {"A": ["A2", "A3", "A5"]}
    # Refusing the tighter rungs (what the GUI spacing filter does) costs the
    # third probe rather than tightening the other two.
    strict = auto_pick_per_gene(rows, n_per_gene=3, relaxed_gap=40, floor_gap=40)
    assert strict == {"A": ["A2", "A5"]}
    # A5 stays reachable at any rung; A1/A4 always lose to the higher-scoring
    # A2 they overlap.
    sel0 = auto_pick_per_gene(rows, n_per_gene=3, min_gap=0, relaxed_gap=0,
                              floor_gap=0)
    assert sel0 == {"A": ["A2", "A3", "A5"]}


def test_find_overlaps_sweep():
    iv = [{"start": 100, "end": 150}, {"start": 140, "end": 190}, {"start": 200, "end": 250}]
    assert sorted(find_overlaps(iv)) == [(0, 1)]
    none = [{"start": 0, "end": 10}, {"start": 20, "end": 30}]
    assert find_overlaps(none) == []
    # With a min-gap, intervals that are merely too close are also flagged.
    assert sorted(find_overlaps(none, min_gap=40)) == [(0, 1)]


def _write_results(out, rows, header_extra=""):
    os.makedirs(out, exist_ok=True)
    cols = ["rank", "probe_id", "gene_id", "start", "end", "mrna_target_50bp",
            "probe_rc_50bp", "target_LHS_25bp", "target_RHS_25bp",
            "lhs_oligo", "rhs_oligo", "composite_score"]
    with open(os.path.join(out, "results.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])


def _manifest(out, params):
    import json
    with open(os.path.join(out, "manifest.json"), "w") as fh:
        json.dump({"params": params}, fh)


def test_export_single_barcode(tmp_path):
    out = str(tmp_path / "run")
    rows = [{"rank": 1, "probe_id": "G_t1", "gene_id": "G", "start": 0, "end": 50,
             "probe_rc_50bp": "A" * 50, "lhs_oligo": "LHS1", "rhs_oligo": "/5Phos/RHS1",
             "composite_score": 0.9}]
    _write_results(out, rows)
    _manifest(out, {"modality": "flex_v2"})
    save_selections(out, {"G": ["G_t1"]})
    paths = export_selected_probes(out, {"G": ["G_t1"]}, overwrite=True)
    with open(paths["idt"]) as fh:
        idt = list(csv.reader(fh))[1:]
    names = [r[1] for r in idt]
    assert names == ["G_t1_LHS", "G_t1_RHS"]
    assert not any("_RHS_RHS" in n for n in names)


def test_export_multi_barcode_expansion(tmp_path):
    out = str(tmp_path / "run")
    rhs = "/5Phos/RHSACGCGGTTAGCACGTANN" + BARCODE_PLACEHOLDER + "CGGTCCTAGCAA"
    rows = [{"rank": 1, "probe_id": "G_t1", "gene_id": "G", "start": 0, "end": 50,
             "probe_rc_50bp": "A" * 50, "lhs_oligo": "LHS1", "rhs_oligo": rhs,
             "composite_score": 0.9}]
    _write_results(out, rows)
    _manifest(out, {"modality": "flex_v1", "probe_barcodes": ["BC001", "BC003"]})
    save_selections(out, {"G": ["G_t1"]})
    paths = export_selected_probes(out, {"G": ["G_t1"]}, overwrite=True)
    with open(paths["idt"]) as fh:
        idt = list(csv.reader(fh))[1:]
    pools = sorted({r[0] for r in idt})
    # Pool name is the barcode/sample name itself (no "_pool" suffix).
    assert pools == ["BC001", "BC003"]
    # No placeholder leaks into the final CSV.
    assert not any(BARCODE_PLACEHOLDER in r[2] for r in idt)
    # BC003 -> AGTAGGCT substituted.
    assert any("AGTAGGCT" in r[2] for r in idt if r[0] == "BC003")


# ---- project bulk export (combine selections across runs) ----

def _make_run(tmp_path, name, gene, params):
    out = str(tmp_path / name)
    rows = [{"rank": 1, "probe_id": f"{gene}_t1", "gene_id": gene,
             "start": 0, "end": 50, "probe_rc_50bp": "A" * 50,
             "lhs_oligo": f"LHS_{name}", "rhs_oligo": f"/5Phos/RHS_{name}",
             "composite_score": 0.9}]
    _write_results(out, rows)
    _manifest(out, params)
    return out


def test_bulk_export_combines_runs_with_run_prefixed_pools(tmp_path):
    out_a = _make_run(tmp_path, "virusA", "G1", {"modality": "flex_v2"})
    out_b = _make_run(tmp_path, "virusB", "G2", {"modality": "flex_v2"})
    dest = str(tmp_path / "bulk")
    res = export_project_bulk(
        [
            {"name": "virusA", "output_dir": out_a, "selections": {"G1": ["G1_t1"]}},
            {"name": "virusB", "output_dir": out_b, "selections": {"G2": ["G2_t1"]}},
        ],
        dest, project_name="myproj", overwrite=True,
    )
    assert res["summary"] == {"virusA": 1, "virusB": 1}
    assert res["skipped"] == []
    with open(res["paths"]["idt"]) as fh:
        idt = list(csv.reader(fh))[1:]
    pools = {r[0] for r in idt}
    # Run name is prefixed so two organisms never share a pool.
    assert pools == {"virusA__G1_pool", "virusB__G2_pool"}
    assert all(r[1].startswith(("virusA__", "virusB__")) for r in idt)
    # Manifest records what went in.
    man = json.load(open(f"{dest}/myproj_bulk_manifest.json"))
    assert man["total_pairs"] == 2 and man["runs_included"] == {"virusA": 1, "virusB": 1}


def test_bulk_export_skips_runs_without_selection(tmp_path):
    out_a = _make_run(tmp_path, "virusA", "G1", {"modality": "flex_v2"})
    dest = str(tmp_path / "bulk")
    res = export_project_bulk(
        [
            {"name": "virusA", "output_dir": out_a, "selections": {"G1": ["G1_t1"]}},
            {"name": "noSel", "output_dir": out_a, "selections": {}},
            {"name": "missing", "output_dir": str(tmp_path / "nope"),
             "selections": {"G": ["G_t1"]}},
        ],
        dest, overwrite=True,
    )
    assert res["summary"] == {"virusA": 1}
    assert set(res["skipped"]) == {"noSel", "missing"}


def test_bulk_export_raises_when_nothing_selected(tmp_path):
    out_a = _make_run(tmp_path, "virusA", "G1", {"modality": "flex_v2"})
    with pytest.raises(WriterError):
        export_project_bulk(
            [{"name": "virusA", "output_dir": out_a, "selections": {}}],
            str(tmp_path / "bulk"), overwrite=True,
        )


def test_bulk_export_multi_barcode_pools_per_sample(tmp_path):
    out = str(tmp_path / "fluA")
    rhs = "/5Phos/RHS" + BARCODE_PLACEHOLDER + "TAIL"
    rows = [{"rank": 1, "probe_id": "G_t1", "gene_id": "G", "start": 0, "end": 50,
             "probe_rc_50bp": "A" * 50, "lhs_oligo": "LHS", "rhs_oligo": rhs,
             "composite_score": 0.9}]
    _write_results(out, rows)
    _manifest(out, {"modality": "flex_v1", "probe_barcodes": ["BC001", "BC003"]})
    res = export_project_bulk(
        [{"name": "fluA", "output_dir": out, "selections": {"G": ["G_t1"]}}],
        str(tmp_path / "bulk"), overwrite=True,
    )
    assert res["summary"] == {"fluA": 2}      # one row per barcode
    with open(res["paths"]["idt"]) as fh:
        idt = list(csv.reader(fh))[1:]
    pools = sorted({r[0] for r in idt})
    assert pools == ["fluA__BC001", "fluA__BC003"]
    assert not any(BARCODE_PLACEHOLDER in r[2] for r in idt)
