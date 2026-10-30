"""End-to-end pipeline runs covering the key correctness guarantees."""
import csv
import json
import os

from flexprobe.models import (
    ScreenResources, ExportSettings, FilterThresholds, ModalityDefinition, RunConfig,
)
from flexprobe.pipeline import run_pipeline

from conftest import SYNTH_200, write_fasta


def _cfg(fa, out, **kw):
    base = dict(
        fasta_path=fa, output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        filters=FilterThresholds(gc_min=0.30, gc_max=0.80, max_n_per_half=2,
                                 max_homopolymer_run=6, require_tn_junction=False),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        tile_length=50, tile_step=5,
    )
    base.update(kw)
    return RunConfig(**base)


def test_pipeline_produces_results_and_specificity_status(tmp_path):
    fa = write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)])
    out = str(tmp_path / "out")
    r = run_pipeline(_cfg(fa, out))
    assert r["stats"]["pairs_output"] > 0
    assert r["stats"]["specificity_screened"] is False
    with open(os.path.join(out, "results.csv")) as fh:
        rows = list(csv.DictReader(fh))
    # specificity_status present and honest about not being screened.
    assert all(x["specificity_status"] == "not_screened" for x in rows)
    needed = {"mrna_target_50bp", "probe_rc_50bp", "target_LHS_25bp",
              "target_RHS_25bp", "lhs_oligo", "rhs_oligo", "specificity_status",
              "tm_target_lhs"}
    assert needed <= set(rows[0])
    # Sub-score columns are exposed for rank auditability (#6/#7).
    assert {"score_gc_proximity", "score_specificity_margin"} <= set(rows[0])


def test_manifest_records_tool_versions_and_stage_timings(tmp_path):
    """Provenance a methods section needs, recorded in the manifest.

    ``tool_versions`` must not be ``{}``, and the runtime must be more than a
    single timestamp. Without per-stage timings, a whole-transcriptome run that
    takes many hours has no recoverable per-stage cost to compare or report.
    """
    fa = write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)])
    out = str(tmp_path / "out")
    run_pipeline(_cfg(fa, out))
    with open(os.path.join(out, "manifest.json")) as fh:
        man = json.load(fh)

    tools = man["tool_versions"]
    assert tools, "tool_versions must not be empty"
    # python is always resolvable; blastn/primer3 depend on the environment, so
    # assert their shape rather than their presence (an absent key honestly
    # means 'not found' — better than a fabricated placeholder).
    assert tools["python"]
    for key in ("blastn", "primer3"):
        if key in tools:
            assert isinstance(tools[key], str) and tools[key]

    elapsed = man["elapsed_seconds"]
    # Every stage this run actually passed through is timed, plus the total.
    assert {"ingest", "tile", "filter", "rank", "chemistry", "export"} <= set(elapsed)
    assert all(v >= 0 for v in elapsed.values())
    # The total covers the whole run, so it is >= the sum of any subset.
    stages = sum(v for k, v in elapsed.items() if k != "total")
    assert elapsed["total"] >= stages - 0.01     # float slack


def test_warns_when_tn_junction_off_for_chemistry_that_wants_it(tmp_path, caplog):
    # Audit #29: Flex v2 chemistry expects the TN junction; running with the
    # filter off must emit a warning so the user knows they relaxed it.
    import logging
    fa = write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)])
    out = str(tmp_path / "out")
    with caplog.at_level(logging.WARNING, logger="flexprobe.pipeline"):
        run_pipeline(_cfg(fa, out))   # _cfg uses require_tn_junction=False
    assert any("require_tn_junction is OFF" in r.message for r in caplog.records)


def test_discovery_defers_order_csvs(tmp_path):
    fa = write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)])
    out = str(tmp_path / "out")
    run_pipeline(_cfg(fa, out, discovery_only=True))
    assert os.path.exists(os.path.join(out, "results.csv"))
    assert not os.path.exists(os.path.join(out, "idt_opools.csv"))


def test_filter_report_written(tmp_path):
    fa = write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)])
    out = str(tmp_path / "out")
    run_pipeline(_cfg(fa, out, filters=FilterThresholds(
        gc_min=0.44, gc_max=0.72, max_n_per_half=0, max_homopolymer_run=4,
        require_tn_junction=True), tile_step=1))
    rep = json.load(open(os.path.join(out, "filter_report.json")))
    assert "counts" in rep and "examples" in rep


def test_multi_strain_writes_consensus_and_strain_columns(tmp_path):
    b = list(SYNTH_200)
    for i in (60, 130):
        b[i] = "C" if b[i] != "C" else "A"
    fa = write_fasta(tmp_path / "OROV_L.fa", [("sA", SYNTH_200), ("sB", "".join(b))])
    out = str(tmp_path / "out")
    r = run_pipeline(_cfg(fa, out, multi_strain=True, consensus_strategy="conserved"))
    assert r["stats"]["multi_strain"] is True
    # Audit #39: assert VALUES, not just column presence.
    # consensus.fa carries a real IUPAC consensus (first record).
    with open(os.path.join(out, "consensus.fa")) as fh:
        cons_lines = [ln.strip() for ln in fh if ln.strip()]
    assert cons_lines[0].startswith(">")
    first_seq = []
    for ln in cons_lines[1:]:
        if ln.startswith(">"):
            break
        first_seq.append(ln)
    consensus_seq = "".join(first_seq)
    assert consensus_seq                                  # non-empty
    assert set(consensus_seq) <= set("ACGTRYSWKMBDHVN")   # valid IUPAC

    with open(os.path.join(out, "results.csv")) as fh:
        rows = list(csv.DictReader(fh))
    assert rows
    assert any(h.startswith("mismatches__") for h in rows[0])
    assert any(h.startswith("est_dTm__") for h in rows[0])
    for row in rows:
        mm = {s: int(row[f"mismatches__{s}"]) for s in ("sA", "sB")}
        assert all(v >= 0 for v in mm.values())
        # summary columns are the true min/max of the per-strain counts
        assert int(row["min_strain_mismatches"]) == min(mm.values())
        assert int(row["max_strain_mismatches"]) == max(mm.values())
        # per-strain ΔTm is a real number; worst_strain_dTm is the minimum
        dtm = {s: float(row[f"est_dTm__{s}"]) for s in ("sA", "sB")}
        assert float(row["worst_strain_dTm"]) == min(dtm.values())


def test_idt_no_double_rhs_suffix(tmp_path):
    fa = write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)])
    out = str(tmp_path / "out")
    run_pipeline(_cfg(fa, out, discovery_only=False,
                      modality=ModalityDefinition(name="flex", version="v1")))
    with open(os.path.join(out, "idt_opools.csv")) as fh:
        names = [r[1] for r in list(csv.reader(fh))[1:]]
    assert names and not any("_RHS_RHS" in n for n in names)


def test_one_shot_flex_v1_multiplex_expands_barcodes(tmp_path):
    # One-shot (non-discovery, e.g. CLI) Flex v1 multiplex must EXPAND the order
    # CSV per barcode — never emit the literal {BC} placeholder.
    fa = write_fasta(tmp_path / "g.fa", [("G", SYNTH_200)])
    out = str(tmp_path / "out")
    run_pipeline(_cfg(
        fa, out, discovery_only=False,
        modality=ModalityDefinition(name="flex", version="v1"),
        multiplex=True, probe_barcodes=["BC001", "BC002", "BC003", "BC004"],
    ))
    with open(os.path.join(out, "idt_opools.csv")) as fh:
        rows = list(csv.reader(fh))[1:]
    pools = {r[0] for r in rows}
    seqs = "".join(r[2] for r in rows)
    assert pools == {"BC001", "BC002", "BC003", "BC004"}
    assert "{BC}" not in seqs                      # placeholder fully substituted
