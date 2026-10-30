"""Auto-export outcomes must reach the run record honestly.

A discovery run with ``--auto-select-per-gene`` skips the review dialog: the
GUI picks the probes and writes the order CSVs itself. A raised exporter or an
empty pick must not be presented to the user as "auto-exported" when no order
CSVs exist on disk.

The failures here are engineered without monkeypatching: a *directory* sitting
where ``idt_opools.csv`` must be written makes the real exporter raise on the
real code path, so the test cannot drift away from what production does.
"""
import csv
import json
import os

import pytest

from flexprobe_gui.run_history import RunHistory, RunRecord
from flexprobe_gui.run_manager import RunManager


RESULT_COLS = ["rank", "probe_id", "gene_id", "start", "end", "mrna_target_50bp",
               "probe_rc_50bp", "target_LHS_25bp", "target_RHS_25bp",
               "lhs_oligo", "rhs_oligo", "composite_score"]


def _run_dir(tmp_path, name="run", *, rows=None):
    """A finished discovery run: results.csv + manifest, nothing exported yet."""
    out = str(tmp_path / name)
    os.makedirs(out, exist_ok=True)
    if rows is None:
        rows = [
            {"rank": 1, "probe_id": "G_t1", "gene_id": "G", "start": 0, "end": 50,
             "probe_rc_50bp": "A" * 50, "lhs_oligo": "LHS1",
             "rhs_oligo": "/5Phos/RHS1", "composite_score": 0.9},
            {"rank": 2, "probe_id": "G_t2", "gene_id": "G", "start": 900,
             "end": 950, "probe_rc_50bp": "C" * 50, "lhs_oligo": "LHS2",
             "rhs_oligo": "/5Phos/RHS2", "composite_score": 0.8},
        ]
    with open(os.path.join(out, "results.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(RESULT_COLS)
        for r in rows:
            w.writerow([r.get(c, "") for c in RESULT_COLS])
    with open(os.path.join(out, "manifest.json"), "w") as fh:
        json.dump({"params": {"modality": "flex_v2"}}, fh)
    return out


def _manager_with_run(tmp_path, out):
    m = RunManager(RunHistory(path=str(tmp_path / "runs.json")))
    rec = RunRecord(name="auto", output_dir=out, status="running")
    m.history.append(rec)
    return m, rec


def _finish(m, rec, n_per_gene=2):
    m._on_finished_ok(rec.run_id, {
        "stats": {"discovery_only": True, "auto_select_per_gene": n_per_gene},
    })
    return m.history.get(rec.run_id)


# --- a failed export must not be reported as a finished, exported run ----

def test_failed_auto_export_is_not_reported_complete(tmp_path):
    out = _run_dir(tmp_path)
    # A directory where the IDT order sheet has to go: the atomic rename at the
    # end of export_selected_probes cannot land, so the exporter raises for real.
    os.makedirs(os.path.join(out, "idt_opools.csv"))
    m, rec = _manager_with_run(tmp_path, out)

    got = _finish(m, rec)

    assert got.status == "ready_for_review", (
        "a run whose order-CSV export raised was presented to the user as "
        "complete and exported"
    )
    assert got.status != "complete"
    assert "failed" in got.stage
    # The pick itself succeeded and is on disk, so the review dialog has
    # something to offer — that is why this is recoverable, not fatal.
    assert os.path.exists(os.path.join(out, "selections.json"))


def test_auto_export_with_nothing_selected_is_not_reported_complete(tmp_path):
    # results.csv exists but holds no usable probe (zero-length window), so the
    # picker returns nothing. Nothing was exported; nothing may claim it was.
    out = _run_dir(tmp_path, rows=[
        {"rank": 1, "probe_id": "G_t1", "gene_id": "G", "start": 0, "end": 0,
         "probe_rc_50bp": "A" * 50, "lhs_oligo": "L", "rhs_oligo": "R",
         "composite_score": 0.9},
    ])
    m, rec = _manager_with_run(tmp_path, out)

    got = _finish(m, rec)

    assert got.status == "ready_for_review"
    assert got.status != "complete"
    assert not os.path.exists(os.path.join(out, "idt_opools.csv"))


def test_missing_results_csv_is_not_reported_complete(tmp_path):
    out = str(tmp_path / "empty")
    os.makedirs(out)
    m, rec = _manager_with_run(tmp_path, out)

    got = _finish(m, rec)

    assert got.status == "ready_for_review"
    assert got.status != "complete"


def test_successful_auto_export_still_completes(tmp_path):
    out = _run_dir(tmp_path)
    m, rec = _manager_with_run(tmp_path, out)

    got = _finish(m, rec)

    assert got.status == "complete"
    assert got.progress == 100
    assert "auto-exported" in got.stage
    assert os.path.exists(os.path.join(out, "idt_opools.csv"))
    assert os.path.exists(os.path.join(out, "cellranger_probe_set.csv"))


def test_manual_review_run_is_untouched(tmp_path):
    """Discovery without auto-select still stops at ready_for_review."""
    out = _run_dir(tmp_path)
    m, rec = _manager_with_run(tmp_path, out)
    m._on_finished_ok(rec.run_id, {
        "stats": {"discovery_only": True, "auto_select_per_gene": 0}})
    got = m.history.get(rec.run_id)
    assert got.status == "ready_for_review"
    assert not os.path.exists(os.path.join(out, "selections.json"))
# --- one implementation, so the CLI and the GUI cannot drift --------------

def test_neither_caller_reimplements_the_selection_rule():
    """Two copies of the selection rule drift apart. Only a test that notices
    a caller picking probes itself keeps there being one."""
    import inspect

    import flexprobe.cli as cli
    import flexprobe_gui.run_manager as rm

    for label, fn in (("cli._auto_export", cli._auto_export),
                      ("RunManager._auto_export", rm.RunManager._auto_export)):
        src = inspect.getsource(fn)
        assert "auto_pick_with_report" not in src, (
            f"{label} calls the picker directly again — the selection rule "
            f"belongs to flexprobe.auto_export.auto_export only"
        )
        assert "auto_export(" in src, f"{label} no longer calls the shared path"


def test_cli_and_gui_auto_export_produce_identical_orders(tmp_path):
    """Same run dir, same gaps -> byte-identical order CSVs from both paths."""
    import flexprobe.cli as cli

    # Same basename under different parents: the Cell Ranger header embeds the
    # run folder's name as the panel name, and this test is about the probes.
    out_cli = _run_dir(tmp_path, "a/run")
    out_gui = _run_dir(tmp_path, "b/run")

    cli._auto_export(out_cli, 2)
    m, rec = _manager_with_run(tmp_path, out_gui)
    got = _finish(m, rec)
    assert got.status == "complete"

    for fname in ("idt_opools.csv", "cellranger_probe_set.csv"):
        a = open(os.path.join(out_cli, fname), "rb").read()
        b = open(os.path.join(out_gui, fname), "rb").read()
        assert a == b, f"{fname} differs between the CLI and GUI auto-export"
    sel_a = json.loads(open(os.path.join(out_cli, "selections.json")).read())
    sel_b = json.loads(open(os.path.join(out_gui, "selections.json")).read())
    assert sel_a["selections"] == sel_b["selections"]


def test_shared_auto_export_reports_a_failure_to_the_cli_too(capsys, tmp_path):
    import flexprobe.cli as cli
    out = _run_dir(tmp_path)
    os.makedirs(os.path.join(out, "idt_opools.csv"))
    cli._auto_export(out, 2)
    err = capsys.readouterr().err
    assert "auto-export failed" in err


@pytest.mark.parametrize("status_name", ["EXPORTED", "NO_RESULTS",
                                         "NOTHING_SELECTED", "EXPORT_FAILED"])
def test_auto_export_status_enum_covers_every_outcome(status_name):
    from flexprobe.auto_export import AutoExportStatus
    assert hasattr(AutoExportStatus, status_name)
