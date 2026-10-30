"""Run-provenance helpers: stage timing and tool-version capture.

Both back manifest fields a methods section depends on. The end-to-end
assertions live in test_pipeline_e2e.py; these cover the semantics the pipeline
relies on but a full run can't pin down deterministically.
"""
import json

from flexprobe.pipeline import _collect_tool_versions, _StageTimer
from flexprobe.writers.manifest import record_stage_elapsed


# ---- stage timing ----

def test_timer_accumulates_interleaved_stages(monkeypatch):
    """Re-entering a stage adds to its total rather than restarting it.

    The streamed paths interleave — tile/filter alternate per gene, and
    blast/specificity alternate per chunk when specificity is scored inside the
    BLAST loop — so a stage's manifest entry must be its true share of wall
    clock, not just its last contiguous span.
    """
    clock = [100.0]
    monkeypatch.setattr("flexprobe.pipeline.time.monotonic", lambda: clock[0])

    t = _StageTimer()
    t.enter("tile");   clock[0] += 3      # tile: 3
    t.enter("filter"); clock[0] += 1      # filter: 1
    t.enter("tile");   clock[0] += 2      # tile: 3 + 2 = 5
    t.enter("filter"); clock[0] += 4      # filter: 1 + 4 = 5

    snap = t.snapshot()
    assert snap["tile"] == 5
    assert snap["filter"] == 5
    assert snap["total"] == 10


def test_snapshot_banks_the_open_stage_without_ending_it(monkeypatch):
    """The manifest is written from inside the export stage, so the snapshot
    has to include the export time elapsed so far and leave the timer running."""
    clock = [0.0]
    monkeypatch.setattr("flexprobe.pipeline.time.monotonic", lambda: clock[0])

    t = _StageTimer()
    t.enter("export"); clock[0] += 2
    assert t.snapshot()["export"] == 2

    clock[0] += 3                          # still inside export
    later = t.snapshot()
    assert later["export"] == 5            # kept accumulating, not reset
    assert later["total"] == 5


def test_timer_before_any_stage_reports_only_total(monkeypatch):
    clock = [0.0]
    monkeypatch.setattr("flexprobe.pipeline.time.monotonic", lambda: clock[0])
    t = _StageTimer()
    clock[0] += 7
    assert t.snapshot() == {"total": 7}


def test_total_is_wall_clock_not_the_sum_of_stages(monkeypatch):
    """Timer construction to snapshot, including any pre-first-stage setup."""
    clock = [0.0]
    monkeypatch.setattr("flexprobe.pipeline.time.monotonic", lambda: clock[0])
    t = _StageTimer()
    clock[0] += 5                          # config validation etc., pre-ingest
    t.enter("ingest"); clock[0] += 4
    snap = t.snapshot()
    assert snap["ingest"] == 4
    assert snap["total"] == 9


# ---- tool versions ----

def test_collect_tool_versions_always_reports_python():
    assert _collect_tool_versions()["python"]


def test_collect_tool_versions_survives_a_missing_blast(monkeypatch):
    """A run with no BLAST on PATH must still produce a manifest — the version
    is simply omitted. An absent key means 'unknown'; a placeholder string
    would be a fabricated provenance claim."""
    monkeypatch.setattr("flexprobe.blast.blastn_version", lambda *a, **k: None)
    versions = _collect_tool_versions()
    assert "blastn" not in versions
    assert versions["python"]


def test_collect_tool_versions_records_blastn_when_present(monkeypatch):
    monkeypatch.setattr("flexprobe.blast.blastn_version", lambda *a, **k: "2.16.0")
    assert _collect_tool_versions()["blastn"] == "2.16.0"


def test_blastn_version_never_raises(monkeypatch):
    """Called while assembling the manifest, so a broken/absent blastn must not
    take down an otherwise successful multi-day run."""
    from flexprobe import blast

    monkeypatch.setattr(blast, "_find_tool", lambda name: "/nonexistent/blastn")
    assert blast.blastn_version() is None

    monkeypatch.setattr(blast, "_find_tool", lambda name: None)
    assert blast.blastn_version() is None


def test_skip_screen_run_does_not_claim_a_blastn_version(monkeypatch):
    """A --skip-screen run never invoked blastn, so recording the version of the
    binary that happens to be on PATH would be a fabricated provenance claim —
    the same failure mode as a placeholder string."""
    monkeypatch.setattr("flexprobe.blast.blastn_version", lambda *a, **k: "2.16.0")
    assert "blastn" not in _collect_tool_versions(screen_used=False)
    assert _collect_tool_versions(screen_used=True)["blastn"] == "2.16.0"


# ---- post-pipeline stage timing (select-coverage) ----

def test_record_stage_elapsed_folds_into_total(tmp_path):
    """select-coverage runs as its own invocation after the pipeline, so it has
    to add its share back. On the July 2026 mouse run this step was 3 h 05 m of
    a 44 h 14 m job — invisible to a manifest the pipeline wrote and closed."""
    man = tmp_path / "manifest.json"
    man.write_text(json.dumps({"elapsed_seconds": {"export": 10.0, "total": 100.0}}))

    assert record_stage_elapsed(str(man), "select_coverage", 25.0) is True
    elapsed = json.loads(man.read_text())["elapsed_seconds"]
    assert elapsed["select_coverage"] == 25.0
    assert elapsed["total"] == 125.0          # total tracks the whole job
    assert elapsed["export"] == 10.0          # existing stages untouched
    # the invariant the e2e test relies on still holds
    assert elapsed["total"] >= sum(v for k, v in elapsed.items() if k != "total")


def test_record_stage_elapsed_accumulates_on_rerun(tmp_path):
    """select-coverage is safe to re-run standalone; a second pass adds to the
    stage rather than silently replacing it."""
    man = tmp_path / "manifest.json"
    man.write_text(json.dumps({"elapsed_seconds": {"total": 5.0}}))
    record_stage_elapsed(str(man), "select_coverage", 2.0)
    record_stage_elapsed(str(man), "select_coverage", 3.0)
    elapsed = json.loads(man.read_text())["elapsed_seconds"]
    assert elapsed["select_coverage"] == 5.0
    assert elapsed["total"] == 10.0


def test_record_stage_elapsed_is_non_fatal(tmp_path):
    """The selection step is explicitly non-fatal, so a provenance update must
    never be the thing that fails it."""
    assert record_stage_elapsed(str(tmp_path / "absent.json"), "s", 1.0) is False

    broken = tmp_path / "broken.json"
    broken.write_text("{not json")
    assert record_stage_elapsed(str(broken), "s", 1.0) is False

    foreign = tmp_path / "foreign.json"
    foreign.write_text(json.dumps({"elapsed_seconds": "not-a-dict"}))
    assert record_stage_elapsed(str(foreign), "s", 1.0) is False


def test_record_stage_elapsed_seeds_a_manifest_without_timings(tmp_path):
    """Older runs predate elapsed_seconds; adding to them must still work."""
    man = tmp_path / "manifest.json"
    man.write_text(json.dumps({"app_version": "0.3.9"}))
    assert record_stage_elapsed(str(man), "select_coverage", 4.0) is True
    data = json.loads(man.read_text())
    assert data["elapsed_seconds"] == {"select_coverage": 4.0, "total": 4.0}
    assert data["app_version"] == "0.3.9"     # untouched


def test_blastn_version_does_not_enforce_the_minimum(monkeypatch):
    """Unlike check_version, this records what actually ran — including a BLAST
    too old for the pipeline's floor. Provenance, not a gate."""
    from flexprobe import blast

    class _R:
        returncode = 0
        stdout = "blastn: 2.2.31+\n Package: blast 2.2.31, build Jan 1 2016"

    monkeypatch.setattr(blast, "_find_tool", lambda name: "/usr/bin/blastn")
    monkeypatch.setattr(blast.subprocess, "run", lambda *a, **k: _R())
    assert blast.blastn_version() == "2.2.31"


def test_blastn_version_ignores_a_failed_invocation(monkeypatch):
    """A blastn that runs but exits non-zero gives empty/garbage stdout; report
    'unknown' rather than whatever the parser makes of it."""
    from flexprobe import blast

    class _R:
        returncode = 127
        stdout = ""

    monkeypatch.setattr(blast, "_find_tool", lambda name: "/usr/bin/blastn")
    monkeypatch.setattr(blast.subprocess, "run", lambda *a, **k: _R())
    assert blast.blastn_version() is None


def test_manifest_records_the_db_it_actually_screened_against(tmp_path):
    """A configured FASTA path is not the BLAST DB prefix.

    blastn reads ``<fasta-minus-extension>.nin`` etc., so globbing the
    configured ``…/foo.fa`` for sidecars finds nothing. Both manifest helpers
    globbed the configured path directly, which is why every shipped manifest
    recorded a 0-byte database and an empty ``db_sha256s`` — the two fields that
    answer "what was this screened against?".
    """
    from flexprobe.pipeline import _blast_db_sizes, _blast_db_sha256s

    fasta = tmp_path / "ref.fa"
    fasta.write_text(">r\nACGTACGTACGTACGTACGTACGT\n")
    # Stand-in sidecars: the helpers only need them to exist and have bytes.
    for ext in (".nhr", ".nin", ".nsq"):
        (tmp_path / f"ref{ext}").write_bytes(b"x" * 128)

    sizes = _blast_db_sizes([str(fasta)])
    assert sizes.get(str(fasta), 0) > 0, "size must come from the resolved prefix"
    digests = _blast_db_sha256s([str(fasta)])
    assert len(digests.get(str(fasta), "")) == 64

    # A DB prefix given directly still works.
    prefix = str(tmp_path / "ref")
    assert _blast_db_sizes([prefix]).get(prefix, 0) > 0
    assert len(_blast_db_sha256s([prefix]).get(prefix, "")) == 64


def test_missing_db_is_absent_not_zero(tmp_path):
    """Unknown must not masquerade as a real measurement of an empty DB."""
    from flexprobe.pipeline import _blast_db_sizes, _blast_db_sha256s
    gone = str(tmp_path / "nope.fa")
    assert gone not in _blast_db_sizes([gone])
    assert gone not in _blast_db_sha256s([gone])


def test_emptied_funnel_writes_a_summary(tmp_path):
    """"0 probes" and "never ran" must be distinguishable from the outputs."""
    import json
    from flexprobe.errors import FlexprobeError
    from flexprobe.summary import write_failed_summary

    out = tmp_path / "run"
    out.mkdir()
    (out / "filter_report.json").write_text(json.dumps(
        {"counts": {"total": 72, "passed": 0, "tn_junction": 60}}))
    write_failed_summary(str(out), FlexprobeError("No candidates passed filters"))

    payload = json.loads((out / "design_summary.json").read_text())
    assert payload["status"] == "no_candidates"
    assert payload["stats"]["tiles_total"] == 72
    assert payload["stats"]["pairs_output"] == 0
    assert payload["stats"]["filter_stats"]["tn_junction"] == 60


def test_gui_run_that_empties_the_funnel_writes_a_summary(tmp_path):
    """End-to-end on the GUI path, not just the CLI's.

    The 0.3.10 mislabelling fix landed in cli.py only, so a run launched from
    the app that emptied its funnel still wrote no design_summary.json and no
    manifest.json -- exactly the state that cannot be told apart from a run
    that never started. Drives the real worker so the wiring is what is tested.
    """
    import importlib.util
    import json
    import pytest as _pytest
    if importlib.util.find_spec("PySide6") is None:
        _pytest.skip("PySide6 not installed")
    from flexprobe.models import (
        ScreenResources, ModalityDefinition, RunConfig,
    )
    from flexprobe_gui.run_manager import _DesignWorker

    fa = tmp_path / "in.fa"
    fa.write_text(">G1\n" + "AAAAAAAAAA" * 12 + "\n")   # homopolymers: all rejected
    out = tmp_path / "out"
    cfg = RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=True),
    )
    w = _DesignWorker("r1", cfg)
    w.run()                     # synchronous; no event loop needed

    summary = out / "design_summary.json"
    assert summary.exists(), "a GUI run that produced 0 probes left no summary"
    payload = json.loads(summary.read_text())
    assert payload["status"] in ("no_candidates", "failed")
    if payload["status"] == "no_candidates":
        assert payload["stats"]["tiles_total"] > 0
        assert payload["stats"]["pairs_output"] == 0


def test_failed_run_summary_uses_the_same_schema_as_a_successful_one(tmp_path):
    """A summary's SHAPE must not depend on whether the run succeeded.

    If the abort path emitted `candidates_passed_filters` and omitted
    `specificity_passed`, a failed run would write a design_summary.json that
    is present, looks ordinary, and is silently keyed differently, and every
    consumer would die with a KeyError instead of reading "this run produced
    nothing". That is what `status` is for.
    """
    import json
    from flexprobe.summary import write_failed_summary
    out = tmp_path / "run"
    out.mkdir()
    (out / "filter_report.json").write_text(json.dumps(
        {"counts": {"total": 218, "passed": 0, "tn_junction": 218}}))
    write_failed_summary(str(out), RuntimeError("no candidates"))
    st = json.loads((out / "design_summary.json").read_text())
    assert st["status"] == "no_candidates"
    for k in ("tiles_total", "filter_passed", "specificity_passed",
              "filter_stats"):
        assert k in st["stats"], f"failed-run summary is missing {k}"
    assert st["stats"]["filter_passed"] == 0
    assert st["stats"]["specificity_passed"] == 0


# ---- a pair half carries the same provenance as its results.csv row ----

def test_pair_halves_carry_the_full_provenance_of_their_results_row(tmp_path):
    """Every provenance field on ``pair.lhs``/``pair.rhs`` matches results.csv.

    Halves built by hand from a curated subset of fields would leave
    everything outside that list at its default: a fully screened run would
    return ``pair.lhs.offtarget_recount_basis == "not_screened"`` while
    results.csv correctly said ``exact``, and the off-target evidence, the
    score components and the degenerate role would be absent on the halves.
    The artefacts are written from the scored candidates, so they stay right,
    but a GUI or script reading ``result["pairs"]`` would be told the run had
    never been screened.

    Asserted against the writer's OWN header and row renderer rather than a
    hand-listed set of fields, so a column added to results.csv is covered here
    the day it is added.
    """
    import csv
    import shutil

    import pytest

    from flexprobe.blast import build_db, discover_blast
    from flexprobe.screen import ENGINE_BLAST
    from flexprobe.models import (
        ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
        SpecificityThresholds,
    )
    from flexprobe.pipeline import run_pipeline
    from flexprobe.writers.results_table import BASE_HEADERS, _row

    if shutil.which("blastn") is None or shutil.which("makeblastdb") is None:
        pytest.skip("ncbi-blast+ (blastn/makeblastdb) not installed")

    target = (
        "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGTTACCGATGCATGC"
        "ATGCATGCATGCATGCATGCAACGTACGTACGTACGTACGTTGCATGCATG"
        "CATGCATGCATG"
    )
    fa = tmp_path / "target.fasta"
    fa.write_text(f">GENE1\n{target}\n")
    # Screen against the transcript itself: a real, fully screened run, so the
    # recount basis is decided rather than defaulted.
    db = str(tmp_path / "db")
    build_db(discover_blast()["makeblastdb"], str(fa), db, title="self")
    out = tmp_path / "out"
    result = run_pipeline(RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        # Pinned to BLAST: this test builds a real BLAST DB and passes the
        # PREFIX, which the k-mer default refuses by design (it indexes
        # sequence, so it wants the FASTA). The provenance contract under test
        # is engine-independent — see the k-mer twin in test_isoform_coverage.
        screen=ScreenResources(skip_screen=False, db_paths=[db],
                             engine=ENGINE_BLAST),
        # The junction-duplex gate is k-mer-only, so a BLAST config now has
        # to disable it explicitly (models.RunConfig.validate refuses the
        # combination rather than reporting a gate that never ran).
        specificity=SpecificityThresholds(off_target_junction_duplex_bp=0),
        export=ExportSettings(overwrite=True), discovery_only=True,
    ))
    pairs = result["pairs"]
    assert pairs, "expected the run to produce pairs"

    with open(out / "results.csv") as fh:
        rows = {r["probe_id"]: r for r in csv.DictReader(fh)}

    # `sequence` is the one field a half deliberately overrides: it holds the
    # oligo that half is ordered as, not the 50-bp mRNA window.
    #
    # `screen_engine` is excluded for a different reason, and it is not an
    # exception to the rule this test enforces: it is not candidate provenance
    # at all. It describes the RUN, is passed to the writer per call, and is
    # identical on every row — so it cannot live on a ProbeCandidate without
    # stamping run state onto every one of them. `_row(half)` therefore renders
    # the writer's default rather than this run's engine, and comparing the two
    # would be comparing a default against a value.
    overridden = {"mrna_target_50bp", "screen_engine"}

    for pair in pairs:
        row = rows[pair.lhs.probe_id]
        for half, oligo_col in ((pair.lhs, "lhs_oligo"), (pair.rhs, "rhs_oligo")):
            rendered = dict(zip(BASE_HEADERS, (str(v) for v in _row(half))))
            assert set(rendered) == set(BASE_HEADERS)
            for column, value in rendered.items():
                if column in overridden:
                    continue
                assert value == row[column], (
                    f"{column} differs between the pair half and its "
                    f"results.csv row for {pair.lhs.probe_id}"
                )
            assert half.sequence == row[oligo_col]

    # The field the bug was found on, named explicitly: a screened run must not
    # read as unscreened, and the off-target evidence must travel with it.
    lhs = pairs[0].lhs
    assert lhs.specificity_status != "not_screened"
    assert lhs.offtarget_recount_basis == rows[lhs.probe_id]["offtarget_recount_basis"]
    assert lhs.offtarget_recount_basis != "not_screened"
    # Both halves report the same evidence — the screen scored a site, not a
    # strand. (`==` is NaN-safe here via repr: offtarget_min_qcov is NaN when
    # there were no non-target hits at all, which is itself the claim.)
    for name in ("offtarget_min_qcov", "offtarget_worst_mismatches",
                 "offtarget_worst_junction_mm", "offtarget_nearmiss_count",
                 "best_offtarget_pident", "best_offtarget_bitscore",
                 "best_offtarget_evalue", "score_components", "junction_base",
                 "specificity_pass", "specificity_status", "rank"):
        assert repr(getattr(lhs, name)) == repr(getattr(pairs[0].rhs, name))


def test_pair_halves_do_not_share_mutable_state_with_the_scored_candidate():
    """A half's list/dict fields are copies, not aliases of the candidate's.

    ``dataclasses.replace`` is shallow. If the halves aliased the scored
    candidate's containers, anything mutating a half — a GUI annotating a
    selection, a downstream filter trimming isoform ids — would silently rewrite
    the object results.csv was written from, and the other half besides.
    """
    from flexprobe.models import ProbeCandidate
    from flexprobe.pipeline import _pair_half

    source = ProbeCandidate(
        gene_id="G", probe_id="G_0001", sequence="ACGT" * 12 + "AC",
        isoform_hit_ids=["ENST1", "ENST2"],
        score_components={"gc_proximity": 0.5},
        filter_fail_reasons=["gc"],
        strain_binding=[{"strain_id": "s1", "total_mismatches": 0}],
    )
    lhs = _pair_half(source, "LHS_OLIGO")
    rhs = _pair_half(source, "RHS_OLIGO")

    lhs.isoform_hit_ids.append("ENST3")
    lhs.score_components["gc_proximity"] = 0.0
    lhs.strain_binding[0]["total_mismatches"] = 99

    assert source.isoform_hit_ids == ["ENST1", "ENST2"]
    assert source.score_components == {"gc_proximity": 0.5}
    assert source.strain_binding[0]["total_mismatches"] == 0
    assert rhs.isoform_hit_ids == ["ENST1", "ENST2"]
    assert rhs.strain_binding[0]["total_mismatches"] == 0
    assert (lhs.sequence, rhs.sequence) == ("LHS_OLIGO", "RHS_OLIGO")
