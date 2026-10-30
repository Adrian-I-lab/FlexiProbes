"""Every step must be recorded AS IT RUNS, not summarised at the end.

The 0.5.2 `candidates_per_gene` regression cost 12-24 h per whole-transcriptome
run and hid for weeks because the only timing artefact was written before the
slow part ran. These tests pin the two properties that would have surfaced it:
the tail of the run is a named, measured stage, and the timeline is on disk
before the run finishes (so a walltime kill still points at the culprit).
"""
import json
import os

from flexprobe.models import (
    ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
)
from flexprobe.pipeline import run_pipeline

GENES = {
    "geneA": ("GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
              "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"),
    "geneB": ("TTGGCCAATTGGCCAATTGGCCAATAGGCCTTAAGGCCTTAAGGCCTTAA"
              "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGG"),
}


def _run(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in GENES.items()))
    out = str(tmp_path / "out")
    run_pipeline(RunConfig(
        fasta_path=str(fa), output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
    ))
    return out


def test_writes_both_timing_artefacts(tmp_path):
    out = _run(tmp_path)
    assert os.path.exists(os.path.join(out, "run_timings.json"))
    assert os.path.exists(os.path.join(out, "run_timeline.jsonl"))


def test_summary_build_is_a_measured_stage(tmp_path):
    # The stats dict (candidates_per_gene lives there) must be built inside a
    # stage boundary, or its cost is attributed to nothing.
    out = _run(tmp_path)
    stages = json.load(open(os.path.join(out, "run_timings.json")))["stages"]
    assert "summary" in stages, (
        "the post-export summary build is unmeasured again — a slow one would "
        "be invisible in both the manifest and run_timings.json"
    )


def test_every_stage_reports_wall_cpu_and_memory(tmp_path):
    out = _run(tmp_path)
    stages = json.load(open(os.path.join(out, "run_timings.json")))["stages"]
    assert "total" in stages
    for name, s in stages.items():
        assert set(s) == {"wall_s", "cpu_s", "cores_mean", "peak_rss_mb",
                          "events"}, f"{name} is missing a metric"
        assert s["wall_s"] >= 0.0
        assert s["cpu_s"] >= 0.0
    assert stages["total"]["peak_rss_mb"] > 0.0


def test_timeline_is_appended_as_the_run_goes(tmp_path):
    # Records must be flushed per boundary, in order, each naming its step —
    # that is what makes a killed run diagnosable.
    out = _run(tmp_path)
    with open(os.path.join(out, "run_timeline.jsonl")) as fh:
        recs = [json.loads(line) for line in fh if line.strip()]
    assert len(recs) > 3
    assert [r["seq"] for r in recs] == list(range(1, len(recs) + 1))
    for r in recs:
        assert r["stage"] and r["msg"]
        assert r["step_s"] >= 0.0 and r["step_cpu_s"] >= 0.0
    # Monotonic elapsed: the file is a timeline, not a bag of records.
    assert all(a["elapsed_s"] <= b["elapsed_s"]
               for a, b in zip(recs, recs[1:]))


def test_resumed_link_appends_rather_than_truncating(tmp_path):
    # A PBS chain runs the same output dir repeatedly; the links must read as
    # one continuous timeline or a resumed run loses its own history.
    out = _run(tmp_path)
    path = os.path.join(out, "run_timeline.jsonl")
    first = len(open(path).readlines())
    _run(tmp_path)
    assert len(open(path).readlines()) > first
