"""``--prepare-only``: split one run into a prepare job and a screen job.

WHY THE SPLIT EXISTS. Python does not return freed pages to the OS, so a
process that tiled a transcriptome carries the tiling high-water mark for the
rest of its life — and the off-target screen then forks ``screen_workers``
children off that inflated parent. On a whole transcriptome, the same k-mer
screen needs several times more memory in a link that tiled first than in a
link that resumed from a prepared checkpoint, and a single job sized for that
peak is easily killed by a memory or walltime limit.

WHAT THESE TESTS PIN. The split is only worth anything if it changes nothing
about the answer, so the load-bearing test here is
:func:`test_resume_after_prepare_matches_single_process`: prepare-then-resume
must produce a byte-identical ``results.csv`` to a single-process run on the
same input. Everything else guards a specific way the split could go wrong
silently:

* a prepared output directory must be unmistakable for a finished one (no
  ``results.csv``, no ``manifest.json``, no ``.flexiprobes_DONE``, and a
  ``design_summary.json`` status of ``prepared``),
* the candidate checkpoint must SURVIVE (a finished run deletes it), and
* the two jobs must agree on the checkpoint fingerprint, or the screen job
  re-tiles the transcriptome inside the smaller allocation the split just
  gave it.
"""
import csv
import json
import os

import pytest

from flexprobe import candidate_cache
from flexprobe.constants import PREPARED_MARKER
from flexprobe.errors import ConfigError
from flexprobe.models import (
    ExportSettings, ModalityDefinition, RunConfig, ScreenResources,
    SpecificityThresholds,
)
from flexprobe.pipeline import PREPARE_MANIFEST_FILE, run_pipeline

# The DONE sentinel is the PBS script's "this run is finished" marker. Named
# here as a literal rather than imported from hpc_export so that a rename over
# there cannot quietly make this assertion vacuous.
DONE_SENTINEL = ".flexiprobes_DONE"

# Two ~200 bp records with balanced GC and short homopolymer runs, so step-1
# tiling yields plenty of filter-passing windows, and each record is its own
# off-target reference (db_paths points at the same FASTA). No external tools:
# the k-mer engine searches a FASTA directly.
GENES = {
    "geneA": ("ACGTACGTGCATGCATACGTACGTGCATGCATTGCATGCAACGTTGCAAC"
              "GTACGTGCATGCATACGTGTCAGCATGCATACGTACGTGCATTGCATGCA"
              "ACGTACGTGCATGCATACGTACGTGCATGCATTGCAACGTTGCATGCATA"
              "CGTACGTGCATGCATACGTGTCAGCATGCATACGTACGTGCATTGCATGC"),
    "geneB": ("TGCATGCAACGTACGTGCATTGCATGCAACGTACGTGCATGCATACGTAC"
              "GTGCATGCATTGCATGCAACGTTGCAACGTACGTGCATGCATACGTGTCA"
              "GCATGCATACGTACGTGCATTGCATGCAACGTACGTGCATGCATACGTAC"
              "GTGCATGCATTGCAACGTTGCATGCATACGTACGTGCATGCATACGTGTC"),
}


def _fasta(tmp_path, name="g.fa"):
    fa = tmp_path / name
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in GENES.items()))
    return fa


def _cfg(fa, out, **kw):
    """A config whose only variable is what `kw` says.

    Everything a resumed run must agree with the original about lives here, so
    the equivalence test cannot pass by accident on a config that differs.
    """
    kw.setdefault("prepare_only", False)
    return RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(db_paths=[str(fa)], query_chunk=5000),
        specificity=SpecificityThresholds(),
        export=ExportSettings(overwrite=True),
        discovery_only=True,
        keep_candidate_cache=False,
        **kw,
    )


def _timeline(out):
    path = os.path.join(str(out), "run_timeline.jsonl")
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


# ---------------------------------------------------------------------------
# where the run stops
# ---------------------------------------------------------------------------

def test_prepare_only_stops_before_rank_and_export(tmp_path):
    fa = _fasta(tmp_path)
    out = tmp_path / "out"
    result = run_pipeline(_cfg(fa, out, prepare_only=True))

    # Nothing a finished run produces.
    for name in ("results.csv", "manifest.json", "idt_opools.csv",
                 "cellranger_probe_set.csv", DONE_SENTINEL):
        assert not os.path.exists(out / name), \
            f"{name} must not exist after a prepare-only run"

    # No stage past the filter was even entered.
    stages = {r["stage"] for r in _timeline(out)}
    assert stages & {"ingest", "tile", "filter"}, stages
    assert not stages & {"kmer", "blast", "specificity", "rank",
                         "chemistry", "export", "done"}, stages

    assert result["pairs"] == []
    assert result["stats"]["prepare_only"] is True
    assert result["stats"]["run_complete"] is False
    assert result["stats"]["pairs_output"] == 0
    assert result["stats"]["prepared_candidates"] > 0


def test_prepare_only_writes_the_marker(tmp_path):
    fa = _fasta(tmp_path)
    out = tmp_path / "out"
    result = run_pipeline(_cfg(fa, out, prepare_only=True))

    marker = out / PREPARED_MARKER
    assert marker.exists()
    payload = json.loads(marker.read_text())
    assert payload["candidates"] == result["stats"]["prepared_candidates"]
    assert payload["resumed_from_existing_checkpoint"] is False
    assert payload["prepared_at"] and payload["app_version"]
    # The marker has to say, in words, that this is not a finished run: a human
    # finding it months later is one of its two readers.
    assert "NOT finished" in payload["note"]


def test_prepare_only_keeps_the_candidate_checkpoint(tmp_path):
    """A finished run deletes .candidate_cache/. A prepared one must not: the
    checkpoint is the entire deliverable of the prepare job."""
    fa = _fasta(tmp_path)
    out = tmp_path / "out"
    run_pipeline(_cfg(fa, out, prepare_only=True))
    assert os.path.isdir(candidate_cache.cache_dir(str(out)))

    # Contrast, on the same input: a full run with keep_candidate_cache off
    # clears it, which is what makes the assertion above meaningful.
    out2 = tmp_path / "out2"
    run_pipeline(_cfg(fa, out2))
    assert not os.path.isdir(candidate_cache.cache_dir(str(out2)))


def test_prepare_only_does_not_write_the_done_sentinel(tmp_path):
    """Separate from the sweep above because this is the specific mistake that
    would strand a bundle: `.flexiprobes_DONE` makes every remaining chain
    link no-op, so a prepared run that wrote it would never be screened."""
    fa = _fasta(tmp_path)
    out = tmp_path / "out"
    run_pipeline(_cfg(fa, out, prepare_only=True))
    assert not os.path.exists(out / DONE_SENTINEL)
    assert os.path.exists(out / PREPARED_MARKER)


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------

def test_prepare_manifest_records_prepare_only(tmp_path):
    fa = _fasta(tmp_path)
    out = tmp_path / "out"
    run_pipeline(_cfg(fa, out, prepare_only=True))

    man = json.loads((out / PREPARE_MANIFEST_FILE).read_text())
    assert man["params"]["prepare_only"] is True
    assert man["params"]["run_complete"] is False
    assert man["params"]["prepared_marker"] == PREPARED_MARKER
    assert man["params"]["prepared_candidates"] > 0
    # It must NOT claim work this process did not do.
    assert man["db_sha256s"] == {}
    assert man["reference_composition"] == {}
    assert "offtarget_recount" not in man["params"]
    assert "junction_duplex_measured" not in man["params"]
    # ...and it must not masquerade as the finished-run manifest.
    assert not os.path.exists(out / "manifest.json")


def test_screen_job_writes_the_real_manifest_over_the_top(tmp_path):
    """The prepare manifest must not block the screen job's own.

    write_manifest refuses to overwrite by default, so a prepare job that had
    written `manifest.json` would make the screen job fail at the last step of
    a multi-day run.
    """
    fa = _fasta(tmp_path)
    out = tmp_path / "out"
    run_pipeline(_cfg(fa, out, prepare_only=True))
    run_pipeline(_cfg(fa, out))

    man = json.loads((out / "manifest.json").read_text())
    assert man["params"].get("prepare_only") in (None, False)
    # Fields only a run that actually screened may carry — the prepare manifest
    # deliberately omits all three (see test_prepare_manifest_records_prepare_only).
    assert man["params"]["junction_duplex_measured"] is True
    assert man["params"]["offtarget_recount"]
    assert "kmer" in man["elapsed_seconds"]
    # The prepare manifest stays as the record of the first half.
    assert (out / PREPARE_MANIFEST_FILE).exists()


def test_prepared_summary_status_is_prepared(tmp_path, capsys):
    """The CLI must report a prepare job as a SUCCESS with its own status.

    write_ok_summary would stamp "no_probes" and cmd_design would exit 1, which
    the PBS chain classifies as a deterministic failure and aborts the whole
    chain on — so this is the difference between a working split and a chain
    that dies at link 1.
    """
    from flexprobe.cli import main

    fa = _fasta(tmp_path)
    out = tmp_path / "out"
    rc = main(["design", "--fasta", str(fa), "--out", str(out),
               "--blast-db", str(fa), "--prepare-only", "--force", "--quiet"])
    assert rc == 0

    summary = json.loads((out / "design_summary.json").read_text())
    assert summary["status"] == "prepared"
    assert summary["stats"]["prepare_only"] is True
    assert summary["stats"]["pairs_output"] == 0

    printed = capsys.readouterr().out
    assert "PREPARE-ONLY" in printed
    assert "INCOMPLETE" in printed


# ---------------------------------------------------------------------------
# the resume must not change the answer
# ---------------------------------------------------------------------------

def test_resume_after_prepare_matches_single_process(tmp_path):
    """THE test. A split run and a single-process run on the same input must
    produce a byte-identical results.csv.

    A speed win that changes the answer is not a win. results.csv is the
    primary deliverable and every downstream artefact is derived from it, so
    byte equality here is the strongest available statement that the split is
    result-neutral.
    """
    fa = _fasta(tmp_path)
    split_out = tmp_path / "split"
    single_out = tmp_path / "single"

    # Arm NEW: job A prepares, job B resumes the same directory.
    run_pipeline(_cfg(fa, split_out, prepare_only=True))
    run_pipeline(_cfg(fa, split_out))
    # Arm OLD: one process, all stages.
    run_pipeline(_cfg(fa, single_out))

    split_csv = (split_out / "results.csv").read_bytes()
    single_csv = (single_out / "results.csv").read_bytes()
    assert split_csv == single_csv, (
        "the split changed the answer: results.csv differs between a "
        "prepare+resume run and a single-process run on the same input"
    )
    # Not a vacuous pass on two empty files.
    assert len(list(csv.DictReader(single_csv.decode().splitlines()))) > 0


def test_resume_after_prepare_reloads_rather_than_retiling(tmp_path):
    """The screen job must READ the checkpoint, not rebuild it.

    Without this the split is worse than no split: two allocations, and the
    expensive half runs in both of them.
    """
    fa = _fasta(tmp_path)
    out = tmp_path / "out"
    run_pipeline(_cfg(fa, out, prepare_only=True))
    run_pipeline(_cfg(fa, out))

    msgs = [r["msg"] for r in _timeline(out)]
    assert any("from checkpoint" in m for m in msgs), msgs[:20]
    assert any("skipped tiling + filtering" in m for m in msgs)


def test_prepare_only_does_not_change_the_checkpoint_fingerprint(tmp_path):
    """Job A and job B share ONE config.json and differ only by the flag.

    So the flag must be outside the candidate-cache fingerprint. Measured
    before the fix: the same config hashed to d774c45d... with the flag and
    4c88fb62... without it, which would have made the screen job miss the
    checkpoint and re-tile the transcriptome in its smaller allocation.
    """
    fa = _fasta(tmp_path)
    with_flag = candidate_cache.fingerprint(
        _cfg(fa, tmp_path / "a", prepare_only=True))
    without = candidate_cache.fingerprint(
        _cfg(fa, tmp_path / "a", prepare_only=False))
    assert with_flag == without
    assert "prepare_only" in candidate_cache.FINGERPRINT_IGNORED

    # Positive control: the fingerprint is not simply insensitive to
    # everything. A field that DOES change the candidate set must move it.
    other = candidate_cache.fingerprint(
        _cfg(fa, tmp_path / "a", prepare_only=True, tile_step=7))
    assert other != with_flag


def test_a_second_prepare_is_idempotent(tmp_path):
    """Re-running prepare over a prepared directory must not corrupt it.

    The PBS chain uses `afterany`, so a prepare job that is walltime-killed
    after writing the checkpoint but before the marker is retried — and the
    retry lands on a directory that already holds a checkpoint.
    """
    fa = _fasta(tmp_path)
    out = tmp_path / "out"
    first = run_pipeline(_cfg(fa, out, prepare_only=True))
    second = run_pipeline(_cfg(fa, out, prepare_only=True))

    assert (second["stats"]["prepared_candidates"]
            == first["stats"]["prepared_candidates"])
    assert second["stats"]["prepared_from_existing_checkpoint"] is True
    marker = json.loads((out / PREPARED_MARKER).read_text())
    assert marker["resumed_from_existing_checkpoint"] is True
    # And the screen job still finishes from it.
    run_pipeline(_cfg(fa, out))
    assert (out / "results.csv").exists()


# ---------------------------------------------------------------------------
# refused combinations
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kw,needle", [
    ({"resume_screen": False}, "resume_screen"),
    ({"multi_strain": True}, "multi_strain"),
])
def test_prepare_only_refuses_configs_that_write_no_checkpoint(
        tmp_path, kw, needle):
    fa = _fasta(tmp_path)
    cfg = _cfg(fa, tmp_path / "out", prepare_only=True, **kw)
    with pytest.raises(ConfigError) as exc:
        cfg.validate()
    assert needle in str(exc.value)


def test_prepare_only_refuses_skip_screen(tmp_path):
    fa = _fasta(tmp_path)
    cfg = _cfg(fa, tmp_path / "out", prepare_only=True)
    cfg.screen.skip_screen = True
    with pytest.raises(ConfigError) as exc:
        cfg.validate()
    assert "skip_screen" in str(exc.value)


def test_prepare_only_refuses_cached_hits(tmp_path):
    fa = _fasta(tmp_path)
    cfg = _cfg(fa, tmp_path / "out", prepare_only=True)
    cfg.screen.use_cached_hits = True
    with pytest.raises(ConfigError) as exc:
        cfg.validate()
    assert "use_cached_hits" in str(exc.value)


def test_plain_run_is_unaffected(tmp_path):
    """The default is off, and a run without the flag behaves exactly as
    before: results, manifest, and no prepare artefacts."""
    fa = _fasta(tmp_path)
    out = tmp_path / "out"
    result = run_pipeline(_cfg(fa, out))
    assert (out / "results.csv").exists()
    assert (out / "manifest.json").exists()
    assert not (out / PREPARED_MARKER).exists()
    assert not (out / PREPARE_MANIFEST_FILE).exists()
    assert "prepare_only" not in result["stats"]
