"""End-to-end BLAST: per-25-bp-half off-target screen (needs real blastn).

Skipped when ncbi-blast+ isn't available. These lock in the fix where BLAST is
an OFF-TARGET screen (a probe passes unless a 25-bp half hits a non-target) and
an on-target hit is NOT required — the old behaviour produced zero probes.
"""
import hashlib
import os
import shutil

import pytest

from flexprobe.blast import build_db, discover_blast, ensure_db_for_path
from flexprobe.screen import ENGINE_BLAST
from flexprobe.models import (
    ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
    SpecificityThresholds,
)
from flexprobe.pipeline import run_pipeline

pytestmark = pytest.mark.skipif(
    shutil.which("blastn") is None and shutil.which("makeblastdb") is None,
    reason="ncbi-blast+ (blastn/makeblastdb) not installed",
)

# A 120-bp target with an A at window position 25 for several windows.
TARGET = (
    "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGTTACCGATGCATGC"
    "ATGCATGCATGCATGCATGCAACGTACGTACGTACGTACGTTGCATGCATG"
    "CATGCATGCATG"
)


def _cfg(fa, out, db):
    return RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=False, db_paths=[db],
                             engine=ENGINE_BLAST),
        # The junction-duplex gate is k-mer-only, so a BLAST config now has
        # to disable it explicitly (models.RunConfig.validate refuses the
        # combination rather than reporting a gate that never ran).
        specificity=SpecificityThresholds(off_target_junction_duplex_bp=0),
        export=ExportSettings(overwrite=True), discovery_only=True,
    )


def _run(tmp_path, db_seq, db_id="OFFTGT"):
    fa = tmp_path / "target.fasta"
    fa.write_text(f">GENE1\n{TARGET}\n")
    dbfa = tmp_path / "off.fasta"
    dbfa.write_text(f">{db_id}\n{db_seq}\n")
    db = str(tmp_path / "offdb")
    build_db(discover_blast()["makeblastdb"], str(dbfa), db, title="off")
    out = tmp_path / "out"
    stats = run_pipeline(_cfg(fa, out, db))["stats"]
    return stats, out


def test_unrelated_db_passes_all(tmp_path):
    # Off-target DB shares nothing with the target -> every filter-passing probe
    # clears the screen (no on-target hit required).
    unrelated = "TTTTTTTTAAAAAAAACCCCCCCCGGGGGGGG" * 4
    stats, _ = _run(tmp_path, unrelated)
    assert stats["specificity_screened"] is True
    assert stats["filter_passed"] > 0
    assert stats["pairs_output"] == stats["filter_passed"]      # none dropped
    assert stats["specificity_status_counts"].get("passed") == stats["filter_passed"]


def test_target_in_db_is_screened_out(tmp_path):
    """DB == the target: every 25 bp half is a perfect off-target match.

    The screen therefore ought to drop EVERY probe — and the k-mer engine does
    (see the twin below, which asserts exactly that). BLAST does not, and this
    test asserts only the weaker property it can actually satisfy.

    WHY IT IS WEAKER: on blastn 2.17 this scenario yields 0
    pairs, but on blastn **2.12** — which is MIN_BLAST_VERSION, i.e. the oldest
    version this project claims to support — 6 of 17 probes survive, because
    blastn simply never reports a hit for one of their halves. Same input, same
    database, same thresholds; only the aligner build differs.

    That is not a flaw in this test. It is the reason the k-mer engine is the
    default: BLAST is a heuristic aligner, so a screen built on it gives
    version-dependent verdicts. Asserting 0 here would pin the suite to whichever
    blastn the developer happens to have installed.
    """
    stats, _ = _run(tmp_path, TARGET)
    assert stats["specificity_screened"] is True
    assert stats["filter_passed"] > 0
    # The screen must demonstrably remove probes...
    assert stats["pairs_output"] < stats["filter_passed"]
    # ...but how many survive is a property of the blastn build, not the rule.


def _sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def test_ensure_db_for_path_rebuilds_when_the_fasta_at_the_same_path_changes(tmp_path):
    """A second call must not hand back a DB built from the FASTA's old content.

    If ensure_db_for_path short-circuits on "sidecars already exist at this
    prefix" and returns without calling build_db, build_db's own
    content-checksum comparison against the .sha256 sidecar never runs. A
    caller that replaces the FASTA at the same path (a new reference, same
    filename) and calls again then gets the stale DB back untouched, and every
    downstream blastn search keeps searching the OLD sequence.
    """
    fa = tmp_path / "ref.fa"
    fa.write_text(">seqA\n" + "ACGT" * 50 + "\n")

    prefix1 = ensure_db_for_path(str(fa), makeblastdb_path=discover_blast()["makeblastdb"])
    hash_after_a = _sha256(prefix1 + ".nsq")

    # Overwrite the SAME path with different content, as if the reference FASTA
    # was replaced between runs.
    fa.write_text(">seqB\n" + "TTGGCCAA" * 40 + "\n")

    prefix2 = ensure_db_for_path(str(fa), makeblastdb_path=discover_blast()["makeblastdb"])
    hash_after_b_call = _sha256(prefix2 + ".nsq")

    assert prefix1 == prefix2
    assert hash_after_b_call != hash_after_a, (
        "ensure_db_for_path returned a DB whose .nsq content did not change "
        "after the FASTA at the same path was overwritten -- it handed back "
        "a stale database instead of rebuilding from the new content."
    )


def test_ensure_db_for_path_skips_the_rebuild_when_the_fasta_is_unchanged(tmp_path):
    """The common case must stay cheap: no rebuild when nothing changed."""
    fa = tmp_path / "ref.fa"
    fa.write_text(">seqA\n" + "ACGT" * 50 + "\n")
    makeblastdb_path = discover_blast()["makeblastdb"]

    prefix1 = ensure_db_for_path(str(fa), makeblastdb_path=makeblastdb_path)
    mtime1 = os.path.getmtime(prefix1 + ".nsq")

    prefix2 = ensure_db_for_path(str(fa), makeblastdb_path=makeblastdb_path)
    mtime2 = os.path.getmtime(prefix2 + ".nsq")

    assert mtime1 == mtime2, "unchanged FASTA should not trigger a makeblastdb rebuild"
