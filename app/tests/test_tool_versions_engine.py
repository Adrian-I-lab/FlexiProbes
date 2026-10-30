"""``tool_versions`` must name only the tools that actually ran.

The k-mer engine has been the default since 0.5.1 and never invokes blastn, yet
the full-run manifest recorded a blastn version whenever any screen ran,
because the call site passed ``screen_used=not skip_screen``. Every k-mer run
therefore cited a BLAST+ version in its provenance. These tests run the real
pipeline under both engines, with the blastn version probe stubbed so the
answer does not depend on what is installed.
"""
import json
import os

import pytest

from flexprobe import blast
from flexprobe.models import (
    ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
    SpecificityThresholds,
)
from flexprobe.pipeline import run_pipeline
from flexprobe.screen import ENGINE_BLAST, ENGINE_KMER

GENES = {
    "geneA": ("GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
              "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"),
    "geneB": ("TTGGCCAATTGGCCAATTGGCCAATAGGCCTTAAGGCCTTAAGGCCTTAA"
              "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGG"),
}


def _manifest_tools(tmp_path, engine, skip_screen=False):
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in GENES.items()))
    out = str(tmp_path / "out")
    run_pipeline(RunConfig(
        fasta_path=str(fa), output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(engine=engine, db_paths=[str(fa)],
                               skip_screen=skip_screen),
        specificity=SpecificityThresholds(off_target_junction_duplex_bp=0),
        export=ExportSettings(overwrite=True),
    ))
    with open(os.path.join(out, "manifest.json")) as fh:
        return json.load(fh)["tool_versions"]


def test_kmer_run_does_not_claim_blastn(tmp_path, monkeypatch):
    monkeypatch.setattr(blast, "blastn_version", lambda *a, **k: "2.16.0")
    tools = _manifest_tools(tmp_path, ENGINE_KMER)
    assert "blastn" not in tools, (
        "a k-mer run never invoked blastn, so recording its version is a "
        "fabricated provenance claim")
    assert tools["python"]


def test_skip_screen_run_does_not_claim_blastn(tmp_path, monkeypatch):
    monkeypatch.setattr(blast, "blastn_version", lambda *a, **k: "2.16.0")
    assert "blastn" not in _manifest_tools(tmp_path, ENGINE_BLAST,
                                           skip_screen=True)


def test_blast_run_still_records_blastn(tmp_path, monkeypatch):
    # The other value of the same check: the BLAST engine did run blastn.
    try:
        blast.discover_blast()
    except blast.BlastError:
        pytest.skip("BLAST+ is not installed, so the BLAST engine cannot run")
    monkeypatch.setattr(blast, "blastn_version", lambda *a, **k: "2.16.0")
    assert _manifest_tools(tmp_path, ENGINE_BLAST)["blastn"] == "2.16.0"
