"""The screen must be named after the engine that actually ran it.

The k-mer engine is the default. A progress event, log line or
`elapsed_seconds` key that says BLAST for a screen that never invoked blastn
leads anyone reading the log to write the wrong methods section. These tests
pin the naming to `config.screen.engine` in both directions, because the BLAST
engine is still selectable and must still say BLAST.
"""
import json
import os

import pytest

from flexprobe.models import (
    ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
    SpecificityThresholds,
)
from flexprobe.pipeline import run_pipeline
from flexprobe.screen import ENGINE_BLAST, ENGINE_KMER, screen_labels

GENES = {
    "geneA": ("GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
              "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"),
    "geneB": ("TTGGCCAATTGGCCAATTGGCCAATAGGCCTTAAGGCCTTAAGGCCTTAA"
              "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGG"),
}


def _run(tmp_path, engine, query_chunk=5000):
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in GENES.items()))
    out = str(tmp_path / "out")
    run_pipeline(RunConfig(
        fasta_path=str(fa), output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(engine=engine, db_paths=[str(fa)],
                             query_chunk=query_chunk),
        # Off for BOTH engines so the two arms of the parametrisation stay
        # comparable: the junction gate is k-mer-only and RunConfig.validate
        # now refuses to pair it with the BLAST engine.
        specificity=SpecificityThresholds(off_target_junction_duplex_bp=0),
        export=ExportSettings(overwrite=True),
    ))
    with open(os.path.join(out, "run_timeline.jsonl")) as fh:
        return out, [json.loads(line) for line in fh if line.strip()]


def test_labels_differ_per_engine():
    assert screen_labels(ENGINE_KMER)["stage"] == ENGINE_KMER
    assert screen_labels(ENGINE_BLAST)["stage"] == ENGINE_BLAST
    assert "k-mer" in screen_labels(ENGINE_KMER)["chunk"]
    assert "BLAST" in screen_labels(ENGINE_BLAST)["chunk"]


def test_unknown_engine_falls_back_rather_than_raising():
    # make_screen reports an unknown engine properly; this must not pre-empt it
    # with a KeyError from the label table.
    assert screen_labels("nonesuch") == screen_labels(ENGINE_BLAST)


@pytest.mark.parametrize("engine", [ENGINE_KMER, ENGINE_BLAST])
def test_stage_key_is_the_engine(tmp_path, engine):
    out, recs = _run(tmp_path, engine)
    assert engine in {r["stage"] for r in recs}
    other = ENGINE_BLAST if engine == ENGINE_KMER else ENGINE_KMER
    assert other not in {r["stage"] for r in recs}
    stages = json.load(open(os.path.join(out, "run_timings.json")))["stages"]
    assert engine in stages, "run_timings.json must name the engine that ran"


def test_kmer_run_never_claims_to_be_blast(tmp_path):
    # The whole point: a k-mer log must not contain the word BLAST in a screen
    # message. (Artefact names like candidates_prescreen.csv.gz are on-disk
    # contract and keep their names, so only screen-stage messages are checked.)
    _, recs = _run(tmp_path, ENGINE_KMER)
    offenders = [r["msg"] for r in recs
                 if r["stage"] == ENGINE_KMER and "blast" in r["msg"].lower()
                 and "prescreen" not in r["msg"].lower()
                 and "pre-blast" not in r["msg"].lower()]
    assert not offenders, f"k-mer run emitted BLAST wording: {offenders}"


@pytest.mark.parametrize("engine,word", [(ENGINE_KMER, "k-mer chunk"),
                                         (ENGINE_BLAST, "BLAST chunk")])
def test_chunk_lines_name_the_engine(tmp_path, engine, word):
    # More than one chunk is what gates the per-chunk progress line. GENES
    # yields 26 candidates, so query_chunk=13 makes exactly two chunks: one
    # blastn call per chunk, not one per candidate as query_chunk=1 did.
    _, recs = _run(tmp_path, engine, query_chunk=13)
    chunks = [r["msg"] for r in recs if "chunk" in r["msg"]]
    assert len(chunks) >= 2, f"expected several chunk lines, got {chunks}"
    assert all(m.startswith(word) for m in chunks), chunks
