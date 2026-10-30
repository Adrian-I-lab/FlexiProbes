"""L1: the opt-in parallel tile+filter path must match the serial path exactly."""
import csv

from flexprobe.models import (
    ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
)
from flexprobe.pipeline import run_pipeline

# Three distinct ~120-bp "genes" so the run is genuinely multi-gene.
GENES = {
    "geneA": ("GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
              "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"),
    "geneB": ("TTGGCCAATTGGCCAATTGGCCAATAGGCCTTAAGGCCTTAAGGCCTTAA"
              "GGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGGTTCCAAGG"),
    "geneC": ("ACACACGTGTGTACACACGTGTGTACGTGTACACGTGTACACGTGTACAC"
              "TGTGACACTGTGACACTGTGACACTGTGACACTGTGACACTGTGACACTG"),
}


def _fasta(tmp_path):
    p = tmp_path / "multi.fasta"
    p.write_text("".join(f">{gid}\n{seq}\n" for gid, seq in GENES.items()))
    return str(p)


def _cfg(fa, out, workers):
    return RunConfig(
        fasta_path=fa, output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        workers=workers,
    )


def _rows(out):
    with open(f"{out}/results.csv") as fh:
        return [(r["probe_id"], r["mrna_target_50bp"], r["rank"])
                for r in csv.DictReader(fh)]


def test_parallel_matches_serial(tmp_path):
    fa = _fasta(tmp_path)
    out1 = str(tmp_path / "serial")
    out2 = str(tmp_path / "parallel")

    r_serial = run_pipeline(_cfg(fa, out1, workers=1))
    r_parallel = run_pipeline(_cfg(fa, out2, workers=2))

    # Identical stats and identical results.csv (same rows, same order).
    assert r_serial["stats"]["filter_passed"] == r_parallel["stats"]["filter_passed"]
    assert r_serial["stats"]["pairs_output"] == r_parallel["stats"]["pairs_output"]
    assert _rows(out1) == _rows(out2)
    assert _rows(out1), "expected non-empty results"
