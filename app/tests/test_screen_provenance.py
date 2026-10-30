# flexiprobes — probe design for 10x Genomics Flex v1 / Flex v2 / Visium.
# Copyright (C) 2026 Adrian Ilich and contributors
#
# This file is part of flexiprobes.
#
# flexiprobes is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# flexiprobes is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Which engine produced a result must be inseparable from the result.

``offtarget_recount_basis``, ``offtarget_min_qcov``,
``best_offtarget_bitscore``/``_evalue`` and ``offtarget_nearmiss_count`` all
mean different things per engine. A results.csv circulates detached from its
manifest, so the engine has to be on the rows too — otherwise two runs get
compared on columns that are not comparable.
"""
from __future__ import annotations

import csv
import json
import random

import pytest

from flexprobe.models import (
    ScreenResources, ExportSettings, ModalityDefinition, ProbeCandidate,
    RunConfig,
)
from flexprobe.pipeline import run_pipeline
from flexprobe.screen import ENGINE_KMER
from flexprobe.writers.results_table import BASE_HEADERS, write_results_table


def _rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


# ---------------------------------------------------------------------------
# results.csv
# ---------------------------------------------------------------------------

def test_results_table_has_a_screen_engine_column():
    assert "screen_engine" in BASE_HEADERS


def test_screen_engine_column_records_the_engine(tmp_path):
    c = ProbeCandidate(gene_id="g", probe_id="p0", sequence="ACGT" * 12 + "AC",
                       start=0, end=50)
    out = str(tmp_path / "r.csv")
    write_results_table([c], out, overwrite=True, screen_engine=ENGINE_KMER)
    with open(out) as fh:
        rows = list(csv.DictReader(fh))
    assert rows[0]["screen_engine"] == ENGINE_KMER


def test_screen_engine_label_tracks_the_run_default(tmp_path):
    """The writer's label default must follow the RUN default, not be frozen.

    If they diverge, a caller that omits the argument gets one engine's
    behaviour under the other engine's label — which is worse than no column,
    because it is confidently wrong.
    """
    c = ProbeCandidate(gene_id="g", probe_id="p0", sequence="A" * 50,
                       start=0, end=50)
    out = str(tmp_path / "r.csv")
    write_results_table([c], out, overwrite=True)
    with open(out) as fh:
        rows = list(csv.DictReader(fh))
    from flexprobe.constants import SCREEN_ENGINE_DEFAULT
    assert rows[0]["screen_engine"] == SCREEN_ENGINE_DEFAULT
    assert rows[0]["screen_engine"] == ScreenResources().engine


def test_column_count_matches_header_count(tmp_path):
    """Inserting a column mid-table is how a CSV silently shifts by one."""
    c = ProbeCandidate(gene_id="g", probe_id="p0", sequence="ACGT" * 12 + "AC",
                       start=0, end=50)
    out = str(tmp_path / "r.csv")
    write_results_table([c], out, overwrite=True, screen_engine=ENGINE_KMER)
    with open(out) as fh:
        rows = list(csv.reader(fh))
    assert len(rows[0]) == len(rows[1]) == len(BASE_HEADERS)


# ---------------------------------------------------------------------------
# End to end under the k-mer engine
# ---------------------------------------------------------------------------

def _run(tmp_path, engine, db_path):
    rng = random.Random(99)
    target = _rand_seq(600, rng)
    fa = tmp_path / "target.fa"
    fa.write_text(f">GENE1\n{target}\n")
    off = tmp_path / "off.fa"
    off.write_text(f">OFFTGT\n{_rand_seq(3000, rng)}\n")
    out = tmp_path / f"out_{engine}"
    cfg = RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=False, db_paths=[db_path],
                             engine=engine),
        export=ExportSettings(overwrite=True), discovery_only=True,
    )
    run_pipeline(cfg)
    return out


def test_kmer_engine_runs_end_to_end_and_labels_its_outputs(tmp_path):
    """A full pipeline run on the k-mer engine — no blastn on the machine
    required, which is itself part of the point."""
    rng = random.Random(99)
    off = tmp_path / "off.fa"
    off.write_text(f">OFFTGT\n{_rand_seq(3000, rng)}\n")
    out = _run(tmp_path, ENGINE_KMER, str(off))

    with open(out / "results.csv") as fh:
        rows = list(csv.DictReader(fh))
    assert rows, "the k-mer run produced no probes"
    assert {r["screen_engine"] for r in rows} == {ENGINE_KMER}
    # The property the engine exists for: exact by construction, full coverage.
    assert {r["offtarget_recount_basis"] for r in rows} == {"exact"}

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["params"]["screen_engine"] == ENGINE_KMER
    assert manifest["params"]["screen"]["engine"] == ENGINE_KMER
    # The manifest's run-level basis must AGREE with the per-probe column.
    # Keying it off SubjectWindows alone reported 'not_screened' for a fully
    # screened k-mer run, which reads as "no specificity was done" and would
    # invite a reader to throw the run away.
    assert manifest["params"]["offtarget_recount"]["basis"] == "exact"


def test_kmer_run_needs_no_blastdbcmd_or_parse_seqids(tmp_path, monkeypatch):
    """The k-mer engine must not reach for BLAST tooling at all.

    If it did, the whole premise fails: the fork storm, the -parse_seqids hard
    error and the conservative fallback all come in through these calls.
    """
    from flexprobe import blast as fp_blast

    def _boom(*_a, **_k):
        raise AssertionError("the k-mer engine invoked BLAST tooling")

    monkeypatch.setattr(fp_blast, "discover_blast", _boom)
    monkeypatch.setattr(fp_blast, "find_blastdbcmd", _boom)
    monkeypatch.setattr(fp_blast, "run_blastn", _boom)

    rng = random.Random(99)
    off = tmp_path / "off.fa"
    off.write_text(f">OFFTGT\n{_rand_seq(3000, rng)}\n")
    out = _run(tmp_path, ENGINE_KMER, str(off))
    assert (out / "results.csv").exists()


def test_cli_exposes_the_engine_choice():
    from flexprobe.cli import main
    with pytest.raises(SystemExit) as e:
        main(["design", "--fasta", "x", "--out", "y",
              "--screen-engine", "nonsense"])
    assert e.value.code != 0


TARGET = (
    "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGTTACCGATGCATGC"
    "ATGCATGCATGCATGCATGCAACGTACGTACGTACGTACGTTGCATGCATG"
    "CATGCATGCATG"
)


def test_target_in_db_is_screened_out_kmer_engine(tmp_path):
    """A perfect off-target match must reject EVERY probe, guaranteed.

    Twin of test_blast_e2e.test_target_in_db_is_screened_out, which had to give
    this assertion up: on blastn 2.12 (MIN_BLAST_VERSION) 6 of 17 probes survive
    the identical scenario because the aligner never reports one of their
    halves.

    The k-mer search is exhaustive over its mismatch band, so a perfect match
    cannot be missed — there is no aligner heuristic to be version-dependent
    about. Deliberately NOT skipped when blastn is absent: needing no BLAST is
    half the point, and this is the assertion the BLAST twin above had to give
    up.
    """
    fa = tmp_path / "target.fasta"
    fa.write_text(f">GENE1\n{TARGET}\n")
    dbfa = tmp_path / "off.fasta"
    dbfa.write_text(f">OFFTGT\n{TARGET}\n")
    out = tmp_path / "out"
    cfg = RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
        # The k-mer engine indexes the FASTA directly — no database to build.
        screen=ScreenResources(skip_screen=False, db_paths=[str(dbfa)],
                             engine=ENGINE_KMER),
        export=ExportSettings(overwrite=True), discovery_only=True,
    )
    stats = run_pipeline(cfg)["stats"]
    assert stats["specificity_screened"] is True
    assert stats["filter_passed"] > 0
    assert stats["pairs_output"] == 0
