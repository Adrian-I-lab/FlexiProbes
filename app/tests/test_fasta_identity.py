"""Input integrity: FASTA record ids must survive parsing intact and be unique.

The bug this pins down: ``io_fasta`` truncated every record id at the first
``|``, and the pipeline keys sequences by id in a dict. A legacy NCBI-style
header set (``>gi|12345|ref|NM_000546.6| TP53``, ...) therefore collapsed onto
the single id ``gi`` — 5 records in, "Loaded 5 genes" reported, ONE gene
actually designed against, no warning anywhere. Silent data loss in the input
stage.
"""
import pytest

from flexprobe.errors import FastaParseError
from flexprobe.io_fasta import fasta_to_dict, iter_fasta, parse_fasta
from flexprobe.models import (
    ScreenResources, ExportSettings, FilterThresholds, ModalityDefinition,
    RunConfig,
)
from flexprobe.pipeline import run_pipeline

from conftest import SYNTH_200


# An NCBI-style pipe-delimited header set: five DISTINCT records whose headers
# all begin with the same first pipe-field. Keyed on that field alone, they
# would collapse to one gene.
NCBI_HEADERS = [
    "gi|12345|ref|NM_000546.6| TP53 transcript variant 1",
    "gi|22345|ref|NM_001126112.3| TP53 transcript variant 2",
    "gi|32345|ref|NM_007294.4| BRCA1 transcript variant 1",
    "gi|42345|ref|NM_000059.4| BRCA2",
    "gi|52345|ref|NM_004333.6| BRAF",
]


def _write(path, headers, seqs=None):
    seqs = seqs or [SYNTH_200] * len(headers)
    with open(path, "w") as fh:
        for h, s in zip(headers, seqs):
            fh.write(f">{h}\n{s}\n")
    return str(path)


# ---------------------------------------------------------------------------
# N records in => N genes out
# ---------------------------------------------------------------------------

def test_ncbi_pipe_headers_yield_one_record_each(tmp_path):
    """The regression: 5 pipe-delimited records must parse as 5 distinct ids."""
    fa = _write(tmp_path / "ncbi.fa", NCBI_HEADERS)
    records = parse_fasta(fa)
    assert len(records) == 5
    ids = [gid for gid, _d, _s in records]
    assert len(set(ids)) == 5, f"records merged: {ids}"
    # The id is the FULL first whitespace token, not the part before '|'.
    assert ids[0] == "gi|12345|ref|NM_000546.6|"
    assert fasta_to_dict(fa).keys() == set(ids)
    # ...and the description is still everything after the first token.
    assert records[0][1] == "TP53 transcript variant 1"


@pytest.mark.parametrize("n", [1, 2, 5, 17])
def test_n_distinct_records_yield_n_genes(tmp_path, n):
    """Plain (non-piped) headers: the count must round-trip exactly."""
    fa = _write(tmp_path / f"n{n}.fa", [f"TX{i:03d} desc {i}" for i in range(n)])
    assert len(parse_fasta(fa)) == n
    assert len(fasta_to_dict(fa)) == n


def test_pipeline_designs_against_every_ncbi_record(tmp_path):
    """End-to-end: the pipeline must see 5 genes, not 1.

    This is the assertion that would have caught the bug in production — the
    old parser reported "Loaded 5 genes" while ``gene_seqs`` held exactly one.
    """
    # Distinct sequences so a merge would be visible in the outputs too.
    seqs = [SYNTH_200[i:] + SYNTH_200[:i] for i in range(len(NCBI_HEADERS))]
    fa = _write(tmp_path / "ncbi.fa", NCBI_HEADERS, seqs)
    cfg = RunConfig(
        fasta_path=fa, output_dir=str(tmp_path / "out"),
        modality=ModalityDefinition(name="flex", version="v2"),
        filters=FilterThresholds(gc_min=0.30, gc_max=0.80, max_n_per_half=2,
                                 max_homopolymer_run=6),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        tile_length=50, tile_step=10,
    )
    result = run_pipeline(cfg)
    assert result["stats"]["genes"] == 5
    designed = {p.lhs.gene_id for p in result["pairs"]}
    assert designed == set(NCBI_HEADERS[i].split()[0] for i in range(5))


# ---------------------------------------------------------------------------
# duplicate ids fail closed
# ---------------------------------------------------------------------------

def test_duplicate_ids_raise_and_name_the_collisions(tmp_path):
    fa = _write(tmp_path / "dup.fa", ["G1 a", "G2 b", "G1 c", "G3 d", "G2 e"])
    with pytest.raises(FastaParseError) as exc:
        parse_fasta(fa)
    msg = str(exc.value)
    assert "'G1'" in msg and "'G2'" in msg      # every collision named
    assert "'G3'" not in msg                    # the innocent record is not
    assert "5 records" in msg and "3 unique" in msg


def test_opt_in_pipe_split_still_detects_the_merge_it_causes(tmp_path):
    """The legacy behaviour is available, but can no longer merge silently."""
    fa = _write(tmp_path / "ncbi.fa", NCBI_HEADERS)
    with pytest.raises(FastaParseError) as exc:
        parse_fasta(fa, split_header_on_pipe=True)
    assert "'gi'" in str(exc.value)
    assert "split_header_on_pipe" in str(exc.value)   # names the cause


def test_pipe_split_opt_in_works_when_the_first_field_is_unique(tmp_path):
    """The committed OROV validation inputs are this shape."""
    fa = _write(tmp_path / "orov.fa",
                ["OROV_BeAn19991|Large|KP052850.1",
                 "OROV_PEIAM4578|Large|PQ073184.1"])
    ids = [g for g, _d, _s in parse_fasta(fa, split_header_on_pipe=True)]
    assert ids == ["OROV_BeAn19991", "OROV_PEIAM4578"]
    # Default (no split) keeps the full token — also unique, just longer.
    ids_full = [g for g, _d, _s in parse_fasta(fa)]
    assert ids_full == ["OROV_BeAn19991|Large|KP052850.1",
                        "OROV_PEIAM4578|Large|PQ073184.1"]


def test_streaming_parser_fails_closed_on_the_second_occurrence(tmp_path):
    """iter_fasta cannot buffer every id's line number, so it stops at the
    first collision rather than streaming a record that will be overwritten."""
    fa = _write(tmp_path / "dup.fa", ["G1", "G2", "G1"])
    with pytest.raises(FastaParseError, match="Duplicate FASTA record id 'G1'"):
        list(iter_fasta(fa))
    # ...and the check can be turned off for callers doing their own pass.
    assert len(list(iter_fasta(fa, check_duplicates=False))) == 3


def test_config_flag_reaches_the_parser(tmp_path):
    """The opt-in is a RunConfig field, so it is recorded in the manifest and
    folded into the resume fingerprint like every other output-affecting knob."""
    assert RunConfig().split_header_on_pipe is False
    fa = _write(tmp_path / "orov.fa",
                ["OROV_A|Large|KP052850.1", "OROV_B|Large|PQ073184.1"])
    cfg = RunConfig(
        fasta_path=fa, output_dir=str(tmp_path / "out"),
        modality=ModalityDefinition(name="flex", version="v2"),
        filters=FilterThresholds(gc_min=0.30, gc_max=0.80, max_n_per_half=2,
                                 max_homopolymer_run=6),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        tile_length=50, tile_step=10,
        split_header_on_pipe=True,
    )
    result = run_pipeline(cfg)
    assert {p.lhs.gene_id for p in result["pairs"]} <= {"OROV_A", "OROV_B"}
