"""Cell Ranger / Space Ranger probe-set CSV FILE FORMAT.

The exported "Cell Ranger probe set" CSV must be a valid Cell Ranger probe
set. That means a leading '#' metadata block, an `included` column, no
`region` value fabricated without a gene model, and a warning when a FASTA
record id is a transcript accession rather than a gene id. Three writers
produce this format, and this file exercises all three:
  - writers/cellranger_probe_set.py (pipeline one-shot export)
  - select_export.export_selected_probes (per-run selection export)
  - select_export.export_project_bulk (project bulk export)

These are FILE FORMAT tests, not probe-selection tests. Nothing here
depends on, or changes, any window, probe_seq, or ranking.
"""
import csv
import json
import logging
import os
import pathlib

import pytest

from flexprobe.models import ProbeCandidate, ProbePair
from flexprobe.select_export import export_project_bulk, export_selected_probes
from flexprobe.writers.cellranger_probe_set import (
    CR_PANEL_TYPE, CR_PROBE_SET_FILE_FORMAT, write_cellranger_probe_set,
)

# ---------------------------------------------------------------------------
# Locating a real 10x vendor probe-set CSV to diff against.
#
# The vendor CSVs are ~6 MB third-party downloads from 10x Genomics and are
# not part of this repository, so this one test needs a copy on disk. Point an
# environment variable at it:
#
#   FLEXIPROBES_TENX_PROBE_SET=<dir or .csv file>
#
# When the variable is unset and no copy is found, the test XFAILs rather than
# skips, so "we did not actually compare against a vendor file" is visible in
# the run summary instead of indistinguishable from a pass.
#
# Point it at a directory and the fixture picks the right file itself: it must
# be probe_set_file_format 1.0. That is not a stylistic preference. Format 1.0
# is the minimal four-column variant (gene_id, probe_seq, probe_id, included);
# 2.0 adds `region` and 3.0 adds `region` + `gene_name`, and flexiprobes
# deliberately emits neither (see test_region_column_is_never_fabricated —
# there is no gene model in this pipeline to populate them honestly). Diffing
# our export against a 3.0 file would therefore fail for a reason that is not
# a defect. Of 10x's published Flex v1 panels, only
# Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv is format 1.0.
# ---------------------------------------------------------------------------
_TENX_ENV = "FLEXIPROBES_TENX_PROBE_SET"
_REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
_DEFAULT_TENX_DIR = _REPO_ROOT / "validation" / "tenx_flex_v1" / "data"

# The probe_set_file_format our export corresponds to.
_WANTED_FORMAT = "1.0"


def _is_format_1_0(path: pathlib.Path) -> bool:
    """True if ``path`` is a vendor probe set declaring format 1.0."""
    try:
        with open(path, encoding="utf-8") as fh:
            first = fh.readline().strip()
    except OSError:
        return False
    return first == f"#probe_set_file_format={_WANTED_FORMAT}"


@pytest.fixture(scope="module")
def vendor_probe_set_csv() -> str:
    """Path to a real 10x format-1.0 probe-set CSV, or xfail if there is none.

    XFAIL, not skip. A skip here is indistinguishable from a pass in the
    summary line, so a run that never compared against a vendor file would
    look like a success.
    """
    raw = os.environ.get(_TENX_ENV)
    where = pathlib.Path(raw).expanduser() if raw else _DEFAULT_TENX_DIR
    source = f"${_TENX_ENV}={raw!r}" if raw else f"default {where}"

    if where.is_file():
        if not _is_format_1_0(where):
            pytest.xfail(
                f"{source} is not a probe_set_file_format={_WANTED_FORMAT} "
                f"file. 2.0 adds a `region` column and 3.0 adds "
                f"`region`+`gene_name`; flexiprobes emits neither by design, "
                f"so only a 1.0 file is a like-for-like structural diff."
            )
        return str(where)

    if where.is_dir():
        matches = sorted(c for c in where.glob("*.csv") if _is_format_1_0(c))
        if matches:
            return str(matches[0])
        pytest.xfail(
            f"no probe_set_file_format={_WANTED_FORMAT} CSV under {where} "
            f"({source}). Found {len(list(where.glob('*.csv')))} CSV(s), none "
            f"of them format {_WANTED_FORMAT}."
        )

    pytest.xfail(
        f"no 10x vendor probe-set CSV at {where} ({source}). These are large "
        f"third-party downloads and are not committed. Fetch 10x's Chromium "
        f"Human Transcriptome Probe Set v1.0 (GRCh38-2020-A) and set "
        f"{_TENX_ENV} to a directory or file holding one."
    )


def _cand(gene_id: str) -> ProbeCandidate:
    c = ProbeCandidate(gene_id=gene_id, probe_id=gene_id + "_p", sequence="A" * 50)
    c.probe_rc = "A" * 50
    return c


def _pair(gene_id: str) -> ProbePair:
    return ProbePair(lhs=_cand(gene_id), rhs=_cand(gene_id))


def _read_header_and_columns(path):
    """Return ``([raw '#' metadata lines], [column header names])``."""
    with open(path) as fh:
        lines = fh.readlines()
    i = 0
    meta = []
    while i < len(lines) and lines[i].startswith("#"):
        meta.append(lines[i].rstrip("\n"))
        i += 1
    columns = lines[i].rstrip("\n").split(",")
    return meta, columns


def _meta_fields(meta_lines):
    return dict(line[1:].split("=", 1) for line in meta_lines)


# ---------------------------------------------------------------------------
# Leading '#' metadata block
# ---------------------------------------------------------------------------

def test_pipeline_writer_emits_the_metadata_header(tmp_path):
    out = str(tmp_path / "cr.csv")
    write_cellranger_probe_set(
        [_pair("GENE1")], out, overwrite=True,
        panel_name="flexiprobes custom probe set: demo",
        reference_genome="demo.fasta",
    )
    meta, columns = _read_header_and_columns(out)
    assert len(meta) == 5
    fields = _meta_fields(meta)
    assert fields["probe_set_file_format"] == CR_PROBE_SET_FILE_FORMAT
    assert fields["panel_type"] == CR_PANEL_TYPE == "custom"
    assert fields["panel_name"] == "flexiprobes custom probe set: demo"
    assert fields["reference_genome"] == "demo.fasta"
    assert fields["reference_version"]     # non-empty (tool-version derived)
    assert columns == ["gene_id", "probe_seq", "probe_id", "included"]


def test_pipeline_writer_header_has_sensible_fallbacks_when_unspecified(tmp_path):
    """Direct/test callers that don't pass panel_name/reference_genome still
    get a non-empty, honest header rather than a crash or a blank field."""
    out = str(tmp_path / "cr.csv")
    write_cellranger_probe_set([_pair("GENE1")], out, overwrite=True)
    fields = _meta_fields(_read_header_and_columns(out)[0])
    assert fields["panel_name"]
    assert fields["reference_genome"]


def test_per_run_export_emits_the_metadata_header(tmp_path):
    out = str(tmp_path / "run")
    _write_results(out, [_result_row("G", "G_t1")])
    _manifest(out, {"modality": "flex_v2", "fasta_path": "/inputs/my_target.fasta"})
    paths = export_selected_probes(out, {"G": ["G_t1"]}, overwrite=True)
    meta, columns = _read_header_and_columns(paths["cellranger"])
    fields = _meta_fields(meta)
    assert len(meta) == 5
    assert fields["probe_set_file_format"] == CR_PROBE_SET_FILE_FORMAT
    assert fields["panel_type"] == "custom"
    # Panel identity is sourced from the run's own FASTA input.
    assert "my_target" in fields["panel_name"]
    assert fields["reference_genome"] == "my_target"
    assert columns == ["gene_id", "probe_seq", "probe_id", "included"]


def test_bulk_export_emits_the_metadata_header(tmp_path):
    out_a = str(tmp_path / "virusA")
    _write_results(out_a, [_result_row("G1", "G1_t1")])
    _manifest(out_a, {"modality": "flex_v2"})
    dest = str(tmp_path / "bulk")
    res = export_project_bulk(
        [{"name": "virusA", "output_dir": out_a, "selections": {"G1": ["G1_t1"]}}],
        dest, project_name="myproj", overwrite=True,
    )
    meta, columns = _read_header_and_columns(res["paths"]["cellranger"])
    fields = _meta_fields(meta)
    assert len(meta) == 5
    assert "myproj" in fields["panel_name"]
    assert fields["reference_genome"] == "myproj"
    assert columns == ["gene_id", "probe_seq", "probe_id", "included"]


# ---------------------------------------------------------------------------
# The included column, always TRUE for a written (= included) probe
# ---------------------------------------------------------------------------

def test_included_column_present_and_true_for_every_row(tmp_path):
    out = str(tmp_path / "cr.csv")
    write_cellranger_probe_set([_pair("G1"), _pair("G2")], out, overwrite=True)
    with open(out) as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    assert len(rows) == 2
    assert all(r["included"] == "TRUE" for r in rows)


def test_per_run_and_bulk_export_included_column(tmp_path):
    out = str(tmp_path / "run")
    _write_results(out, [_result_row("G", "G_t1")])
    _manifest(out, {"modality": "flex_v1"})
    paths = export_selected_probes(out, {"G": ["G_t1"]}, overwrite=True)
    with open(paths["cellranger"]) as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    assert rows and all(r["included"] == "TRUE" for r in rows)

    dest = str(tmp_path / "bulk")
    res = export_project_bulk(
        [{"name": "run1", "output_dir": out, "selections": {"G": ["G_t1"]}}],
        dest, overwrite=True,
    )
    with open(res["paths"]["cellranger"]) as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    assert rows and all(r["included"] == "TRUE" for r in rows)


# ---------------------------------------------------------------------------
# The region column is never fabricated
# ---------------------------------------------------------------------------

def test_region_column_is_never_fabricated(tmp_path):
    """No exon/intron/gene model exists anywhere in this pipeline (confirmed
    by inspection: no GTF or annotation input path in app/src/flexprobe/), so
    'region' (spliced/unspliced) cannot be honestly populated for any window.
    It is omitted -- matching real vendor format 1.0, which has no region
    column at all -- rather than hard-coded to a value that would be wrong
    for effectively all of this tool's actual use cases (viral targets have
    no introns at all; gene-panel inputs are unknown without a GTF)."""
    for spatial in (False, True):
        out = str(tmp_path / f"cr_{spatial}.csv")
        write_cellranger_probe_set(
            [_pair("G1")], out, spatial=spatial, overwrite=True,
        )
        _, columns = _read_header_and_columns(out)
        assert "region" not in columns, f"spatial={spatial}"


# ---------------------------------------------------------------------------
# gene_id: warn on a transcript accession rather than pass it through
# silently as though it were a gene id.
# ---------------------------------------------------------------------------

def test_transcript_accession_gene_id_triggers_a_warning(tmp_path, caplog):
    out = str(tmp_path / "cr.csv")
    with caplog.at_level(
        logging.WARNING, logger="flexprobe.writers.cellranger_probe_set",
    ):
        write_cellranger_probe_set([_pair("ENST00000019317.8")], out, overwrite=True)
    assert any("transcript accession" in r.message for r in caplog.records)
    # The id is still written through unchanged -- a warning, not a rewrite,
    # since the pipeline has no gene model to translate it correctly.
    with open(out) as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    assert rows[0]["gene_id"] == "ENST00000019317.8"


@pytest.mark.parametrize("gene_id", ["ENSG00000074966", "OROV_BeAn19991", "TSPAN6"])
def test_plausible_gene_ids_do_not_warn(tmp_path, caplog, gene_id):
    out = str(tmp_path / "cr.csv")
    with caplog.at_level(
        logging.WARNING, logger="flexprobe.writers.cellranger_probe_set",
    ):
        write_cellranger_probe_set([_pair(gene_id)], out, overwrite=True)
    assert not any("transcript accession" in r.message for r in caplog.records)


def test_select_export_writers_also_warn_on_transcript_ids(tmp_path, caplog):
    out = str(tmp_path / "run")
    _write_results(out, [_result_row("ENST00000019317.8", "t1")])
    _manifest(out, {"modality": "flex_v1"})
    with caplog.at_level(
        logging.WARNING, logger="flexprobe.writers.cellranger_probe_set",
    ):
        export_selected_probes(out, {"ENST00000019317.8": ["t1"]}, overwrite=True)
    assert any("transcript accession" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# Structural diff against a real 10x vendor probe-set CSV, when available.
# ---------------------------------------------------------------------------

def test_structural_diff_against_real_vendor_csv(tmp_path, vendor_probe_set_csv):
    vendor_meta, vendor_columns = _read_header_and_columns(vendor_probe_set_csv)
    vendor_keys = set(_meta_fields(vendor_meta).keys())

    out = str(tmp_path / "cr.csv")
    write_cellranger_probe_set([_pair("G1")], out, overwrite=True)
    ours_meta, ours_columns = _read_header_and_columns(out)
    ours_keys = set(_meta_fields(ours_meta).keys())

    # Same metadata KEYS as the vendor file (values legitimately differ:
    # flexiprobes has no genome build to report).
    assert ours_keys == vendor_keys
    # Our column set matches vendor format 1.0 (the minimal variant, which is
    # what our column set -- no region/gene_name -- corresponds to) exactly.
    assert ours_columns == vendor_columns == [
        "gene_id", "probe_seq", "probe_id", "included",
    ]


# ---------------------------------------------------------------------------
# Shared fixtures (mirrors test_auto_select_and_export.py's helpers)
# ---------------------------------------------------------------------------

def _result_row(gene_id, probe_id):
    return {
        "rank": 1, "probe_id": probe_id, "gene_id": gene_id, "start": 0, "end": 50,
        "probe_rc_50bp": "A" * 50, "lhs_oligo": "LHS", "rhs_oligo": "/5Phos/RHS",
        "composite_score": 0.9,
    }


def _write_results(out, rows):
    os.makedirs(out, exist_ok=True)
    cols = ["rank", "probe_id", "gene_id", "start", "end", "mrna_target_50bp",
            "probe_rc_50bp", "target_LHS_25bp", "target_RHS_25bp",
            "lhs_oligo", "rhs_oligo", "composite_score"]
    with open(os.path.join(out, "results.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])


def _manifest(out, params):
    with open(os.path.join(out, "manifest.json"), "w") as fh:
        json.dump({"params": params}, fh)
