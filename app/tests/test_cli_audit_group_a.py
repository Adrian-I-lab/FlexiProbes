"""CLI help text, group A: each flag's help states what the code does.
Each check pins the current fact and, where one existed, the absence of an
inaccurate older wording. The spacing-ladder flags are covered by
test_cli_audit_group_b.py.
"""
import pytest

from flexprobe import cli


def _sub(command):
    return cli.build_parser()._subparsers._group_actions[0].choices[command]


def _help(dest, command="design"):
    for a in _sub(command)._actions:
        if a.dest == dest and a.help:
            return " ".join(a.help.split())
    raise AssertionError(f"no help for {dest} in {command}")


def _help_of_option(option, command="design"):
    for a in _sub(command)._actions:
        if option in a.option_strings:
            return " ".join((a.help or "").split())
    raise AssertionError(option)


CASES = [
    # (id, option, command, must contain, must not contain)
    ("A2", "--no-off-target-same-subject", "design", "17.6%", "17.7%"),
    ("A3", "--no-resume", "design", ".candidate_cache", "resumable BLAST"),
    ("A4", "--screen-top-per-gene", "design", "either engine", "Cuts BLAST"),
    ("A5", "--screen-cache-mode", "design", "screen checkpoint", "each BLAST checkpoint"),
    ("A6", "--null-model", "design", "--null-model-shuffles", "extra BLAST queries"),
    ("A7", "--emit-prescreen-candidates", "design", "prescreen_index.json", "BLAST actually screened"),
    ("A8", "--keep-candidate-cache", "design", "--skip-screen", None),
    ("A9", "--group-by-gene", "design", "Ignored with --multi-strain", "BLAST DB"),
    ("A10", "--min-isoform-coverage", "design", "Ignored with --multi-strain", None),
    ("A11", "--strategy", "design", "--conserved-max-degenerate-per-half", "conserved (all records)"),
    ("A14", "--prepare-only", "design", "Refused with --no-resume", None),
    ("A15", "--enforce-gc-range", "design", "no effect without this flag", None),
    ("A16", "--no-primer3", "design", "--tm-min/--tm-max gate still applies", None),
    ("A17", "--screen-workers", "design", "usable cores", None),
    ("A18", "--singleplex", "design", "Flex v2 only", "Relevant for Flex v1/v2"),
    ("A19", "--cached-hits", "design", "nothing in this package writes", "previously saved BLAST-hits"),
    ("A20", "--off-target-ref", "design", "skipped with a warning", "Omit it together with --skip-screen"),
    ("A21", "--off-target-same-subject", "design", "--off-target-either-half", None),
    ("A22", "--no-require-tn-junction", "design", "explicit --tn-junction-mode", None),
    ("A23", "--dest", "export-hpc", "empty directory", "an existing bundle is overwritten"),
    ("A24", "--remote-db-path", "export-hpc", "FASTA", "BLAST DB path(s)"),
    ("A29", "--image", "export-hpc", "relative path", None),
    ("A29b", "--cores", "export-hpc", "--screen-workers", None),
    ("A30", "--modality", "design", "flex_v2", None),
    ("A30b", "--tile-step", "design", "bp", None),
    ("A30c", "--tile-size", "design", "50", None),
]


@pytest.mark.parametrize("cid,option,command,present,absent", CASES,
                         ids=[c[0] for c in CASES])
def test_help_is_corrected(cid, option, command, present, absent):
    text = _help_of_option(option, command)
    assert present in text, (cid, text)
    if absent:
        assert absent not in text, (cid, text)


def test_a1_top_level_description():
    d = cli.build_parser().description
    assert "BLAST" not in d.split("(")[1] or "k-mer" in d
    assert "off-target screen" in d


def test_a28_export_hpc_description():
    d = " ".join(_sub("export-hpc").description.split())
    assert "README.txt" in d and "--split-prepare-screen" in d


def test_a20_the_runtime_warning_names_the_real_flag():
    import inspect
    import flexprobe.pipeline as pl
    src = inspect.getsource(pl)
    assert "or pass --blast-db to enable" not in src
    assert '"--off-target-ref to enable real off-target screening."' in src
