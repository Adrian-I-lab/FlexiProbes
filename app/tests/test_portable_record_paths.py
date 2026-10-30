"""The two engine record writers emit portable paths.

Absolute input, output and reference paths in run_config.json and
manifest.json would carry home-directory and cluster layouts into any shared
or published copy of a run. So absolute paths are written relative to the
directory holding the file, and the file names exactly which fields it
converted, so load_config can resolve them back: a written config still
reloads to the same files and reruns identically. The manifest also keeps the
input's file name beside its sha256, so the record stays identifiable when the
tree is moved. The non-engine writers follow the same rule.
"""
import dataclasses
import json
import os

import golden_scenarios as g
from flexprobe.config_io import config_to_dict, load_config, save_config
from flexprobe.pipeline import run_pipeline


def _strings(x):
    if isinstance(x, dict):
        for k, v in x.items():
            yield str(k)
            yield from _strings(v)
    elif isinstance(x, list):
        for v in x:
            yield from _strings(v)
    elif isinstance(x, str):
        yield x


def _run(tmp_path, name="kmer_screened"):
    work = tmp_path / "work"
    work.mkdir(exist_ok=True)
    cfg = g.SCENARIOS[name](str(work))
    run_pipeline(cfg)
    return cfg, str(work)


def test_no_absolute_path_in_either_record(tmp_path):
    cfg, work = _run(tmp_path)
    real = os.path.realpath(work)
    for f in ("manifest.json", "run_config.json"):
        with open(os.path.join(cfg.output_dir, f)) as fh:
            doc = json.load(fh)
        leaked = [s for s in _strings(doc) if work in s or real in s]
        assert not leaked, (f, leaked)


def test_the_written_config_reloads_and_reruns_identically(tmp_path):
    cfg, _work = _run(tmp_path)
    from flexprobe.config_io import RELATIVE_MARKER
    written = os.path.join(cfg.output_dir, "run_config.json")
    with open(written) as fh:
        assert RELATIVE_MARKER in json.load(fh)
    again = load_config(written)
    assert again == cfg
    rerun = dataclasses.replace(again, output_dir=str(tmp_path / "rerun"))
    run_pipeline(rerun)
    with open(os.path.join(cfg.output_dir, "results.csv"), "rb") as a, \
            open(os.path.join(rerun.output_dir, "results.csv"), "rb") as b:
        assert a.read() == b.read()


def test_the_manifest_keeps_name_and_hash(tmp_path):
    cfg, _work = _run(tmp_path)
    with open(os.path.join(cfg.output_dir, "manifest.json")) as fh:
        man = json.load(fh)
    assert man["inputs"]["fasta_name"] == os.path.basename(cfg.fasta_path)
    assert len(man["inputs"]["fasta_sha256"]) == 64
    # Reference identity is kept, keyed by the portable path.
    assert len(man["db_sha256s"]) == 1
    assert list(man["db_sha256s"].values())[0]
    assert list(man["db_sha256s"]) == list(man["reference_composition"])


def test_an_ordinary_save_is_unchanged(tmp_path):
    """save_config without portable=True keeps what it was given, no marker:
    the HPC bundle's config.json holds CLUSTER paths, written on another
    machine, and must not be rewritten against the local directory."""
    from flexprobe.config_io import RELATIVE_MARKER
    cfg = g.SCENARIOS["single_gene_v2"](str(tmp_path))
    p = save_config(cfg, str(tmp_path / "plain.json"))
    with open(p) as fh:
        data = json.load(fh)
    assert RELATIVE_MARKER not in data
    assert data == json.loads(json.dumps(config_to_dict(cfg)))


def test_relative_paths_keep_their_meaning(tmp_path, monkeypatch):
    """A hand-written relative path is cwd-relative; only fields the writer
    converted are resolved against the file's directory."""
    cfg = g.SCENARIOS["single_gene_v2"](str(tmp_path))
    cfg = dataclasses.replace(cfg, fasta_path="in.fa", output_dir="out")
    (tmp_path / "sub").mkdir()
    p = save_config(cfg, str(tmp_path / "sub" / "c.json"), portable=True)
    loaded = load_config(p)
    assert loaded.fasta_path == "in.fa" and loaded.output_dir == "out"


def test_the_prepare_only_records_are_portable_too(tmp_path):
    from flexprobe.pipeline import PREPARE_MANIFEST_FILE
    work = tmp_path / "work"
    work.mkdir()
    cfg = dataclasses.replace(g.SCENARIOS["kmer_screened"](str(work)),
                              prepare_only=True)
    run_pipeline(cfg)
    real = os.path.realpath(str(work))
    for f in (PREPARE_MANIFEST_FILE, "run_config.json"):
        with open(os.path.join(cfg.output_dir, f)) as fh:
            doc = json.load(fh)
        leaked = [s for s in _strings(doc) if str(work) in s or real in s]
        assert not leaked, (f, leaked)
