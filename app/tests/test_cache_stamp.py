"""Resume checkpoints carry an engine stamp, and a mismatch stops the run.

Resume only when the engine that wrote a checkpoint is the engine running now.
A missing or different stamp is a named error that prints both stamps. The
checkpoint is never deleted and never silently recomputed. Only engine files
count (see ENGINE_FREEZE.json), so a writer or logging edit after the freeze
must not block a resume.

A cache keyed on inputs alone is not enough. An engine change can alter what
the screen returns while leaving the version string and ``CACHE_SCHEMA``
alone, and a checkpoint written before it would then be read back as if
nothing had changed.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import flexprobe
from flexprobe import candidate_cache, engine_identity, screen_cache
from flexprobe.errors import CacheStampMismatchError, FlexprobeError
from flexprobe.models import ProbeCandidate

PKG = Path(flexprobe.__file__).resolve().parent

STAMP_A = {"engine_version": "0.7.0", "engine_sha256": "a" * 64,
           "container_build_id": None, "tool_versions": {"python": "3.11.15"}}
STAMP_B = dict(STAMP_A, engine_sha256="b" * 64)


# --- which files are the engine -------------------------------------------

def _package_files():
    return {
        p.relative_to(PKG).as_posix()
        for p in PKG.rglob("*.py") if "__pycache__" not in p.parts
    }


def test_every_module_is_classified_exactly_once():
    files = _package_files()
    eng, non = engine_identity.ENGINE_FILES, engine_identity.NON_ENGINE_FILES
    assert not (eng & non), f"classified as both: {sorted(eng & non)}"
    unclassified = files - eng - non
    assert not unclassified, (
        f"new module(s) {sorted(unclassified)} are in neither ENGINE_FILES nor "
        "NON_ENGINE_FILES in flexprobe/engine_identity.py. Classify each one: an "
        "unlisted engine module would escape the resume stamp."
    )
    missing = (eng | non) - files
    assert not missing, f"listed but absent: {sorted(missing)}"


# --- engine_sha256 moves with engine code only -----------------------------

_PROBE = ("import json; from flexprobe import engine_identity as e; "
          "print(json.dumps(e.engine_sha256()))")


def _copy_pkg(dest: Path) -> Path:
    root = dest / "src"
    shutil.copytree(PKG, root / "flexprobe",
                    ignore=shutil.ignore_patterns("__pycache__"))
    return root


def _sha(src_root: Path) -> str:
    env = dict(os.environ, PYTHONPATH=str(src_root),
               PYTHONDONTWRITEBYTECODE="1")
    out = subprocess.run([sys.executable, "-c", _PROBE], env=env,
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out)


def _edited(tmp_path: Path, name: str, rel: str) -> str:
    root = _copy_pkg(tmp_path / name)
    f = root / "flexprobe" / rel
    f.write_text(f.read_text(encoding="utf-8") + "\n# edit\n", encoding="utf-8")
    return _sha(root)


def test_engine_sha256_is_stable_across_identical_copies(tmp_path):
    a = _sha(_copy_pkg(tmp_path / "a"))
    b = _sha(_copy_pkg(tmp_path / "b"))
    assert a == b and len(a) == 64


def test_engine_edit_changes_engine_sha256(tmp_path):
    base = _sha(_copy_pkg(tmp_path / "base"))
    assert _edited(tmp_path, "eng", "specificity.py") != base


def test_non_engine_edit_leaves_engine_sha256_alone(tmp_path):
    base = _sha(_copy_pkg(tmp_path / "base"))
    assert _edited(tmp_path, "non", "writers/manifest.py") == base


def test_current_stamp_fields():
    s = engine_identity.current_stamp()
    assert set(s) == {"engine_version", "engine_sha256",
                      "container_build_id", "tool_versions"}
    assert s["engine_version"] == flexprobe.__version__
    assert s["engine_sha256"] == engine_identity.engine_sha256()
    assert "blastn" not in s["tool_versions"]


def test_container_build_id_is_read_from_env(monkeypatch):
    monkeypatch.setenv(engine_identity.CONTAINER_BUILD_ID_ENV, "img-123")
    assert engine_identity.current_stamp()["container_build_id"] == "img-123"
    monkeypatch.delenv(engine_identity.CONTAINER_BUILD_ID_ENV)
    assert engine_identity.current_stamp()["container_build_id"] is None


# --- screen cache ------------------------------------------------------------

def _meta(cdir):
    with open(os.path.join(cdir, "_meta.json"), encoding="utf-8") as fh:
        return json.load(fh)


def test_screen_cache_writes_stamp_and_resumes_on_match(tmp_path):
    cdir = screen_cache.open_cache(str(tmp_path), "fp1", stamp=STAMP_A)
    assert _meta(cdir)["stamp"] == STAMP_A
    screen_cache.save_chunk(cdir, 0, {"p": [{"sseqid": "x"}]})
    again = screen_cache.open_cache(str(tmp_path), "fp1", stamp=STAMP_A)
    assert screen_cache.load_chunk(again, 0) == {"p": [{"sseqid": "x"}]}


def test_screen_cache_refuses_other_engine_and_keeps_it(tmp_path):
    cdir = screen_cache.open_cache(str(tmp_path), "fp1", stamp=STAMP_A)
    screen_cache.save_chunk(cdir, 0, {"p": [{"sseqid": "x"}]})
    with pytest.raises(CacheStampMismatchError) as ei:
        screen_cache.open_cache(str(tmp_path), "fp1", stamp=STAMP_B)
    msg = str(ei.value)
    assert "a" * 64 in msg and "b" * 64 in msg, "both stamps must be printed"
    assert isinstance(ei.value, FlexprobeError)   # CLI exits 2 on these
    # Not deleted, not rewritten.
    assert _meta(cdir)["stamp"] == STAMP_A
    assert screen_cache.load_chunk(cdir, 0) == {"p": [{"sseqid": "x"}]}


def test_screen_cache_refuses_unstamped_legacy_cache(tmp_path):
    cdir = screen_cache.cache_dir(str(tmp_path))
    os.makedirs(cdir)
    with open(os.path.join(cdir, "_meta.json"), "w") as fh:
        json.dump({"fingerprint": "fp1"}, fh)          # schema-5 shape
    with pytest.raises(CacheStampMismatchError) as ei:
        screen_cache.open_cache(str(tmp_path), "fp1", stamp=STAMP_A)
    assert ei.value.cached is None
    assert "carries no engine stamp" in str(ei.value)
    assert os.path.isdir(cdir)


def test_screen_cache_same_engine_new_inputs_starts_fresh(tmp_path):
    cdir = screen_cache.open_cache(str(tmp_path), "fp1", stamp=STAMP_A)
    screen_cache.save_chunk(cdir, 0, {"p": []})
    cdir2 = screen_cache.open_cache(str(tmp_path), "fp2", stamp=STAMP_A)
    assert screen_cache.load_chunk(cdir2, 0) is None
    assert _meta(cdir2) == {"fingerprint": "fp2", "stamp": STAMP_A}


# --- candidate cache ---------------------------------------------------------

def _cands():
    return [ProbeCandidate(gene_id="g1", probe_id="g1_0000",
                           sequence="ACGT" * 12 + "AC", start=0, end=50)]


def test_candidate_cache_round_trip_on_match(tmp_path):
    candidate_cache.save(str(tmp_path), "fp1", _cands(), stamp=STAMP_A)
    hit = candidate_cache.load(str(tmp_path), "fp1", stamp=STAMP_A)
    assert hit is not None and [c.probe_id for c in hit[0]] == ["g1_0000"]


def test_candidate_cache_refuses_other_engine_and_keeps_it(tmp_path):
    candidate_cache.save(str(tmp_path), "fp1", _cands(), stamp=STAMP_A)
    with pytest.raises(CacheStampMismatchError) as ei:
        candidate_cache.load(str(tmp_path), "fp1", stamp=STAMP_B)
    assert ei.value.cached == STAMP_A and ei.value.current == STAMP_B
    assert candidate_cache.load(str(tmp_path), "fp1", stamp=STAMP_A) is not None


def test_candidate_cache_refuses_unstamped_legacy_checkpoint(tmp_path):
    candidate_cache.save(str(tmp_path), "fp1", _cands(), stamp=STAMP_A)
    meta_path = os.path.join(candidate_cache.cache_dir(str(tmp_path)),
                             "_meta.json")
    meta = json.load(open(meta_path))
    del meta["stamp"]                                  # schema-2 shape
    json.dump(meta, open(meta_path, "w"))
    with pytest.raises(CacheStampMismatchError):
        candidate_cache.load(str(tmp_path), "fp1", stamp=STAMP_A)


# --- through the pipeline ----------------------------------------------------

def _candidate(i):
    return ProbeCandidate(gene_id="g1", probe_id=f"g1_{i:04d}",
                          sequence="ACGT" * 12 + "AC", start=i, end=i + 50)


def test_pipeline_screen_resume_under_other_engine_raises(tmp_path, monkeypatch):
    import flexprobe.pipeline as pl
    from flexprobe.screen import ENGINE_KMER

    ref = tmp_path / "reference.fa"
    ref.write_text(">refA\n" + ("ACGTTGCA" * 60) + "\n")
    out = tmp_path / "out"
    out.mkdir()
    args = dict(min_pident=80.0, min_aln=20, resume=True, engine=ENGINE_KMER)

    pl._run_screen([_candidate(0)], [str(ref)], str(out), **args)
    cdir = screen_cache.cache_dir(str(out))
    written = _meta(cdir)["stamp"]
    assert written == engine_identity.current_stamp()   # k-mer: no blastn

    # Same inputs, same engine: resumes.
    pl._run_screen([_candidate(0)], [str(ref)], str(out), **args)

    # Same inputs, different engine code: stops, and keeps the checkpoint.
    real = engine_identity.current_stamp
    monkeypatch.setattr(engine_identity, "current_stamp",
                        lambda **kw: dict(real(**kw), engine_sha256="0" * 64))
    with pytest.raises(CacheStampMismatchError):
        pl._run_screen([_candidate(0)], [str(ref)], str(out), **args)
    assert _meta(cdir)["stamp"] == written
