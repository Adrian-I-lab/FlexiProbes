"""A k-mer index records the builder code that wrote it.

``k``, ``format``, ``primary_only`` and ``source_sha256`` alone do not identify
an index: ``INDEX_FORMAT`` need not change when ``kmer_index.py`` does, so an
index built by older builder code would be reused without a word. So
``.meta.json`` also carries ``builder_sha256`` (a hash of ``kmer_index.py``
only, so a ranking change does not force a rebuild) plus the engine stamp, for
the record.

Three paths: no index builds; a matching stamp reuses; a missing or different
stamp raises ``KmerIndexStampMismatchError`` naming both stamps and deletes
nothing. ``rebuild=True`` (CLI ``--rebuild-kmer-index``) rebuilds on purpose.
"""
import json
import os

import pytest

from flexprobe import kmer_index
from flexprobe.errors import FlexprobeError, KmerIndexStampMismatchError
from flexprobe.kmer_index import KmerIndex

_SEQ = ("ACGTTGCAAGCTTGACCTAGGATCCGATTACAGGCATGCTTAAGCGGCCGCATCGATGCAT" * 20)


@pytest.fixture()
def ref(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">t1\n{_SEQ}\n>t2\n{_SEQ[::-1]}\n")
    return str(fasta)


def _meta(fasta, k=12):
    prefix = KmerIndex.cache_prefix(fasta, k)
    with open(prefix + ".meta.json") as fh:
        return prefix, json.load(fh)


def _no_build(monkeypatch):
    def boom(self, *a, **k):
        raise AssertionError("the index was rebuilt")
    monkeypatch.setattr(KmerIndex, "__init__", boom)


def test_builder_sha256_is_a_hash_of_the_builder_source():
    import hashlib
    src = open(kmer_index.__file__, "rb").read()
    assert kmer_index.builder_sha256() == hashlib.sha256(src).hexdigest()


def test_a_new_index_records_the_builder_and_the_engine_stamp(ref):
    from flexprobe.engine_identity import current_stamp
    KmerIndex.build_or_load(ref, 12)
    _, meta = _meta(ref)
    assert meta["builder_sha256"] == kmer_index.builder_sha256()
    assert meta["engine_stamp"] == current_stamp()


def test_a_matching_stamp_reuses_the_index(ref, monkeypatch):
    KmerIndex.build_or_load(ref, 12)
    _no_build(monkeypatch)
    idx, _ = KmerIndex.build_or_load(ref, 12)
    assert idx.k == 12


@pytest.mark.parametrize("stale", ["different", "missing"])
def test_a_stale_stamp_is_a_named_loud_error_and_deletes_nothing(
        ref, monkeypatch, stale):
    KmerIndex.build_or_load(ref, 12)
    prefix, meta = _meta(ref)
    if stale == "different":
        meta["builder_sha256"] = "0" * 64
    else:
        del meta["builder_sha256"]
    with open(prefix + ".meta.json", "w") as fh:
        json.dump(meta, fh)
    before = {p: os.path.getmtime(prefix + p)
              for p in (".pos.npy", ".off.npy", ".meta.json")}
    _no_build(monkeypatch)

    with pytest.raises(KmerIndexStampMismatchError) as ei:
        KmerIndex.build_or_load(ref, 12)

    msg = str(ei.value)
    assert isinstance(ei.value, FlexprobeError)
    assert "builder_sha256" in msg and "--rebuild-kmer-index" in msg
    assert kmer_index.builder_sha256() in msg              # current stamp
    if stale == "different":
        assert "0" * 64 in msg                             # the index's stamp
    assert before == {p: os.path.getmtime(prefix + p) for p in before}


def test_the_rebuild_flag_rebuilds_a_stale_index(ref):
    KmerIndex.build_or_load(ref, 12)
    prefix, meta = _meta(ref)
    meta["builder_sha256"] = "0" * 64
    with open(prefix + ".meta.json", "w") as fh:
        json.dump(meta, fh)

    KmerIndex.build_or_load(ref, 12, rebuild=True)

    _, meta = _meta(ref)
    assert meta["builder_sha256"] == kmer_index.builder_sha256()


def test_the_kmer_screen_passes_the_rebuild_flag_through(ref, monkeypatch):
    from flexprobe.screen import ENGINE_KMER, make_screen
    seen = {}
    real = KmerIndex.build_or_load.__func__

    def spy(cls, *a, **k):
        seen["rebuild"] = k.get("rebuild")
        return real(cls, *a, **k)
    monkeypatch.setattr(KmerIndex, "build_or_load", classmethod(spy))
    make_screen(ENGINE_KMER, rebuild_kmer_index=True).prepare([ref])
    assert seen["rebuild"] is True


def test_the_cli_flag_reaches_the_run_config(tmp_path, ref):
    from flexprobe.cli import CONFIG_OVERRIDES, _build_run_config, build_parser
    p = build_parser()
    base = ["design", "--fasta", ref, "--out", str(tmp_path / "o")]
    assert _build_run_config(p.parse_args(base)).screen.rebuild_kmer_index \
        is False
    args = p.parse_args(base + ["--rebuild-kmer-index"])
    assert _build_run_config(args).screen.rebuild_kmer_index is True
    assert "rebuild_kmer_index" in CONFIG_OVERRIDES
