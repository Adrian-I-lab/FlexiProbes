"""Pytest fixtures + path setup for the flexprobe test suite."""
import pathlib
import sys

import pytest

from _interpreter_guard import check as _check_interpreter

# Make ``src`` importable even when the package isn't pip-installed.
_SRC = pathlib.Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


def pytest_collection_modifyitems(config, items):
    """Stop a Qt-less interpreter from aborting the run mid-flight.

    See ``_interpreter_guard`` for why this is a behavioural probe and not a
    check on the interpreter's path or version.
    """
    _check_interpreter(
        items,
        pathlib.Path(__file__).resolve().parent,
        raise_usage_error=pytest.UsageError,
    )


@pytest.fixture(autouse=True)
def _isolated_flexprobe_home(tmp_path, monkeypatch):
    """Never let the suite touch a developer's real ``~/.flexprobe``.

    Every ``*_path()``/``_imported_root()`` helper across ``flexprobe_gui``
    (settings, project/run history, query drafts) and ``flexprobe.hpc_export``
    shares one fallback: ``os.environ.get("FLEXPROBE_HOME") or
    os.path.join(os.path.expanduser("~"), ".flexprobe")`` — and several of them
    call ``os.makedirs(...)`` as soon as they're invoked, not only when they're
    actually written to. A test that constructs a GUI store or calls
    ``import_hpc_output()`` with the default ``dest_root`` — even a test that
    only means to check that it *raises* — creates (or writes into) that real
    directory unless ``FLEXPROBE_HOME`` is already redirected. That happened in
    practice: ``tests/test_hpc_export.py`` did exactly this before the fix that
    added this fixture, and left a real ``~/.flexprobe/imported/`` behind on
    the machine this suite runs on.

    Autouse + function-scoped ``tmp_path``/``monkeypatch`` means every test
    gets its own throwaway redirect and it is undone automatically afterwards
    — including undoing a test that sets ``FLEXPROBE_HOME`` directly via
    ``os.environ[...]`` rather than through ``monkeypatch`` (monkeypatch still
    restores the value that was present when THIS fixture set it). A test that
    needs a specific store location keeps working unchanged: it either passes
    an explicit path (most of the suite already does) or overrides
    ``FLEXPROBE_HOME`` itself, same as before.
    """
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / ".flexprobe_home"))


@pytest.fixture
def fixtures_dir() -> pathlib.Path:
    return pathlib.Path(__file__).parent / "fixtures"


@pytest.fixture
def tmp_run_dir(tmp_path) -> pathlib.Path:
    d = tmp_path / "run"
    d.mkdir()
    return d


# A 200-bp synthetic sequence used by several pipeline tests.
SYNTH_200 = ("AGCTGTAGCTAGCTAGCATGCATGCATGCAGTCAGTCGATGCATCGTAGCATCGATCGTACTAGCTAGCAT"
             * 4)[:200]


@pytest.fixture
def synth_fasta(tmp_path) -> str:
    p = tmp_path / "synth.fa"
    p.write_text(f">SYN\n{SYNTH_200}\n")
    return str(p)


def write_fasta(path, records):
    """records: list of (id, seq)."""
    with open(path, "w") as fh:
        for sid, seq in records:
            fh.write(f">{sid}\n{seq}\n")
    return str(path)


# ---------------------------------------------------------------------------
# Helpers for the full-half off-target mismatch RECOUNT
# (specificity.full_half_mismatches). A hit dict alone is not enough to test it:
# the recount also needs the probe half itself and the subject bases opposite
# the positions BLAST never aligned.
# ---------------------------------------------------------------------------

class FakeSubjects:
    """An in-memory stand-in for ``blast.SubjectWindows``.

    Exposes the same ``window(sseqid, start, end)`` contract and the same
    ``n_exact`` / ``n_fallback`` counters, so the real counting and strand
    arithmetic are exercised with deterministic data and no ``blastdbcmd``
    subprocess. Returns None for an out-of-range request exactly as the real
    class does, so the conservative fallback is reachable from a test.

    Not a subclass: inheriting would drag in the subprocess plumbing this exists
    to avoid, and duck typing is all ``evaluate_specificity`` asks for.
    """

    def __init__(self, sequences=None, *, available=True):
        self.sequences = dict(sequences or {})
        self.available = available
        self.n_exact = 0
        self.n_fallback = 0
        self.n_calls = 0

    def window(self, sseqid, start, end):
        self.n_calls += 1
        seq = self.sequences.get(sseqid)
        if seq is None or start < 1 or end < start or end > len(seq):
            return None
        return seq[start - 1:end]


def aligned_hit(
    sseqid, query, subject, hsp_start, hsp_end,
    *, subject_offset=1, half=None, probe="p1",
):
    """A BLAST-shaped hit for an UNGAPPED plus-strand HSP.

    ``query`` and ``subject`` are equal-length strings in the probe's reading
    frame (the subject being what the DB really holds opposite the probe).
    ``hsp_start``/``hsp_end`` are the 1-based inclusive query positions BLAST
    chose to report — which is the whole point, since BLAST reports the
    maximum-scoring LOCAL segment and routinely trims the ends.
    ``subject_offset`` is the 1-based subject coordinate of query position 1.

    Every derived field (pident, length, mismatch, sstart, send, qseq, sseq) is
    computed from the two sequences, so a test cannot accidentally assert against
    a hit that is internally inconsistent.
    """
    qseq = query[hsp_start - 1:hsp_end]
    sseq = subject[hsp_start - 1:hsp_end]
    length = len(qseq)
    mismatch = sum(1 for a, b in zip(qseq, sseq) if a != b)
    return {
        "qseqid": f"{probe}__{half}" if half else probe,
        "sseqid": sseqid,
        "pident": round(100.0 * (length - mismatch) / length, 3) if length else 0.0,
        "length": length, "mismatch": mismatch, "gapopen": 0,
        "qstart": hsp_start, "qend": hsp_end,
        "sstart": subject_offset + hsp_start - 1,
        "send": subject_offset + hsp_end - 1,
        "evalue": 1e-5, "bitscore": 50.0,
        "qseq": qseq, "sseq": sseq,
    }
