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

"""BLAST integration: discovery, DB build, search, and TSV parsing.

NO JUNCTION-DUPLEX RECORDS COME OUT OF THIS MODULE, AND NONE CAN.
-----------------------------------------------------------------
The junction-duplex gate (``specificity.junction_clean_duplex``, default 9 bp)
rejects a pair when some non-target gene presents that many perfectly matched
bases on BOTH sides of the ligation nick in contiguous reference sequence. The
k-mer engine answers it by seeding on the nick itself
(``screen.KmerScreen._junction_scan``); this module cannot, for a reason that is
about BLAST's output rather than about effort:

* THE SITES ARE ABSENT FROM THE HIT LIST, not under-described in it. The shape
  the gate exists to catch is a non-target gene matching perfectly across the
  nick with every mismatch pushed out to the oligo ends — 6+ mismatches in each
  25 nt half. blastn-short reports no alignment there at all, at any
  ``-outfmt``. The alignment strings (``qseq``/``sseq``) this module already
  requests give exact per-position mismatches for the hits that ARE reported,
  so post-processing them would return "clean" for precisely the probes that
  are not.
* A SEPARATE 12-mer SEARCH WOULD NOT RESTORE COMPLETENESS. blastn's word
  extension, ``-evalue`` filtering and ``-max_target_seqs`` cap are heuristic;
  the gate's guarantee is that no qualifying site can be missed, and a
  heuristic search cannot supply it. It would also cost a second full pass over
  the reference.

So the combination is REFUSED rather than silently degraded:
``models.RunConfig.validate`` raises when ``--screen-engine blast`` is paired
with ``off_target_junction_duplex_bp > 0``, telling the caller to switch engines
or set the option to 0. A run that reports the gate at 9 while nothing measured
it is the failure being prevented — the evidence columns
(``offtarget_clean_left_bp`` / ``_right_bp``) would sit at -1 and the gate would
never fire.
"""
from __future__ import annotations
import csv
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
from collections import OrderedDict
from typing import List, Dict, Any, Iterable, Optional, Sequence, Tuple

from .constants import (BLAST_MAX_TARGET_SEQS, BLASTN_PENALTY,
                        BLASTN_REWARD, MIN_BLAST_VERSION)
from .errors import BlastError, BlastVersionError, BlastDatabaseError

logger = logging.getLogger("flexprobe.blast")

MIN_VERSION = MIN_BLAST_VERSION

BLASTN_OUTFMT_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore",
    # The aligned segments themselves, for the full-half mismatch RECOUNT
    # (specificity.full_half_mismatches): mismatches inside the HSP are counted
    # column-by-column from these strings rather than derived from a rounded
    # ``pident``, and ``qstart``/``qend`` then say exactly which query positions
    # were never compared at all. Deriving the count from pident*length would
    # round silently and, far worse, treat every unaligned query position as a
    # mismatch, certifying probes clean against a close relative they can
    # actually ligate on.
    # NOTE: these columns set the on-disk TSV column count. Any change here
    # must bump screen_cache.CACHE_SCHEMA in step, so a cache written in another
    # column format is never re-read (the fingerprint won't match and the stale
    # dir is wiped, rather than the columns being silently mis-parsed).
    "qseq", "sseq",
]

_COL_TYPES = {
    "pident": float, "length": int, "mismatch": int,
    "gapopen": int, "qstart": int, "qend": int,
    "sstart": int, "send": int, "evalue": float, "bitscore": float,
    # qseq/sseq are deliberately absent: they stay strings (csv already yields
    # str), and adding them here would only invite a cast that can't fail
    # informatively.
}


def _blast_env() -> Dict[str, str]:
    """Environment for BLAST subprocesses, with NCBI usage reporting disabled.

    NCBI BLAST+ phones home to ``*.ncbi.nlm.nih.gov`` with anonymous usage
    statistics **by default**. Because flexiprobes *bundles* BLAST+ and
    invokes it automatically, we suppress that network call so the app makes no
    silent outbound connection. ``setdefault`` is used so a user/institution that
    deliberately exports ``BLAST_USAGE_REPORT=true`` keeps their choice.
    """
    env = dict(os.environ)
    env.setdefault("BLAST_USAGE_REPORT", "false")
    env.setdefault("NCBI_USAGE_REPORT_ENABLED", "false")
    return env


def _parse_version(version_str: str) -> Tuple[int, ...]:
    m = re.search(r'(\d+\.\d+\.\d+)', version_str)
    if not m:
        raise BlastVersionError(f"Cannot parse BLAST version from: {version_str!r}")
    return tuple(int(x) for x in m.group(1).split('.'))


def _is_windows() -> bool:
    """True on Windows. A function, not a constant, so tests can simulate it."""
    return os.name == "nt"


# Common install locations to search when a tool isn't on PATH. This matters
# for the packaged GUI: a macOS .app launched from Finder/Dock inherits only a
# minimal PATH (/usr/bin:/bin:/usr/sbin:/sbin) — NOT the user's shell PATH — so
# Homebrew's /opt/homebrew/bin (and similar) are invisible and blastn looks
# "missing" even though it's installed.
_COMMON_BLAST_DIRS = (
    "/opt/homebrew/bin",     # Homebrew, Apple Silicon
    "/usr/local/bin",        # Homebrew, Intel / manual installs
    "/opt/local/bin",        # MacPorts
    "/usr/bin", "/bin",
    "/opt/conda/bin",
)


def _bundled_blast_dir() -> Optional[str]:
    """Return the dir of BLAST binaries shipped *inside* the packaged app.

    The desktop build vendors ``blastn`` + ``makeblastdb`` (see
    ``packaging/vendor/blast`` and the PyInstaller spec) so an end user needs no
    separate BLAST+ install. Only meaningful in a frozen (PyInstaller) app;
    returns ``None`` when running from source. Checks ``sys._MEIPASS`` plus the
    macOS ``.app`` layouts relative to the launcher so it works regardless of
    where PyInstaller places bundled data.
    """
    if not getattr(sys, "frozen", False):
        return None
    candidates: List[str] = []
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        candidates.append(os.path.join(meipass, "blast", "bin"))
    exe_dir = os.path.dirname(os.path.realpath(sys.executable))
    candidates += [
        os.path.join(exe_dir, "blast", "bin"),                       # onedir
        os.path.join(exe_dir, "..", "Frameworks", "blast", "bin"),   # .app
        os.path.join(exe_dir, "..", "Resources", "blast", "bin"),    # .app
    ]
    for d in candidates:
        d = os.path.normpath(d)
        if os.path.isdir(d):
            return d
    return None


def _blast_search_dirs() -> List[str]:
    """Directories to look in for BLAST binaries, in priority order."""
    dirs: List[str] = []
    env = os.environ.get("FLEXPROBE_BLAST_DIR")
    if env:
        dirs.append(env)
    # The active Python's / conda env's bin (covers `pip`/conda blast installs).
    for base in (sys.prefix, os.environ.get("CONDA_PREFIX")):
        if base:
            dirs.append(os.path.join(base, "bin"))
    dirs.extend(_COMMON_BLAST_DIRS)
    # De-dup, preserving order.
    seen: set = set()
    out: List[str] = []
    for d in dirs:
        if d and d not in seen:
            seen.add(d)
            out.append(d)
    return out


def _executable_names(name: str) -> List[str]:
    """The file names ``name`` can have on disk on this platform.

    On Windows an executable carries an extension: the vendored tree ships
    ``blastn.exe``, never a bare ``blastn``. Joining the bare name onto a
    directory finds nothing there, which would make the bundled BLAST+,
    ``$FLEXPROBE_BLAST_DIR`` and the fallback dirs all invisible on Windows. ``shutil.which`` already applies ``PATHEXT`` to PATH;
    this does the same for an explicit directory. ``PATHEXT`` is split on ``;``
    rather than ``os.pathsep`` so the Windows branch can be exercised off
    Windows, and lower-cased so it also matches on a case-sensitive filesystem.
    """
    if not _is_windows() or os.path.splitext(name)[1]:
        return [name]
    exts = os.environ.get("PATHEXT") or ".COM;.EXE;.BAT;.CMD"
    return [name + e.lower() for e in exts.split(";") if e]


def _tool_in_dir(directory: str, name: str) -> Optional[str]:
    """Path of executable ``name`` inside ``directory``, or None."""
    for fname in _executable_names(name):
        cand = os.path.join(directory, fname)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None


def _find_tool(name: str) -> Optional[str]:
    """Locate a BLAST executable.

    Priority: ``$FLEXPROBE_BLAST_DIR`` (explicit override) → the copy bundled
    inside the packaged app → PATH → common install dirs. Bundled-before-PATH so
    the packaged app uses its own pinned, known-good BLAST+ rather than whatever
    happens to be on the machine.
    """
    env = os.environ.get("FLEXPROBE_BLAST_DIR")
    if env:
        cand = _tool_in_dir(env, name)
        if cand:
            return cand
    bundled = _bundled_blast_dir()
    if bundled:
        cand = _tool_in_dir(bundled, name)
        if cand:
            return cand
    path = shutil.which(name)
    if path:
        return path
    for d in _blast_search_dirs():
        cand = _tool_in_dir(d, name)
        if cand:
            return cand
    return None


def find_blastdbcmd() -> Optional[str]:
    """Locate ``blastdbcmd``, or None if it isn't installed.

    Deliberately NOT part of :func:`discover_blast`, which raises when a tool is
    missing. ``blastdbcmd`` is needed only to retrieve the subject bases opposite
    the *unaligned* part of a probe half (see ``specificity.full_half_mismatches``
    and :class:`SubjectWindows`). Its absence must degrade the screen to the
    conservative count, never abort a run: the fallback treats unaligned query
    positions as MATCHES, which can only over-reject, whereas failing the run
    would deny the user a design entirely.
    """
    return _find_tool("blastdbcmd")


def _missing_blast_message(missing: Sequence[str]) -> str:
    """The "not found" error, with install advice for THIS platform.

    A Windows user told to run ``brew`` or to look in ``/opt/homebrew/bin`` has
    been given nothing they can act on.
    """
    head = (f"{', '.join(missing)} not found. The packaged desktop app ships "
            "its own BLAST+ (blastn + makeblastdb) and finds it automatically, "
            "so this normally only happens when running from source. ")
    if _is_windows():
        return head + (
            "Install NCBI BLAST+ for Windows (the ncbi-blast-*-win64 installer "
            "from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
            "). If it IS installed but isn't on PATH, set the "
            "FLEXPROBE_BLAST_DIR environment variable to the folder containing "
            "blastn.exe (e.g. C:\\Program Files\\NCBI\\blast-2.17.0+\\bin).")
    return head + (
        "Install ncbi-blast+ (macOS: `brew install blast`; Debian/Ubuntu: "
        "`apt-get install ncbi-blast+`). If it IS installed but isn't on "
        "PATH, set the FLEXPROBE_BLAST_DIR environment variable to the "
        "folder containing blastn (e.g. /opt/homebrew/bin).")


def discover_blast() -> Dict[str, str]:
    """Find blastn and makeblastdb; return {tool: path}.

    Searches PATH, then ``$FLEXPROBE_BLAST_DIR``, the active env's bin, and
    common install locations — so the binaries are found even when the (Finder-
    launched) GUI doesn't inherit the user's shell PATH.

    ``blastdbcmd`` is looked up too but is OPTIONAL — it is added to the returned
    dict only when present, because a missing copy costs exactness (conservative
    recount), not correctness. Use :func:`find_blastdbcmd` when that is all you
    need.
    """
    tools = {}
    missing = []
    for name in ("blastn", "makeblastdb"):
        path = _find_tool(name)
        if path is None:
            missing.append(name)
        else:
            tools[name] = path
    if missing:
        raise BlastError(_missing_blast_message(missing))
    dbcmd = find_blastdbcmd()
    if dbcmd:
        tools["blastdbcmd"] = dbcmd
    return tools


def check_version(blastn_path: str, min_version: Tuple[int, ...] = MIN_VERSION) -> str:
    """Verify blastn version >= min_version. Returns version string."""
    result = subprocess.run(
        [blastn_path, "-version"],
        capture_output=True, text=True, timeout=30, env=_blast_env(),
    )
    ver_str = result.stdout.strip()
    ver = _parse_version(ver_str)
    if ver < min_version:
        raise BlastVersionError(
            f"BLAST {'.'.join(map(str, ver))} < required {'.'.join(map(str, min_version))}"
        )
    logger.info("BLAST version %s OK", '.'.join(map(str, ver)))
    return '.'.join(map(str, ver))


def blastn_version(blastn_path: Optional[str] = None) -> Optional[str]:
    """Best-effort ``blastn`` version for provenance, or None if unavailable.

    Unlike :func:`check_version` this never raises and never enforces a minimum
    — it is called when assembling the run manifest, where the goal is to record
    what actually ran (including a too-old BLAST) rather than to gate the run.
    Discovers ``blastn`` itself when no path is given.
    """
    try:
        if blastn_path is None:
            blastn_path = _find_tool("blastn")
        if not blastn_path:
            return None
        result = subprocess.run(
            [blastn_path, "-version"],
            capture_output=True, text=True, timeout=30, env=_blast_env(),
        )
        if result.returncode != 0:
            # A blastn that runs but fails would otherwise feed empty/garbage
            # stdout to the parser; better to report "unknown" than a misparse.
            logger.debug("blastn -version exited %s", result.returncode)
            return None
        return '.'.join(map(str, _parse_version(result.stdout)))
    except Exception:            # missing / unparseable / timed out
        logger.debug("Could not determine blastn version", exc_info=True)
        return None


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


# Bumped whenever the makeblastdb COMMAND LINE changes, and folded into the
# .sha256 sidecar so a DB built by an older recipe is rebuilt rather than reused.
# The sidecar is compared verbatim, so changing what we write is itself the
# invalidation.
#   1 = plain `makeblastdb -dbtype nucl`
#   2 = adds -parse_seqids, so blastdbcmd can retrieve by accession
_DB_BUILD_RECIPE = 3   # 3: -parse_seqids is mandatory (no fallback build)


def build_db(
    makeblastdb_path: str,
    fasta_path: str,
    db_prefix: str,
    title: str = "flexprobe_db",
) -> str:
    """Build a BLAST nucleotide DB. Returns db_prefix.

    Uses checksum caching: if a .sha256 sidecar matches the input FASTA *and*
    the build recipe, the existing DB is reused.

    ``-parse_seqids`` is requested because without it the DB carries no accession
    index and ``blastdbcmd -entry <sseqid>`` fails outright with "DB contains no
    accession info".
    That lookup is how the full-half mismatch recount retrieves the subject bases
    opposite the unaligned part of a probe half, so without it every off-target
    hit falls back to the conservative count and the screen over-rejects.

    There is NO fallback build without the flag. makeblastdb refuses to parse
    seqids when the FASTA contains duplicate ids. Retrying without the flag
    would produce a database that looks fine and silently makes every
    off-target count conservative, which discards a large share of usable
    probes. A DB that cannot be recounted against is therefore not a degraded
    option but an unusable one, and a duplicate id is a hard error naming the
    fix.
    """
    sha_path = db_prefix + ".sha256"
    current_sha = f"{_sha256_file(fasta_path)}\trecipe={_DB_BUILD_RECIPE}"
    if os.path.exists(sha_path):
        with open(sha_path) as f:
            cached_sha = f.read().strip()
        if cached_sha == current_sha:
            logger.info("BLAST DB cache hit for %s", db_prefix)
            return db_prefix

    base_cmd = [
        makeblastdb_path,
        "-in", fasta_path,
        "-dbtype", "nucl",
        "-out", db_prefix,
        "-title", title,
    ]
    cmd = base_cmd + ["-parse_seqids"]
    logger.info("Building BLAST DB: %s", ' '.join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                            env=_blast_env())
    if result.returncode != 0:
        # There is deliberately NO retry without -parse_seqids. A database
        # without an accession index cannot answer blastdbcmd -entry, and the
        # off-target recount reads the subject bases opposite the unaligned part
        # of a probe half through exactly that call. Without it the screen has no
        # way to count and falls back to treating unexamined positions as
        # matches, which discards a large share of usable probes for no
        # reason. A DB that cannot be recounted against is not a
        # cheaper option; it is an unusable one, so this fails loudly here rather
        # than silently degrading every verdict downstream.
        raise BlastDatabaseError(
            f"makeblastdb could not build {db_prefix} with -parse_seqids.\n"
            f"{result.stderr.strip()}\n\n"
            "The usual cause is duplicate sequence ids in the FASTA — "
            "makeblastdb refuses to index accessions when ids repeat. Make the "
            "ids unique and retry; for example:\n"
            "  awk '/^>/{n[$1]++; if(n[$1]>1) $1=$1\"_\"n[$1]} 1' in.fa > uniq.fa\n"
            "The accession index is required: the off-target recount retrieves "
            "subject bases with `blastdbcmd -entry`, and a DB without it cannot "
            "be screened against."
        )
    with open(sha_path, 'w') as f:
        f.write(current_sha)
    return db_prefix


def run_blastn(
    blastn_path: str,
    query_fasta: str,
    db_prefix: str,
    evalue: float = 1e9,
    word_size: int = 5,
    num_threads: int = 1,
    perc_identity: Optional[float] = None,
    qcov_hsp_perc: Optional[float] = None,
    max_hsps: Optional[int] = None,
    max_target_seqs: int = BLAST_MAX_TARGET_SEQS,
    strand: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
) -> str:
    """Run blastn-short and return path to TSV output.

    ``perc_identity`` / ``qcov_hsp_perc`` / ``max_hsps`` push the downstream
    identity / length / redundancy filters into BLAST itself so it never emits
    the (potentially enormous) volume of sub-threshold short seeds the caller
    would only discard. They must be set no stricter than the Python post-filter
    so BLAST never drops a hit that would have counted.

    ``max_target_seqs`` is set explicitly (not left at BLAST's default of 500)
    so the off-target screen doesn't silently truncate a query's hit list and
    let a non-specific probe pass — see :data:`BLAST_MAX_TARGET_SEQS`.

    ``strand`` is passed straight through to blastn's own ``-strand`` argument
    (``"both"``/``"plus"``/``"minus"``); ``None`` omits the flag, which is
    blastn's ``-strand both`` default. See
    :data:`flexprobe.constants.REFERENCE_STRAND_DEFAULT` for why the caller
    must not default this to ``"plus"`` on its own — a probe genuinely cannot
    bind a reference's minus strand ONLY when that reference is sense-oriented
    cDNA, which is not true of every reference this tool accepts (a genome, or
    a negative-sense viral FASTA).
    """
    out_tsv = query_fasta + ".blast_hits.tsv"
    cmd = [
        blastn_path,
        "-task", "blastn-short",
        "-query", query_fasta,
        "-db", db_prefix,
        "-evalue", str(evalue),
        "-word_size", str(word_size),
        "-num_threads", str(num_threads),
        "-max_target_seqs", str(max_target_seqs),
        # THE COMPLETE-COVERAGE CONFIGURATION. Every flag below, together with
        # the permissive -evalue and small -word_size defaults above, exists to
        # stop BLAST reporting a TRIMMED local alignment or none at all,
        # because the specificity criterion is defined across the whole 25 nt
        # half and a trimmed alignment cannot answer it. With blastn-short's
        # stock settings, disqualifying off-targets that an exhaustive search
        # finds can go unreported.
        #
        # The cost is real: this search is far slower than the k-mer engine.
        # That trade is taken deliberately. BLAST is the fallback for
        # references the k-mer engine cannot process, so its user has no
        # faster correct option, and a fast wrong answer is worse than a slow
        # right one.
        #
        # -reward/-penalty: 1 and -1 rather than blastn-short's 1 and -3. A
        # 25-mer four mismatches from a subject then scores 21 rather than 9,
        # so the full-length alignment outscores short exact cores and gets
        # reported instead of them. constants.BLASTN_SHORT_MISMATCH_COST is
        # reward + |penalty| and MUST be changed with these two.
        "-reward", str(BLASTN_REWARD),
        "-penalty", str(BLASTN_PENALTY),
        # -ungapped: the criterion counts substitutions across a fixed-length
        # half. A gapped alignment answers a different question, and gaps
        # only ever widen the pushdown margin, so removing them is safe.
        "-ungapped",
        # -dust / -soft_masking: low-complexity masking removes exactly the
        # repetitive sites an off-target screen most needs to see.
        "-dust", "no",
        "-soft_masking", "false",
        "-outfmt", "6 " + " ".join(BLASTN_OUTFMT_COLS),
        "-out", out_tsv,
    ]
    if perc_identity is not None:
        cmd += ["-perc_identity", str(perc_identity)]
    if qcov_hsp_perc is not None:
        cmd += ["-qcov_hsp_perc", str(qcov_hsp_perc)]
    if max_hsps is not None:
        cmd += ["-max_hsps", str(max_hsps)]
    if strand is not None:
        cmd += ["-strand", strand]
    if extra_args:
        cmd.extend(extra_args)
    logger.info("Running blastn: %s", ' '.join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600,
                            env=_blast_env())
    if result.returncode != 0:
        raise BlastError(f"blastn failed: {result.stderr}")
    return out_tsv


def parse_hits_tsv(tsv_path: str) -> List[Dict[str, Any]]:
    """Parse blastn tab-separated output with explicit dtype casting."""
    hits: List[Dict[str, Any]] = []
    if not os.path.exists(tsv_path):
        return hits
    with open(tsv_path, 'r') as fh:
        reader = csv.DictReader(fh, fieldnames=BLASTN_OUTFMT_COLS, delimiter='\t')
        for row in reader:
            typed_row: Dict[str, Any] = {}
            for col, val in row.items():
                if col in _COL_TYPES:
                    # A truncated or malformed BLAST row (e.g. an interrupted
                    # run, or a non-numeric value) must surface as a clear,
                    # located error rather than a raw ValueError/TypeError.
                    try:
                        typed_row[col] = _COL_TYPES[col](val)
                    except (TypeError, ValueError) as exc:
                        raise BlastError(
                            f"Malformed BLAST output in {tsv_path} at line "
                            f"{reader.line_num}: column {col!r} has a "
                            f"non-numeric value {val!r}. The BLAST run may have "
                            f"been interrupted; delete the file and re-run."
                        ) from exc
                else:
                    typed_row[col] = val
            hits.append(typed_row)
    return hits


def load_cached_hits(path: str) -> Any:
    """Load cached BLAST hits from a JSON file.

    The expected structure is the same dict the pipeline consumes:
    ``{probe_id: [hit, ...]}``. Typed ``Any`` because it's ``json.load``
    and the pipeline keys it by probe_id.
    """
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Convenience: accept a raw FASTA where a DB prefix is expected.
# ---------------------------------------------------------------------------

_FASTA_EXTS = (".fa", ".fasta", ".fna", ".ffn", ".frn")


def is_db_prefix(path: str) -> bool:
    """Return True if ``path`` looks like an existing BLAST DB prefix."""
    return any(os.path.exists(path + ext) for ext in (".nhr", ".nin", ".nsq"))


def resolve_db_prefix(path: str) -> Optional[str]:
    """The BLAST DB prefix ``path`` refers to, or None if no DB exists there.

    Read-only counterpart of :func:`ensure_db_for_path`, using the same prefix
    derivation (a FASTA's sidecars sit beside it with the extension stripped)
    but never invoking makeblastdb. Callers that need to *describe* the DB — the
    manifest's size and checksum fields — must resolve exactly the way the
    search did, or they silently describe nothing: globbing ``<prefix>.*`` on a
    configured path like ``…cdna.all.fa`` finds no sidecars, because they are
    named ``…cdna.all.nin``. That mismatch would make every manifest record a
    0-byte database with no checksum.
    """
    if is_db_prefix(path):
        return path
    lower = path.lower()
    for ext in _FASTA_EXTS:
        if lower.endswith(ext):
            cand = path[: -len(ext)]
            return cand if is_db_prefix(cand) else None
    return None


def ensure_db_for_path(
    path: str,
    *,
    cache_dir: Optional[str] = None,
    makeblastdb_path: Optional[str] = None,
) -> str:
    """Return a usable BLAST DB prefix for ``path``.

    Behaviour:
      • If ``path`` is already a DB prefix (its sidecars exist), it's
        returned unchanged.
      • If ``path`` is a FASTA, a DB is built next to it (or under
        ``cache_dir`` if given) and the prefix is returned. The build
        uses SHA-256 caching so re-runs with the same FASTA are instant.

    Lets the GUI accept either DB prefixes or raw FASTAs in the BLAST
    field. Raises BlastError if neither form matches.
    """
    if is_db_prefix(path):
        return path
    lower = path.lower()
    if any(lower.endswith(ext) for ext in _FASTA_EXTS):
        if not os.path.exists(path):
            raise BlastError(f"FASTA does not exist: {path}")
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
            base = os.path.basename(path)
            for ext in _FASTA_EXTS:
                if base.lower().endswith(ext):
                    base = base[: -len(ext)]
                    break
            db_prefix = os.path.join(cache_dir, base)
        else:
            # Strip the FASTA extension so the DB sidecars sit next to the
            # FASTA without colliding with it.
            db_prefix = path
            for ext in _FASTA_EXTS:
                if path.lower().endswith(ext):
                    db_prefix = path[: -len(ext)]
                    break
        # Always route through build_db, even if a DB already sits at this
        # prefix: build_db is what compares the FASTA's current content hash
        # against the .sha256 sidecar, so it's the only thing that notices
        # when the FASTA at `path` was overwritten since the DB was built.
        # Trusting "the sidecar files exist" here would silently hand back a
        # stale DB built from the FASTA's old content.
        if makeblastdb_path is None:
            makeblastdb_path = discover_blast()["makeblastdb"]
        return build_db(makeblastdb_path, path, db_prefix, title=os.path.basename(db_prefix))
    raise BlastError(
        f"{path!r} is neither a BLAST DB prefix (no .nhr/.nin/.nsq) nor a "
        f"FASTA file (no {'/'.join(_FASTA_EXTS)} extension)."
    )


def save_cached_hits(hits: List[Dict[str, Any]], path: str) -> None:
    """Save BLAST hits to a JSON file for caching."""
    with open(path, 'w') as f:
        json.dump(hits, f)


# ---------------------------------------------------------------------------
# Subject retrieval for the full-half mismatch recount
# ---------------------------------------------------------------------------

# How many distinct subjects to hold in memory, and the total bases across them.
# BOTH bounds matter: an entry cap alone is unbounded in bytes (one chromosome
# blows the budget), a byte cap alone is unbounded in dict overhead. Both are
# fixed regardless of DB size.
# Prefetch must not evict its own results before scoring reads them. A BLAST
# chunk is 5,000 probes = 10,000 halves, and on a transcriptome self-screen those
# touch tens of thousands of distinct subjects. A cache too small to hold them
# lets the batch fetch overflow it, and every evicted subject is then re-fetched
# one subprocess at a time during scoring: the same fork-per-subject cost the
# batching exists to remove, arriving by a different door, with blastdbcmd
# running continuously at near-zero CPU.
#
# Sized so a whole chunk's subjects fit. 512 MB is above the entire human cDNA
# set (390 Mbases), so in practice the cap is never the binding constraint for a
# single species; it remains a cap, and only what is actually touched is held.
SUBJECT_CACHE_MAX_ENTRIES = 500000
SUBJECT_CACHE_MAX_BASES = 512 * 1024 * 1024

# A subject longer than this is never cached whole; its windows are fetched
# per-range instead (and those windows are cached, keyed by range). Transcriptome
# subjects are a few kb, so this only trips on genome-scale records, where
# holding the sequence would be the memory bug rather than the optimisation.
SUBJECT_CACHE_MAX_SEQ_BASES = 4 * 1024 * 1024

_COMPLEMENT = str.maketrans("ACGTUNacgtunRYKMBVDHrykmbvdh",
                            "TGCAANtgcaanYRMKVBHDyrmkvbhd")


def reverse_complement(seq: str) -> str:
    """Reverse complement, preserving IUPAC ambiguity codes.

    Local to blast.py rather than imported from chemistry: this one is used to
    put a *subject* window into the probe's reading frame for a minus-strand HSP,
    which is a BLAST-coordinate concern, not an oligo-chemistry one.
    """
    return seq.translate(_COMPLEMENT)[::-1]


class SubjectWindows:
    """Retrieve subject bases opposite the UNALIGNED part of a probe half.

    ``blastn`` reports only the aligned segment. The 10x rule is stated over the
    whole 25 nt half, so the positions BLAST never compared have to be compared
    here — and that needs the subject's actual bases, which only the database
    has. This class wraps ``blastdbcmd`` with the caching that makes that
    affordable: a whole-transcriptome screen hits the same abundant subject
    thousands of times, and one ``blastdbcmd`` process per hit would cost more
    than the BLAST search itself.

    Every failure mode returns ``None`` — no blastdbcmd on the machine, the DB
    built without ``-parse_seqids``, a subject absent from the DB, a range that
    runs off the end of the record, a non-zero exit, a timeout. The caller must
    treat ``None`` as "count the unaligned positions as MATCHES" (see
    ``specificity.full_half_mismatches``), which under-counts mismatches and so
    over-rejects. Never fall back to inferring the count from HSP geometry: that
    is the error this class exists to prevent.

    Negative results are cached too. Without that, a DB with no accession index
    would spawn one doomed subprocess per hit — the exact cost the cache exists
    to avoid, incurred precisely in the case where it buys nothing.
    """

    def __init__(
        self,
        db_prefixes: Sequence[str],
        blastdbcmd_path: Optional[str] = None,
        *,
        max_entries: int = SUBJECT_CACHE_MAX_ENTRIES,
        max_bases: int = SUBJECT_CACHE_MAX_BASES,
        max_seq_bases: int = SUBJECT_CACHE_MAX_SEQ_BASES,
    ) -> None:
        self.db_prefixes = [p for p in db_prefixes if p]
        self.blastdbcmd_path = (
            blastdbcmd_path if blastdbcmd_path is not None else find_blastdbcmd()
        )
        self.max_entries = max_entries
        self.max_bases = max_bases
        self.max_seq_bases = max_seq_bases
        # (db, sseqid) -> sequence, or None for "unavailable", or "" for
        # "present but too long to cache; fetch ranges instead".
        self._seqs: "OrderedDict[Tuple[str, str], Optional[str]]" = OrderedDict()
        self._bases = 0
        # (db, sseqid, start, end) -> window, for the oversized-subject path.
        self._ranges: "OrderedDict[Tuple[str, str, int, int], Optional[str]]" = \
            OrderedDict()
        # Aggregate counters so the caller can emit ONE warning per run rather
        # than one log line per hit (a whole-transcriptome run has tens of
        # millions of hits; per-hit logging is itself an outage).
        self.n_exact = 0
        self.n_fallback = 0
        self.n_calls = 0
        # DBs that can NEVER serve retrieval, learned from the first failure.
        # A database built without -parse_seqids has no accession index at all
        # (its deflines are BL_ORD_ID:N), so every -entry lookup fails the same
        # way. Without this the screen spawns one doomed blastdbcmd per distinct
        # subject to re-learn it, which on a transcriptome can take many hours
        # while doing no useful work. One failure is enough; after it, this DB
        # is skipped and counts go conservative.
        self._dead_dbs: set = set()

    @property
    def available(self) -> bool:
        """True if an exact recount is even possible right now.

        False once every configured DB has proved unable to serve retrieval —
        which is a property of how the DB was BUILT, not of this run, so it can
        only be discovered by asking. Callers use this to stop asking.
        """
        if not self.blastdbcmd_path:
            return False
        return any(p not in self._dead_dbs for p in self.db_prefixes)

    @property
    def live_dbs(self) -> List[str]:
        """Configured DBs not yet proved incapable of retrieval."""
        return [p for p in self.db_prefixes if p not in self._dead_dbs]

    # An accession that cannot exist, used to ask a DB whether it is capable of
    # -entry lookups at all without needing to know a real accession first.
    _CAPABILITY_PROBE = "__flexiprobes_capability_probe__"

    def verify(self) -> List[str]:
        """Probe each DB once for -entry capability; return the usable ones.

        A DB built without -parse_seqids answers EVERY lookup with "DB contains
        no accession info", while a capable DB answers this deliberately absent
        accession with "Entry not found". The two stderrs are what distinguish
        "cannot answer" from "does not contain", and one call per DB settles it
        before any hit is scored.

        This exists because the failure it detects is silent and total. At the
        exact pushdown floor BLAST reports every chance 9 nt core; if the subject
        cannot be retrieved, the conservative fallback counts the other 16
        positions as MATCHES, so a meaningless core reads as a 0-mismatch
        off-target and EVERY probe can be rejected. Discovering that after a
        full run, or worse after a cluster run, is not acceptable when one call
        up front settles it.
        """
        for db in list(self.db_prefixes):
            if db in self._dead_dbs:
                continue
            self._run(["-db", db, "-entry", self._CAPABILITY_PROBE,
                       "-outfmt", "%a"], db=db)
        return self.live_dbs

    def prefetch(self, sseqids: Iterable[str]) -> int:
        """Warm the cache for many subjects in ONE blastdbcmd call per DB.

        Returns the number of subjects newly cached.

        This is the difference between the exact recount being usable and not.
        ``_sequence`` spawns one blastdbcmd process per cache miss, and each
        invocation re-opens the database index; a human-vs-human screen touches
        tens of thousands of distinct subjects, so the per-hit path is dominated
        entirely by process startup. Even a small panel against GRCh38 cdna.all
        can spend hours forking and almost no time computing.

        ``-entry_batch`` takes a file of accessions and streams them all back
        from one process, which collapses that to a single call per chunk.

        Failures are silent by design: a subject that does not come back simply
        stays uncached and falls through to the per-entry path, which already
        degrades to a conservative count. Prefetching is an optimisation, and an
        optimisation must never be able to change a verdict.
        """
        if not self.available:
            return 0
        wanted = [s for s in dict.fromkeys(sseqids) if s]
        if not wanted:
            return 0
        cached = 0
        for db in self.live_dbs:
            missing = [s for s in wanted if (db, s) not in self._seqs]
            if not missing:
                continue
            batch = self._run_batch(db, missing)
            if not batch:
                continue
            for sid, seq in batch.items():
                key = (db, sid)
                if key in self._seqs:
                    continue
                # Same oversize policy as _sequence: too long to hold whole is
                # recorded as "" so window() knows to fetch ranges instead.
                self._store_seq(key, "" if len(seq) > self.max_seq_bases else seq)
                cached += 1
            # Only subjects this DB actually returned are resolved; the rest fall
            # through to the next DB exactly as the per-entry path would.
            wanted = [s for s in wanted if (db, s) in self._seqs
                      and self._seqs[(db, s)] is None]
            if not wanted:
                break
        return cached

    def _run_batch(self, db: str, sseqids: Sequence[str]) -> Dict[str, str]:
        """One ``-entry_batch`` call; ``{accession: sequence}``, empty on failure."""
        if not self.blastdbcmd_path:
            return {}
        import tempfile
        self.n_calls += 1
        out: Dict[str, str] = {}
        tmp_path = ""
        try:
            with tempfile.NamedTemporaryFile("w", suffix=".ids", delete=False) as fh:
                fh.write("\n".join(sseqids) + "\n")
                tmp_path = fh.name
            result = subprocess.run(
                [self.blastdbcmd_path, "-db", db, "-entry_batch", tmp_path,
                 "-outfmt", "%a\t%s"],
                capture_output=True, text=True, timeout=1800, env=_blast_env(),
            )
            if result.returncode != 0:
                self._note_failure(db, result.stderr)
                logger.debug("blastdbcmd -entry_batch exited %s: %s",
                             result.returncode, result.stderr.strip()[:200])
            for line in result.stdout.splitlines():
                acc, _, seq = line.partition("\t")
                if acc and seq:
                    out[acc.strip()] = "".join(seq.split()).upper()
        except (OSError, subprocess.SubprocessError):
            logger.debug("blastdbcmd -entry_batch failed", exc_info=True)
            return {}
        finally:
            if tmp_path:
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass
        return out

    def window(self, sseqid: str, start: int, end: int) -> Optional[str]:
        """Subject bases for 1-based inclusive ``start..end``, or None.

        ``start``/``end`` are PLUS-strand subject coordinates and must satisfy
        ``1 <= start <= end``. Minus-strand handling is the caller's job (see
        ``specificity``): the caller knows which orientation the probe reads in,
        this method only knows the database.
        """
        if start < 1 or end < start or not sseqid:
            return None
        seq = self._sequence(sseqid)
        if seq is not None:
            # Whole record cached: an out-of-range request means our coordinate
            # arithmetic was wrong, so return None rather than a silent short
            # slice that would read as extra matches.
            if end > len(seq):
                return None
            return seq[start - 1:end]
        return self._range(sseqid, start, end)

    def window_clamped(
        self, sseqid: str, start: int, end: int,
    ) -> Optional[Tuple[str, int]]:
        """``(bases available, positions that lie off the record)``, or None.

        None means the SUBJECT ITSELF could not be retrieved — genuinely unknown.
        A request running past either end of a *retrievable* record is neither an
        error nor unknown: there is no sequence there, so the probe has nothing
        to base-pair with at those positions. Reporting them separately lets the
        caller count them as mismatches, which is what they physically are.

        Conflating the two is what :func:`window` does (it returns None for both),
        and the cost is large: every hit near a transcript end would degrade to
        the conservative count, and conservative counts over-reject, so a human
        screen can lose almost every probe.
        """
        if not sseqid or end < start:
            return None
        seq = self._sequence(sseqid)
        if seq is None:
            lo = max(1, start)
            win = self._range(sseqid, lo, end) if end >= lo else ""
            if win is None:
                return None      # the record is the thing we could not get
            return win, (end - start + 1) - len(win)
        lo, hi = max(1, start), min(len(seq), end)
        win = seq[lo - 1:hi] if hi >= lo else ""
        return win, (end - start + 1) - len(win)

    # -- internals ---------------------------------------------------------

    def _sequence(self, sseqid: str) -> Optional[str]:
        """The whole subject, or None if unavailable / too long to cache."""
        for db in self.live_dbs:
            key = (db, sseqid)
            if key in self._seqs:
                self._seqs.move_to_end(key)
                cached = self._seqs[key]
                if cached:
                    return cached
                continue            # None (absent here) or "" (oversized)
            seq = self._run(["-db", db, "-entry", sseqid, "-outfmt", "%s"], db=db)
            if seq is None:
                self._store_seq(key, None)
                continue
            if len(seq) > self.max_seq_bases:
                self._store_seq(key, "")     # present, but fetch by range
                return None
            self._store_seq(key, seq)
            return seq
        return None

    def _range(self, sseqid: str, start: int, end: int) -> Optional[str]:
        """Ranged fetch, for subjects too long to hold whole."""
        for db in self.live_dbs:
            if self._seqs.get((db, sseqid), None) != "":
                continue           # not the oversized-in-this-DB case
            key = (db, sseqid, start, end)
            if key in self._ranges:
                self._ranges.move_to_end(key)
                return self._ranges[key]
            got = self._run([
                "-db", db, "-entry", sseqid,
                "-range", f"{start}-{end}", "-outfmt", "%s",
            ], db=db)
            self._ranges[key] = got
            while len(self._ranges) > self.max_entries:
                self._ranges.popitem(last=False)
            if got is not None:
                return got
        return None

    # stderr fragments meaning "this database can never answer -entry", as
    # opposed to "this particular accession is not in it". The first is a build
    # property and disables the DB; the second is normal and must not.
    _NO_ACCESSION_MARKERS = ("no accession info", "not a valid database")

    def _note_failure(self, db: str, stderr: str) -> None:
        low = (stderr or "").lower()
        if any(m in low for m in self._NO_ACCESSION_MARKERS):
            if db not in self._dead_dbs:
                self._dead_dbs.add(db)
                logger.warning(
                    "BLAST DB %s has no accession index (built without "
                    "-parse_seqids), so off-target mismatch counts cannot be "
                    "recounted against it. Every count from this DB will be "
                    "CONSERVATIVE — it over-rejects, never over-passes. Rebuild "
                    "with `makeblastdb -parse_seqids` to enable the exact "
                    "recount. Not retrying this DB.",
                    os.path.basename(db),
                )

    def _store_seq(self, key: Tuple[str, str], value: Optional[str]) -> None:
        self._seqs[key] = value
        self._bases += len(value or "")
        while self._seqs and (
            len(self._seqs) > self.max_entries or self._bases > self.max_bases
        ):
            _, evicted = self._seqs.popitem(last=False)
            self._bases -= len(evicted or "")

    def _run(self, args: List[str], db: str = "") -> Optional[str]:
        """One blastdbcmd invocation; None on ANY failure."""
        if not self.blastdbcmd_path:
            return None
        self.n_calls += 1
        try:
            result = subprocess.run(
                [self.blastdbcmd_path] + args,
                capture_output=True, text=True, timeout=120, env=_blast_env(),
            )
        except (OSError, subprocess.SubprocessError):
            logger.debug("blastdbcmd invocation failed", exc_info=True)
            return None
        if result.returncode != 0:
            # Expected and common: "DB contains no accession info" (built
            # without -parse_seqids), "Entry not found". debug, not warning —
            # the caller aggregates and reports once.
            self._note_failure(db, result.stderr)
            logger.debug("blastdbcmd exited %s: %s",
                         result.returncode, result.stderr.strip()[:200])
            return None
        seq = "".join(result.stdout.split()).upper()
        return seq or None
