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

"""Multi-strain consensus + per-strain probe-binding quality scoring.

This module supports two run modes that go beyond the single-sequence
pipeline:

  • **Conserved consensus** — every input FASTA record is treated as a
    different strain of the same target fragment (e.g. OROV L-segment
    strain A + strain B). The module aligns them, builds a per-position
    consensus (majority + IUPAC ambiguity for variant positions), and
    exposes the windows that are identical across all strains. Probes
    designed from those windows bind every strain perfectly.

  • **Strain-unique** — the user picks one target strain; the module
    returns windows where the target differs from EVERY other strain by
    a configurable mismatch count (default 5, matching the 10X
    off-target ≥ 5 mismatches rule). Useful for diagnostic panels.

For every candidate probe in either mode, ``score_probe_against_all_strains()``
returns a per-strain breakdown of mismatch counts + an estimated ΔTm so
the review CSV can show how well a given consensus probe still binds
each individual strain.

Alignment policy:
  • Same length within a small tolerance (≤ 5 bp): the strains are
    trimmed to the shortest and matched position-by-position.
  • Larger length differences: a Biopython pairwise aligner is used to
    handle insertions/deletions, then aligned positions of the reference
    strain (the first record) are used as the consensus coordinate
    system. If Biopython is unavailable, raises ``ValueError`` with a
    pointer to install or to pre-align with MAFFT/MUSCLE.
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple


logger = logging.getLogger("flexprobe.consensus")


# IUPAC ambiguity codes (DNA) — used to summarize positions where strains
# disagree but the bases are still in {A,C,G,T}.
IUPAC = {
    frozenset("A"): "A",
    frozenset("C"): "C",
    frozenset("G"): "G",
    frozenset("T"): "T",
    frozenset("AC"): "M",
    frozenset("AG"): "R",
    frozenset("AT"): "W",
    frozenset("CG"): "S",
    frozenset("CT"): "Y",
    frozenset("GT"): "K",
    frozenset("ACG"): "V",
    frozenset("ACT"): "H",
    frozenset("AGT"): "D",
    frozenset("CGT"): "B",
    frozenset("ACGT"): "N",
}

# Reverse map: IUPAC code -> set of constituent bases. Used by the binding-
# quality scorer to treat an IUPAC code in a consensus probe as "any of
# these bases" rather than as a hard mismatch.
IUPAC_BASES = {
    "A": frozenset("A"), "C": frozenset("C"),
    "G": frozenset("G"), "T": frozenset("T"),
    "M": frozenset("AC"), "R": frozenset("AG"), "W": frozenset("AT"),
    "S": frozenset("CG"), "Y": frozenset("CT"), "K": frozenset("GT"),
    "V": frozenset("ACG"), "H": frozenset("ACT"),
    "D": frozenset("AGT"), "B": frozenset("CGT"),
    "N": frozenset("ACGT"),
}


# Mismatch / ΔTm heuristics. These are intentionally conservative and
# back-of-the-envelope; the goal is comparative (probe X binds strain A
# slightly better than strain B), not absolute thermodynamic prediction.
#
# Sources for the rule of thumb:
#   • CG000621 Rev D: "SNPs and mismatches should be at least four bp
#     away from the ligation junction" → we weight a mismatch by how close
#     it sits to the nick, on a linear ramp from 3x at the junction down to
#     1x at either probe end (constants.mismatch_position_weight). The ramp
#     gives the same cost as a binary "double-weight anything within
#     JUNCTION_WINDOW" step at the nick and at the probe ends, and grades
#     the positions in between.
#   • Owczarzy et al. nearest-neighbour estimates put a single internal
#     mismatch at roughly –4 to –6 °C of Tm depending on context; junction
#     mismatches in ligation chemistry are functionally worse (the ligase
#     simply won't seal).
# Canonical home is flexprobe.constants; re-exported here for back-compat.
from .constants import (  # noqa: E402
    DEGENERATE_JUNCTION_CLEARANCE, DTM_JUNCTION_BONUS, DTM_PER_MISMATCH,
    JUNCTION_WINDOW, OFF_TARGET_JUNCTION_DUPLEX_BP,
    STRAIN_UNIQUE_PER_HALF_DEFAULT,
    mismatch_position_weight,
)
# The strain-unique gate's within-half positional check reuses the OFF-target
# junction-duplex scan verbatim (see docstring on strain_unique_start_positions
# below for why importing it, rather than re-deriving the same arithmetic, is
# the point). specificity.py imports only from .blast / .models / .constants,
# none of which import this module, so this is not circular.
from .specificity import junction_clean_duplex  # noqa: E402
# ``DTM_JUNCTION_BONUS`` does not price the ON-target ΔTm: the graded ramp
# does, and gives the same endpoint (-6.0 C at the nick). It is imported
# because this module re-exports it for back-compat
# and specificity.py refers to it by name when explaining why the very same
# mismatch is PROTECTIVE off-target. The binding below is what marks that as a
# deliberate re-export rather than a leftover import.
_ = DTM_JUNCTION_BONUS


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class StrainAlignment:
    """Aligned strains + consensus for one fragment."""
    fragment_name: str
    strain_ids: List[str]
    aligned_seqs: List[str]
    consensus: str                       # may contain IUPAC ambiguity codes
    variant_positions: List[int] = field(default_factory=list)
    # Set to ``True`` if a pairwise aligner was used (i.e. there were
    # indels). When ``False`` we used the cheap trim-to-shortest path.
    used_pairwise: bool = False
    # Which aligner produced ``aligned_seqs``, e.g. "MAFFT v7.526 (bundled)",
    # "MUSCLE", "built-in progressive", or "trim-to-shortest" when the records
    # were close enough in length that no alignment was performed. Recorded so
    # the choice survives into a run's outputs and a finished result can say
    # which aligner produced it. This matters because the built-in fallback is
    # NOT a true multiple alignment.
    aligner: str = ""
    # Why an external aligner that WAS found could not be used, when the
    # built-in fallback ran in its place; "" otherwise, including when no
    # aligner was found at all. Recorded so a finished run can tell "no MAFFT
    # available" from "MAFFT shipped and broke"; a log warning alone does not
    # survive into the outputs.
    aligner_fallback: str = ""

    @property
    def n_strains(self) -> int:
        return len(self.strain_ids)

    @property
    def length(self) -> int:
        return len(self.consensus)


# ---------------------------------------------------------------------------
# Consensus building
# ---------------------------------------------------------------------------

def build_consensus(
    records: Sequence[Tuple[str, str]],
    *,
    fragment_name: str = "",
    length_tolerance: int = 5,
) -> StrainAlignment:
    """Build an alignment + consensus from ``(strain_id, sequence)`` records.

    Same-length strains (within ``length_tolerance``) are trimmed to the
    shortest length and matched position-by-position. For longer length
    differences, Biopython's pairwise aligner is invoked and the
    reference strain (the first record) is used as the coordinate system
    for the consensus.

    Raises ``ValueError`` if no records are given or if a non-ACGT base
    (other than N) is encountered after uppercasing.
    """
    if not records:
        raise ValueError("No strain records provided")
    cleaned = [(sid, (seq or "").upper()) for sid, seq in records]

    if len(cleaned) == 1:
        sid, seq = cleaned[0]
        return StrainAlignment(
            fragment_name=fragment_name or sid,
            strain_ids=[sid],
            aligned_seqs=[seq],
            consensus=seq,
            aligner="none (single record)",
        )

    lens = [len(seq) for _, seq in cleaned]
    if max(lens) - min(lens) <= length_tolerance:
        # Closely related, same length → cheap trim path.
        min_len = min(lens)
        aligned = [seq[:min_len] for _, seq in cleaned]
        used_pairwise = False
        aligner = "trim-to-shortest"
        fallback = ""
    else:
        # Prefer a real MSA tool (MAFFT/MUSCLE) when one is on PATH — the
        # built-in progressive-to-reference merge is only a fallback.
        msa, fallback = _external_msa_status(cleaned)
        used_pairwise = msa is None
        if msa is not None:
            aligned, aligner = msa
        else:
            aligned = _pairwise_align(cleaned)
            aligner = "built-in progressive"

    consensus, variants = _consensus_from_aligned(aligned)
    return StrainAlignment(
        fragment_name=fragment_name or cleaned[0][0],
        strain_ids=[sid for sid, _ in cleaned],
        aligned_seqs=aligned,
        consensus=consensus,
        variant_positions=variants,
        used_pairwise=used_pairwise,
        aligner=aligner,
        aligner_fallback=fallback,
    )


def _consensus_from_aligned(
    aligned: Sequence[str],
) -> Tuple[str, List[int]]:
    """Per-position majority-or-IUPAC consensus + list of variant positions."""
    out: List[str] = []
    variants: List[int] = []
    n = len(aligned[0])
    for i in range(n):
        bases = {seq[i] for seq in aligned}
        # Gap positions (from pairwise alignment) collapse to '-'; treat
        # them as variants and emit 'N' so downstream tiling skips them.
        if "-" in bases:
            out.append("N")
            variants.append(i)
            continue
        if len(bases) == 1:
            out.append(next(iter(bases)))
        else:
            valid = bases & set("ACGT")
            code = IUPAC.get(frozenset(valid), "N") if valid == bases else "N"
            out.append(code)
            variants.append(i)
    return "".join(out), variants


# Common install locations to search when MAFFT/MUSCLE is not on PATH. This is
# the same problem the BLAST resolver solves (see blast._COMMON_BLAST_DIRS): a
# macOS .app launched from Finder/Dock inherits only a minimal PATH
# (/usr/bin:/bin:/usr/sbin:/sbin), NOT the user's shell PATH, so a Homebrew
# MAFFT is invisible to a bare shutil.which() and the run drops silently to the
# built-in aligner.
_COMMON_MSA_DIRS = (
    "/opt/homebrew/bin",     # Homebrew, Apple Silicon
    "/usr/local/bin",        # Homebrew, Intel / manual installs
    "/opt/local/bin",        # MacPorts
    "/usr/bin", "/bin",
    "/opt/conda/bin",
)


def _msa_search_dirs() -> List[str]:
    """Directories to look in for an MSA executable, in priority order."""
    dirs: List[str] = []
    env = os.environ.get("FLEXPROBE_MAFFT_DIR")
    if env:
        dirs.append(env)
    # The active Python's / conda env's bin (covers conda-installed MAFFT).
    for base in (sys.prefix, os.environ.get("CONDA_PREFIX")):
        if base:
            dirs.append(os.path.join(base, "bin"))
    dirs.extend(_COMMON_MSA_DIRS)
    seen: set = set()
    out: List[str] = []
    for d in dirs:
        if d and d not in seen:
            seen.add(d)
            out.append(d)
    return out


def _find_msa_tool(name: str) -> Optional[str]:
    """Locate an MSA executable on PATH, else in the common install dirs."""
    path = shutil.which(name)
    if path:
        return path
    for d in _msa_search_dirs():
        cand = os.path.join(d, name)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None


def _is_windows() -> bool:
    """True on Windows. A function, not a constant, so tests can simulate it."""
    return os.name == "nt"


def _bundled_mafft() -> Optional[Tuple[str, str]]:
    """Return ``(driver, root)`` for the MAFFT vendored inside the app.

    The desktop build vendors MAFFT's official **all-in-one / portable**
    package, which upstream recommends for exactly this purpose: the Linux and
    macOS pages both describe it as intended for redistributing MAFFT inside
    another program package rather than for standalone use. Only meaningful in
    a frozen (PyInstaller) app; returns ``None`` when running from source.

    Every platform's package puts a driver named ``mafft.bat`` at the top of a
    relocatable folder, and upstream states no environment variable is needed
    because the driver locates its helpers relative to itself. What sits beside
    the driver differs, so ONLY ``mafft.bat`` is probed for. Layouts of the
    MAFFT 7.526 portable packages:

      * ``mafft-mac``      ``mafft.bat`` + ``mafftdir/{bin,libexec}``
      * ``mafft-linux64``  ``mafft.bat`` + ``mafftdir/{bin,libexec}``
      * ``mafft-win``      ``mafft.bat`` + ``mafft-signed.ps1`` + ``usr/`` + ``tmp/``

    Requiring ``mafftdir`` as well would silently fail on Windows, where the
    helpers live under ``usr/`` inside a partial MSYS tree instead.

    Note the layout deliberately does NOT match a distribution install such as
    Homebrew's (``bin/mafft`` plus a separate ``libexec/mafft``). A packaged
    build must ship the portable tree, not a distribution build, which is
    compiled against its build host and is not a redistribution artefact.

    ``mafft.bat`` is a Bourne shell script on macOS and Linux despite the
    extension, so it is invoked through an explicit interpreter there. On
    Windows it is a real batch file; see :func:`_bundled_mafft_invocation`.
    """
    if not getattr(sys, "frozen", False):
        return None
    roots: List[str] = []
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        roots.append(os.path.join(meipass, "mafft"))
    exe_dir = os.path.dirname(os.path.realpath(sys.executable))
    roots += [
        os.path.join(exe_dir, "mafft"),                       # onedir
        os.path.join(exe_dir, "..", "Frameworks", "mafft"),   # .app
        os.path.join(exe_dir, "..", "Resources", "mafft"),    # .app
    ]
    for r in roots:
        r = os.path.normpath(r)
        driver = os.path.join(r, "mafft.bat")
        if os.path.isfile(driver):
            return driver, r
    return None


def _no_console_window() -> Dict[str, Any]:
    """``subprocess`` keywords that stop a console window opening on Windows.

    ``bash.exe`` is a console program. Launched from the windowed desktop app
    without this, each call opens a console window of its own. Empty elsewhere,
    so macOS and Linux calls are unchanged.
    """
    return {"creationflags": 0x08000000} if _is_windows() else {}  # CREATE_NO_WINDOW


def _tool_version(cmd: List[str], env: Optional[dict] = None) -> str:
    """First line of a tool's version output, or "" if it cannot be read.

    Never raises: an unreadable version degrades the provenance label, and must
    not fail a run that is otherwise fine.
    """
    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=30, env=env, **_no_console_window())
    except (subprocess.SubprocessError, OSError):
        return ""
    text = ((p.stdout or "") + (p.stderr or "")).strip()
    return text.splitlines()[0].strip() if text else ""


def _bundled_mafft_invocation(
    driver: str, env: Dict[str, str],
) -> Tuple[List[str], Dict[str, str]]:
    """Command prefix and environment that run the vendored MAFFT ``driver``.

    macOS and Linux: ``mafft.bat`` is a Bourne shell script, and PyInstaller
    strips the executable bit off files shipped as data, so it is run through
    ``bash`` rather than exec'd. MAFFT_BINARIES tells a driver where its helper
    executables live; this package locates its own ``mafftdir`` relative to the
    driver and needs no such variable, but an inherited one from the user's
    shell would point at a DIFFERENT MAFFT install and silently mix the two, so
    it is removed rather than set.

    Windows: ``mafft.bat`` is a real batch file and there is no ``bash`` on
    PATH, so ``["bash", mafft.bat]`` would raise OSError and the run would fall
    back silently to the built-in aligner. The package's two launchers,
    ``mafft.bat`` and ``mafft-signed.ps1``, both set MAFFT_BINARIES, TMPDIR,
    MAFFT_TMPDIR and a ``/usr/bin``-first PATH, then run
    ``usr\\bin\\bash.exe /usr/bin/mafft``. This does the same directly, with
    the batch file's own values, instead of going through ``cmd.exe``: a batch
    file re-parses its arguments, and with delayed expansion on (as this one
    sets) a ``!``, ``%`` or ``^`` in the temp path would be mangled.
    """
    import tempfile

    env = dict(env)
    if _is_windows():
        root = os.path.dirname(driver)
        tmp = env.get("TMP") or env.get("TEMP") or tempfile.gettempdir()
        env["PATH"] = "/usr/bin/:" + env.get("PATH", "")
        env["MAFFT_BINARIES"] = "/usr/lib/mafft"
        env["TMPDIR"] = tmp
        env["MAFFT_TMPDIR"] = tmp
        bash = os.path.join(root, "usr", "bin", "bash.exe")
        return [bash, "/usr/bin/mafft"], env
    env.pop("MAFFT_BINARIES", None)
    return ["bash", driver], env


def _external_msa(
    records: Sequence[Tuple[str, str]],
) -> Optional[Tuple[List[str], str]]:
    """Align with MAFFT or MUSCLE. Returns ``(aligned, tool_label)`` or None.

    MAFFT is resolved bundled-first, for the same reason ``blast._find_tool``
    prefers its bundled BLAST+: the packaged app should align with its own
    pinned, known-good MAFFT rather than whatever version happens to be on the
    machine, so the same input gives the same alignment on two machines. After
    the bundle it falls back to ``$FLEXPROBE_MAFFT_DIR``, PATH, and the common
    install dirs.

    MUSCLE is consulted only when no MAFFT is found anywhere. It is never
    bundled: MAFFT is already preferred here, so a bundled MAFFT would make a
    bundled MUSCLE unreachable, and shipping two aligners whose selection
    depends on the host machine is a reproducibility hazard, not a feature.

    Returns the aligned sequences in the SAME order as ``records`` (equal
    length, gapped), or None when no MSA tool is available or it fails - the
    caller then falls back to the built-in progressive aligner.
    :func:`_external_msa_status` also says why a found tool failed.
    """
    return _external_msa_status(records)[0]


def _external_msa_status(
    records: Sequence[Tuple[str, str]],
) -> Tuple[Optional[Tuple[List[str], str]], str]:
    """:func:`_external_msa`, plus why an aligner that WAS found went unused.

    The second element is "" when the alignment succeeded or when no aligner
    was found at all, and otherwise a one-line reason naming the tool. It is
    what reaches ``StrainAlignment.aligner_fallback`` and the manifest.
    """
    import tempfile

    bundled = _bundled_mafft()
    mafft = bundled[0] if bundled else _find_msa_tool("mafft")
    muscle = None if mafft is not None else _find_msa_tool("muscle")

    if mafft is None and muscle is None:
        logger.warning(
            "No MAFFT or MUSCLE found: falling back to the built-in "
            "progressive aligner, which aligns each record against the FIRST "
            "record rather than performing a true multiple alignment. Install "
            "MAFFT, or set FLEXPROBE_MAFFT_DIR to the directory containing it.")
        return None, ""

    env = dict(os.environ)
    run_kw: Dict[str, Any] = {}
    mafft_cmd: List[str] = []
    if bundled is not None:
        mafft_cmd, env = _bundled_mafft_invocation(bundled[0], env)
        run_kw = _no_console_window()
        label = f"bundled MAFFT ({bundled[0]})"
    elif mafft is not None:
        mafft_cmd = [mafft]
        label = f"MAFFT ({mafft})"
    else:
        label = f"MUSCLE ({muscle})"
    instead = "; the built-in progressive aligner was used instead"

    with tempfile.TemporaryDirectory() as td:
        inp = os.path.join(td, "in.fa")
        outp = os.path.join(td, "out.fa")
        # Use index-based ids so we can restore input order regardless of how
        # the tool reorders/renames records.
        with open(inp, "w") as fh:
            for i, (_sid, seq) in enumerate(records):
                fh.write(f">{i}\n{seq}\n")
        try:
            if mafft is not None:
                with open(outp, "w") as out:
                    subprocess.run(mafft_cmd + ["--auto", "--quiet", inp],
                                   stdout=out, check=True, timeout=600,
                                   env=env, **run_kw)
                ver = _tool_version(mafft_cmd + ["--version"], env=env)
                tool = f"MAFFT {ver}".strip()
                if bundled is not None:
                    tool += " (bundled)"
            else:
                assert muscle is not None       # guaranteed by the early return
                try:                                    # MUSCLE v5
                    subprocess.run([muscle, "-align", inp, "-output", outp],
                                   check=True, timeout=600,
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)
                except subprocess.CalledProcessError:   # MUSCLE v3 syntax
                    subprocess.run([muscle, "-in", inp, "-out", outp],
                                   check=True, timeout=600,
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)
                tool = f"MUSCLE {_tool_version([muscle, '-version'])}".strip()
            by_idx = {int(sid): seq for sid, seq in _read_fasta_pairs(outp)}
        except (subprocess.SubprocessError, OSError, ValueError) as exc:
            logger.warning("External MSA failed (%s); using built-in aligner", exc)
            return None, (f"{label} failed ({type(exc).__name__}: {exc})"
                          + instead)

    opt = [by_idx.get(i) for i in range(len(records))]
    if any(a is None for a in opt):
        logger.warning("External MSA output malformed; using built-in aligner")
        return None, f"{label} output was missing records" + instead
    aligned = [a for a in opt if a is not None]    # all non-None past the guard
    if len({len(a) for a in aligned}) != 1:
        logger.warning("External MSA output malformed; using built-in aligner")
        return None, f"{label} output rows differed in length" + instead
    logger.info("Aligned %d strains with %s", len(records), tool)
    return ([a.upper() for a in aligned], tool), ""


def _read_fasta_pairs(path: str) -> List[Tuple[str, str]]:
    """Minimal (id, seq) FASTA reader for MSA tool output."""
    pairs: List[Tuple[str, str]] = []
    sid: Optional[str] = None
    chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if sid is not None:
                    pairs.append((sid, "".join(chunks)))
                sid, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if sid is not None:
        pairs.append((sid, "".join(chunks)))
    return pairs


def _pairwise_align(records: Sequence[Tuple[str, str]]) -> List[str]:
    """Align multiple strains to the first record using Biopython pairwise.

    For 2 strains: one alignment. For N strains: each non-reference strain
    is aligned against the reference and the reference's gapped form is
    consolidated so all aligned sequences share the same coordinate
    system. This is not true multiple-sequence alignment, but for ≤ 4
    closely related strains it produces sensible per-position consensus.
    """
    try:
        from Bio import Align
    except ImportError as exc:
        raise ValueError(
            "Strain sequences differ in length by more than the tolerance "
            "and Biopython is not installed. Either pre-align with MAFFT/"
            "MUSCLE and re-export equal-length FASTAs, or `pip install "
            "biopython` to enable the built-in pairwise aligner."
        ) from exc

    # The progressive-to-reference merge is only an approximation (each strain is
    # aligned to the reference, not to the growing MSA); beyond a handful of
    # closely related strains a true MSA is preferable. Warn and point at one.
    if len(records) > 4:
        logger.warning(
            "Building a consensus from %d strains with the built-in "
            "progressive pairwise aligner (not a true MSA). For this many "
            "strains, or divergent inputs, install MAFFT/MUSCLE (auto-used when "
            "on PATH) or pre-align and re-export equal-length FASTAs.",
            len(records),
        )

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    # Score scheme tuned for closely related strains (high match reward).
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    ref_id, ref_seq = records[0]
    aligned_strings: List[str] = []
    ref_gapped: Optional[str] = None
    for sid, seq in records:
        if sid == ref_id and ref_gapped is None:
            # Will be filled in by the first alignment.
            aligned_strings.append("")
            continue
        alns = aligner.align(ref_seq, seq)
        best = alns[0]
        r_gapped, s_gapped = _aligned_to_strings(best)
        if ref_gapped is None:
            ref_gapped = r_gapped
            aligned_strings[0] = ref_gapped
            aligned_strings.append(s_gapped)
        else:
            # Re-align s_gapped to the existing ref_gapped coordinate system
            # by inserting gaps at the union positions. new_col marks the
            # freshly-inserted columns so prior strains re-lift correctly.
            ref_gapped, s_gapped, new_col = _merge_alignment(
                ref_gapped, r_gapped, s_gapped)
            # Update existing aligned strings to the new ref coordinates.
            new = []
            for prev in aligned_strings:
                new.append(_relift_with_gaps(prev, new_col))
            aligned_strings = new
            aligned_strings.append(s_gapped)
    return aligned_strings


def _aligned_to_strings(alignment) -> Tuple[str, str]:
    """Extract the two gapped strings from a Biopython PairwiseAligner hit."""
    # Biopython 1.81+: alignment[0] is target gapped, alignment[1] is query
    # gapped. Convert to plain Python strings.
    target = str(alignment[0])
    query = str(alignment[1])
    return target, query


def _merge_alignment(
    ref_existing: str,
    ref_new: str,
    strain_new: str,
) -> Tuple[str, str, List[bool]]:
    """Merge a new (ref, strain) pairwise alignment into the existing ref
    coordinate system so all strains can share columns.

    Walks both gap-strings in lockstep, inserting matching gaps as needed.
    Also returns a per-output-column ``new_col`` mask: True where a brand-new
    column was inserted (the existing reference — and therefore every already-
    aligned strain — had no character there). That mask is what re-lifts the
    earlier strains without assuming they were gap-free (audit M1).
    """
    i = j = 0
    merged_ref = []
    merged_strain = []
    new_col: List[bool] = []
    while i < len(ref_existing) or j < len(ref_new):
        a = ref_existing[i] if i < len(ref_existing) else "-"
        b = ref_new[j] if j < len(ref_new) else "-"
        if a == "-" and b == "-":
            merged_ref.append("-")
            merged_strain.append("-")
            new_col.append(False)   # existing column consumed
            i += 1
            j += 1
        elif a == "-":
            # Existing has a gap the new alignment doesn't — keep the gap,
            # insert a gap in the new strain to stay aligned.
            merged_ref.append("-")
            merged_strain.append("-")
            new_col.append(False)   # existing column consumed
            i += 1
        elif b == "-":
            # New alignment has a gap; insert it everywhere — a brand-new column
            # the existing ref (and every prior strain) doesn't have.
            merged_ref.append("-")
            merged_strain.append(strain_new[j])
            new_col.append(True)
            j += 1
        else:
            merged_ref.append(a)
            merged_strain.append(strain_new[j])
            new_col.append(False)   # existing column consumed
            i += 1
            j += 1
    return "".join(merged_ref), "".join(merged_strain), new_col


def _relift_with_gaps(prev: str, new_col: List[bool]) -> str:
    """Re-lift an already-aligned strain into the widened coordinate system.

    ``new_col[k]`` is True where a brand-new column was inserted by the latest
    merge; there ``prev`` gets a fresh gap. Everywhere else ``prev`` supplies its
    next character *unchanged* — including its own gaps — so a strain that was
    already gapped is re-lifted correctly (the old approach re-derived positions
    from the ref and corrupted gapped strains; audit M1).
    """
    out = []
    pi = 0
    for is_new in new_col:
        if is_new:
            out.append("-")
        else:
            out.append(prev[pi] if pi < len(prev) else "-")
            pi += 1
    return "".join(out)


def resolve_to_reference(window: str, reference: str) -> str:
    """Return ``window`` with every ambiguous base replaced by ``reference``.

    A consensus built from several input sequences carries an ``N`` (or an
    IUPAC ambiguity code) wherever the sequences disagree. Those codes can't
    be ordered as a concrete probe, so at output time each ambiguous position
    defaults to the **first input sequence** (the reference / default strain):
    any base in ``window`` that is not a plain ``A/C/G/T`` is taken from the
    same position of ``reference`` instead. The result is a clean, orderable
    ACGT sequence with no N's.

    If ``reference`` is itself ambiguous or too short at a position there is
    nothing better to default to, so the original base is kept.
    """
    wu = window.upper()
    if all(b in "ACGT" for b in wu):
        return wu
    ru = reference.upper()
    out: List[str] = []
    for i, b in enumerate(wu):
        if b in "ACGT":
            out.append(b)
            continue
        r = ru[i] if i < len(ru) else b
        out.append(r if r in "ACGT" else b)
    return "".join(out)


# ---------------------------------------------------------------------------
# Region queries
# ---------------------------------------------------------------------------

def is_window_conserved(
    alignment: StrainAlignment, start: int, end: int,
) -> bool:
    """True iff every strain has the same bases in ``[start, end)``."""
    if not alignment.aligned_seqs:
        return False
    if end > len(alignment.aligned_seqs[0]):
        return False
    ref = alignment.aligned_seqs[0][start:end]
    if "-" in ref or "N" in ref:
        return False
    return all(s[start:end] == ref for s in alignment.aligned_seqs[1:])


def conserved_start_positions(
    alignment: StrainAlignment,
    *,
    tile_length: int = 50,
    tile_step: int = 1,
) -> List[int]:
    """Start positions (0-indexed) of windows that are conserved across
    all strains."""
    out: List[int] = []
    L = alignment.length
    for start in range(0, L - tile_length + 1, tile_step):
        if is_window_conserved(alignment, start, start + tile_length):
            out.append(start)
    return out


_ACGT = frozenset("ACGT")


@dataclass
class DegenerateWindow:
    """One window that qualifies under the degenerate tolerance rule.

    ``sequence`` is the IUPAC consensus over the window — the *parent* probe
    string. ``degenerate_positions`` are 0-based offsets INTO THE WINDOW (not
    into the alignment) where the records disagree, and ``iupac_codes`` is the
    consensus code at each, in the same order. Both are empty for a window that
    is identical across every record, which is what makes the exact-identity
    path (``max_degenerate_per_half=0``) fall out of the same call.
    """
    start: int
    sequence: str
    degenerate_positions: List[int] = field(default_factory=list)
    iupac_codes: List[str] = field(default_factory=list)

    @property
    def n_degenerate(self) -> int:
        return len(self.degenerate_positions)


def degenerate_conserved_windows(
    alignment: StrainAlignment,
    *,
    tile_length: int = 50,
    tile_step: int = 1,
    max_degenerate_per_half: int = 0,
    junction_clearance: int = DEGENERATE_JUNCTION_CLEARANCE,
) -> List[DegenerateWindow]:
    """Pan-record windows admitting at most N variant positions per probe half.

    This is :func:`conserved_start_positions` with a tolerance dial.
    ``max_degenerate_per_half=0`` (the default) is exact identity across every
    record and reproduces that function's result; higher values admit windows
    that carry IUPAC ambiguity where the records disagree. A window qualifies
    when all three hold:

    1. **No gap and no ambiguous input.** Every record must supply a plain
       A/C/G/T at every column of the window. This is the important one, and it
       is not the same as "the consensus has few non-ACGT bases":
       :func:`_consensus_from_aligned` emits ``'N'`` both for a *gap* column and
       for a four-way disagreement, and ``'N'`` is in
       ``constants.DEGENERATE_BASES``. A per-half budget that simply counted
       non-ACGT bases in the consensus would therefore silently admit windows
       spanning an indel — where one record genuinely has a base the other does
       not. A 50-mer designed there does not mispair against the other record,
       it *bulges*, and everything downstream of the indel is frame-shifted; it
       is not designable as a pan-record probe at all. On a real alignment with
       indels this can admit many extra, undesignable windows.
       ``is_window_conserved`` and ``constitutive_start_positions``
       already take this position (both reject any window containing ``'-'`` or
       ``'N'``), and nothing downstream will catch it either — ``filters``'
       gap guard tests for ``'-'``/``'.'``, so a consensus ``'N'`` arriving from
       a gap column sails straight past it.
    2. **Per-half budget.** Each ``tile_length // 2`` half carries at most
       ``max_degenerate_per_half`` variant positions. Per HALF, not per window:
       the two halves are separate oligos that hybridise independently, so two
       mismatches in one half is a different molecule from one in each.
    3. **Junction clearance.** No variant position falls within
       ``junction_clearance`` bp either side of the ligation junction — for a
       50-mer at the default clearance of 4 that is offsets 21..28 inclusive.
       CG000621 Rev D asks for SNPs and mismatches at least four bp away from
       the junction, and the ligase will not seal across a junction-proximal
       mispair. The band is derived from ``tile_length``, never from a hard-
       coded 25.

    Because rule 1 excludes gap columns and ambiguous input alike, the only
    non-ACGT bases that can reach a consensus parent built from these windows
    are TRUE IUPAC ambiguity codes — never a gap-derived ``'N'``. That is what
    makes the per-record fan-out well defined: every degenerate position has a
    real, concrete base in every record, so resolving the parent against a
    record's own aligned sequence yields the exact oligo that binds it.

    Returns the qualifying windows in ascending ``start`` order.
    """
    seqs = alignment.aligned_seqs
    if not seqs:
        return []
    length = min(len(s) for s in seqs)
    if tile_length <= 0 or length < tile_length:
        return []
    step = max(1, int(tile_step))
    budget = max(0, int(max_degenerate_per_half))
    clearance = max(0, int(junction_clearance))
    half = tile_length // 2

    # Per-column pass, once, so the window scan is prefix-sum arithmetic rather
    # than rebuilding a set of bases for all 50 columns of every window.
    #   bad[i]  — some record has a gap or a non-ACGT base here (rule 1)
    #   var[i]  — the records disagree here (all in ACGT)
    #   cons[i] — the consensus character: the shared base, or the IUPAC code
    upper = [s.upper() for s in seqs]
    bad = [False] * length
    var = [False] * length
    cons: List[str] = []
    for i in range(length):
        bases = {s[i] for s in upper}
        if not bases <= _ACGT:
            bad[i] = True
            cons.append("N")
            continue
        if len(bases) == 1:
            cons.append(next(iter(bases)))
        else:
            var[i] = True
            cons.append(IUPAC.get(frozenset(bases), "N"))

    # Prefix sums: n_bad[k] / n_var[k] = count over columns [0, k).
    n_bad = [0] * (length + 1)
    n_var = [0] * (length + 1)
    for i in range(length):
        n_bad[i + 1] = n_bad[i] + (1 if bad[i] else 0)
        n_var[i + 1] = n_var[i] + (1 if var[i] else 0)

    # The junction band as window-relative offsets: [half - clearance,
    # half + clearance). Clearance 0 leaves it empty (no band check); use
    # FilterThresholds.reject_junction_base_n to turn the gate off explicitly.
    j_lo = max(0, half - clearance)
    j_hi = min(tile_length, half + clearance)

    out: List[DegenerateWindow] = []
    for start in range(0, length - tile_length + 1, step):
        end = start + tile_length
        if n_bad[end] - n_bad[start] > 0:
            continue                                    # rule 1
        n_left = n_var[start + half] - n_var[start]
        n_right = n_var[end] - n_var[start + half]
        if n_left > budget or n_right > budget:
            continue                                    # rule 2
        if j_hi > j_lo and (n_var[start + j_hi] - n_var[start + j_lo]) > 0:
            continue                                    # rule 3
        positions = [i for i in range(tile_length) if var[start + i]]
        out.append(DegenerateWindow(
            start=start,
            sequence="".join(cons[start:end]),
            degenerate_positions=positions,
            iupac_codes=[cons[start + i] for i in positions],
        ))
    return out


def resolve_window_per_record(
    window: str,
    aligned_seqs: Sequence[str],
    strain_ids: Sequence[str],
    start: int,
) -> List[Tuple[str, List[str]]]:
    """Fan a degenerate ``window`` out into one concrete oligo per record.

    Returns ``[(resolved_sequence, [record_id, ...]), ...]`` in first-seen
    record order. Records whose sequences resolve identically COLLAPSE into a
    single entry listing all of them — two strains that happen to agree at every
    variant position of this window need one oligo between them, not two.

    Resolutions are taken per record, never as the cartesian product of the
    IUPAC codes. The product enumerates sequences no input organism carries: a
    window with an 'R' at offset 3 and a 'Y' at offset 40 has four product
    members but only ever two real strains, and the two invented ones would be
    synthesised, screened and ordered for nothing.
    """
    end = start + len(window)
    out: List[Tuple[str, List[str]]] = []
    index: dict = {}
    for sid, aseq in zip(strain_ids, aligned_seqs):
        resolved = resolve_to_reference(window, aseq[start:end])
        slot = index.get(resolved)
        if slot is None:
            index[resolved] = len(out)
            out.append((resolved, [sid]))
        else:
            out[slot][1].append(sid)
    return out


def strain_unique_start_positions(
    alignment: StrainAlignment,
    *,
    target_strain_idx: int,
    tile_length: int = 50,
    tile_step: int = 1,
    min_mismatches_vs_others: int = 5,
    per_half: bool = STRAIN_UNIQUE_PER_HALF_DEFAULT,
    junction_duplex_bp: int = OFF_TARGET_JUNCTION_DUPLEX_BP,
) -> List[Tuple[int, int]]:
    """Start positions where the target strain differs from every other
    strain by at least ``min_mismatches_vs_others`` bp.

    ``per_half`` (default True) requires those mismatches to fall within a
    SINGLE probe half — i.e. the target must differ from every other strain by
    ``min_mismatches_vs_others`` in its LHS *or* in its RHS. This is the same
    unit the 10x threshold is written in ("at least five mismatches in at least
    one of the LHS or RHS probes", CG000621 Rev D / CG000839 Rev A) and the same
    unit ``specificity.evaluate_specificity`` scores off-targets in.

    It matters because RTL only makes a library molecule when BOTH halves
    hybridise and the ligase seals the nick between them. Mismatches split
    across the junction (say 3 LHS + 2 RHS) total 5 but leave both halves
    perfectly capable of binding the non-target strain, so such a probe is not
    strain-discriminating however the whole-window arithmetic reads.

    Per-half mode also counts each gap column in the other strain as a single
    mismatch, matching :func:`score_probe_against_all_strains`. The legacy
    whole-window mode instead treats any gap as automatically discriminating,
    which turns a 2 nt indel into "25 mismatches" and lets through windows whose
    real discrimination is 3 or 4 bp.

    Set ``per_half=False`` for the legacy whole-window count, retained so panels
    designed with it can be reproduced exactly.

    WITHIN-HALF POSITION — ``junction_duplex_bp``. Per-half settles "which half"
    but not "where in the half": a raw >= 5 count says nothing about whether
    those mismatches sit next to the nick or clear out at the probe's terminal
    end, and ligation cares only about the nick. Two windows can carry the
    identical per-half count and behave oppositely:

      * 5 mismatches bunched at the OUTER end of the LHS leave 20+ contiguous
        matched bases spanning the nick — the exact condition under which the
        ligase seals and the pair still makes a library molecule on the
        "excluded" strain, satisfying the count rule while being wrong about
        the biology.
      * The same 5 mismatches sitting against the nick break that run
        immediately, so no clean duplex can form there regardless of how the
        rest of the half reads.

    This is the mirror image of a gate this codebase already applies on the
    off-target side (see ``constants.py`` around ``OFF_TARGET_JUNCTION_DUPLEX_BP``,
    and ``specificity.py``): :func:`specificity.junction_clean_duplex` measures
    the unbroken run of matched bases immediately either side of the nick in
    CONTIGUOUS reference sequence, and a run of >= ``junction_duplex_bp`` on
    BOTH sides means the duplex can seal there. Off-target uses that to reject
    a pair the count rule alone would pass; here it is reused, unmodified, to
    withhold "strain-unique" from a window the count rule alone would grant —
    a half only certifies discrimination against a given other strain when it
    BOTH reaches ``min_mismatches_vs_others`` AND that other strain's window
    does not present a sealing duplex around the nick. Because this can only
    remove windows the raw count already admitted, it never contradicts the
    10x count rule — same relationship the off-target gate has to it — and
    reusing the same function/threshold (rather than inventing a new one) means
    the two "does this half actually ligate here" questions are answered by one
    piece of arithmetic, not two that could silently drift apart. Without it,
    the per-half count alone can over-certify a distal-mismatch window as
    strain-unique when the duplex would still seal.

    Only engaged when ``per_half`` is True: the legacy whole-window mode stays
    bit-for-bit reproducible (see above), and a duplex check keyed to a
    two-half nick has no meaning without the per-half split. Set
    ``junction_duplex_bp=0`` to disable it and use the per-half count alone,
    the same "0 disables" convention as the off-target twin of this knob.

    Returns ``[(start, min_mm_observed)]`` so the caller can rank windows by how
    distinctive they are. ``min_mm_observed`` is reported in whichever unit the
    gate used: the whole-window mismatch count when ``per_half`` is False, and
    the best-half count (the larger of the two halves, minimised over the other
    strains) when it is True — so in both modes the value is the quantity that
    was actually compared against ``min_mismatches_vs_others``. The duplex
    check is a pass/fail veto layered on top and does not change what this
    number means.
    """
    out: List[Tuple[int, int]] = []
    target = alignment.aligned_seqs[target_strain_idx]
    others = [
        s for i, s in enumerate(alignment.aligned_seqs)
        if i != target_strain_idx
    ]
    # Follow the tile size rather than a hard-coded 25, matching
    # score_probe_against_all_strains.
    half_length = tile_length // 2
    L = len(target)
    for start in range(0, L - tile_length + 1, tile_step):
        t_win = target[start:start + tile_length]
        if "-" in t_win or "N" in t_win:
            continue
        min_mm = None
        # Flips True the moment ANY other strain seals a junction-spanning
        # duplex against this window: one strain doing so is enough to make
        # the window not-discriminating against it, so the whole position
        # fails regardless of how the mismatch counts read.
        duplex_seals_any = False
        for o in others:
            o_win = o[start:start + tile_length]
            if per_half:
                # Each gap column is one position the probe cannot base-pair,
                # so it counts as one mismatch — exactly what the '-' -> 'X'
                # substitution in score_probe_against_all_strains does, which is
                # what fills the lhs_mm__* / rhs_mm__* columns of results.csv.
                # The legacy branch below instead treats ANY gap as fully
                # discriminating; that overstates a short indel badly (a 2 nt
                # deletion is a small bulge, not 25 mismatches) but is kept
                # as-is so whole-window runs reproduce bit-for-bit.
                o_cmp = o_win.replace("-", "X")
                # The discriminating power of the better half: a probe only
                # fails to ligate if ONE half is too mismatched to bind.
                lhs = sum(1 for a, b in zip(t_win[:half_length],
                                            o_cmp[:half_length]) if a != b)
                rhs = sum(1 for a, b in zip(t_win[half_length:],
                                            o_cmp[half_length:]) if a != b)
                mm = max(lhs, rhs)
                if junction_duplex_bp > 0:
                    clean_left, clean_right = junction_clean_duplex(
                        t_win, o_cmp, nick=half_length)
                    if min(clean_left, clean_right) >= junction_duplex_bp:
                        duplex_seals_any = True
            elif "-" in o_win:
                # gap in other strain — by definition very different here;
                # we'll count as len(window) so it always passes
                mm = tile_length
            else:
                mm = sum(1 for a, b in zip(t_win, o_win) if a != b)
            if min_mm is None or mm < min_mm:
                min_mm = mm
        if (min_mm is not None and min_mm >= min_mismatches_vs_others
                and not duplex_seals_any):
            out.append((start, min_mm))
    return out


def constitutive_start_positions(
    alignment: StrainAlignment,
    *,
    tile_length: int = 50,
    tile_step: int = 1,
    min_isoforms: Optional[int] = None,
) -> List[Tuple[int, str, Tuple[int, ...]]]:
    """Isoform-collapse: windows shared by a QUORUM of isoforms.

    For each aligned window, group the isoforms by their (ungapped, N-free)
    window sequence and take the LARGEST group. If that group has at least
    ``min_isoforms`` members the window qualifies: one probe designed from the
    shared sequence binds every isoform in the group.

    ``min_isoforms`` defaults to all records (strict constitutive == the same as
    ``conserved`` mode). Lower it to cover, say, the dominant isoforms when a
    rare isoform has no shared region.

    Returns ``[(start, shared_sequence, covered_isoform_indices)]``, sorted by
    start, with the covered group preferring more isoforms then lower indices so
    the result is deterministic.
    """
    seqs = alignment.aligned_seqs
    n = len(seqs)
    if n == 0:
        return []
    need = n if min_isoforms is None else max(1, min(int(min_isoforms), n))
    out: List[Tuple[int, str, Tuple[int, ...]]] = []
    L = alignment.length
    for start in range(0, L - tile_length + 1, tile_step):
        groups: dict = {}
        for i, s in enumerate(seqs):
            w = s[start:start + tile_length]
            if "-" in w or "N" in w:
                continue
            groups.setdefault(w, []).append(i)
        if not groups:
            continue
        # Largest group wins; ties broken by the lowest member index for
        # determinism.
        best_w, best_idx = max(
            groups.items(), key=lambda kv: (len(kv[1]), -min(kv[1]))
        )
        if len(best_idx) >= need:
            out.append((start, best_w, tuple(sorted(best_idx))))
    return out


# ---------------------------------------------------------------------------
# Per-strain binding quality
# ---------------------------------------------------------------------------

@dataclass
class StrainBindingScore:
    strain_id: str
    total_mismatches: int
    lhs_mismatches: int          # bp 0..24 (LHS-bound)
    rhs_mismatches: int          # bp 25..49 (RHS-bound)
    junction_mismatches: int     # bp 21..28 (within JUNCTION_WINDOW of pos 25)
    mismatch_positions: List[int]
    # Rough ΔTm vs perfect-match, with every mismatch weighted by its distance
    # from the ligation nick (constants.mismatch_position_weight).
    estimated_dTm: float
    likely_binds: bool           # True iff total <= 2 AND junction == 0


def _binding_score_for_window(
    consensus_window: str,
    strain_window: str,
    strain_id: str,
    half_length: int = 25,
    junction_window: int = JUNCTION_WINDOW,
) -> StrainBindingScore:
    # Treat IUPAC codes in the consensus as "any of these bases", per IDT's
    # standard handling of mixed-base oligos. A consensus 'M' at position i
    # is NOT a mismatch against a strain's 'A' or 'C' — it would be a
    # mismatch against 'G' or 'T'.
    positions: List[int] = []
    for i, (c, s) in enumerate(zip(consensus_window, strain_window)):
        cu = c.upper()
        su = s.upper()
        allowed = IUPAC_BASES.get(cu)
        if allowed is None:
            # Unknown character (e.g. '-' or 'X' for gaps) — count as mismatch.
            positions.append(i)
            continue
        if su not in allowed:
            positions.append(i)
    n_total = len(positions)
    n_lhs = sum(1 for p in positions if p < half_length)
    n_rhs = sum(1 for p in positions if p >= half_length)
    junction_range = range(half_length - junction_window, half_length + junction_window)
    # Reported unweighted, and still the basis of ``likely_binds``: the graded
    # ΔTm below changes how a mismatch is PRICED, never how many there are.
    n_junction = sum(1 for p in positions if p in junction_range)
    # Every mismatch is priced by how far it sits from the ligation nick — the
    # ramp in constants.mismatch_position_weight, which costs -6.0 C at the nick
    # and -2.0 C at either probe end, the same as the step form
    # ``DTM_PER_MISMATCH * n_total + DTM_JUNCTION_BONUS * n_junction`` at those
    # two extremes. Unlike a step at ``junction_window``, the ramp has a
    # gradient across the rest of the half, so a mid-half mismatch costs more
    # than an end-of-probe one.
    est_dTm = DTM_PER_MISMATCH * sum(
        mismatch_position_weight(p, half_length) for p in positions
    )
    return StrainBindingScore(
        strain_id=strain_id,
        total_mismatches=n_total,
        lhs_mismatches=n_lhs,
        rhs_mismatches=n_rhs,
        junction_mismatches=n_junction,
        mismatch_positions=positions,
        estimated_dTm=est_dTm,
        likely_binds=(n_total <= 2 and n_junction == 0),
    )


def score_probe_against_all_strains(
    alignment: StrainAlignment,
    start: int,
    tile_length: int = 50,
    *,
    probe_seq: Optional[str] = None,
    junction_window: int = JUNCTION_WINDOW,
) -> List[StrainBindingScore]:
    """Score one 50-bp window against every strain in the alignment.

    By default the consensus is treated as the "probe sequence" — fine
    for conserved-mode runs (consensus == every strain in the window).
    For strain-unique mode the caller should pass ``probe_seq`` equal to
    the target strain's actual sequence; otherwise IUPAC codes in the
    consensus would silently match every strain and the per-strain
    mismatch columns in results.csv would be misleading zeros.
    """
    probe_window = (
        probe_seq
        if probe_seq is not None
        else alignment.consensus[start:start + tile_length]
    )
    scores: List[StrainBindingScore] = []
    for sid, seq in zip(alignment.strain_ids, alignment.aligned_seqs):
        strain_window = seq[start:start + tile_length]
        # If the strain has gaps in this region, treat each '-' as a
        # mismatch by substituting an X that won't match anything.
        if "-" in strain_window:
            strain_window = strain_window.replace("-", "X")
        scores.append(_binding_score_for_window(
            probe_window, strain_window, sid,
            half_length=tile_length // 2,   # follow the tile size, not a fixed 25
            junction_window=junction_window,
        ))
    return scores
