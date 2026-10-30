"""No personal or site-specific path may sit in a tracked text file.

A home directory, a scratch directory or a cluster layout in a tracked file
publishes infrastructure detail, and in a script it also makes the script run
only on one machine. So every file ``git ls-files`` reports is scanned, code,
docs, docstrings, help text and generated artefacts alike.

The built-in patterns are generic. A maintainer adds site-specific fragments
(a cluster project root, a username, a group name) through the
``FLEXIPROBES_PRIVATE_PATTERNS`` environment variable, separated by ``|``, so
those fragments are never themselves committed.

The generic patterns are assembled from fragments so that this file does not
match itself. A guard that exempts its own source by name would also exempt any
path written into it later.
"""
import os
import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

_J = "".join
FORBIDDEN = [
    _J(["/Users/", r"(?!me/|you/|<)[a-z][a-z0-9._-]+/"]),  # macOS home dir
    _J(["/home/", r"(?!me/|you/|<)[a-z][a-z0-9._-]+/"]),   # Linux home dir
    _J(["/private", "/tmp"]),                               # macOS scratch dir
    _J(["/tmp/", "claude-"]),                               # agent scratchpad
    _J(["/mnt/", "lustre"]),                                # cluster mount
]
_EXTRA = [re.escape(p) for p in
          os.environ.get("FLEXIPROBES_PRIVATE_PATTERNS", "").split("|") if p]
PATTERN = re.compile("|".join(FORBIDDEN + _EXTRA))


def _tracked_files():
    try:
        out = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "ls-files", "-z"],
            check=True, capture_output=True).stdout
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("not a git checkout (installed or sdist tree)")
    return [p for p in out.decode().split("\0") if p]


def _offenders():
    """Map each offending tracked text file to its first offending line."""
    files = _tracked_files()
    hits, scanned = {}, 0
    for rel in files:
        path = REPO_ROOT / rel
        if not path.is_file():
            continue
        data = path.read_bytes()
        if b"\0" in data[:8192]:
            continue
        scanned += 1
        text = data.decode("utf-8", errors="replace")
        m = PATTERN.search(text)
        if m:
            line = text.count("\n", 0, m.start()) + 1
            hits[rel] = f"{rel}:{line}: {m.group(0)}"
    # A scan that read nothing would pass. Refuse that.
    assert scanned > 200, f"only {scanned} tracked text files scanned"
    return hits


def test_pattern_matches_each_forbidden_fragment():
    # Positive control: the pattern fires on a realistic line for each
    # generic pattern, and not on a neutral or placeholder path.
    for line in ('cd "/Users/' + 'jsmith/projects/x"', "ls /home/" + "jsmith/data/",
                 "/private" + "/tmp/x", "/tmp/" + "claude-1/x",
                 "/mnt/" + "lustre/x"):
        assert PATTERN.search(line), line
    for line in ("/path/to/project/cellranger", "/Users/me/blast_dbs/human_tx",
                 "/home/you/refs/"):
        assert not PATTERN.search(line), line


def test_no_private_paths_in_tracked_files():
    bad = sorted(_offenders().values())
    assert not bad, (
        f"{len(bad)} tracked file(s) carry a site-specific or personal path. "
        "Use an environment variable in scripts and /path/to/... in docs:\n"
        + "\n".join(bad))
