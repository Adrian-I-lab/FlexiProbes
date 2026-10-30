"""HEAD must not claim to be a version that has already been RELEASED.

WHY THIS EXISTS. test_version_coherence.py checks that every declaration inside
the tree agrees with every other: pyproject, the two packages and CITATION.cff.
All of them can agree perfectly and still be wrong
in one specific way it cannot see — the tree can claim a version that was cut as
a tag and then diverge from it.

If a version is tagged and further commits land under app/src while
`__version__` stays the same, the released artefact and the working tree are
DIFFERENT ENGINES carrying ONE version string, and the difference can be
observable (for example, a manifest key one writes and the other lacks).
Every internal check stays green.

A version string is a promise that two things built from it are the same thing.
This test enforces that promise across the tag boundary.
"""
from __future__ import annotations

import subprocess

import pytest

from flexprobe import __version__ as engine_version


def _git(*args: str) -> tuple[int, str]:
    try:
        p = subprocess.run(("git",) + args, capture_output=True, text=True,
                           timeout=30)
    except (OSError, subprocess.SubprocessError):          # pragma: no cover
        return 127, ""
    return p.returncode, p.stdout.strip()


def test_head_version_is_not_an_already_released_version():
    """If a tag exists for this version, the tree must match that tag.

    Skips rather than fails when the source is not a git checkout (an installed
    wheel, the packaged app, a release tarball), because there is nothing to
    compare against there and a hard failure would break a legitimate install.
    """
    rc, _ = _git("rev-parse", "--git-dir")
    if rc != 0:
        pytest.skip("not a git checkout; nothing to compare against")

    tag = f"v{engine_version}"
    rc, _ = _git("rev-parse", "--verify", "--quiet", f"{tag}^{{tag}}")
    if rc != 0:
        rc, _ = _git("rev-parse", "--verify", "--quiet", tag)
    if rc != 0:
        # No tag for this version yet: this is an UNRELEASED version, which is
        # the normal state during development. Nothing to enforce.
        pytest.skip(f"{tag} does not exist yet; version is unreleased")

    # ":/app/src" is a ROOT-RELATIVE magic pathspec. A plain "app/src" is
    # resolved against the CURRENT directory, and pytest runs this suite from
    # app/, so it became app/app/src, matched nothing, and the assertion below
    # passed while five files were in fact diverged. A guard that reports green
    # on the defect it was written for is worse than no guard.
    rc, changed = _git("diff", "--name-only", f"{tag}..HEAD", "--", ":/app/src")
    assert rc == 0, "git diff against the tag failed"
    files = [f for f in changed.splitlines() if f.strip()]

    rc, behind = _git("rev-list", "--count", f"{tag}..HEAD")
    n_commits = behind if rc == 0 else "?"

    assert not files, (
        f"__version__ is {engine_version!r} and tag {tag} EXISTS, but the tree "
        f"has diverged from it under app/src:\n"
        + "".join(f"    {f}\n" for f in files)
        + f"  ({n_commits} commit(s) since the tag.)\n"
        f"  Two different engines would ship under one version string, and an "
        f"artefact built from either would be indistinguishable by version. "
        f"BUMP __version__ (and pyproject, CITATION.cff and the bundle gates "
        f"with it) rather than moving the published tag."
    )
