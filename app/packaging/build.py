#!/usr/bin/env python3
"""build.py – CLI helper for flexiprobes packaging.

Usage:
    python packaging/build.py [--clean] [--spec SPEC_PATH]   # from app/

⚠️ BUILD FROM A CLEAN VENV with the app and PyInstaller installed
(``pip install -e ".[gui]" pyinstaller -c constraints.txt``), NEVER the
system or anaconda Python.

PyInstaller bundles whatever Qt it finds in the interpreter that runs it. The
anaconda Python carries its own **Qt5**, which then gets collected alongside
PySide6's Qt6. ``libqcocoa.dylib`` is built against Qt6, resolves the Qt5
libraries first, fails to load, and Qt calls qFatal — so the .app dies at
launch with "flexiprobes quit unexpectedly" (abort inside
QGuiApplicationPrivate::createPlatformIntegration). The same anaconda/Qt
conflict makes the GUI tests segfault.

The build REPORTS SUCCESS either way — only launching reveals it. Verify:

    find dist/flexiprobes.app -iname 'libQt5*' | wc -l    # must be 0
    open dist/flexiprobes.app                             # must stay up

This script wraps PyInstaller invocation with sensible defaults and
post-build checks so that CI and local developer workflows stay
consistent.
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_SPEC = Path(__file__).resolve().parent / "flexprobe.spec"
# Always emit the built artifact to the repository-root dist/ (never app/dist/),
# regardless of the current working directory the build is launched from.
# build.py lives at <repo>/app/packaging/build.py, so parents[2] is the repo root.
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DIST = REPO_ROOT / "dist"


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Build flexiprobes desktop binary.")
    parser.add_argument(
        "--spec",
        type=Path,
        default=DEFAULT_SPEC,
        help="Path to the PyInstaller .spec file (default: packaging/flexprobe.spec)",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        default=False,
        help="Remove build/ and dist/ directories before building.",
    )
    parser.add_argument(
        "--dist-dir",
        type=Path,
        default=DEFAULT_DIST,
        help="Output directory for the built artifact "
             "(default: the repository-root dist/, not app/dist/).",
    )
    parser.add_argument(
        "--sign",
        default=os.environ.get("FLEXIPROBES_CODESIGN_IDENTITY"),
        help="macOS only: codesign the built .app with this identity "
             "(e.g. 'Developer ID Application: Name (TEAMID)'), or '-' for "
             "ad-hoc signing. Defaults to $FLEXIPROBES_CODESIGN_IDENTITY. Notarization "
             "(notarytool) is a separate manual step — see app/README.md.",
    )
    return parser.parse_args(argv)


def codesign_macos(app_path: Path, identity: str) -> None:
    """Deep-sign a macOS .app bundle.

    ``identity='-'`` does ad-hoc signing for local use: NO hardened runtime, so
    the PyInstaller/Python bundle still launches (hardened runtime blocks
    Python's JIT/unsigned executable memory unless granted entitlements).

    A real Developer ID identity additionally enables the hardened runtime +
    secure timestamp (required for notarization) and applies the entitlements
    plist that grants Python the JIT/unsigned-memory/library-validation
    exceptions it needs to run under hardened runtime.
    """
    adhoc = identity.strip() == "-"
    cmd = ["codesign", "--force", "--deep", "--sign", identity]
    if not adhoc:
        ents = Path(__file__).resolve().parent / "entitlements.plist"
        cmd += ["--timestamp", "--options", "runtime"]
        if ents.exists():
            cmd += ["--entitlements", str(ents)]
    cmd.append(str(app_path))
    print(f"Signing ({'ad-hoc' if adhoc else 'Developer ID'}): {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        sys.exit(f"codesign failed with code {result.returncode}")
    verify = subprocess.run(
        ["codesign", "--verify", "--verbose=2", str(app_path)], check=False,
    )
    if verify.returncode != 0:
        sys.exit("codesign verification failed")
    print("OK  Signature verified.")


def clean(dist_dir: Path):
    """Remove previous build artefacts.

    Removes the PyInstaller ``build/`` work dir and the previously-built app
    artifacts inside ``dist_dir`` — sibling release files in ``dist_dir`` (the
    distributable ``.zip`` + checksum + source offer) are left untouched, since
    ``dist_dir`` is now the shared repository-root ``dist/``.
    """
    for d in [Path("build"),
              dist_dir / "flexiprobes.app",
              dist_dir / "flexiprobes"]:
        if d.exists():
            print(f"Removing {d} ...")
            shutil.rmtree(d)


def build(spec: Path, dist_dir: Path):
    """Run PyInstaller with the given spec and output directory."""
    if not spec.exists():
        sys.exit(f"ERROR: spec file not found: {spec}")

    cmd = [
        sys.executable,
        "-m",
        "PyInstaller",
        str(spec),
        "--noconfirm",
        "--clean",
        "--distpath",
        str(dist_dir),
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        sys.exit(f"PyInstaller exited with code {result.returncode}")


def verify(dist_dir: Path):
    """Basic post-build sanity checks."""
    system = platform.system()
    # Windows/Linux now build "onedir" (a flexiprobes/ folder) for LGPL
    # relinkability; the launcher lives inside that folder.
    if system == "Darwin":
        expected = dist_dir / "flexiprobes.app"
    elif system == "Windows":
        expected = dist_dir / "flexiprobes" / "flexiprobes.exe"
    else:
        expected = dist_dir / "flexiprobes" / "flexiprobes"

    if expected.exists():
        size_mb = expected.stat().st_size / (1024 * 1024) if expected.is_file() else 0
        print(f"OK  Artifact found: {expected}")
        if size_mb:
            print(f"    Size: {size_mb:.1f} MB")
    else:
        sys.exit(f"FAIL  Expected artifact not found: {expected}")

    _check_no_qt5_contamination(expected)


def _check_no_qt5_contamination(bundle: Path) -> None:
    """Fail the build if a foreign Qt5 got collected alongside PySide6's Qt6.

    Building with the anaconda Python pulls in ITS Qt5. libqcocoa.dylib is built
    against Qt6, resolves the Qt5 libs first, fails to load, and Qt qFatal()s —
    the .app aborts at launch with "quit unexpectedly" while the build itself
    reports success. This check makes that failure loud. Build from a clean
    venv, never the anaconda Python.
    """
    strays = sorted(p.name for p in bundle.rglob("libQt5*"))
    if not strays:
        print("OK  No foreign Qt5 libraries in the bundle.")
        return
    shown = ", ".join(strays[:4]) + (" …" if len(strays) > 4 else "")
    sys.exit(
        f"FAIL  {len(strays)} Qt5 librar{'y' if len(strays) == 1 else 'ies'} "
        f"collected alongside PySide6's Qt6: {shown}\n"
        f"      This bundle WILL abort at launch (qFatal: no 'cocoa' platform "
        f"plugin).\n"
        f"      You are building with the wrong interpreter — almost certainly "
        f"anaconda.\n"
        f"      Rebuild with:  ../.venv-build/bin/python build.py --clean"
    )


def main(argv=None):
    args = _parse_args(argv)
    if args.clean:
        clean(args.dist_dir)
    build(args.spec, args.dist_dir)
    verify(args.dist_dir)
    if args.sign:
        if platform.system() != "Darwin":
            print("--sign is macOS-only; ignoring on this platform.")
        else:
            app = args.dist_dir / "flexiprobes.app"
            if app.exists():
                codesign_macos(app, args.sign)
            else:
                sys.exit(f"Cannot sign: {app} not found")
    print("Build complete.")


if __name__ == "__main__":
    main()
