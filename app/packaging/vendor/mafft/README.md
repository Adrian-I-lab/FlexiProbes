# Vendored MAFFT

MAFFT is bundled into the desktop app so multi-strain runs get a real multiple
sequence alignment with **no separate install**, on all three platforms. The
packaged app finds it automatically (`flexprobe.consensus._bundled_mafft`);
`FLEXPROBE_MAFFT_DIR` overrides, and PATH plus the common install directories
are searched after that.

## What is vendored, and what must not be

MAFFT's official **all-in-one / portable** packages. Upstream recommends them
for precisely this use: the Linux and macOS download pages both call the
package "Recommended only when redistributing MAFFT within another program
package".

Do **not** vendor a Homebrew, apt, or conda build. Those are compiled against
their build host, are not redistribution artefacts, and use a different layout
(`bin/mafft` plus a separate `libexec/mafft`) that the runtime resolver
deliberately does not look for.

## Layout

Every package puts a driver named `mafft.bat` at the top of a relocatable
folder, and the driver locates its helpers relative to itself. **What sits
beside the driver differs by platform**, so the resolver and the spec key on
`mafft.bat` alone. Requiring `mafftdir` as well, as an earlier revision did,
would have failed silently on Windows.

```
vendor/mafft/
├── LICENSE.mafft
├── macos/          mafft.bat + mafftdir/{bin,libexec}      (from mafft-mac)
├── linux-x64/      mafft.bat + mafftdir/{bin,libexec}      (from mafft-linux64)
└── windows-x64/    mafft.bat + mafft-signed.ps1 + usr/ + tmp/   (from mafft-win)
```

Unpack each archive and move the **contents** of its folder into the platform
directory, so `mafft.bat` sits directly under it.

On macOS and Linux `mafft.bat` is a Bourne shell script despite the extension,
so `consensus._external_msa` invokes it through an explicit interpreter. On
Windows the helpers live under `usr/` inside a partial MSYS tree, not Cygwin.

`packaging/flexprobe.spec` picks the subdirectory matching the **build host**,
walks it, and sends executables to PyInstaller's `binaries` (which preserves
the executable bit) and scripts and data to `datas`, keeping the relative
layout under `mafft/`. macOS `codesign --deep` signs the executables with the
rest of the bundle. A missing subdirectory does not fail the build, it warns
and ships without MAFFT.

`windows-x64/tmp/` ships empty and the driver writes into it. It carries a
`.gitkeep` because Git cannot track an empty directory.

## Line endings

`.gitattributes` marks this whole tree `-text`. Upstream's Windows page says
"Make sure to keep LF", and Git's end-of-line conversion would rewrite
`mafft.bat` to CRLF on a Windows checkout and break the driver. Do not remove
that rule.

## Provenance

Version **7.526 on every platform**, deliberately. Upstream's macOS page marks
the 7.526 all-in-one as alpha and offers 7.490 as stable, but a per-platform
version split would reintroduce the cross-platform reproducibility problem
this bundle exists to remove. One version everywhere is worth more than a
stable label on one of the three.

Downloaded from <https://mafft.cbrc.jp/alignment/software/> on 2026-09-10.

| Platform | Archive | SHA-256 | Bytes | Unpacks to |
|---|---|---|--:|---|
| macOS | `mafft-7.526-mac.zip` | `704e0d2284f87b1b40e7f065359880967649136d4d328809b00288e2370be71f` | 21,442,050 | `mafft-mac/` |
| Linux | `mafft-7.526-linux.tgz` | `bf59d016f1b2030bc7fc83b7715ae1f0823570de3bd059c26eb522a72f9c2952` | 16,021,747 | `mafft-linux64/` |
| Windows | `mafft-7.526-win64-signed.zip` | `7c3a508c8512af04007d68019c6fbd4ff00d891244b58f87ec9aa0e53fbae8a7` | 27,476,513 | `mafft-win/` |

Verified on 2026-09-10:

* **macOS binaries are universal** (`Mach-O universal binary with 2
  architectures: x86_64, arm64`), so no Rosetta dependency on Apple Silicon.
* The macOS package produces a **byte-identical** alignment to a Homebrew
  MAFFT v7.526 on the OROV L segment (6,852 columns).
* Resolved and executed from a simulated frozen bundle, reporting
  `MAFFT v7.526 (2024/Apr/26) (bundled)`.
* The spec collects 36 executables and 5 data files for macOS.

Unpacked footprint: 45 MB (macOS), 34 MB (Linux), 70 MB (Windows), 149 MB
total. For scale, the existing `vendor/blast/` tree is 137 MB. The packages are
shipped **unmodified**. `dash_client` accounts for roughly 23 MB across the
three and is not reachable from `mafft --auto`, but trimming an upstream
package costs a clean provenance statement for a 15% saving and was not done.

## Licence

MAFFT core is BSD, Copyright (c) 2009 Kazutaka Katoh. Compatible with this
project's GPL-3.0-or-later.

Note the all-in-one archives ship **no licence file of their own**. The text in
`LICENSE.mafft` here, and its copy in `app/THIRD-PARTY-LICENSES/mafft-BSD.txt`,
was taken from the MAFFT distribution's own `LICENSE` (sha256
`5b88dda69a361f21c241a3938c8b454d327c87137499acd963ffa13ad7f9426d`). Because
the redistributable archive carries none, including it here is not optional.
