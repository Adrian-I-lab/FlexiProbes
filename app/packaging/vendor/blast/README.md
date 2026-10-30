# Vendored BLAST+ binaries

`blastn` and `makeblastdb` are bundled into the desktop app so end users need
**no separate BLAST+ install**, on every platform. The packaged app finds them
automatically (`flexprobe.blast._bundled_blast_dir`); `FLEXPROBE_BLAST_DIR` still
overrides.

## Layout

One subdirectory per OS/arch (the macOS and Linux `blastn` share a filename, so
they can't live together):

```
vendor/blast/
├── LICENSE.ncbi-blast        # shared — BLAST+ is public domain
├── macos-arm64/bin/          # blastn, makeblastdb, blastdbcmd (Apple Silicon .app)
├── windows-x64/bin/          # the same three, .exe        (Windows .exe)
└── linux-x64/bin/            # blastn, makeblastdb, blastdbcmd (Linux onedir)
```

`packaging/flexprobe.spec` picks the subdir matching the **build host** (by
`sys.platform` + `platform.machine()`) and adds the three tools to PyInstaller's
`binaries` with destination `blast/bin`, so they land at
`<app>/.../blast/bin/`. The build's `codesign --deep` (macOS) signs them with the
rest of the bundle. If a host's subdir is missing, the build still succeeds — it
just won't ship BLAST+ (the spec prints a warning).

## Provenance

All from the NCBI BLAST+ **2.17.0** official builds (build Jul 1 2025),
`ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/`:

| Subdir         | NCBI tarball                              | Binary type            |
|----------------|-------------------------------------------|------------------------|
| `macos-arm64`  | `ncbi-blast-2.17.0+-aarch64-macosx.tar.gz`| Mach-O arm64           |
| `windows-x64`  | `ncbi-blast-2.17.0+-x64-win64.tar.gz`     | PE32+ x86-64           |
| `linux-x64`    | `ncbi-blast-2.17.0+-x64-linux.tar.gz`     | ELF x86-64 (glibc)     |

Three tools per platform: `blastn`, `makeblastdb` and `blastdbcmd`
(`.exe` on Windows).
The BLAST engine uses it to read the subject bases opposite the unaligned part
of each probe half, and without it every off-target count falls back to the
conservative estimate, which over-rejects. Before then no platform shipped it,
so the desktop app's BLAST engine always ran in that fallback.

It came from the same three NCBI tarballs, checked against NCBI's `.md5`
files, and the `blastn` and `makeblastdb` in those tarballs are byte-identical
to the ones already here, so all three tools are one release. sha256:

| Subdir        | blastdbcmd sha256                                                  |
|---------------|--------------------------------------------------------------------|
| `linux-x64`   | `7b7985134ac696b2f28428baf8a210554c4e2d6c38cc1dfdd46d207dc4cc762f` |
| `macos-arm64` | `5f0736c6464bd914c4cf8bf7fcea873bdc9545d6f296de3bdf2f9a4ffeae4808` |
| `windows-x64` | `dc21b96926657847ac2da8583132d65a0390831b442acee261d1095604db39a5` |

No external `lib/` dir is needed, but the Linux binaries are dynamically
linked: they need glibc 2.28 or later plus `libgomp`, `libbz2` and `zlib`
(`libgomp1` is installed explicitly by the container recipes). **License:**
US-Government **public domain** (see `LICENSE.ncbi-blast`); free to
redistribute, no attribution obligation. Per-build metadata is in each
subdir's `ncbi_package_info`, which is NCBI's own file, unchanged, and
identical to the copy inside each tarball.

> Not vendored: `macos-x64` (Intel Mac) and `linux-aarch64`, because no CI build targets
> them. The spec already knows their subdir names; drop the matching NCBI build's
> `blastn`, `makeblastdb` and `blastdbcmd` into `vendor/blast/<plat>/bin/` to enable them.

## Refreshing

```bash
# e.g. macOS arm64
curl -O https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-aarch64-macosx.tar.gz
tar xzf ncbi-blast-2.17.0+-aarch64-macosx.tar.gz
md5 -q ncbi-blast-2.17.0+-aarch64-macosx.tar.gz   # compare with NCBI's .md5
cp ncbi-blast-2.17.0+/bin/{blastn,makeblastdb,blastdbcmd} macos-arm64/bin/
```

After a refresh, update the sha256 values the container recipes check
(`packaging/container/Dockerfile` and `flexiprobe.def`), or their builds fail.
`app/tests/test_container_pinned.py` fails first if they disagree.
