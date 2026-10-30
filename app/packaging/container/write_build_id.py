"""Write /opt/flexiprobe/BUILD_ID.json inside a container build.

Run by BOTH Dockerfile and flexiprobe.def, so the two images describe
themselves identically. It names the image's inputs, so a run and its resume
caches can say which build made them.

``content_id`` is a sha256 over every field except ``build_time_utc`` and
``recipe``: two images built from the same base, source, BLAST and MAFFT carry
the same content ID whichever recipe built them, and differ only in build time.

It is also what ``FLEXIPROBES_CONTAINER_BUILD_ID`` is set to, the variable the
resume-cache stamps read (``engine_identity.CONTAINER_BUILD_ID_ENV``). A Docker
``ENV`` line cannot hold it, because it is computed during the build, and an
entrypoint wrapper would be skipped by ``apptainer exec``, which is how the
cluster runs the image. So a ``.pth`` file in site-packages sets it (if unset)
at every interpreter start. That covers docker run, apptainer exec and run,
and both recipes identically. The plain value is also written to
``/opt/flexiprobe/BUILD_ID`` for shells.

    python write_build_id.py <recipe name> <base image reference>
"""
import hashlib
import json
import os
import platform
import subprocess
import sys
import time

OUT = "/opt/flexiprobe/BUILD_ID.json"
OUT_PLAIN = "/opt/flexiprobe/BUILD_ID"
ENV_NAME = "FLEXIPROBES_CONTAINER_BUILD_ID"
PTH_NAME = "flexiprobes_container_build_id.pth"


def tree_sha256(root):
    """sha256 over every file under ``root``: relative path, NUL, bytes."""
    h = hashlib.sha256()
    for d, dirs, files in os.walk(root):
        dirs.sort()
        for f in sorted(files):
            p = os.path.join(d, f)
            h.update(os.path.relpath(p, root).encode() + b"\0")
            with open(p, "rb") as fh:
                h.update(fh.read())
    return h.hexdigest()


def main(recipe, base_image):
    ident = {
        "build_time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "recipe": recipe,
        "base_image": base_image,
        "flexiprobes_version": subprocess.run(
            ["flexiprobes", "version"], capture_output=True, text=True,
            check=True).stdout.strip(),
        "python": platform.python_version(),
        "source_sha256": tree_sha256("/opt/flexiprobe-src/src"),
        "blast_sha256": tree_sha256("/opt/blast/bin"),
        "mafft_sha256": tree_sha256("/opt/mafft"),
    }
    stable = {k: v for k, v in ident.items() if k not in ("build_time_utc", "recipe")}
    ident["content_id"] = hashlib.sha256(
        json.dumps(stable, sort_keys=True).encode()).hexdigest()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(ident, fh, indent=2, sort_keys=True)
    with open(OUT_PLAIN, "w") as fh:
        fh.write(ident["content_id"] + "\n")
    write_pth(ident["content_id"], site_dir())
    print(json.dumps(ident, indent=2, sort_keys=True))


def site_dir():
    import sysconfig
    return sysconfig.get_paths()["purelib"]


def pth_line(content_id):
    """One ``import`` line: site.py executes it at every interpreter start."""
    return (f"import os; os.environ.setdefault({ENV_NAME!r}, {content_id!r})\n")


def write_pth(content_id, directory):
    path = os.path.join(directory, PTH_NAME)
    with open(path, "w") as fh:
        fh.write(pth_line(content_id))
    return path


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
