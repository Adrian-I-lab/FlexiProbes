# flexiprobes HPC container (Docker / OCI).
#
# An alternative to flexiprobe.def for sites that build OCI images and convert
# them to Apptainer. Build from the app/ directory (the build context):
#
#     cd app
#     docker build --platform linux/amd64 -f packaging/container/Dockerfile -t flexiprobe:latest .
#
# Convert to a Singularity/Apptainer image for the cluster:
#
#     apptainer build flexiprobe.sif docker-daemon://flexiprobe:latest
#
# BASE IMAGE PINNED BY DIGEST. A tag such as
# python:3.11-slim moves: two builds a month apart get different Debian
# packages and a different Python patch release under one name. The digest is
# the MULTI-PLATFORM INDEX digest, looked up on 2026-09-22 with
#     docker buildx imagetools inspect python:3.11-slim
# so the platform is fixed separately, here and on the build command line.
# flexiprobe.def carries the same digest; test_container_pinned checks it.
FROM --platform=linux/amd64 python:3.11-slim@sha256:da047cb8f9d1d98e5c070f5300ba9f7274e33b8fc0e5be5ed88740aed1b95ba9

# BLAST+ (blastn, makeblastdb, blastdbcmd) is the VENDORED 2.17.0 (build Jul 1 2025) the desktop app ships, not
# apt's ncbi-blast+, whose version is whatever the base image's Debian release
# carries on build day. The binaries are dynamically linked against glibc up
# to 2.28 and need libgomp1, libbz2 and zlib at run time; libgomp1 is the one
# the slim base may lack, so it is installed explicitly.
# The sha256 check fails the build if the copied bytes are not the vendored
# ones, and the version check fails it if they are not 2.17.0.
#
# MAFFT (below) is in the image as well. Multi-strain consensus designs align their
# input records, and without a real MSA tool the pipeline falls back to a
# built-in progressive merge that consensus.py itself documents as NOT a true
# multiple alignment. The desktop build vendors MAFFT (consensus._bundled_mafft),
# so a container without it would align cluster runs differently from desktop
# runs of the same input. Installing it here closes that gap for every bundle,
# present and future, at the cost of a few MB in an image that already carries
# BLAST+. Single-species runs never call it; they simply never reach that path.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY packaging/vendor/blast/linux-x64/bin /opt/blast/bin
RUN printf '%s  %s\n' \
        33b64bc67d3149cee2459b2f7766b363323df632cf12c099546de00aea9698b5 /opt/blast/bin/blastn \
        c1ffdcf6f15d1d8d75377cc9f37afa42bc9fa06678903bad77c2641e60778fce /opt/blast/bin/makeblastdb \
        7b7985134ac696b2f28428baf8a210554c4e2d6c38cc1dfdd46d207dc4cc762f /opt/blast/bin/blastdbcmd \
      | sha256sum -c - \
    && ln -s /opt/blast/bin/blastn /usr/local/bin/blastn \
    && ln -s /opt/blast/bin/makeblastdb /usr/local/bin/makeblastdb \
    && blastn -version | grep -q 'blastn: 2.17.0+'

# blastdbcmd is REQUIRED. The BLAST engine uses it to read the subject bases
# opposite the unaligned part of each probe half; without it every count falls
# back to CONSERVATIVE (unaligned positions counted as matches), which
# over-rejects, silently. apt's ncbi-blast+ carried it, so the vendored tree
# now does too (same 2.17.0 tarballs), and the build fails without it.
RUN ln -s /opt/blast/bin/blastdbcmd /usr/local/bin/blastdbcmd \
    && blastdbcmd -version | grep -q 'blastdbcmd: 2.17.0+'

# MAFFT comes from the SAME vendored portable 7.526 the desktop app ships,
# NOT from apt. apt would close the "no MSA tool" gap and open a worse one: a
# cluster run would align with a different MAFFT version from a desktop run of
# the same input. flexiprobe.def does the identical thing, so the two recipes
# cannot drift apart from each other either.
COPY packaging/vendor/mafft/linux-x64 /opt/mafft
RUN printf '#!/bin/sh\nexec bash /opt/mafft/mafft.bat "$@"\n' > /usr/local/bin/mafft \
    && chmod 0755 /usr/local/bin/mafft

WORKDIR /opt/flexiprobe-src
COPY pyproject.toml README.md constraints.txt ./
COPY src ./src
# Each check FAILS THE BUILD if its tool is missing, so an image can never ship
# silently without one. `command -v` rather than `mafft --version` because that
# writes to stderr and its exit status is not dependable across versions.
# -c constraints.txt PINS the runtime deps to the exact versions this release
# was tested against. Without it pip resolves pyproject.toml's deliberately
# loose ranges (e.g. primer3-py>=2.0,<3) and takes whatever is newest on the
# day the image is built, so two bundles built a week apart could screen with
# different thermodynamics. A different primer3-py release changes the golden
# output baselines. Constraints only
# CAP versions, they never add packages, so the GUI pins in the file are
# inert here (this install omits the [gui] extra).
RUN pip install --no-cache-dir -c constraints.txt . \
    && flexiprobes version \
    && blastn -version \
    && (mafft --version 2>&1 | grep -q v7.526)

# BUILD ID: /opt/flexiprobe/BUILD_ID.json names this
# image's inputs, and FLEXIPROBES_CONTAINER_BUILD_ID is set to its content_id
# in every Python process the image starts (see the script's docstring for why
# that is not an ENV line). flexiprobe.def runs the same script.
COPY packaging/container/write_build_id.py /opt/flexiprobe/write_build_id.py
RUN python /opt/flexiprobe/write_build_id.py Dockerfile \
        python:3.11-slim@sha256:da047cb8f9d1d98e5c070f5300ba9f7274e33b8fc0e5be5ed88740aed1b95ba9

ENTRYPOINT ["flexiprobes"]
CMD ["version"]
