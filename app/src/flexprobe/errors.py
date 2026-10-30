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

"""Custom exception hierarchy for flexprobe."""


class FlexprobeError(Exception):
    """Base exception for all flexprobe errors."""


class FastaParseError(FlexprobeError):
    """Raised when a FASTA file cannot be parsed."""


class BlastError(FlexprobeError):
    """Raised on BLAST invocation or parsing failures."""


class BlastVersionError(BlastError):
    """Raised when BLAST version is below required minimum."""


class BlastDatabaseError(BlastError):
    """Raised when a BLAST database cannot be built or found."""


class WriterError(FlexprobeError):
    """Raised when an output writer fails."""


class PipelineError(FlexprobeError):
    """Raised for top-level pipeline orchestration failures."""


class CancellationError(PipelineError):
    """Raised when a pipeline run is cancelled."""


class ConfigError(FlexprobeError):
    """Raised on invalid configuration."""


class KmerIndexStampMismatchError(FlexprobeError):
    """A k-mer index on disk was built by different builder code, or has no stamp.

    Raised instead of reusing or silently rebuilding it. ``cached`` is the
    ``builder_sha256`` in the index's ``.meta.json`` (``None`` when absent) and
    ``current`` is the hash of this build's ``kmer_index.py``. Both are in the
    message. Nothing is deleted.
    """

    def __init__(self, prefix: str, cached, current) -> None:
        self.prefix = prefix
        self.cached = cached
        self.current = current
        what = ("carries no builder stamp" if cached is None
                else "was built by different k-mer index code")
        super().__init__(
            f"k-mer index {prefix}.* {what}. It has NOT been deleted or "
            "reused. Re-run with --rebuild-kmer-index to rebuild it with this "
            "build's code, or use the build that wrote it.\n"
            f"  index builder_sha256:   {cached}\n"
            f"  current builder_sha256: {current}"
        )


class CacheStampMismatchError(FlexprobeError):
    """A resume checkpoint was written by different engine code, or carries no stamp.

    Raised instead of discarding the checkpoint. ``cached`` is the stamp read
    from the checkpoint (``None`` when it has none) and ``current`` is the stamp
    this build would write. Both are in the message.
    """

    def __init__(self, cache_path: str, cached, current) -> None:
        import json

        self.cache_path = cache_path
        self.cached = cached
        self.current = current
        what = ("carries no engine stamp" if cached is None
                else "was written by a different engine")
        super().__init__(
            f"resume checkpoint {cache_path} {what}. It has NOT been deleted. "
            "Resume with the engine that wrote it, or remove that directory "
            "(or pass --no-resume) to recompute the stage.\n"
            f"  checkpoint stamp: {json.dumps(cached, sort_keys=True)}\n"
            f"  current stamp:    {json.dumps(current, sort_keys=True)}"
        )
