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

"""What is actually IN the off-target reference, and what that does to a run.

The screen answers "does this probe half match anything else?" — a question
whose answer is a property of the FASTA as much as of the probe. A reference's
SHA-256 says two runs used the same file, not what that file contained. This
module classifies each record as it is loaded, puts the tally in the manifest,
and excludes the one class of record known to corrupt results.

THE TWO AXES
------------
Both are visible when 10x's public Human Chromium Transcriptome Probe Set
v1.1.0 is screened against Ensembl ``Homo_sapiens.GRCh38.cdna.all.fa``, and
they push in OPPOSITE directions.

**Axis 1 — ALT haplotypes and assembly PATCH scaffolds.** Ensembl's
``cdna.all`` ships them; 10x's Cell Ranger references do not. Counted from
the public GRCh38 file:

    all of cdna.all      : 207,249 records, 390,035,142 bases
    primary chromosomes  : 187,501 records, 360,342,996 bases (92.4%)
    ALT + PATCH          :  19,713 records,  29,635,255 bases (7.6%)
    unplaced scaffolds   :      35 records,       56,891 bases (0.01%)

They are 7.6% of bases, yet they account for the spurious perfect-50 bp
off-target calls that survive a locus-overlap filter. That is expected, not
surprising: an ALT or PATCH record is a near-verbatim COPY of a primary region,
so it manufactures perfect matches at a rate nothing else in the file
approaches. GRK1, for instance, is ENSG00000185974 on chromosome 13 *and*
ENSG00000288263 on HG1524_PATCH; a probe against the real gene "matches two
genes" and is thrown away. These are FALSE NEGATIVES (a usable probe called
bad), and they are why :data:`EXCLUDE_ALT_SCAFFOLDS_DEFAULT` is True.

**Axis 1b - unplaced scaffolds are NOT that class, and are retained.** An
unplaced or unlocalised scaffold is sequence that could not be assigned to a
chromosome. It is UNIQUE, not a copy, so removing it cannot suppress a spurious
self-match and can only hide a real off-target, which is a bad probe that
PASSES. The public file bears this out: of GRCh38 cdna.all's 35 unplaced
records, **0** carry a gene symbol that also appears on a chromosome, against
**23 of the 27** ``HG*_NOVEL_TEST`` patch records that do.

It matters most away from human. Ensembl's mouse GRCm39 and rat mRatBN7.2
cdna.all carry **no ALT or PATCH record at all**, so a rule that dropped every
non-primary record would remove only unique sequence from them (54 and 69
records), including 5 rat pseudogenes that axis 2 below says must never be
excluded. Under this module those references have nothing excluded and are
searched whole.

**Axis 2 — gene biotype.** Cell Ranger references apply a biotype allowlist;
Ensembl's ``cdna.all`` does not. Pseudogenes are 3.6% of GRCh38 primary bases
(13,077,569 bp) yet are over-represented among real off-target findings. Here
the arithmetic reverses. A TRANSCRIBED pseudogene that is ABSENT from the
reference is an off-target the screen cannot see, so a bad probe passes — a
FALSE POSITIVE. **Pseudogenes are therefore reported and never excluded.**
That is a deliberate position, not an omission. It holds because ALT/PATCH is
the only excluded class, and a pseudogene there is a copy of one that is kept.

GRACEFUL DEGRADATION IS THE LOAD-BEARING REQUIREMENT
----------------------------------------------------
Most runs screen against a small custom FASTA — a viral panel, a handful of
transcripts, a user's own file — whose headers carry none of these fields. Such
a reference classifies as :data:`UNKNOWN` throughout, the exclusion becomes a
no-op, and NOTHING is warned about. A viral panel run must not start lecturing
about scaffolds it does not have. The rule that makes this safe: a record is
only ever dropped when the same file also proved it has parseable primary
records, so the filter can never empty a reference it did not understand.

WHAT AN ENSEMBL HEADER LOOKS LIKE
---------------------------------
    >ENST00000415118.1 cdna chromosome:GRCh38:14:22438547:22438554:1 \
gene:ENSG00000223997.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene \
gene_symbol:TRDD1 description:...
    >ENST00000632585.1 cdna scaffold:GRCh38:HG2290_PATCH:6058:6610:1 \
gene:ENSG00000282172.1 gene_biotype:IG_V_gene ...

    >ENSRNOT00000094517.1 cdna primary_assembly:mRatBN7.2:15:26770499:26770801:1 \
gene:ENSRNOG00000064161.1 gene_biotype:TR_V_gene transcript_biotype:TR_V_gene

The location field is ``<region_type>:<build>:<name>:<start>:<end>:<strand>``.
PRIMARY means ``region_type`` is in :data:`_PRIMARY_REGION_TYPES`
(``chromosome`` or ``primary_assembly``) AND ``name`` matches
``^(\\d{1,2}|X|Y|MT)$`` — which is what makes this work unchanged for mouse
(1-19), rat (1-20) or any other Ensembl species, and what puts both
``scaffold:...PATCH`` and the ALT contigs
(``scaffold:GRCh38:HSCHR6_MHC_APD_CTG1:...``) on the non-primary side.

BOTH halves are load-bearing, and rat is why. Ensembl labels human and mouse
``chromosome`` but rat ``primary_assembly``, and it uses ``primary_assembly``
for rat's UNPLACED SCAFFOLDS too, so the region type alone cannot decide. The
name check is what separates them. Checked against the real files:

    Homo_sapiens.GRCh38.cdna.all.fa      187,501 primary /  19,748 non-primary
    Mus_musculus.GRCm39.cdna.all.fa      116,819 primary /      54 non-primary
    Rattus_norvegicus.mRatBN7.2.cdna.all  46,787 primary /      69 non-primary

which is exactly what the region types and names in those files say. Neither
the human nor the mouse file contains a single record of region type
``primary_assembly``.

COST
----
Header parsing happens in the SAME single pass the loader already makes and
touches only the header branch, so per-base work is unchanged. The extra cost
is :func:`classify_header` over every header plus deferring the keep/drop
decision to a second pass over the per-record arrays: a fraction of a second
on GRCh38 cdna.all, with no change in peak memory.

In context that is free twice over. The index build that follows a cache miss
takes minutes and a whole-transcriptome screen takes far longer, and with the
exclusion ON the index is built over 7.6% fewer bases on GRCh38, so the
DEFAULT path is faster end to end, not slower.
"""
from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger("flexprobe.reference_qc")

# Assembly classes. Small ints rather than strings because the loader keeps one
# per record for a 200k-record reference; :data:`ASSEMBLY_LABELS` is what
# reaches a human or a manifest.
PRIMARY = 0
ALT_PATCH = 1
UNPLACED = 2
UNKNOWN = 3
N_ASSEMBLY_CLASSES = 4
ASSEMBLY_LABELS = ("primary", "alt_patch", "unplaced", "unknown")

#: What :func:`should_exclude` may drop. Named separately from the class so the
#: exclusion's SCOPE is a stated fact rather than an inference from a number.
EXCLUDABLE = (ALT_PATCH,)

#: Compatibility alias. It names ONLY the class that is dropped (ALT/PATCH),
#: not "everything off a primary chromosome": unplaced scaffolds are retained.
#: The union of ALT/PATCH and unplaced is reported in the manifest as
#: ``non_primary_*``.
NON_PRIMARY = ALT_PATCH

#: Recorded for a record whose header carries no ``gene_biotype:`` field. Named
#: rather than "" so a manifest reader can tell "no biotype in the header" from
#: a biotype that happens to be empty.
BIOTYPE_UNKNOWN = "unknown"

#: Exclude ALT/PATCH records from the off-target reference by default.
#: ON because the exclusion removes spurious self-matches (module docstring,
#: axis 1), and OFF is available for reproducing a run made without it.
EXCLUDE_ALT_SCAFFOLDS_DEFAULT = True

#: Warn when non-primary records are being SEARCHED and they are at least this
#: fraction of the reference's bases. 1% rather than 0 so a reference with a
#: stray unplaced scaffold does not produce a warning nobody can act on; GRCh38
#: cdna.all is at 7.6%, i.e. well clear of it.
ALT_PATCH_WARN_BASE_FRACTION = 0.01

# Ensembl location field: <region_type>:<build>:<name>:<start>:<end>:<strand>.
# Matched structurally (6 colon-separated parts, numeric coordinates) rather
# than by a list of region-type words, because Ensembl uses several
# ("chromosome", "scaffold", "supercontig", "primary_assembly", ...) and a
# missing one would silently classify a whole species as unknown.
_LOCATION_PARTS = 6
_BIOTYPE_PREFIX = "gene_biotype:"
_PRIMARY_NAMES = frozenset({"X", "Y", "MT"})
#: Region types that denote the primary assembly. Ensembl uses "chromosome" for
#: human and mouse and "primary_assembly" for others — rat's entire
#: ``mRatBN7.2.cdna.all`` (46,856 records) is the latter, with no "chromosome"
#: record at all. Recognising only "chromosome" would put every rat record on
#: the non-primary side. The NAME check still applies on top, so rat's 69
#: unplaced ``primary_assembly:...:JACYVU...`` scaffolds stay non-primary.
_PRIMARY_REGION_TYPES = frozenset({"chromosome", "primary_assembly"})


def _is_alt_or_patch(name: str) -> bool:
    """Does this sequence name identify an ALT haplotype or an assembly PATCH?

    THE POINT OF THE WHOLE MODULE IS THIS ONE DISTINCTION, so the rule is
    positive: a name is excluded only when it is RECOGNISED as an alternate
    representation, and anything unrecognised is retained and searched. The
    failure mode of a wrong "yes" is a good probe rejected; the failure mode of
    a wrong "no" is nothing at all, because the record simply stays in the
    reference where it always was.

    The four shapes, all seen in Ensembl ``Homo_sapiens.GRCh38.cdna.all.fa``,
    which of the human, mouse and rat references is the only one that carries
    any of them:

        ``HSCHR6_MHC_APD_CTG1``      ALT haplotype
        ``HG1524_PATCH``             GRC FIX patch
        ``HG151_NOVEL_TEST``         GRC NOVEL patch
        ``*_alt``                    Ensembl's other ALT suffix

    ``HG<digits>_`` is deliberately matched as a prefix rather than only
    ``_PATCH``, and the 27 ``HG142_HG150_NOVEL_TEST`` / ``HG151_NOVEL_TEST``
    records are why. They carry no ``PATCH`` in the name, but **23 of the 27
    are a second copy of a gene that is already on a chromosome** (an olfactory
    receptor cluster: OR8J2, OR5T1, OR9G4, ...), so they behave exactly like a
    PATCH and must be excluded with them. The true unplaced scaffolds are the
    control: of their 35 records, **0** duplicate a primary gene symbol.
    """
    up = name.upper()
    if up.startswith("HSCHR") or "_PATCH" in up or up.endswith("_ALT"):
        return True
    # GRC patch scaffolds are HG<digits>_<something>. Spelled out rather than
    # regexed so "HG" alone, or "HGSOMETHING", cannot match.
    if up.startswith("HG"):
        i = 2
        while i < len(up) and up[i].isdigit():
            i += 1
        return i > 2 and i < len(up) and up[i] == "_"
    return False


def _is_primary_name(name: str) -> bool:
    """``^(\\d{1,2}|X|Y|MT)$`` — the primary assembly's own sequence names.

    ``isascii`` guards ``isdigit``, which is True for superscripts and other
    unicode digit forms that are not chromosome names.
    """
    if name in _PRIMARY_NAMES:
        return True
    return len(name) <= 2 and name.isascii() and name.isdigit()


def classify_header(header: str) -> Tuple[int, str]:
    """``(assembly_class, gene_biotype)`` for one FASTA header (no leading ``>``).

    Returns ``(UNKNOWN, BIOTYPE_UNKNOWN)`` for anything that is not an
    Ensembl-style header — which is the common case, and must stay silent.

    Written as a truncate-then-split rather than a regex because it is called
    once per record on the load hot path, where it runs about twice as fast as
    the equivalent regex pair. The truncation at
    ``description:`` matters — that field is most of the line's length and none
    of its information.
    """
    cut = header.find(" description:")
    if cut > 0:
        header = header[:cut]
    assembly = UNKNOWN
    biotype = BIOTYPE_UNKNOWN
    for tok in header.split(" "):
        if tok.startswith(_BIOTYPE_PREFIX):
            biotype = tok[len(_BIOTYPE_PREFIX):]
        elif assembly == UNKNOWN and ":" in tok:
            parts = tok.split(":")
            # Coordinates are what identify this as a location field; a
            # ``gene:ENSG...`` or ``description:`` token cannot impersonate it.
            if (len(parts) == _LOCATION_PARTS
                    and parts[3].isdigit() and parts[4].isdigit()):
                if (parts[0] in _PRIMARY_REGION_TYPES
                        and _is_primary_name(parts[2])):
                    assembly = PRIMARY
                elif _is_alt_or_patch(parts[2]):
                    assembly = ALT_PATCH
                else:
                    # Unplaced or unlocalised scaffold, or a name shape nothing
                    # here recognises. Searchable either way: see
                    # :func:`_is_alt_or_patch` on why the default is retain.
                    assembly = UNPLACED
    return assembly, biotype


def location_of(header: str) -> Tuple[str, str]:
    """``(region_type, sequence_name)``, or ``("", "")`` for a non-Ensembl header.

    Only used for reporting (how many distinct scaffolds of each class a
    reference carries, and which region types it used), so it re-parses rather
    than complicating the hot-path return of :func:`classify_header`. The
    loader already calls it once per NON-primary record and no others, so on a
    healthy reference it runs over a few percent of the file.
    """
    cut = header.find(" description:")
    if cut > 0:
        header = header[:cut]
    for tok in header.split(" "):
        parts = tok.split(":")
        if (len(parts) == _LOCATION_PARTS
                and parts[3].isdigit() and parts[4].isdigit()):
            return parts[0], parts[2]
    return "", ""


def seq_region_of(header: str) -> str:
    """The location field's sequence name, or ``""`` for a non-Ensembl header."""
    return location_of(header)[1]


def is_pseudogene(biotype: str) -> bool:
    """Ensembl spells ~10 biotypes with ``pseudogene`` in the name.

    ``processed_pseudogene``, ``unprocessed_pseudogene``,
    ``transcribed_*_pseudogene``, ``polymorphic_pseudogene``,
    ``IG_V_pseudogene``, ``rRNA_pseudogene``, ... Matching the substring is
    deliberate: the list changes between Ensembl releases and a hard-coded set
    would quietly under-count on the next one.
    """
    return "pseudogene" in biotype


class ReferenceComposition:
    """Running tally of what a reference FASTA contains, by record class.

    Accumulated during the loader's existing single pass (:meth:`add` per
    record), then finalised once the keep/drop decision is known
    (:meth:`finalise`). Sizes are small — three biotype dicts of a few dozen
    entries each — so this costs nothing next to the sequence itself.
    """

    def __init__(self, path: str) -> None:
        self.path = path
        self.records = [0] * N_ASSEMBLY_CLASSES   # indexed by assembly class
        self.bases = [0] * N_ASSEMBLY_CLASSES
        # {assembly class: {gene_biotype: [records, bases]}}. Kept per class,
        # not pooled, so the biotype view can be reported for the records
        # actually SEARCHED — which is the set a user's off-target verdicts
        # came from — rather than for a file the run partly ignored.
        self._biotypes: List[Dict[str, List[int]]] = [
            {} for _ in range(N_ASSEMBLY_CLASSES)]
        # Distinct sequence names per class. Kept per class, not pooled, so the
        # report can say how many ALT/PATCH scaffolds were dropped WITHOUT
        # counting the unplaced ones that were kept.
        self._regions: List[Set[str]] = [
            set() for _ in range(N_ASSEMBLY_CLASSES)]
        # Region-type words actually seen ("chromosome", "primary_assembly",
        # "scaffold", ...). Diagnostic only, and the thing that makes
        # :func:`warn_if_no_primary` say WHY rather than just THAT.
        self._region_types: Set[str] = set()
        # Set by :meth:`finalise`; None means "the loader never decided", which
        # is the state a header-only scan for the manifest stays in.
        self.excluded_non_primary: Optional[bool] = None
        # SHA-256 of the source FASTA, when whoever loaded it already had one.
        # This is the only place a k-mer run's reference gets HASHED in the
        # manifest: ``RunManifest.db_sha256s`` digests a BLAST database's
        # on-disk index files, and a k-mer run has none — it reads the FASTA
        # directly — so that field is ``{}`` for every such run.
        # None (and an absent key) honestly means "not computed here".
        self.source_sha256: Optional[str] = None

    # -- accumulation ------------------------------------------------------

    def add(self, assembly: int, biotype: str, bases: int,
            seq_region: str = "", region_type: str = "") -> None:
        self.records[assembly] += 1
        self.bases[assembly] += bases
        slot = self._biotypes[assembly].setdefault(biotype, [0, 0])
        slot[0] += 1
        slot[1] += bases
        if seq_region and assembly in (ALT_PATCH, UNPLACED):
            self._regions[assembly].add(seq_region)
        if region_type:
            self._region_types.add(region_type)

    def finalise(self, *, excluded_non_primary: bool) -> "ReferenceComposition":
        self.excluded_non_primary = excluded_non_primary
        return self

    # -- derived views -----------------------------------------------------

    @property
    def total_records(self) -> int:
        return sum(self.records)

    @property
    def total_bases(self) -> int:
        return sum(self.bases)

    @property
    def ensembl_style(self) -> bool:
        """True when the headers carried parseable Ensembl location fields.

        The gate on every behaviour in this module. False for a viral panel, a
        user's own FASTA, or anything else without the fields — and then the
        exclusion is a no-op and no warning is emitted.
        """
        return (self.records[PRIMARY] + self.records[ALT_PATCH]
                + self.records[UNPLACED]) > 0

    @property
    def non_primary_seq_regions(self) -> int:
        """How many distinct ALT/PATCH/scaffold sequence names are present.

        A count of records says how much there is; a count of scaffolds says
        how spread out it is, which is what tells a reader this is the
        assembly's alternate-locus machinery rather than one odd contig.
        """
        return len(self._regions[ALT_PATCH] | self._regions[UNPLACED])

    @property
    def region_types_seen(self) -> Set[str]:
        """Ensembl region-type words observed, for diagnostics."""
        return set(self._region_types)

    @property
    def alt_patch_seq_regions(self) -> int:
        """Distinct ALT/PATCH scaffolds, i.e. the ones actually excluded."""
        return len(self._regions[ALT_PATCH])

    @property
    def unplaced_seq_regions(self) -> int:
        """Distinct unplaced/unlocalised scaffolds, which are RETAINED."""
        return len(self._regions[UNPLACED])

    @property
    def non_primary_base_fraction(self) -> float:
        """ALT/PATCH plus unplaced, as a fraction of all bases.

        Reported as a summary of everything off the primary chromosomes.
        :attr:`alt_patch_base_fraction` is the one that describes the exclusion.
        """
        total = self.total_bases
        if not total:
            return 0.0
        return (self.bases[ALT_PATCH] + self.bases[UNPLACED]) / total

    @property
    def alt_patch_base_fraction(self) -> float:
        total = self.total_bases
        return (self.bases[ALT_PATCH] / total) if total else 0.0

    def searched_classes(self) -> Tuple[int, ...]:
        """Which record classes the screen will actually see.

        UNPLACED is on BOTH sides. An unplaced scaffold is unique sequence that
        could not be placed on a chromosome; it duplicates nothing, so removing
        it cannot suppress a spurious self-match and can only hide a real
        off-target (module docstring, axis 1b).
        """
        if self.excluded_non_primary:
            return (PRIMARY, UNPLACED, UNKNOWN)
        return (PRIMARY, ALT_PATCH, UNPLACED, UNKNOWN)

    @property
    def records_searched(self) -> int:
        return sum(self.records[c] for c in self.searched_classes())

    @property
    def bases_searched(self) -> int:
        return sum(self.bases[c] for c in self.searched_classes())

    def biotype_bases(self) -> Dict[str, int]:
        """``{gene_biotype: bases}`` over the records the screen searches."""
        out: Dict[str, int] = {}
        for cls in self.searched_classes():
            for bio, (_n, b) in self._biotypes[cls].items():
                out[bio] = out.get(bio, 0) + b
        return out

    def biotype_records(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for cls in self.searched_classes():
            for bio, (n, _b) in self._biotypes[cls].items():
                out[bio] = out.get(bio, 0) + n
        return out

    def pseudogene_bases(self) -> int:
        return sum(b for bio, b in self.biotype_bases().items()
                   if is_pseudogene(bio))

    # -- output ------------------------------------------------------------

    def to_manifest(self) -> Dict[str, object]:
        """The manifest's record of this reference's composition.

        Deliberately absolute counts AND fractions: a fraction alone cannot be
        recomputed against a different reference, and a count alone makes the
        reader do arithmetic to see that 7.6% of GRCh38 cdna.all is copies of
        itself.
        """
        searched_bases = self.bases_searched
        pseudo = self.pseudogene_bases()
        out: Dict[str, object] = {
            "records": self.total_records,
            "bases": self.total_bases,
            "primary_records": self.records[PRIMARY],
            "primary_bases": self.bases[PRIMARY],
            # PER EXCLUSION REASON. One class is dropped and one is kept, and
            # a reader must be able to see which is which without knowing the
            # release: a single "non_primary" number cannot distinguish "19,713
            # duplicate copies removed" from "69 unique scaffolds removed".
            "alt_patch_records": self.records[ALT_PATCH],
            "alt_patch_bases": self.bases[ALT_PATCH],
            "alt_patch_base_fraction": round(self.alt_patch_base_fraction, 6),
            "alt_patch_seq_regions": self.alt_patch_seq_regions,
            "unplaced_records": self.records[UNPLACED],
            "unplaced_bases": self.bases[UNPLACED],
            "unplaced_seq_regions": self.unplaced_seq_regions,
            # ALT/PATCH + unplaced: everything off the primary chromosomes,
            # reported as one quantity beside the per-reason fields above.
            "non_primary_records": (self.records[ALT_PATCH]
                                    + self.records[UNPLACED]),
            "non_primary_bases": (self.bases[ALT_PATCH]
                                  + self.bases[UNPLACED]),
            "non_primary_base_fraction": round(
                self.non_primary_base_fraction, 6),
            "non_primary_seq_regions": self.non_primary_seq_regions,
            # Which classes the exclusion was allowed to remove, by name, so
            # the record does not depend on knowing this release's rule.
            "excluded_classes": (
                [ASSEMBLY_LABELS[c] for c in EXCLUDABLE]
                if self.excluded_non_primary else []),
            "retained_classes": [ASSEMBLY_LABELS[c]
                                 for c in self.searched_classes()],
            # Records whose headers are not Ensembl-style. For a viral panel or
            # a user FASTA this is everything, and that is the normal case, not
            # a defect.
            "unclassified_records": self.records[UNKNOWN],
            "unclassified_bases": self.bases[UNKNOWN],
            "ensembl_style_headers": self.ensembl_style,
            # None when nothing decided (a manifest-time scan of a reference the
            # k-mer loader never opened, e.g. a BLAST-engine run).
            "non_primary_excluded": self.excluded_non_primary,
            "records_searched": self.records_searched,
            "bases_searched": searched_bases,
            # Pseudogenes are REPORTED, never excluded — a transcribed
            # pseudogene missing from the reference is an off-target the screen
            # cannot see, i.e. a bad probe that passes. See the module
            # docstring, axis 2.
            "pseudogene_bases": pseudo,
            "pseudogene_base_fraction_of_searched": (
                round(pseudo / searched_bases, 6) if searched_bases else 0.0),
            "gene_biotype_bases": dict(sorted(self.biotype_bases().items())),
            "gene_biotype_records": dict(sorted(self.biotype_records().items())),
        }
        if self.source_sha256:
            # Absent rather than null when unknown — the same convention
            # tool_versions uses, and the opposite of a fabricated value.
            out["sha256"] = self.source_sha256
        return out


def should_exclude(comp: ReferenceComposition, requested: bool) -> bool:
    """Whether non-primary records may actually be dropped from this reference.

    The second term is the whole safety argument, and it is deliberately
    "has PRIMARY records" rather than "has Ensembl-style headers". A reference
    with no parseable location fields has every record in :data:`UNKNOWN`, so
    "drop the non-primary ones" would either do nothing or — read the other way
    — delete the file. And a reference that is ENTIRELY ALT contigs (someone
    pulled out the MHC haplotypes on purpose) would be emptied by a rule that
    only checked for parseability. Requiring proof that the file HAS a primary
    assembly to fall back on is what makes the default safe to ship ON.
    """
    return (bool(requested)
            and comp.records[PRIMARY] > 0
            and comp.records[ALT_PATCH] > 0)


def warn_if_no_primary(comp: ReferenceComposition) -> Optional[str]:
    """Shout when a reference parses as Ensembl but yields NO primary records.

    THE GUARD FOR AN UNRECOGNISED NAMING SCHEME. If a region type such as
    rat's ``primary_assembly`` were not recognised, every one of the 46,856
    records in ``Rattus_norvegicus.mRatBN7.2.cdna.all.fa`` would classify
    non-primary and none primary. Nothing would fail: :func:`should_exclude`
    requires PRIMARY > 0, so it would correctly drop nothing, and the only
    symptom would be a warning telling the user the opposite of the truth. A
    whole reference silently changing class is the failure mode a new Ensembl
    release or species produces, and it must be loud rather than inferred
    from a tally nobody reads.

    Returns the message, or None when there is nothing to say.
    """
    if not comp.ensembl_style or comp.records[PRIMARY] > 0:
        return None
    msg = (
        f"{os.path.basename(comp.path)}: {comp.total_records:,} records carry "
        f"parseable Ensembl location fields but NONE classified as primary "
        f"assembly. That is almost always this build failing to recognise the "
        f"file's region type or sequence naming, not a reference with no "
        f"chromosomes. Nothing has been excluded (the exclusion requires "
        f"primary records), but any primary/non-primary figure from this run "
        f"is unreliable. Region types seen: "
        f"{', '.join(sorted(comp.region_types_seen)) or 'none'}."
    )
    logger.warning("%s", msg)
    return msg


def warn_if_searching_alt(comp: ReferenceComposition, *,
                          disable_flag: str = "--no-exclude-alt-scaffolds",
                          ) -> Optional[str]:
    """One WARNING per reference when ALT/PATCH records are being searched.

    ONE LINE PER RUN, and that is a hard requirement rather than a style
    preference: a whole-transcriptome screen produces tens of millions of hits,
    and per-record logging of any kind is itself an outage.

    Returns the message (for tests and callers that want to surface it) or None.
    """
    if comp.excluded_non_primary or not comp.ensembl_style:
        return None
    if comp.alt_patch_base_fraction < ALT_PATCH_WARN_BASE_FRACTION:
        return None
    msg = (
        f"{os.path.basename(comp.path)}: {comp.records[ALT_PATCH]:,} ALT/PATCH "
        f"records ({comp.alt_patch_base_fraction:.1%} of bases, "
        f"{comp.alt_patch_seq_regions:,} scaffolds) are being SEARCHED as "
        f"off-target subjects. They are near-verbatim COPIES of primary regions "
        f"— on GRCh38 cdna.all they are the source of spurious perfect-50bp "
        f"off-target calls, i.e. good probes rejected for matching a second copy "
        f"of their own gene. Remove {disable_flag} to exclude them."
    )
    logger.warning("%s", msg)
    return msg


# ---------------------------------------------------------------------------
# Registry: composition measured at load time, read back at manifest time
# ---------------------------------------------------------------------------
#
# The reference is parsed by the screen during ``prepare()``; the manifest is
# assembled by the pipeline several stages later. Rather than thread the object
# through every layer between them (or re-read a 450 MB file), the loader
# publishes what it measured and the manifest asks for it by path. Keyed on the
# absolute path, which is what both ends have.
_MEASURED: Dict[str, ReferenceComposition] = {}


def register(comp: ReferenceComposition) -> None:
    _MEASURED[os.path.abspath(comp.path)] = comp


def measured(path: str) -> Optional[ReferenceComposition]:
    """What the loader measured for ``path``, or None if it never loaded it."""
    return _MEASURED.get(os.path.abspath(path))


def clear_registry() -> None:
    """Drop measured compositions (tests; a GUI session running many designs)."""
    _MEASURED.clear()


def scan_fasta(path: str) -> ReferenceComposition:
    """Measure composition from headers alone, without loading the sequence.

    The fallback for the manifest when the k-mer loader never ran — a
    BLAST-engine run, where the search is done by ``blastn`` against a database
    and no :class:`~flexprobe.kmer_index.Reference` is ever built. It still
    reads the whole file (base counts come from the sequence lines' lengths)
    but allocates nothing per base, so it costs I/O and not memory.

    Never raises: an unreadable reference yields an empty tally, because a
    provenance record must not be the thing that fails a finished run.
    """
    comp = ReferenceComposition(path)
    assembly = UNKNOWN
    biotype = BIOTYPE_UNKNOWN
    region = ""
    region_type = ""
    bases = 0
    started = False
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    if started:
                        comp.add(assembly, biotype, bases, region, region_type)
                    started = True
                    header = line[1:].strip()
                    assembly, biotype = classify_header(header)
                    if assembly == PRIMARY:
                        region_type, region = "", ""
                    else:
                        region_type, region = location_of(header)
                    bases = 0
                else:
                    bases += len(line.strip())
        if started:
            comp.add(assembly, biotype, bases, region, region_type)
    except OSError:
        logger.warning("Could not scan reference composition for %s", path,
                       exc_info=True)
    return comp
