"""Fixed pipeline scenarios + a canonical output digest.

WHY THIS EXISTS. ``pipeline.run_pipeline`` is being broken up into per-stage
functions. Every step of that work is a pure refactor: the outputs must not
move by a single byte. Unit tests assert *properties* ("the manifest has
tool_versions", "min_strain_mismatches is the minimum") and a refactor can
satisfy every one of them while silently reordering a selection or dropping a
column. A whole-output digest cannot be satisfied by accident.

The scenarios are chosen to reach every branch the refactor touches:

===========================  ======================================================
scenario                     what it is the only cover for
===========================  ======================================================
``single_gene_v2``           the plain path: ingest -> tile -> filter -> rank -> export
``flex_v1_order_csvs``       chemistry/pair + the order-CSV writers (discovery off,
                             v1 multiplex, so barcodes actually expand)
``multi_strain_conserved``   the 327-line ``config.multi_strain`` tiling branch,
                             conserved strategy + consensus.fa
``strain_unique``            the strain-unique consensus strategy, screened (the
                             strain rule is only observable via specificity)
``isoform_constitutive``     constitutive/isoform collapse, screened -- the only
                             cover for ``consensus_strain_rule = PAN_STRAIN``
``kmer_screened``            the off-target screen actually running: stage 4/5,
                             the junction-duplex gate, and every screen-derived
                             manifest field
``isoform_gene_aware``       ``gene_aware`` mode (group_by_gene ON, multi_strain
                             OFF) -- eight branches including the isoform-coverage
                             ranking weight, the min_isoform_coverage gate and
                             three manifest fields. No other scenario sets it,
                             because gene_aware is defined as
                             ``group_by_gene and not multi_strain``.
``multi_strain_capped``      stage 3b, the BATCH pre-screen cap. Reachable only
                             when multi-strain and screening are BOTH on (the
                             normal path streams the cap instead), so no other
                             scenario touches it -- verified by deleting the
                             post-cap ``_reconcile_degenerate_families`` call
                             and watching the other six stay green.
===========================  ======================================================

Kept in its own module (not the test file) so the same scenarios can be driven
from a throwaway script when bisecting a digest change by hand.
"""
import hashlib
import json
import os
from typing import Any, Dict, List, Tuple

from flexprobe.models import (
    ExportSettings,
    FilterThresholds,
    ModalityDefinition,
    RunConfig,
    ScreenResources,
)

# ---------------------------------------------------------------------------
# Sequences. Literal and self-contained: a golden baseline whose *input* can
# drift is worthless, so these are never imported from a fixture that another
# test could edit.
# ---------------------------------------------------------------------------

_BASE_200 = (
    "AGCTGTAGCTAGCTAGCATGCATGCATGCAGTCAGTCGATGCATCGTAGCATCGATCGTACTAGCTAGCAT"
    * 4
)[:200]


def _mutate(seq: str, positions: Tuple[int, ...]) -> str:
    """Flip the given positions to a different base, deterministically."""
    b = list(seq)
    for i in positions:
        b[i] = "C" if b[i] != "C" else "A"
    return "".join(b)


_STRAIN_A = _BASE_200
_STRAIN_B = _mutate(_BASE_200, (60, 130))
_STRAIN_C = _mutate(_BASE_200, (60, 90, 175))

# A 120-bp target and an unrelated decoy, same shape as the engine fixtures in
# test_junction_duplex.py, so the k-mer screen has something to reject against.
_SCREEN_TARGET = (
    "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGTTACCGATGCATGC"
    "ATGCATGCATGCATGCATGCAACGTACGTACGTACGTACGTTGCATGCATG"
    "CATGCATGCATG"
)
_SCREEN_DECOY = "TTTTTTTTAAAAAAAACCCCCCCCGGGGGGGG" * 4


def _write_fasta(path: str, records: List[Tuple[str, str]],
                 descriptions: Dict[str, str] = None) -> str:
    """``descriptions`` adds text after the id, e.g. Ensembl's ``gene:ENSG...``.

    Needed because ``_build_gene_groups`` reads the gene key out of the FASTA
    *description*, so a headerless record is always its own group and the whole
    gene-aware branch stays unreachable.
    """
    descriptions = descriptions or {}
    with open(path, "w") as fh:
        for sid, seq in records:
            desc = descriptions.get(sid, "")
            header = f"{sid} {desc}".rstrip()
            fh.write(f">{header}\n{seq}\n")
    return path


def _relaxed() -> FilterThresholds:
    """Permissive thresholds so the synthetic sequences yield probes at all.

    A scenario that outputs zero probes digests identically no matter what the
    refactor breaks, which is the one failure mode a golden test must not have.
    """
    return FilterThresholds(
        gc_min=0.30,
        gc_max=0.80,
        max_n_per_half=2,
        max_homopolymer_run=6,
        require_tn_junction=False,
    )


# ---------------------------------------------------------------------------
# Scenario construction
# ---------------------------------------------------------------------------

def _cfg(fasta: str, out: str, **kw: Any) -> RunConfig:
    base: Dict[str, Any] = dict(
        fasta_path=fasta,
        output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        filters=_relaxed(),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        tile_length=50,
        tile_step=5,
    )
    base.update(kw)
    return RunConfig(**base)


def _single_gene_v2(work: str) -> RunConfig:
    fa = _write_fasta(os.path.join(work, "single.fa"), [("G", _BASE_200)])
    return _cfg(fa, os.path.join(work, "out"))


def _flex_v1_order_csvs(work: str) -> RunConfig:
    fa = _write_fasta(os.path.join(work, "v1.fa"),
                      [("G1", _BASE_200), ("G2", _STRAIN_C)])
    return _cfg(
        fa, os.path.join(work, "out"),
        modality=ModalityDefinition(name="flex", version="v1"),
        discovery_only=False,
        auto_select_per_gene=3,
        probe_barcodes=["BC001", "BC002"],
        multiplex=True,
    )


def _multi_strain_conserved(work: str) -> RunConfig:
    # Basename is the fragment label the consensus path uses, so name it.
    fa = _write_fasta(os.path.join(work, "OROV_L.fa"),
                      [("sA", _STRAIN_A), ("sB", _STRAIN_B), ("sC", _STRAIN_C)])
    return _cfg(fa, os.path.join(work, "out"),
                multi_strain=True, consensus_strategy="conserved",
                # Non-zero so the record ROTATION runs. build_consensus treats
                # the first record as both the coordinate system and the strain
                # ambiguous bases resolve to, so rotating changes the consensus
                # itself -- and at the default 0 those three lines never
                # execute, which line coverage of this file confirmed.
                # multi_strain_capped keeps the unrotated conserved path.
                consensus_reference=1)


def _strain_unique(work: str) -> RunConfig:
    """The strain-unique consensus strategy, WITH the screen on.

    Screening is not incidental here: ``consensus_strain_rule`` is only ever
    READ by the specificity stage, so a skip-screen run leaves the strategy's
    choice of StrainRule completely unobserved -- and that assignment sits
    inside the multi-strain branch being extracted.
    """
    fa = _write_fasta(os.path.join(work, "OROV_M.fa"),
                      [("sA", _STRAIN_A), ("sB", _STRAIN_B), ("sC", _STRAIN_C)])
    db = _write_fasta(os.path.join(work, "off.fa"), [("UNRELATED", _SCREEN_DECOY)])
    return _cfg(fa, os.path.join(work, "out"),
                multi_strain=True, consensus_strategy="strain_unique",
                strain_unique_target=0, strain_unique_min_mismatches=1,
                screen=ScreenResources(skip_screen=False, db_paths=[db]))


def _isoform_constitutive(work: str) -> RunConfig:
    """Constitutive (isoform-collapse) strategy, WITH the screen on.

    Same reason as _strain_unique: this is the branch that sets
    ``consensus_strain_rule = StrainRule.PAN_STRAIN``, and that assignment has
    no observable effect unless specificity actually runs. Proved by replacing
    the line with ``pass`` and watching an unscreened version stay green.
    """
    fa = _write_fasta(os.path.join(work, "GENE.fa"),
                      [("GENE-201", _STRAIN_A), ("GENE-202", _STRAIN_B),
                       ("GENE-203", _STRAIN_C)])
    db = _write_fasta(os.path.join(work, "off.fa"), [("UNRELATED", _SCREEN_DECOY)])
    return _cfg(fa, os.path.join(work, "out"),
                multi_strain=True, consensus_strategy="constitutive",
                group_by_gene=True,
                # THREE isoforms and a quorum of two, so the kept windows
                # include some that bind only 2 of 3. With every window binding
                # every isoform (the 2-isoform, quorum-of-all version) EXACT_SET
                # and PAN_STRAIN agree on every candidate and the rule is
                # unobservable no matter what the refactor does to it.
                isoform_min_count=2,
                screen=ScreenResources(skip_screen=False, db_paths=[db]))


def _kmer_screened(work: str) -> RunConfig:
    fa = _write_fasta(os.path.join(work, "target.fa"), [("GENE1", _SCREEN_TARGET)])
    db = _write_fasta(os.path.join(work, "off.fa"), [("UNRELATED", _SCREEN_DECOY)])
    return _cfg(
        fa, os.path.join(work, "out"),
        # Default thresholds here, not the relaxed set: the screened scenario is
        # the one that must exercise the real gates.
        filters=FilterThresholds(require_tn_junction=False),
        screen=ScreenResources(skip_screen=False, db_paths=[db]),
        tile_step=1,
    )


def _multi_strain_capped(work: str) -> RunConfig:
    """Multi-strain WITH screening, so the batch cap (stage 3b) actually runs.

    ``_cap_active`` requires ``will_blast``, and the streamed cap is skipped on
    the multi-strain path — so this combination, and only this combination,
    reaches ``select_top_per_gene`` + the post-cap family reconcile in the
    batch fallback. ``screen_top_per_gene`` is set low enough that the cap
    genuinely discards candidates rather than being a no-op.
    """
    fa = _write_fasta(os.path.join(work, "OROV_S.fa"),
                      [("sA", _STRAIN_A), ("sB", _STRAIN_B), ("sC", _STRAIN_C)])
    db = _write_fasta(os.path.join(work, "off.fa"), [("UNRELATED", _SCREEN_DECOY)])
    return _cfg(fa, os.path.join(work, "out"),
                multi_strain=True, consensus_strategy="conserved",
                screen=ScreenResources(skip_screen=False, db_paths=[db]),
                screen_top_per_gene=8)


def _isoform_gene_aware(work: str) -> RunConfig:
    """``gene_aware`` mode: group_by_gene ON and multi_strain OFF.

    ``pipeline.py`` computes ``gene_aware = config.group_by_gene and not
    config.multi_strain``, so ``isoform_constitutive`` -- which sets both --
    leaves it FALSE. That left eight ``gene_aware`` branches unreached,
    including the isoform-coverage ranking weight, the ``min_isoform_coverage``
    hard gate and three manifest fields. Proved by deleting the
    ``isoform_coverage=DEFAULT_ISOFORM_COVERAGE_WEIGHT`` line in stage 6 and
    watching every other scenario stay green.

    ``min_isoform_coverage`` is set below 1.0 but above 0 so the gate is both
    ACTIVE and non-trivial: the two isoforms differ, so some windows bind one
    and some bind both.
    """
    fa = _write_fasta(
        os.path.join(work, "isoforms.fa"),
        [("TX-201", _STRAIN_A), ("TX-202", _STRAIN_B), ("OTHER-201", _STRAIN_C)],
        descriptions={
            "TX-201": "gene:GENEA",
            "TX-202": "gene:GENEA",
            "OTHER-201": "gene:GENEB",
        },
    )
    return _cfg(fa, os.path.join(work, "out"),
                group_by_gene=True, min_isoform_coverage=0.5)


#: name -> factory(work_dir) -> RunConfig. Ordered; the baseline file is keyed
#: by name, so adding a scenario never invalidates the existing entries.
SCENARIOS = {
    "single_gene_v2": _single_gene_v2,
    "flex_v1_order_csvs": _flex_v1_order_csvs,
    "multi_strain_conserved": _multi_strain_conserved,
    "strain_unique": _strain_unique,
    "isoform_constitutive": _isoform_constitutive,
    "kmer_screened": _kmer_screened,
    "multi_strain_capped": _multi_strain_capped,
    "isoform_gene_aware": _isoform_gene_aware,
}


# ---------------------------------------------------------------------------
# Digesting
# ---------------------------------------------------------------------------

#: Files whose content is wall-clock timing and nothing else. Excluded whole:
#: they change on every run by construction, and _StageTimer has its own tests
#: (test_run_timings.py).
_TIMING_FILES = {"run_timings.json", "run_timeline.jsonl"}

#: Manifest/run_config keys that legitimately differ between two runs of the
#: same code. Dropped by name at any depth.
#:   timestamp/elapsed_seconds — wall clock
#:   code_version              — git sha; changes on every commit of this work
#:   tool_versions             — python/primer3 build, differs per machine
_VOLATILE_KEYS = {"timestamp", "elapsed_seconds", "code_version", "tool_versions"}

_JSON_FILES = {"manifest.json", "run_config.json", "filter_report.json"}


def _scrub(obj: Any, subs: List[Tuple[str, str]]) -> Any:
    """Drop volatile keys and normalise absolute paths, recursively.

    Path normalisation is done on *values and keys alike* because
    ``db_sha256s`` and ``reference_composition`` are keyed BY path — a
    tmp-dir-keyed dict would otherwise make every run digest differently.
    """
    if isinstance(obj, dict):
        return {
            _scrub(k, subs): _scrub(v, subs)
            for k, v in obj.items()
            if k not in _VOLATILE_KEYS
        }
    if isinstance(obj, list):
        return [_scrub(v, subs) for v in obj]
    if isinstance(obj, str):
        for needle, token in subs:
            obj = obj.replace(needle, token)
        return obj
    return obj


def _canonical_json(text: str, subs: List[Tuple[str, str]]) -> bytes:
    return json.dumps(
        _scrub(json.loads(text), subs),
        sort_keys=True, separators=(",", ":"), default=str,
    ).encode()


def digest_run(output_dir: str, subs: List[Tuple[str, str]]) -> Dict[str, str]:
    """Return ``{relative_path: sha256}`` for every stable output file.

    ``subs`` is a list of (absolute_prefix, token) replacements applied to text
    content, so a digest taken under ``/var/folders/...`` matches one taken
    under ``/tmp/...``.
    """
    out: Dict[str, str] = {}
    for root, dirs, files in os.walk(output_dir):
        dirs.sort()
        for name in sorted(files):
            if name in _TIMING_FILES or name.endswith((".log", ".pyc")):
                continue
            path = os.path.join(root, name)
            rel = os.path.relpath(path, output_dir)
            if name in _JSON_FILES:
                with open(path) as fh:
                    payload = _canonical_json(fh.read(), subs)
            else:
                with open(path, "rb") as fh:
                    raw = fh.read()
                try:
                    text = raw.decode()
                except UnicodeDecodeError:
                    payload = raw
                else:
                    for needle, token in subs:
                        text = text.replace(needle, token)
                    payload = text.encode()
            out[rel] = hashlib.sha256(payload).hexdigest()
    return out


def run_scenario(name: str, work_dir: str) -> Dict[str, str]:
    """Execute one scenario under ``work_dir`` and return its output digest."""
    from flexprobe.pipeline import run_pipeline

    os.makedirs(work_dir, exist_ok=True)
    cfg = SCENARIOS[name](work_dir)
    run_pipeline(cfg)
    # Both the work dir and its realpath: macOS hands out /var/folders/... which
    # is a symlink to /private/var/folders/..., and the two forms appear in
    # different fields of the same manifest.
    subs = [(os.path.realpath(work_dir), "<WORK>"), (work_dir, "<WORK>")]
    return digest_run(cfg.output_dir, subs)
