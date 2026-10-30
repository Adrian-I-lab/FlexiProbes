#!/usr/bin/env python3
"""flexprobe headless CLI — runs the real probe-design pipeline."""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any, Dict, Optional

from . import __version__
from .constants import (
    MIN_PROBE_GAP_BP, MIN_PROBE_GAP_FLOOR_BP,
    OFF_TARGET_JUNCTION_DUPLEX_BP, OFF_TARGET_JUNCTION_DUPLEX_MIN_BP,
    OFF_TARGET_SAME_SUBJECT_DEFAULT,
    PREPARED_MARKER,
    PREFERRED_PROBE_GAP_BP, PROBE_GAP_RELAX_MIN_PROBES, RELAXED_PROBE_GAP_BP,
    REFERENCE_STRAND_DEFAULT, REFERENCE_STRAND_PLUS, REFERENCE_STRANDS,
    SCREEN_ENGINE_DEFAULT,
    TM_MIN_DEFAULT, TN_MODE_AUTO, TN_MODE_SOFT, TN_MODES,
)
from .errors import BlastError, FlexprobeError
from .reference_qc import EXCLUDE_ALT_SCAFFOLDS_DEFAULT
from .screen import ENGINE_BLAST, ENGINES

from .null_model import (
    DEFAULT_MAX_WINDOWS as NULL_MODEL_DEFAULT_MAX_WINDOWS,
    DEFAULT_SHUFFLES as NULL_MODEL_DEFAULT_SHUFFLES,
)
from .hpc_export import (
    DEFAULT_CHAIN_LINKS, DEFAULT_CORES, DEFAULT_MEM_GB, DEFAULT_WALLTIME,
    DEFAULT_SCREEN_MEM_GB,
    DEFAULT_PREPARE_CORES, DEFAULT_PREPARE_LINKS, DEFAULT_PREPARE_MEM_GB,
    DEFAULT_PREPARE_WALLTIME,
    PREPARE_SCRIPT_STEM, SCREEN_SCRIPT_STEM,
    SCHEDULERS,
)
from .auto_export import AutoExportResult, auto_export, spacing_tier_lines
from .logging_setup import setup_logging
from .models import (
    ScreenResources,
    ExportSettings,
    FilterThresholds,
    ModalityDefinition,
    RunConfig,
    SpecificityThresholds,
)
from .pipeline import ProgressEvent, run_pipeline
from .summary import (
    write_failed_summary, write_ok_summary, write_prepared_summary,
)


MODALITY_PRESETS = {
    "flex_v1":  ModalityDefinition(name="flex",   version="v1", spatial=False),
    "flex_v2":  ModalityDefinition(name="flex",   version="v2", spatial=False),
    "spatial":  ModalityDefinition(name="visium", version="v1", spatial=True),
}


def _ensure_outdir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def _make_progress_printer():
    def _print(ev: ProgressEvent) -> None:
        # pct is None for messages that carry no progress (a pre-stage cost
        # warning, say). Same width as "[100.0%]" so the log stays columnar.
        head = "[  --  ]" if ev.pct is None else f"[{ev.pct:5.1f}%]"
        print(f"{head} {ev.stage:<12} {ev.message}", flush=True)
    return _print


def _build_run_config(args: argparse.Namespace) -> RunConfig:
    # --config carries a full RunConfig serialized by the GUI (HPC export). Load
    # it verbatim so the cluster run is identical to the one configured in the
    # app, then let the path-bearing flags override the values baked at export
    # time — the bundle lives at a different path on the cluster.
    if getattr(args, "config", None):
        from .config_io import load_config
        cfg = load_config(args.config)
        if args.fasta:
            cfg.fasta_path = args.fasta
        if args.out:
            cfg.output_dir = args.out
        if args.off_target_ref:
            cfg.screen.db_paths = list(args.off_target_ref)
        if args.cached_hits:
            cfg.screen.use_cached_hits = True
            cfg.screen.cached_hits_path = args.cached_hits
        if args.skip_screen:
            cfg.screen.skip_screen = True
        # Turn-on-only, and that is the whole mechanism the split chain rests
        # on: ONE config.json is shared by both jobs, job A adds
        # --prepare-only on its command line and job B omits it. Baking the
        # flag into the config instead would leave job B preparing again and
        # never screening. See candidate_cache.FINGERPRINT_IGNORED for why the
        # two jobs still agree on the checkpoint digest.
        if getattr(args, "prepare_only", False):
            cfg.prepare_only = True
        # is not None (not truthiness): --no-exclude-alt-scaffolds must be able
        # to turn OFF an exclusion the baked config turned on, which a
        # turn-on-only override could not express.
        if getattr(args, "exclude_alt_scaffolds", None) is not None:
            cfg.screen.exclude_alt_scaffolds = args.exclude_alt_scaffolds
        if getattr(args, "reference_strand", None) is not None:
            cfg.screen.reference_strand = args.reference_strand
        # Turn-on only: an operator action for this invocation, never baked.
        if getattr(args, "rebuild_kmer_index", False):
            cfg.screen.rebuild_kmer_index = True
        # --workers defaults to None so an explicit value (including 1, to fit a
        # smaller allocation) overrides the baked config; unset keeps it.
        if args.workers is not None:
            cfg.workers = args.workers
        if getattr(args, "screen_workers", None) is not None:
            cfg.screen_workers = args.screen_workers
        if args.force:
            cfg.export.overwrite = True
        if args.no_resume:
            cfg.resume_screen = False
        # is not None (not truthiness) so an explicit `--screen-top-per-gene 0`
        # can reset a baked non-zero cap back to "screen all".
        if args.screen_top_per_gene is not None:
            cfg.screen_top_per_gene = args.screen_top_per_gene
        # Turn-on-only overrides, like --force/--skip-screen above: a baked
        # config that already emits the artefacts is never silently turned off
        # by an invocation that omits the flag.
        if args.emit_prescreen_candidates:
            cfg.emit_prescreen_candidates = True
        if args.emit_prescreen_precap:
            cfg.emit_prescreen_candidates = True
            cfg.emit_prescreen_precap_rows = True
        if args.keep_candidate_cache:
            cfg.keep_candidate_cache = True
        return cfg

    modality = MODALITY_PRESETS[args.modality]
    filters = FilterThresholds(
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        enforce_gc_range=args.enforce_gc_range,
        # A negative --tm-min is the documented way to disable the gate;
        # argparse has no "off" sentinel for a float option.
        tm_min=(None if args.tm_min is not None and args.tm_min < 0
                else args.tm_min),
        tm_max=args.tm_max,
        max_n_per_half=args.max_n_per_half,
        max_homopolymer_run=args.max_homopolymer_run,
        # None == neither flag passed; fall through to the dataclass default
        # (True) rather than pinning it here, so there is one source of truth.
        **({} if args.require_tn_junction is None
           else {"require_tn_junction": args.require_tn_junction}),
        tn_junction_mode=args.tn_junction_mode,
        **({"junction_window": args.junction_window}
           if args.junction_window is not None else {}),
    )
    specificity = SpecificityThresholds(
        off_target_full_half=not args.off_target_local_alignment,
        off_target_require_both_halves=not args.off_target_either_half,
        off_target_same_subject=args.off_target_same_subject,
        off_target_junction_duplex_bp=args.off_target_junction_duplex,
    )
    screen = ScreenResources(
        **({"query_chunk": args.query_chunk}
           if args.query_chunk is not None else {}),
        **({"max_target_seqs": args.max_target_seqs}
           if args.max_target_seqs is not None else {}),
        db_paths=list(args.off_target_ref) if args.off_target_ref else [],
        engine=args.screen_engine,
        exclude_alt_scaffolds=(
            EXCLUDE_ALT_SCAFFOLDS_DEFAULT
            if getattr(args, "exclude_alt_scaffolds", None) is None
            else args.exclude_alt_scaffolds),
        reference_strand=(
            REFERENCE_STRAND_DEFAULT
            if getattr(args, "reference_strand", None) is None
            else args.reference_strand),
        rebuild_kmer_index=bool(getattr(args, "rebuild_kmer_index", False)),
        skip_screen=args.skip_screen,
        use_cached_hits=bool(args.cached_hits),
        cached_hits_path=args.cached_hits,
    )
    export = ExportSettings(
        write_idt_opools=not args.no_idt,
        write_cellranger_probe_set=not args.no_cellranger,
        write_results_table=True,
        write_manifest=True,
        overwrite=args.force,
    )
    return RunConfig(
        fasta_path=args.fasta or "",
        output_dir=args.out or "flexprobe_out",
        modality=modality,
        filters=filters,
        specificity=specificity,
        screen=screen,
        export=export,
        tile_length=args.tile_size,
        tile_step=args.tile_step,
        multiplex=not args.singleplex,
        # --probe-barcode is repeatable so users can specify multiple Flex v1
        # Probe Barcodes for multi-sample multiplex. Falls back to BC001 when
        # nothing was passed.
        probe_barcodes=list(args.probe_barcode) or ["BC001"],
        primer3_scoring=not args.no_primer3,
        **({"min_probe_gap_bp": args.min_probe_gap}
           if args.min_probe_gap is not None else {}),
        **({"relaxed_probe_gap_bp": args.relaxed_probe_gap}
           if getattr(args, "relaxed_probe_gap", None) is not None else {}),
        **({"floor_probe_gap_bp": args.floor_probe_gap}
           if getattr(args, "floor_probe_gap", None) is not None else {}),
        split_header_on_pipe=getattr(args, "split_header_on_pipe", False),
        workers=args.workers if args.workers is not None else 1,
        screen_workers=(args.screen_workers
                        if getattr(args, "screen_workers", None) is not None
                        else 1),
        extra_on_target_ids=list(args.extra_on_target),
        **({"group_by_gene": args.group_by_gene}
           if args.group_by_gene is not None else {}),
        min_isoform_coverage=args.min_isoform_coverage,
        multi_strain=args.multi_strain,
        consensus_strategy=args.strategy,
        isoform_min_count=args.isoform_min,
        **({"conserved_max_degenerate_per_half":
            args.conserved_max_degenerate_per_half}
           if args.conserved_max_degenerate_per_half is not None else {}),
        strain_unique_target=args.strain_unique_target,
        strain_unique_min_mismatches=args.strain_unique_min_mm,
        strain_unique_per_half=not args.strain_unique_whole_window,
        strain_unique_junction_duplex_bp=args.strain_unique_junction_duplex,
        discovery_only=args.discovery_only or args.auto_select_per_gene > 0,
        auto_select_per_gene=args.auto_select_per_gene,
        resume_screen=not args.no_resume,
        prepare_only=getattr(args, "prepare_only", False),
        screen_cache_mode=getattr(args, "screen_cache_mode", "verdicts"),
        screen_top_per_gene=args.screen_top_per_gene or 0,
        # --emit-prescreen-precap implies the base artefacts: the pre-cap dump is
        # an addition to them, not a replacement (RunConfig.validate enforces it).
        emit_prescreen_candidates=(args.emit_prescreen_candidates
                                  or args.emit_prescreen_precap),
        emit_prescreen_precap_rows=args.emit_prescreen_precap,
        keep_candidate_cache=args.keep_candidate_cache,
        # --null-model is shorthand for the recommended 20 shuffles; an explicit
        # --null-model-shuffles N wins so the count can be tuned either way.
        null_model_shuffles=(
            args.null_model_shuffles
            if args.null_model_shuffles is not None
            else (NULL_MODEL_DEFAULT_SHUFFLES if args.null_model else 0)
        ),
        null_model_max_windows=args.null_model_max_windows,
        null_model_seed=args.null_model_seed,
    )


# --------------- subcommand handlers ---------------

def cmd_design(args: argparse.Namespace) -> int:
    cfg = _build_run_config(args)
    if not cfg.fasta_path:
        print(
            "[flexiprobes] ERROR: no input FASTA — pass --fasta "
            "(or a --config whose fasta_path is set).",
            file=sys.stderr,
        )
        return 2
    outdir = _ensure_outdir(cfg.output_dir)

    cached_hits = None
    if cfg.screen.use_cached_hits and cfg.screen.cached_hits_path:
        try:
            from .blast import load_cached_hits
            cached_hits = load_cached_hits(cfg.screen.cached_hits_path)
        except Exception as exc:
            print(f"[flexiprobes] WARNING: could not load cached hits: {exc}", file=sys.stderr)

    mod = f"{cfg.modality.name}/{cfg.modality.version}"
    print(f"[flexiprobes] design: fasta={cfg.fasta_path} out={outdir} modality={mod}")
    try:
        result = run_pipeline(
            cfg,
            progress_callback=None if args.quiet else _make_progress_printer(),
            screen_hits_by_probe=cached_hits,
        )
    except FlexprobeError as exc:
        print(f"[flexiprobes] ERROR: {exc}", file=sys.stderr)
        # Write a summary anyway. A run whose funnel emptied is a RESULT (the
        # tiles existed and the filters rejected all of them), and without a
        # design_summary.json a consumer that keys off that file could not
        # tell it apart from a run that never started, and could record real
        # candidate sets as "0 tiles". The counts below come from
        # filter_report.json, which the pipeline does write before it aborts.
        write_failed_summary(outdir, exc)
        return 2

    # ---- --prepare-only: report and stop, BEFORE any of the finished-run tail
    # Placed here, above the auto-export and above write_ok_summary, because
    # every step below assumes a completed funnel. In particular
    # write_ok_summary would stamp status "no_probes" on a run that designed
    # nothing because it was TOLD to stop, and the pairs_output==0 branch at
    # the end would then return exit 1 — turning a successful prepare job into
    # a failure the PBS chain would classify as a deterministic crash and
    # abort the whole chain on.
    if cfg.prepare_only:
        stats = result["stats"]
        summary_path = write_prepared_summary(outdir, stats)
        print(f"[flexiprobes] PREPARE-ONLY: {stats['prepared_candidates']} "
              f"candidates checkpointed under {outdir}")
        print(f"[flexiprobes] marker: "
              f"{os.path.join(outdir, PREPARED_MARKER)}")
        print(f"[flexiprobes] summary: {summary_path} (status=prepared)")
        print("[flexiprobes] This run is INCOMPLETE: no off-target screen, no "
              "ranking, no results.csv. Re-run the same config without "
              "--prepare-only to finish it.")
        return 0

    # Whole-transcriptome convenience: in discovery mode with an auto-select
    # count, pick the top-N spaced probes per gene and write the order CSVs
    # without a manual review step — the same post-pipeline path the GUI runs in
    # RunManager._auto_export. Non-fatal: a failure leaves results.csv for review.
    if cfg.discovery_only and cfg.auto_select_per_gene > 0:
        _auto_export(outdir, cfg.auto_select_per_gene,
                    min_gap=cfg.min_probe_gap_bp,
                    relaxed_gap=cfg.relaxed_probe_gap_bp,
                    floor_gap=cfg.floor_probe_gap_bp,
                    write_idt=cfg.export.write_idt_opools,
                    write_cellranger=cfg.export.write_cellranger_probe_set)

    summary_path = write_ok_summary(outdir, result["stats"])
    pairs_output = result["stats"]["pairs_output"]
    print(f"[flexiprobes] {pairs_output} pairs written under {outdir}")
    print(f"[flexiprobes] summary: {summary_path}")
    if pairs_output == 0:
        # The pipeline ran clean (nothing raised — every candidate had ALREADY
        # passed the input filters) but the funnel still emptied downstream,
        # typically because specificity screening rejected every remaining
        # candidate. That is a real, reportable result — design_summary.json
        # says so via status="no_probes" (write_ok_summary) — but a caller
        # that only checks the exit code must be able to tell it apart from a
        # run that actually designed something, so this is NOT exit 0.
        print(
            "[flexiprobes] ERROR: 0 probe pairs produced — every filter-"
            "passing candidate was rejected downstream (see "
            "specificity_status_counts in design_summary.json). Nothing to "
            "order.", file=sys.stderr,
        )
        return 1
    return 0


def _print_spacing_tiers(counts, min_gap: int, relaxed_gap: int,
                         floor_gap: int, stream=None) -> None:
    """Print the green/orange/red spacing banner.

    The lines themselves come from ``auto_export.spacing_tier_lines`` so that
    the GUI, which logs them rather than printing them, cannot end up wording
    the same three numbers differently.
    """
    out = stream if stream is not None else sys.stdout
    for line in spacing_tier_lines(counts, min_gap, relaxed_gap, floor_gap):
        print(f"[flexiprobes] {line}", file=out)


def _auto_export(outdir: str, n_per_gene: int, *,
                 min_gap: int = MIN_PROBE_GAP_BP,
                 relaxed_gap: int = RELAXED_PROBE_GAP_BP,
                 floor_gap: int = MIN_PROBE_GAP_FLOOR_BP,
                 write_idt: bool = True,
                 write_cellranger: bool = True) -> AutoExportResult:
    """Auto-pick top-N probes per gene, write order CSVs, and say what happened.

    The selection, the export and the wording of every message live in
    ``flexprobe.auto_export`` — shared with the GUI, which runs the identical
    operation after a whole-transcriptome run. All that is left here is the
    channel: results on stdout, anything the user has to act on to stderr.
    """
    result = auto_export(outdir, n_per_gene, min_gap=min_gap,
                         relaxed_gap=relaxed_gap, floor_gap=floor_gap,
                         write_idt=write_idt,
                         write_cellranger=write_cellranger)
    for line in result.report_lines():
        print(f"[flexiprobes] {line}")
    for line in result.warning_lines():
        print(f"[flexiprobes] {line}", file=sys.stderr)
    return result


def _blast_db_layout(prefix: str) -> Optional[str]:
    """``"single-volume"``, ``"multi-volume"`` or None when no database is there.

    A large makeblastdb output is split into ``<prefix>.00.*``, ``.01.*``, ...
    volumes tied together by a ``<prefix>.nal`` alias, and has no
    ``<prefix>.nhr`` at all, so checking only the single-volume sidecars
    would report such a database as missing.
    """
    if all(os.path.exists(prefix + ext) for ext in (".nhr", ".nin", ".nsq")):
        return "single-volume"
    if os.path.exists(prefix + ".nal"):
        return "multi-volume"
    return None


def cmd_validate_blast(args: argparse.Namespace) -> int:
    """Check what a BLAST-engine run needs, and say which part is missing.

    Tools: blastn and makeblastdb, and blastdbcmd, without which every
    off-target count falls back to the conservative estimate. With --db: the
    database files (single- or multi-volume) and its accession index, i.e.
    whether it was built with -parse_seqids, asked of the database itself
    exactly as the pipeline asks before it screens. Sidecar files alone are
    not enough, since none of them is the pipeline's hard requirement.
    """
    outdir = _ensure_outdir(args.out)
    report: Dict[str, Any] = {"status": "ok", "db": args.db}

    def fail(key: str, value: Any, message: str) -> None:
        report["status"] = "error"
        report[key] = value
        print(f"[flexiprobes] ERROR: {message}", file=sys.stderr)

    try:
        from .blast import check_version, discover_blast
        tools = discover_blast()
        report["tools"] = tools
        report["version"] = check_version(tools["blastn"])
    except BlastError as exc:
        fail("error", str(exc), str(exc))
        tools = {}
    if tools and not tools.get("blastdbcmd"):
        fail("missing_tools", ["blastdbcmd"],
             "blastdbcmd not found: BLAST-engine off-target counts would fall "
             "back to the conservative estimate, which over-rejects.")

    if args.db:
        layout = _blast_db_layout(args.db)
        report["db_layout"] = layout
        if layout is None:
            fail("missing_db_files", [".nhr", ".nin", ".nsq", "or .nal"],
                 f"no BLAST database at {args.db}: neither the single-volume "
                 f".nhr/.nin/.nsq sidecars nor a multi-volume .nal alias.")
        elif tools.get("blastdbcmd"):
            from .blast import SubjectWindows
            live = SubjectWindows([args.db]).verify()
            report["retrieval"] = bool(live)
            if not live:
                fail("retrieval", False,
                     f"{args.db} cannot serve blastdbcmd -entry: it was built "
                     f"without -parse_seqids, and a BLAST-engine run refuses "
                     f"it. Rebuild with: makeblastdb -in <fasta> -dbtype nucl "
                     f"-parse_seqids -out {args.db}")

    rpath = os.path.join(outdir, "blast_validation.json")
    with open(rpath, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"[flexiprobes] blast validation written to {rpath}")
    return 0 if report["status"] == "ok" else 1


def cmd_build_db(args: argparse.Namespace) -> int:
    outdir = _ensure_outdir(args.out)
    db_prefix = args.db_prefix or os.path.join(outdir, "flexprobe_db")
    try:
        from .blast import build_db, discover_blast
        tools = discover_blast()
        result = build_db(tools["makeblastdb"], args.input, db_prefix, title=args.title)
    except BlastError as exc:
        print(f"[flexiprobes] ERROR: {exc}", file=sys.stderr)
        return 2

    report = {"status": "ok", "input": args.input, "db_prefix": result}
    rpath = os.path.join(outdir, "build_db_result.json")
    with open(rpath, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"[flexiprobes] built DB at {result}")
    return 0


def cmd_export(args: argparse.Namespace) -> int:
    """Re-export a previous run by reading its results.csv and writing other formats."""
    src = args.results_csv
    if not os.path.exists(src):
        print(f"[flexiprobes] ERROR: {src} not found", file=sys.stderr)
        return 2
    target = os.path.join(_ensure_outdir(args.out), f"export.{args.format}")
    if os.path.exists(target) and not args.force:
        print(f"[flexiprobes] ERROR: {target} exists. Use --force to overwrite.", file=sys.stderr)
        return 1

    import csv as _csv
    if args.format in ("csv", "tsv"):
        sep = "," if args.format == "csv" else "\t"
        with open(src) as fin, open(target, "w", newline="") as fout:
            rd = _csv.reader(fin)
            wr = _csv.writer(fout, delimiter=sep)
            for row in rd:
                wr.writerow(row)
    elif args.format == "json":
        with open(src) as fin:
            rows = list(_csv.DictReader(fin))
        with open(target, "w") as fout:
            json.dump(rows, fout, indent=2)
    elif args.format == "bed":
        with open(src) as fin, open(target, "w") as fout:
            for drow in _csv.DictReader(fin):
                fout.write(
                    f"{drow.get('gene_id', '')}\t{drow.get('start', 0)}\t{drow.get('end', 0)}"
                    f"\t{drow.get('probe_id', '')}\t{drow.get('composite_score', 0)}\n"
                )
    else:
        print(f"[flexiprobes] ERROR: unsupported format {args.format}", file=sys.stderr)
        return 1
    print(f"[flexiprobes] exported {src} -> {target}")
    return 0


def _load_gene_isoforms(out_dir: str) -> dict:
    """Read ``manifest.json``'s gene -> isoform-ids map (empty if absent)."""
    path = os.path.join(out_dir, "manifest.json")
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as fh:
            man = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}
    params = (man.get("params") or {}) if isinstance(man, dict) else {}
    raw = params.get("gene_isoforms") or {}
    return {str(g): list(ts) for g, ts in raw.items()} if isinstance(raw, dict) else {}


def cmd_select_coverage(args: argparse.Namespace) -> int:
    """Pick up to N spaced probes per gene that maximise isoform coverage.

    Reads a finished run's ``results.csv`` plus the input cDNA FASTA and writes
    two files: the per-probe selection and a flagged-genes report (genes that
    can't reach full isoform coverage or N spaced probes). Used as the final
    step of the HPC job so the bundle ships both the full candidate list and the
    top-N-per-gene picks. Safe to run standalone on any discovery-mode run.
    """
    import csv as _csv
    import time as _time
    from .isoform_select import (
        build_gene_index, coverage_pick, gene_index_from_manifest,
        results_transcripts, write_selection_csvs,
    )
    from .writers.manifest import record_stage_elapsed

    if not os.path.exists(args.results):
        print(f"[flexiprobes] ERROR: {args.results} not found", file=sys.stderr)
        return 2

    # This step runs as its own invocation after the pipeline, so it has to time
    # itself and fold the result back into the manifest — otherwise its cost is
    # invisible to anyone reading provenance from the results alone.
    _t0 = _time.monotonic()

    _csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    out_dir = os.path.dirname(os.path.abspath(args.results))
    out_path = args.out or os.path.join(out_dir, "selected_probes_per_gene.csv")
    flagged_path = args.flagged or os.path.join(out_dir, "flagged_genes.csv")

    with open(args.results, newline="") as fh:
        rows = list(_csv.DictReader(fh))

    # Prefer the FASTA-free path: a current run persists each probe's bound
    # isoforms (isoform_hit_ids) and the gene->isoform map (manifest), so no
    # input FASTA is needed. Fall back to the FASTA only for older runs.
    has_col = any((r.get("isoform_hit_ids") or "").strip() for r in rows)
    gene_isoforms = _load_gene_isoforms(out_dir)
    if has_col and gene_isoforms:
        print(f"[flexiprobes] select-coverage: {len(rows):,} candidate rows; "
              f"using persisted isoform_hit_ids (no FASTA needed)")
        index = gene_index_from_manifest(gene_isoforms)
    else:
        if not args.fasta or not os.path.exists(args.fasta):
            print("[flexiprobes] ERROR: this run predates the isoform_hit_ids "
                  "column, so --fasta (the input cDNA FASTA) is required.",
                  file=sys.stderr)
            return 2
        print(f"[flexiprobes] select-coverage: {len(rows):,} candidate rows; "
              f"indexing {args.fasta}")
        index = build_gene_index(args.fasta, results_transcripts(rows))
    pick = coverage_pick(
        rows, index, n_per_gene=args.n, min_gap=args.min_gap,
        preferred_gap=args.preferred_gap,
        relaxed_gap=args.relaxed_gap, floor_gap=args.floor_gap,
        relax_below=args.relax_below,
        max_mismatches=args.max_mismatches,
        top_k=args.max_candidates_per_gene,
    )
    write_selection_csvs(pick, index, rows, out_path, flagged_path)
    s = pick.summary
    elapsed = _time.monotonic() - _t0
    print(f"[flexiprobes] selected {s['probes']:,} probes across {s['genes']:,} "
          f"genes; {s['flagged']:,} gene(s) flagged")
    _print_spacing_tiers({n: s.get(f"probes_{n}", 0)
                          for n in ("green", "orange", "red")},
                         args.min_gap, args.relaxed_gap, args.floor_gap)
    if s.get("probes_red"):
        print(f"[flexiprobes] WARNING: {s['probes_red']:,} probe(s) across "
              f"{s.get('genes_red', 0):,} gene(s) sit closer than "
              f"{args.relaxed_gap} bp — the last rung of the ladder. See the "
              f"spacing_tier column in the flagged report.")
    elif s.get("genes_below_min_gap"):
        print(f"[flexiprobes] {s['genes_below_min_gap']:,} gene(s) fell below "
              f"the {args.min_gap} bp target; see spacing_relaxed in the "
              f"flagged report.")
    print(f"[flexiprobes]   -> {out_path}")
    print(f"[flexiprobes]   -> {flagged_path}")
    print(f"[flexiprobes] select-coverage took {elapsed:,.1f}s")
    if record_stage_elapsed(
        os.path.join(out_dir, "manifest.json"), "select_coverage", elapsed
    ):
        print("[flexiprobes]   -> folded into manifest.json elapsed_seconds")
    return 0


def cmd_export_hpc(args: argparse.Namespace) -> int:
    """Write a PBS/SLURM job bundle for the configured run.

    The same :func:`flexprobe.hpc_export.write_hpc_bundle` the GUI's
    "Export job bundle for HPC" dialog calls, driven from the command line so a
    bundle can be regenerated headlessly, scripted, and cited in a methods
    section. Everything that shapes the run comes from the shared ``design``
    flags (or a ``--config`` written earlier), so a bundle is byte-identical to
    what ``flexiprobes design`` with the same flags would have executed.
    """
    from .hpc_export import HpcResources, write_hpc_bundle

    cfg = _build_run_config(args)
    if not cfg.fasta_path and not args.remote_input_path:
        print(
            "[flexiprobes] ERROR: no input FASTA — pass --fasta, a --config "
            "whose fasta_path is set, or --remote-input-path.",
            file=sys.stderr,
        )
        return 2
    # The scheduler allocation is the real parallelism budget, so it also sets
    # the config's `workers` unless the caller pinned one explicitly. Without
    # this the bundle ships a config saying (say) 1 worker while the job script
    # passes --workers 16, and the two disagree in the manifest.
    if args.workers is None:
        cfg.workers = args.cores
    # Same argument for the screen, which on a whole-transcriptome run is where
    # essentially all the walltime goes. Left at the config's value when the
    # caller pinned one explicitly.
    if getattr(args, "screen_workers", None) is None:
        cfg.screen_workers = args.cores

    remote_dbs = list(args.remote_db_path) or None
    if not args.stage_inputs and not args.remote_input_path:
        print(
            "[flexiprobes] WARNING: --stage-inputs not set and no "
            "--remote-input-path given; config.json will keep "
            f"fasta_path={cfg.fasta_path!r} verbatim. That path must exist on "
            "the cluster.",
            file=sys.stderr,
        )
    if not args.stage_inputs and remote_dbs is None and cfg.screen.db_paths:
        print(
            "[flexiprobes] WARNING: BLAST DB path(s) "
            f"{cfg.screen.db_paths} are carried over verbatim. Pass "
            "--remote-db-path to point them at the cluster.",
            file=sys.stderr,
        )

    job_name = args.job_name or os.path.splitext(
        os.path.basename(cfg.fasta_path or args.remote_input_path or "")
    )[0] or "flexiprobes"
    split = bool(getattr(args, "split_prepare_screen", False))
    # --mem-gb sizes the single job, or the SCREEN job of a split chain. Their
    # defaults differ: the screen job resumes from the prepare checkpoint and
    # was measured and sized separately (DEFAULT_SCREEN_MEM_GB), so it must
    # not inherit the single-job default.
    mem_gb = (args.mem_gb if args.mem_gb is not None
              else (DEFAULT_SCREEN_MEM_GB if split else DEFAULT_MEM_GB))
    resources = HpcResources(
        job_name=job_name,
        walltime=args.walltime,
        mem_gb=mem_gb,
        cores=args.cores,
        links=args.links,
        partition=args.queue,
        account=args.account,
    )
    prepare_resources = HpcResources(
        job_name=f"{job_name}-prepare",
        walltime=args.prepare_walltime,
        mem_gb=args.prepare_mem_gb,
        cores=args.prepare_cores,
        links=args.prepare_links,
        partition=args.queue,
        account=args.account,
    ) if split else None
    try:
        dest = write_hpc_bundle(
            cfg,
            scheduler=args.scheduler,
            resources=resources,
            image=args.image,
            dest_dir=args.dest,
            stage_inputs=args.stage_inputs,
            remote_input_path=args.remote_input_path,
            remote_db_paths=remote_dbs,
            expected_version=args.expected_version,
            apptainer_module=args.apptainer_module,
            split_prepare_screen=split,
            prepare_resources=prepare_resources,
            job_log_project=args.job_log_project,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"[flexiprobes] ERROR: {exc}", file=sys.stderr)
        return 2

    ext = "slurm" if args.scheduler == "slurm" else "pbs"
    submit_cmd = "sbatch" if args.scheduler == "slurm" else "qsub"
    print(f"[flexiprobes] wrote HPC bundle: {dest}")
    if split:
        # Both requests printed, because the whole point of the split is that
        # they differ — and a reader checking whether the bundle asks for what
        # they intended must not have to open two files to find out.
        assert prepare_resources is not None
        print(f"[flexiprobes]   scheduler={args.scheduler} job={job_name} "
              f"(SPLIT: prepare + screen)")
        print(f"[flexiprobes]   prepare: walltime={prepare_resources.walltime} "
              f"mem={prepare_resources.mem_gb}gb cores={prepare_resources.cores} "
              f"links={prepare_resources.links}")
        print(f"[flexiprobes]   screen:  walltime={args.walltime} "
              f"mem={mem_gb}gb cores={args.cores} links={args.links}")
    else:
        print(f"[flexiprobes]   scheduler={args.scheduler} job={job_name} "
              f"walltime={args.walltime} mem={mem_gb}gb cores={args.cores}")
    print(f"[flexiprobes]   image={args.image} (asserted at run time as "
          f"flexiprobes {args.expected_version})")
    if split:
        print(f"[flexiprobes]   submit with: ./submit_chain.sh   (queues "
              f"{PREPARE_SCRIPT_STEM}.{ext} then {SCREEN_SCRIPT_STEM}.{ext})")
    else:
        print(f"[flexiprobes]   submit with: ./submit_chain.sh   (or "
              f"{submit_cmd} run_flexiprobes.{ext})")
    return 0


def cmd_version(args: argparse.Namespace) -> int:
    print(f"flexiprobes {__version__}")
    return 0


# --------------- main parser ---------------

def _add_design_config_args(p: argparse.ArgumentParser) -> None:
    """Register every flag that shapes a :class:`RunConfig`.

    Shared by ``design`` and ``export-hpc`` so a bundle can be generated
    headlessly from the exact same flags that would have run the job locally.
    Splitting these out is what makes the two commands provably agree: both
    feed the same namespace to :func:`_build_run_config`, so an option added
    to one is an option on the other. Only ``--quiet`` (progress printing)
    stays on ``design``, because it shapes the run, not the config.
    """
    p.add_argument("--config", default=None,
                          help="Load a full RunConfig from a JSON file (as written "
                               "by the GUI's HPC export) and run its values. Only "
                               "these override it, so the bundle can be relocated "
                               "and resized on a cluster: --fasta, --out, "
                               "--off-target-ref, --cached-hits, --skip-screen, "
                               "--prepare-only, --[no-]exclude-alt-scaffolds, "
                               "--reference-strand, --rebuild-kmer-index, "
                               "--workers, --screen-workers, "
                               "--force, --no-resume, --screen-top-per-gene, "
                               "--emit-prescreen-candidates, --emit-prescreen-precap "
                               "and --keep-candidate-cache. Any other design flag "
                               "is refused with an error, not ignored.")
    p.add_argument("--fasta", default=None,
                          help="Input FASTA file (required unless --config sets it).")
    p.add_argument("--out", default=None,
                          help="Output directory (default flexprobe_out, or the "
                               "value in --config).")
    p.add_argument("--modality", choices=list(MODALITY_PRESETS), default="flex_v2",
            help=('Chemistry preset: flex_v1, flex_v2 (default) or spatial '
             '(Visium v1). Sets the probe handles and the junction rule.'))
    p.add_argument("--singleplex", action="store_true",
                          help=('Use the singleplex RHS handle (default is multiplex). '
             "Changes the output for Flex v2 only: Flex v1's two handles "
             'are identical, and Visium has none.'))
    p.add_argument("--probe-barcode", action="append", default=[],
                          help="Probe Barcode label (BC001-BC016 for Flex v1; "
                               "ignored otherwise). Repeatable: pass once for "
                               "singleplex, multiple times for Flex v1 "
                               "multi-sample multiplex.")
    # default=None distinguishes "not passed" from an explicit True/False, which
    # is what lets us reject --require-tn-junction --tn-junction-mode soft as a
    # real contradiction without making the soft mode unreachable.
    g_tn = p.add_mutually_exclusive_group()
    g_tn.add_argument("--require-tn-junction", dest="require_tn_junction",
                      action="store_true", default=None,
                      help="Enforce the TN junction: reject tiles whose mRNA "
                           "position 25 is not A. ON BY DEFAULT (required for "
                           "Flex v1 / Visium per CG000621; recommended for "
                           "Flex v2 per CG000839).")
    g_tn.add_argument("--no-require-tn-junction", dest="require_tn_junction",
                      action="store_false",
                      help=('Disable the ligation-junction gate entirely. Note only ~51%% '
             "of 10x's own probes carry the junction adenine, so this is "
             'defensible, but state it in your methods. Overridden by an '
             'explicit --tn-junction-mode.'))
    p.add_argument("--tn-junction-mode", choices=list(TN_MODES),
                          default=TN_MODE_AUTO,
                          help="Ligation-junction gate. 'auto' (default) "
                               "follows --require-tn-junction. 'strict_a' is "
                               "the hard adenine gate. 'soft' relaxes it to "
                               "'anything but C' — Flex v1/v2 only, per "
                               "CG000839 ('the TN motif is not always strictly "
                               "required'); 10x's own panels carry the A only "
                               "~51%% of the time but avoid C ~99%% of the "
                               "time. 'off' disables the gate.")
    p.add_argument("--tile-size", type=int, default=50,
            help=('Window length, bp (default 50). Every current chemistry '
             'needs 50; any other value is refused, because chemistry '
             'assembly can build no probe from it.'))
    p.add_argument("--tile-step", type=int, default=1,
            help=('Step between candidate windows along each record, bp '
             '(default 1 = every position).'))
    p.add_argument("--gc-min", type=float, default=0.44,
                          help="Minimum GC fraction per 25-mer half (default "
                               "0.44). Only a hard filter with --enforce-gc-range.")
    p.add_argument("--gc-max", type=float, default=0.72,
                          help="Maximum GC fraction per 25-mer half (default "
                               "0.72). Only a hard filter with --enforce-gc-range.")
    p.add_argument("--enforce-gc-range", action="store_true",
                          help=('Restore the hard GC 44-72%%/half band as a filter (pre-0.4 '
             "behaviour). Off by default: 8.4%% of 10x's own probe halves "
             'violate that band, while none falls below the Tm floor, so '
             '--tm-min is the duplex-stability gate. --gc-min/--gc-max '
             'have no effect without this flag.'))
    p.add_argument("--tm-min", type=float, default=TM_MIN_DEFAULT,
                          help="Per-half Tm floor in degrees C (default "
                               f"{TM_MIN_DEFAULT}). RELATIVE (primer3 under "
                               "PCR-primer salt; 10x does not publish the Flex "
                               "buffer), derived from the minimum per-half Tm "
                               "of 10x's shipped panels. Pass a negative value "
                               "to disable.")
    p.add_argument("--tm-max", type=float, default=None,
                          help="Per-half Tm ceiling in degrees C. Off by "
                               "default: 10x's panels show no upper boundary.")
    p.add_argument("--off-target-either-half", action="store_true",
                          help="Reject a pair when EITHER half has an off-target "
                               "match. Pre-0.4 behaviour, contrary to CG000621/"
                               "CG000839 ('at least five mismatches in at least "
                               "one of the LHS or RHS probes'); kept for "
                               "reproducing earlier runs.")
    p.add_argument("--screen-engine", choices=list(ENGINES),
                          default=SCREEN_ENGINE_DEFAULT,
                          help="Search engine for the off-target screen. "
                               "'kmer' (default) searches an exhaustive "
                               "k-mer index end-to-end, so every mismatch "
                               "count is exact by construction and no "
                               "blastdbcmd retrieval is needed; it takes the "
                               "reference FASTA, not a BLAST database, and its "
                               "bitscore/e-value columns are undefined. "
                               "'blast' runs blastn-short and reconstructs "
                               "each probe half's unaligned positions with "
                               "blastdbcmd — needed for a reference the k-mer "
                               "index cannot address. 'blast' CANNOT run the "
                               "junction-duplex gate (--off-target-junction-"
                               "duplex), which is on by default, so 'blast' is "
                               "refused unless you also pass "
                               "--no-off-target-junction-duplex (or "
                               "--off-target-junction-duplex 0), or "
                               "--skip-screen.")
    p.add_argument("--off-target-same-subject",
                          dest="off_target_same_subject",
                          default=OFF_TARGET_SAME_SUBJECT_DEFAULT,
                          action="store_true",
                          help=('Reject a pair only when both halves match the SAME '
             'off-target gene. ON BY DEFAULT — an LHS on gene X and an RHS '
             'on gene Y are never adjacent on one molecule, so the ligase '
             'has nothing to seal and no false signal is produced. No '
             'effect with --off-target-either-half.'))
    p.add_argument("--no-off-target-same-subject",
                          dest="off_target_same_subject",
                          action="store_false",
                          help=('Score the off-target disjunction over the whole probe '
             'instead of per off-target gene, so a pair is rejected when '
             'its two halves are close to DIFFERENT genes. On 10x\'s full '
             'public Flex v1 human panel it rejects '
             '17.6%% of shipped probes vs 6.8%% under the default '
             '(constants.py). Kept for reproducing earlier runs.'))
    # default=None, not the constant: a --config run must be able to tell "the
    # user said nothing" from "the user asked for the default", so a baked HPC
    # bundle keeps whatever it recorded unless a flag actually overrides it.
    # Same reason --workers defaults to None.
    p.add_argument("--exclude-alt-scaffolds",
                          dest="exclude_alt_scaffolds",
                          default=None,
                          action="store_true",
                          help="Drop ALT haplotype and assembly PATCH scaffold "
                               "records from the off-target reference. ON BY "
                               "DEFAULT — they are 7.6%% of GRCh38 cdna.all's "
                               "bases and produced 100%% of its spurious "
                               "perfect-50bp off-target calls, because an "
                               "ALT/PATCH record is a near-verbatim copy of a "
                               "primary region. No-op for a reference without "
                               "Ensembl location fields (viral panels, custom "
                               "FASTAs); k-mer engine only.")
    p.add_argument("--no-exclude-alt-scaffolds",
                          dest="exclude_alt_scaffolds",
                          action="store_false",
                          help="Search ALT/PATCH scaffolds as off-target "
                               "subjects too. Kept for reproducing runs made "
                               "without the exclusion. Expect good probes "
                               "to be rejected for matching a second copy of "
                               "their own gene.")
    # default=None, not the constant: same reason as --exclude-alt-scaffolds
    # above — a --config run must distinguish "the user said nothing" from
    # "the user asked for the default" so a baked HPC bundle's value survives.
    p.add_argument("--rebuild-kmer-index", dest="rebuild_kmer_index",
                   action="store_true", default=False,
                   help="Rebuild a k-mer index whose builder_sha256 does not "
                        "match this build's k-mer index code. Without it such "
                        "an index stops the run with a named error that "
                        "prints both stamps, and nothing is deleted. k-mer "
                        "engine only.")
    p.add_argument("--reference-strand", choices=list(REFERENCE_STRANDS),
                          dest="reference_strand",
                          default=None,
                          help="Which strand(s) of the off-target reference "
                               "to search. 'both' (default) — a match on "
                               "either strand counts, which is "
                               "safe for any reference. "
                               f"{REFERENCE_STRAND_PLUS!r} skips minus-strand/"
                               "reverse-complement matches: a probe is "
                               "single-stranded DNA and physically cannot "
                               "hybridise to one, so counting it is a phantom "
                               "off-target. ONLY correct when the reference "
                               "is sense-oriented cDNA/transcriptome (Ensembl "
                               "cdna.all, GENCODE transcripts, a Cell Ranger "
                               "reference); using it against a genome or a "
                               "negative-sense viral FASTA would silently "
                               "miss real off-targets on the reference's "
                               "other strand, so this is opt-in, not the "
                               "default.")
    p.add_argument("--off-target-junction-duplex", type=int,
                          dest="off_target_junction_duplex",
                          default=OFF_TARGET_JUNCTION_DUPLEX_BP,
                          metavar="N",
                          help="Reject a pair when a NON-TARGET gene presents "
                               "at least N perfectly matched bases on BOTH "
                               "sides of the ligation nick, in contiguous "
                               "reference sequence. ON BY DEFAULT at "
                               f"{OFF_TARGET_JUNCTION_DUPLEX_BP} bp — but ONLY "
                               "under the k-mer engine (the default); "
                               "--screen-engine blast has no equivalent scan, "
                               "so a non-zero N with --screen-engine blast is "
                               "refused (unless --skip-screen). THIS IS "
                               "OURS, NOT 10x's: their rule is a pure count "
                               "and says nothing about where the mismatches "
                               "sit, so a pair whose mismatches are all "
                               "terminal can still hybridise AND ligate across "
                               "the junction on the wrong gene. Measured on "
                               "10x's own public Flex v1 panel, 9 catches 87%% of the "
                               "probes 10x itself rejected while firing on "
                               "0.6%% of the ones it ships — the maximum "
                               "separation of the two (7 fires on 53%% of "
                               "shipped probes; 11 catches only 80%%). Per "
                               "SIDE, not total: 9+9 is rejected, 15+3 is not. "
                               f"0 disables; below {OFF_TARGET_JUNCTION_DUPLEX_MIN_BP} "
                               "is refused, since the scan cannot see that "
                               "short. Rejections are reported as "
                               "specificity_status=failed_junction_duplex, "
                               "never as failed_off_target.")
    p.add_argument("--no-off-target-junction-duplex",
                          dest="off_target_junction_duplex",
                          action="store_const", const=0,
                          help="Disable the junction-duplex gate and apply only "
                               "the published 10x count rule. Equivalent to "
                               "--off-target-junction-duplex 0. Use it to "
                               "reproduce a run made without the gate, or to "
                               "screen strictly to the vendor spec with "
                               "nothing added.")
    p.add_argument("--off-target-local-alignment", action="store_true",
                          help="Score off-targets on a local alignment "
                               "(>=80%% identity over >=20 bp) instead of "
                               "counting mismatches across the full 25 nt half. "
                               "Legacy rule; kept for reproducing "
                               "earlier runs.")
    p.add_argument("--max-n-per-half", type=int, default=1,
                          help="Allowed degenerate bases (N or any IUPAC code) per 25 nt "
                               "half (default 1).")
    p.add_argument("--max-homopolymer-run", type=int, default=4,
                          help="Reject homopolymer runs longer than this (default 4).")
    p.add_argument("--off-target-ref", "--blast-db", action="append", default=[],
                          dest="off_target_ref", metavar="PATH",
                          help=('Off-target reference to screen against, repeatable. Normally '
             'a FASTA. Under --screen-engine blast it may instead be a '
             'BLAST DB prefix. If omitted, the screen is skipped with a '
             'warning and probes are not specificity-gated; pass '
             '--skip-screen to skip it on purpose. `--blast-db` is a '
             'deprecated alias kept for existing scripts: the default '
             'engine is the k-mer screen, not BLAST, and the argument was '
             'never required to be a BLAST database.'))
    p.add_argument("--skip-screen", action="store_true",
                          help="Skip the off-target screen entirely (useful for "
                               "smoke tests). The run then makes no specificity claim.")
    p.add_argument("--prepare-only", action="store_true",
                          help=('Run ingest + tile + filter, write the candidate checkpoint '
             'and the pre-screen artefacts, then STOP before the '
             'off-target screen. Produces no probes and no results.csv: '
             'the output directory is INCOMPLETE by design, marked with '
             'the .flexiprobes_PREPARED file and a design_summary.json '
             "status of 'prepared'. Re-run the same config WITHOUT this "
             'flag to resume from the checkpoint and finish. This is NOT '
             '--skip-screen, which designs probes with no screen and '
             'produces final output. Use it to give the memory-hungry '
             'tiling half and the parallel screen half their own scheduler '
             'allocations (see export-hpc --split-prepare-screen). Refused with '
             '--no-resume, --skip-screen, --cached-hits or --multi-strain.'))
    p.add_argument("--cached-hits", default=None,
                          help=('JSON of pre-computed off-target screen hits, in the '
             'structure flexprobe.blast.load_cached_hits reads; nothing in '
             'this package writes one. If it cannot be loaded, a warning '
             'is printed and the screen runs normally.'))
    # primer3 scoring is ON by default (365 us/candidate — noise beside BLAST).
    # --enable-primer3 is kept as an accepted no-op so every committed script,
    # PBS bundle and run_config.json that names it keeps working rather than
    # dying on an unrecognised argument.
    p.add_argument("--enable-primer3", action="store_true",
                          help="No-op, kept for compatibility: primer3 "
                               "Tm/hairpin/dimer scoring is on by default.")
    p.add_argument("--no-primer3", dest="no_primer3", action="store_true",
                          help=('Disable primer3 Tm / hairpin / dimer scoring. Structure '
             'terms drop out of the composite score, and a degenerate '
             "consensus parent's aggregate columns (hairpin_dg_worst, "
             'tm_balance_max, the Tm/GC ranges) are left blank, since a '
             "parent can only ever report its variants' worst case. The "
             '--tm-min/--tm-max gate still applies, and the Tm/GC ranges '
             'of ordinary probes are still reported.'))
    p.add_argument("--query-chunk", "--blast-query-chunk", type=int, default=None,
                   dest="query_chunk", metavar="N",
                   help="Probe queries per screen chunk (default 5000). "
                        "Bounds peak memory, sets progress granularity, and is "
                        "the unit a killed run resumes at.")
    p.add_argument("--max-target-seqs", type=int, default=None,
                   help="--screen-engine blast only; the k-mer engine has "
                        "its own per-query cap and ignores this. Max subjects "
                        "blastn may return per query half (default 100000). "
                        "BLAST's own default of 500 is applied heuristically "
                        "DURING the search and can silently drop qualifying "
                        "off-targets, so lower this only if you know why.")
    p.add_argument("--min-probe-gap", type=int, default=None,
                   help=f"Target spacing between chosen probes, bp (default "
                        f"{MIN_PROBE_GAP_BP}). Used only by "
                        f"--auto-select-per-gene and the export-hpc selection "
                        f"step. A gene that cannot yield "
                        f"{PROBE_GAP_RELAX_MIN_PROBES} probes at this spacing "
                        f"steps down to --relaxed-probe-gap, then "
                        f"--floor-probe-gap, and is flagged.")
    p.add_argument("--floor-probe-gap", type=int, default=None,
                   help=f"Last-rung spacing, bp (default "
                        f"{MIN_PROBE_GAP_FLOOR_BP}). A rung is used only when it "
                        f"is smaller than the rung above it, so a value >= "
                        f"--relaxed-probe-gap drops this rung and the closest "
                        f"allowed pair is then --relaxed-probe-gap.")
    p.add_argument("--relaxed-probe-gap", type=int, default=None,
                   help=f"Fallback spacing for genes starved by "
                        f"--min-probe-gap, bp (default {RELAXED_PROBE_GAP_BP}). "
                        f"A value >= --min-probe-gap skips this rung only: the "
                        f"--floor-probe-gap rung still runs. To allow no "
                        f"fallback at all, set both to >= --min-probe-gap.")
    p.add_argument("--junction-window", type=int, default=None,
                   help="Half-width in bp of the ligation-junction zone used "
                        "for SCORING (default 4, CG000621): the per-strain dTm "
                        "penalty doubles inside it and off-target hits report "
                        "their mismatches in it. Placement of degenerate bases "
                        "is governed separately.")
    p.add_argument("--split-header-on-pipe", action="store_true",
                          help="Legacy: truncate each FASTA record id at the "
                               "first '|'. Off by default because it merges "
                               "records that share a first pipe-field (an "
                               "NCBI-style '>gi|123|ref|...' set collapses to a "
                               "single id). Only use when the first field really "
                               "is the unique id.")
    p.add_argument("--workers", type=int, default=None,
                          help="Worker processes for tiling+filtering on "
                               "multi-gene runs (default 1 = serial). With "
                               "--config, an explicit value overrides the baked "
                               "config (pass 1 to fit a smaller allocation).")
    p.add_argument("--screen-workers", type=int, default=None,
                          help=('Worker processes for the OFF-TARGET SCREEN (default 1 = '
             'serial). This is the stage that dominates a '
             'whole-transcriptome run; --workers only covers '
             'tiling+filtering. k-mer engine only. Not part of any cache '
             'fingerprint, so raising it resumes an in-flight run rather '
             'than invalidating its checkpoints. Capped at the usable '
             "cores, and serial where the 'fork' start method is "
             'unavailable.'))
    p.add_argument("--extra-on-target", action="append", default=[],
                          metavar="ACCESSION",
                          help="Accession that must NOT count as off-target "
                               "(e.g. another isoform of the target gene). "
                               "Repeatable.")
    # --- Gene-aware / whole-transcriptome isoform mode ---
    g_gbg = p.add_mutually_exclusive_group()
    g_gbg.add_argument("--group-by-gene", dest="group_by_gene",
                       action="store_true", default=None,
                       help="Isoform-aware mode, ON BY DEFAULT: group input "
                            "transcripts by gene (Ensembl 'gene:' header field) "
                            "so sibling isoforms count as on-target, and score "
                            "each probe on isoform coverage. Pair the input and "
                            "--off-target-ref with the same transcriptome. On "
                            "input with no 'gene:' field it has nothing to group "
                            "and the run says so and proceeds ungrouped. Ignored "
                            "with --multi-strain.")
    g_gbg.add_argument("--no-group-by-gene", dest="group_by_gene",
                       action="store_false",
                       help="Turn gene grouping off: every input record counts "
                            "as on-target for every probe.")
    p.add_argument("--screen-top-per-gene", type=int, default=None,
                          metavar="N",
                          help=('Whole-transcriptome speedup: off-target screen only the top '
             'N candidates per gene (by the screen-free pre-score), '
             'dropping the rest before the screen, under either engine. A '
             "'gene' is the record id unless --group-by-gene groups "
             'isoforms. Ignored when no screen runs. 0 = screen all '
             '(default).'))
    p.add_argument("--min-isoform-coverage", type=float, default=0.0,
                          metavar="FRAC",
                          help=('With --group-by-gene: drop probes that bind less than FRAC '
             "(0-1) of their gene's isoforms (0 = no gate, just soft "
             'scoring). Ignored with --multi-strain.'))
    p.add_argument("--discovery-only", action="store_true",
                          help="Write only results.csv + manifest.json (rank every "
                               "candidate) and stop — pick the final probes later in "
                               "the GUI review dialog. Implied by --auto-select-per-gene.")
    p.add_argument("--auto-select-per-gene", type=int, default=0,
                          metavar="N",
                          help="Whole-transcriptome convenience: after ranking, "
                               "auto-pick the top N spaced probes per gene and write "
                               "the order CSVs without manual review (forces "
                               "discovery mode). 0 = off.")
    # --- Multi-strain / isoform consensus modes ---
    p.add_argument("--multi-strain", action="store_true",
                          help="Treat every FASTA record as a related sequence "
                               "(strain or isoform) of one fragment.")
    p.add_argument("--strategy",
                          choices=["conserved", "strain_unique", "constitutive"],
                          default="conserved",
                          help=('Consensus strategy when --multi-strain: conserved '
             '(pan-record; by default tolerates one variant per half, see '
             '--conserved-max-degenerate-per-half), strain_unique (one vs '
             'rest), or constitutive (isoform-collapse: shared by >= '
             '--isoform-min records).'))
    p.add_argument("--isoform-min", type=int, default=0,
                          help="constitutive strategy: min records a window must "
                               "be shared by (0 = all). ")
    p.add_argument("--conserved-max-degenerate-per-half", type=int, default=None,
                          help="conserved strategy: variant positions tolerated "
                               "per 25 nt half (default 1 = a window may carry "
                               "one variant per half provided none falls within "
                               "the junction clearance; each such window is "
                               "emitted as a degenerate consensus parent (never "
                               "orderable) plus one concrete variant per "
                               "record). 0 requires the window to be identical "
                               "across every record (the strict rule). Exact "
                               "identity is severe on a real strain pair: it "
                               "can leave a whole segment with no designable "
                               "window once the chemistry filters run, where "
                               "1 per half yields many times more windows.")
    p.add_argument("--strain-unique-target", type=int, default=0,
                          help="strain_unique strategy: 0-based index of the "
                               "target record.")
    p.add_argument("--strain-unique-min-mm", type=int, default=5,
                          help="strain_unique strategy: min mismatches vs every "
                               "other record (default 5), counted within a "
                               "single 25 nt half unless "
                               "--strain-unique-whole-window is given.")
    p.add_argument("--strain-unique-whole-window", action="store_true",
                          help="strain_unique strategy: count "
                               "--strain-unique-min-mm over the whole 50 nt "
                               "window instead of within a single 25 nt half. "
                               "Pre-0.4 behaviour; contrary to CG000621/"
                               "CG000839 ('at least five mismatches in at least "
                               "one of the LHS or RHS probes'), which is where "
                               "the 5 comes from. Mismatches split across the "
                               "junction leave both halves able to bind the "
                               "non-target record, so the probe still ligates "
                               "there. Kept for reproducing earlier runs.")
    p.add_argument("--strain-unique-junction-duplex", type=int,
                          dest="strain_unique_junction_duplex",
                          default=OFF_TARGET_JUNCTION_DUPLEX_BP,
                          metavar="N",
                          help="strain_unique strategy (per-half mode only): "
                               "withhold 'strain-unique' from a window when a "
                               "NON-TARGET record still presents at least N "
                               "perfectly matched bases on BOTH sides of the "
                               "ligation nick, in contiguous sequence — the "
                               "mirror image of --off-target-junction-duplex, "
                               "reusing the same scan and the same default "
                               f"{OFF_TARGET_JUNCTION_DUPLEX_BP} bp. A window can "
                               "reach --strain-unique-min-mm with every "
                               "mismatch out at the probe's terminal end, which "
                               "leaves the nick-spanning region untouched — that "
                               "still hybridises and ligates on the 'excluded' "
                               "record. THIS IS OURS, NOT 10x's, for the same "
                               "reason the off-target twin is. 0 disables; "
                               f"below {OFF_TARGET_JUNCTION_DUPLEX_MIN_BP} is "
                               "refused, since the scan cannot see that short.")
    # --- pre-screen candidate set as an analysis artefact ---
    p.add_argument("--emit-prescreen-candidates", action="store_true",
                   help=('Write the candidate set the off-target screen actually '
             'screened to <out>/candidates_prescreen.csv.gz, per-gene '
             'pre-cap vs screened counts to '
             '<out>/prescreen_cap_summary.csv.gz, and '
             '<out>/prescreen_index.json, which states which candidate set '
             'the files hold and their schemas. results.csv holds only the '
             'candidates that SURVIVED the screen, so any null model built '
             'from it is conditioned on the outcome; these files are the '
             'denominator. Off by default (large extra file).'))
    p.add_argument("--emit-prescreen-precap", action="store_true",
                   help="Additionally dump EVERY filter-passing window, before "
                        "the --screen-top-per-gene cap, to "
                        "<out>/candidates_precap.csv.gz. The cap is itself an "
                        "outcome-correlated selection step, so this is the "
                        "stricter denominator — and 20-50x the rows (hundreds "
                        "of millions on a whole transcriptome under a soft "
                        "junction gate). Implies "
                        "--emit-prescreen-candidates.")
    # --- null-model calibration of the off-target gate ---
    p.add_argument("--null-model", action="store_true",
                   help=("Calibrate the off-target gate against noise: screen "
                        "dinucleotide-shuffled copies of a bounded sample of "
                        "this run's windows through the SAME gate and write "
                        "the resulting false-positive rate to "
                        "<out>/null_model.json. A gate whose rejection rate on "
                        "composition-matched random sequence is unknown cannot "
                        "be defended in a methods section, and on test "
                        "panels that rate varies widely between runs. "
                        "Costs at most --null-model-max-windows x "
                        "--null-model-shuffles (default "
                        f"{NULL_MODEL_DEFAULT_MAX_WINDOWS} x "
                        f"{NULL_MODEL_DEFAULT_SHUFFLES}) extra screen queries "
                        "regardless of run size. Skipped when no screen runs."))
    p.add_argument("--null-model-shuffles", type=int, default=None,
                   help="Shuffles per sampled window (implies --null-model). "
                        f"Default {NULL_MODEL_DEFAULT_SHUFFLES} when "
                        "--null-model is given; 0 disables. Below ~20 the "
                        "per-window rate is too coarse to quote.")
    p.add_argument("--null-model-max-windows", type=int,
                   default=NULL_MODEL_DEFAULT_MAX_WINDOWS,
                   help="Cap on how many of the run's windows are sampled for "
                        f"calibration (default {NULL_MODEL_DEFAULT_MAX_WINDOWS}). "
                        "This is what bounds the cost on a whole-transcriptome "
                        "run. Sampled at random, not truncated: the first N "
                        "windows all come from one gene and their composition "
                        "is not the run's composition.")
    p.add_argument("--null-model-seed", type=int, default=42,
                   help="Seed for the shuffle RNG (default 42), so a quoted "
                        "false-positive rate is reproducible.")
    p.add_argument("--keep-candidate-cache", action="store_true",
                   help=('Keep <out>/.candidate_cache/ after a successful run instead '
             'of deleting it. A second, independent copy of the pre-screen '
             'candidate set, insurance for runs whose only other copy is '
             'the artefact above. The cache is written only when resume is '
             'on and the run screens by itself: not with --no-resume, '
             '--skip-screen, --cached-hits or --multi-strain.'))
    p.add_argument("--no-idt", action="store_true",
                   help="Skip the IDT oPools order file, including the one "
                        "--auto-select-per-gene writes.")
    p.add_argument("--no-cellranger", action="store_true",
                          help="Skip the Cell Ranger probe-set file, including "
                               "the one --auto-select-per-gene writes.")
    p.add_argument("--force", action="store_true",
                          help="Overwrite existing output files. With "
                               "--auto-select-per-gene, the order files this "
                               "run writes are then replaced by the selected "
                               "subset.")
    p.add_argument("--no-resume", action="store_true",
                          help=('Disable resumable screening, for either engine: no '
             '<out>/.screen_cache/ or <out>/.candidate_cache/ checkpoints '
             'are written, so a killed run starts again from the '
             'beginning. By default each screen chunk is checkpointed so a '
             'run killed by an HPC walltime can be resubmitted and '
             'continue. Refused with --prepare-only, whose checkpoint is '
             'its whole output.'))
    p.add_argument("--screen-cache-mode", choices=("verdicts", "hits"),
                          default="verdicts",
                          help=('What each screen checkpoint stores, for a screen streamed by '
             "this run (not with --cached-hits). 'verdicts' (default) "
             'keeps the scored result per probe: its size is independent '
             'of hit volume, the only safe choice at whole-transcriptome '
             "scale. 'hits' keeps the raw hits, so you can re-score at a "
             'new --off-target-* threshold without re-running the screen, '
             'but it grows with hits per probe (2,884/probe on human tx = '
             '~1.2 TB over a full run, which will fill your disk). Use '
             "'hits' only on small runs while tuning thresholds."))


#: The design flags ``_build_run_config`` applies on top of a ``--config``.
#: Any other design flag would otherwise be dropped without a word, so a
#: cluster run that looked configured by its command line would silently run
#: on the file's values. Keep in step with the override block above;
#: ``test_cli_audit_group_b`` checks the two agree.
CONFIG_OVERRIDES = frozenset({
    "config", "fasta", "out", "off_target_ref", "cached_hits", "skip_screen",
    "prepare_only", "exclude_alt_scaffolds", "reference_strand",
    "rebuild_kmer_index", "workers",
    "screen_workers", "force", "no_resume", "screen_top_per_gene",
    "emit_prescreen_candidates", "emit_prescreen_precap",
    "keep_candidate_cache",
})


def _flags_ignored_by_config(args: argparse.Namespace) -> list:
    """Design flags given on the command line that ``--config`` would ignore.

    A flag counts as given when its value differs from the design parser's
    default, the only record argparse keeps. Returned as ``--flag`` strings,
    sorted, so the error can name them.
    """
    if not getattr(args, "config", None):
        return []
    probe = argparse.ArgumentParser(add_help=False)
    _add_design_config_args(probe)
    defaults = vars(probe.parse_args([]))
    given = vars(args)
    names = []
    for dest, default in defaults.items():
        if dest in CONFIG_OVERRIDES or dest not in given:
            continue
        if given[dest] != default:
            names.append("--" + dest.replace("_", "-"))
    return sorted(names)


def _unused_flag_warnings(args: argparse.Namespace) -> list:
    """Warnings for flags that were given but that this run will not use.

    --max-target-seqs only reaches blastn; the k-mer engine (the default) has
    its own per-query cap.
    """
    out = []
    if (getattr(args, "max_target_seqs", None) is not None
            and getattr(args, "screen_engine", None) != ENGINE_BLAST):
        out.append("--max-target-seqs has no effect under --screen-engine "
                   f"{args.screen_engine}; it applies to --screen-engine blast only.")
    return out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="flexiprobes",
        description="flexiprobes — probe-design CLI (tile, filter, off-target screen (k-mer by "
                    "default, or BLAST), rank, export)",
    )
    parser.add_argument("--verbose", "-v", action="count", default=0,
                        help="Increase log verbosity (-v INFO, -vv DEBUG).")
    sub = parser.add_subparsers(dest="command", required=True)

    p_design = sub.add_parser("design", help="Run the full probe-design pipeline")
    _add_design_config_args(p_design)
    p_design.add_argument("--quiet", action="store_true",
                          help="Suppress per-stage progress output.")
    p_design.set_defaults(func=cmd_design)

    p_vb = sub.add_parser(
        "validate-blast",
        help="Verify the BLAST tools a BLAST-engine run needs (blastn, "
             "makeblastdb, blastdbcmd) and, with --db, that the database "
             "exists (single- or multi-volume) and carries the -parse_seqids "
             "accession index the pipeline requires")
    p_vb.add_argument("--db", default=None,
                      help="Optional BLAST DB prefix to inspect: files present, "
                           "and whether blastdbcmd can retrieve by accession")
    p_vb.add_argument("--out", default="flexprobe_out")
    p_vb.set_defaults(func=cmd_validate_blast)

    p_bd = sub.add_parser("build-db", help="Build a BLAST nucleotide DB from a FASTA")
    p_bd.add_argument("--input", required=True, help="Input FASTA")
    p_bd.add_argument("--db-prefix", default=None,
                      help="Output DB prefix (default: <out>/flexprobe_db)")
    p_bd.add_argument("--title", default="flexprobe_db")
    p_bd.add_argument("--out", default="flexprobe_out")
    p_bd.set_defaults(func=cmd_build_db)

    p_ex = sub.add_parser("export", help="Re-export an existing results.csv")
    p_ex.add_argument("--results-csv", required=True,
                      help="Path to results.csv from a previous run")
    p_ex.add_argument("--format", default="csv", choices=["csv", "tsv", "json", "bed"])
    p_ex.add_argument("--out", default="flexprobe_out")
    p_ex.add_argument("--force", action="store_true", help="Overwrite existing files")
    p_ex.set_defaults(func=cmd_export)

    p_sc = sub.add_parser(
        "select-coverage",
        help="Pick top-N spaced probes per gene maximising isoform coverage")
    p_sc.add_argument("--results", required=True,
                      help="results.csv from a finished discovery run")
    p_sc.add_argument("--fasta", default=None,
                      help="input cDNA FASTA (only needed for runs that predate "
                           "the isoform_hit_ids column; current runs need no FASTA)")
    p_sc.add_argument("--out", default=None,
                      help="per-probe selection CSV (default: alongside results.csv)")
    p_sc.add_argument("--flagged", default=None,
                      help="flagged-genes report CSV (default: alongside results.csv)")
    p_sc.add_argument("--n", type=int, default=3, help="probes per gene (default 3)")
    p_sc.add_argument("--min-gap", type=int, default=MIN_PROBE_GAP_BP,
                      help=f"target bp between chosen probes (default "
                           f"{MIN_PROBE_GAP_BP}; also forbids overlap). A gene "
                           f"that cannot yield min(--n, --relax-below) probes "
                           f"at this spacing is retried at --relaxed-gap, then "
                           f"--floor-gap, and flagged")
    p_sc.add_argument("--relaxed-gap", type=int, default=RELAXED_PROBE_GAP_BP,
                      help=f"fallback spacing for genes starved by --min-gap "
                           f"(default {RELAXED_PROBE_GAP_BP}). A value >= "
                           f"--min-gap skips this rung only: the --floor-gap "
                           f"rung still runs. To allow no fallback, set both "
                           f">= --min-gap")
    p_sc.add_argument("--floor-gap", type=int, default=MIN_PROBE_GAP_FLOOR_BP,
                      help=f"last-rung spacing, bp (default "
                           f"{MIN_PROBE_GAP_FLOOR_BP}). A rung is used only when "
                           f"it is smaller than the rung above it, so a value "
                           f">= --relaxed-gap drops it and red-tier probes "
                           f"become unselectable")
    p_sc.add_argument("--relax-below", type=int, default=PROBE_GAP_RELAX_MIN_PROBES,
                      help=f"relax spacing only for genes yielding fewer than "
                           f"min(--n, this) probes at --min-gap (default "
                           f"{PROBE_GAP_RELAX_MIN_PROBES})")
    p_sc.add_argument("--preferred-gap", type=int, default=PREFERRED_PROBE_GAP_BP,
                      help=f"soft tie-break, bp (default "
                           f"{PREFERRED_PROBE_GAP_BP}): among candidates adding "
                           f"equal new isoform coverage, prefer one this far "
                           f"from those already chosen. Never overrides "
                           f"--min-gap")
    p_sc.add_argument("--max-mismatches", type=int, default=0,
                      help="mismatches allowed when matching a probe to an "
                           "isoform by sequence (0 = exact). Ignored for runs "
                           "whose results.csv carries isoform_hit_ids, which "
                           "record the membership directly")
    p_sc.add_argument("--max-candidates-per-gene", type=int, default=600,
                      help="cap on candidates considered per gene (top by score)")
    p_sc.set_defaults(func=cmd_select_coverage)

    p_hpc = sub.add_parser(
        "export-hpc",
        help="Write a PBS/SLURM job bundle for this run (the GUI's HPC export, "
             "headless)",
        description="Write a self-contained PBS/SLURM job bundle (config.json + "
                    "one batch script, or two with --split-prepare-screen, + "
                    "submit_chain.sh + README.txt) for the run "
                    "described by the same flags 'flexiprobes design' takes. "
                    "Regenerating a bundle from a recorded command line is what "
                    "makes a cluster run reproducible.",
    )
    _add_design_config_args(p_hpc)
    p_hpc.add_argument("--dest", required=True, metavar="DIR",
                       help=('Directory to write the bundle into (created if absent). '
             'Bundle files are rewritten in place and nothing else is '
             'removed, so a leftover output/.flexiprobes_DONE from an '
             'earlier run makes the new chain skip its work: export into '
             'an empty directory.'))
    p_hpc.add_argument("--scheduler", choices=list(SCHEDULERS), default="slurm",
                       help="Batch system (default slurm).")
    p_hpc.add_argument("--image", required=True, metavar="PATH",
                       help=('Path to the flexiprobes .sif ON THE CLUSTER. A relative path '
             'resolves from the bundle directory the job runs in. The job '
             'script asserts the engine version inside it before doing any '
             'work.'))
    p_hpc.add_argument("--job-name", default=None,
                       help="Scheduler job name (default: the input FASTA's "
                            "basename without extension).")
    p_hpc.add_argument("--walltime", default=DEFAULT_WALLTIME,
                       metavar="HH:MM:SS",
                       help=f"Walltime PER CHAIN LINK, not for the whole run "
                            f"(default {DEFAULT_WALLTIME}). submit_chain.sh "
                            f"queues --links copies of the job and each resumes "
                            f"from the last checkpoint, so the total budget is "
                            f"walltime x links. Raise --links rather than this: "
                            f"a short job is backfilled into gaps a long one "
                            f"cannot use, and a job whose walltime would run "
                            f"past a scheduled maintenance window is held until "
                            f"after the window instead of starting.")
    p_hpc.add_argument("--links", type=int, default=DEFAULT_CHAIN_LINKS,
                       metavar="K",
                       help=f"How many links submit_chain.sh queues by default "
                            f"(default {DEFAULT_CHAIN_LINKS}). Over-provisioning "
                            f"is free: once the run finishes, every remaining "
                            f"link sees the DONE sentinel and exits in about a "
                            f"second. Can still be overridden at submit time "
                            f"with ./submit_chain.sh <K>.")
    p_hpc.add_argument("--mem-gb", type=int, default=None,
                       metavar="GB",
                       help=f"Memory request in GB for the job, or for the "
                            f"SCREEN job with --split-prepare-screen (default "
                            f"{DEFAULT_MEM_GB}, or {DEFAULT_SCREEN_MEM_GB} for "
                            f"the split screen job). Sized for an ordinary "
                            f"multi-gene design; a whole transcriptome needs "
                            f"far more, of the order of 300 GB and 24 cores "
                            f"for the screen.")
    p_hpc.add_argument("--cores", type=int, default=DEFAULT_CORES, metavar="N",
                       help=f"CPU cores to request (default {DEFAULT_CORES}). "
                            f"Also passed to the pipeline as --workers and "
                            f"--screen-workers, and used as the config's "
                            f"`workers` unless --workers is given. Each screen "
                            f"worker holds its own share of the reference, so "
                            f"cores and memory have to be raised together.")
    p_hpc.add_argument("--apptainer-module", default="", metavar="NAME",
                       help="Environment module that puts `apptainer` on PATH "
                            "(e.g. apptainer/1.4.1), baked into the job script. "
                            "A BATCH JOB DOES NOT INHERIT THE MODULES YOU "
                            "LOADED INTERACTIVELY, so without this the chain "
                            "aborts in about a second on any site where "
                            "apptainer needs a module. Find the name with "
                            "`module avail apptainer` on a login node. "
                            "Overridable per submission with "
                            "APPTAINER_MODULE=<name> ./submit_chain.sh.")
    p_hpc.add_argument("--job-log-project", default="", metavar="NAME",
                       help="PBS only. Log every job IN REAL TIME to "
                            "$WORKDIR/projects/NAME/hpc_job_logs/"
                            "<YYYYMMDD_HHMMSS>_<job_name>.<jobid>.log (PROJ in "
                            "the job's environment overrides the root), and "
                            "discard PBS's own spooled log with "
                            "#PBS -o /dev/null. PBS copies that spooled log "
                            "back only when a job ends, so without this a "
                            "running job has nothing readable. NAME may use "
                            "letters, digits, '.', '_' and '-'. Default: off.")
    p_hpc.add_argument("--queue", "--partition", dest="queue", default="",
                       help="PBS queue / SLURM partition (optional).")
    p_hpc.add_argument("--account", default="",
                       help="Allocation / account to charge (optional).")
    p_hpc.add_argument("--stage-inputs", action="store_true",
                       help="Copy the input FASTA into the bundle's inputs/. "
                            "Makes the bundle self-contained, but a "
                            "transcriptome FASTA is hundreds of MB per bundle. "
                            "Prefer --remote-input-path for large inputs "
                            "already on the cluster.")
    p_hpc.add_argument("--remote-input-path", default=None, metavar="PATH",
                       help="Path to the input FASTA AS IT EXISTS ON THE "
                            "CLUSTER. Written into config.json in place of the "
                            "local path. Ignored when --stage-inputs is set.")
    p_hpc.add_argument("--remote-db-path", action="append", default=[],
                       metavar="PATH",
                       help=('Off-target reference path(s) as they exist on the cluster '
             '(repeatable): a FASTA for the k-mer engine (the default), a '
             'BLAST DB prefix for --screen-engine blast. For a '
             'whole-transcriptome run this is the same file as '
             '--remote-input-path: the transcriptome is its own off-target '
             'reference.'))
    p_hpc.add_argument("--split-prepare-screen", action="store_true",
                       help="Emit TWO chained job scripts instead of one: a "
                            "prepare job (ingest + tile + filter, "
                            "--prepare-only, memory-hungry and short) and a "
                            "screen job (resume from the checkpoint, screen, "
                            "rank, export). The same screening work needs several "
                            "times more memory in a process that tiled first "
                            "than in one that resumed from the checkpoint, so "
                            "a single job has to request the larger amount "
                            "for its whole walltime. --mem-gb/--cores/--walltime/"
                            "--links then describe the SCREEN job; the prepare "
                            "job is sized by the --prepare-* options below.")
    p_hpc.add_argument("--prepare-mem-gb", type=int,
                       default=DEFAULT_PREPARE_MEM_GB, metavar="GB",
                       help=f"Memory for the prepare job (default "
                            f"{DEFAULT_PREPARE_MEM_GB}). Sized for the largest "
                            f"tile+filter peak seen on a whole transcriptome, "
                            f"plus headroom; do not trim it. Only used with --split-prepare-screen.")
    p_hpc.add_argument("--prepare-cores", type=int,
                       default=DEFAULT_PREPARE_CORES, metavar="N",
                       help=f"Cores for the prepare job, also passed as "
                            f"--workers (default {DEFAULT_PREPARE_CORES}). The "
                            f"tile stage is single-threaded but costs under 90 "
                            f"s; the FILTER stage is parallel and uses several "
                            f"cores on average, so 1 or 2 here can make the "
                            f"job hit its walltime. Only used with "
                            f"--split-prepare-screen.")
    p_hpc.add_argument("--prepare-walltime", default=DEFAULT_PREPARE_WALLTIME,
                       metavar="HH:MM:SS",
                       help=f"Walltime per prepare link (default "
                            f"{DEFAULT_PREPARE_WALLTIME}), about twice the "
                            f"slowest whole-transcriptome tile+filter. Only used with "
                            f"--split-prepare-screen.")
    p_hpc.add_argument("--prepare-links", type=int,
                       default=DEFAULT_PREPARE_LINKS, metavar="K",
                       help=f"Prepare links to queue (default "
                            f"{DEFAULT_PREPARE_LINKS}): one free retry for a "
                            f"node fault. A prepared bundle's spare link "
                            f"no-ops on the marker. Only used with "
                            f"--split-prepare-screen.")
    p_hpc.add_argument("--expected-version", default=__version__,
                       metavar="VERSION",
                       help=f"Engine version the bundle demands of the "
                            f"container (default {__version__}, this build). "
                            f"The job script runs 'flexiprobes version' inside "
                            f"the .sif and refuses to start on a mismatch — a "
                            f"stale image silently reproduces the old engine's "
                            f"results.")
    p_hpc.set_defaults(func=cmd_export_hpc)

    p_ver = sub.add_parser("version", help="Print version")
    p_ver.set_defaults(func=cmd_version)

    return parser


def main(argv: Optional[list] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    # Genuinely contradictory junction request. Only detectable here, where an
    # explicitly-passed --require-tn-junction (True) is distinguishable from the
    # default (None); by the time resolve_tn_mode sees the config, both look
    # like True and it lets the explicit mode win.
    if (getattr(args, "require_tn_junction", None) is True
            and getattr(args, "tn_junction_mode", None) == TN_MODE_SOFT):
        parser.error(
            "--require-tn-junction demands the junction adenine while "
            "--tn-junction-mode soft tolerates G and T. Pass one or the other: "
            "--tn-junction-mode soft alone is enough (the strict gate is the "
            "default, so the flag is redundant unless you want strict)."
        )
    ignored = _flags_ignored_by_config(args)
    if ignored:
        parser.error(
            f"{', '.join(ignored)} {'is' if len(ignored) == 1 else 'are'} "
            f"ignored with --config, which runs the file's values. Only the "
            f"relocation and resource flags override it (see --config in "
            f"--help). Change the setting in the config file, or drop "
            f"--config and pass every flag."
        )
    levels = ["WARNING", "INFO", "DEBUG"]
    setup_logging(level=levels[min(args.verbose, 2)])
    for warning in _unused_flag_warnings(args):
        print(f"[flexiprobes] WARNING: {warning}", file=sys.stderr)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
