"""ModBAM analysis and result-reading CLI, with explicit input paths."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import shutil
import tempfile
import time

from . import __version__


class _HelpFormatter(argparse.RawDescriptionHelpFormatter):
    def __init__(self, prog):
        super().__init__(prog, max_help_position=36, width=100)

    def _format_usage(self, usage, actions, groups, prefix):
        return super()._format_usage(usage, actions, groups, "Usage: " if prefix is None else prefix)

    def _format_action(self, action):
        if isinstance(action, argparse._SubParsersAction):
            return "".join(self._format_action(item) for item in action._get_subactions())
        return super()._format_action(action)


class _HelpParser(argparse.ArgumentParser):
    """Samtools-style banner, concise usage, and grouped option descriptions."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("formatter_class", _HelpFormatter)
        super().__init__(*args, **kwargs)
        self._optionals.title = "General options"

    def format_help(self):
        formatter = self._get_formatter()
        formatter.add_text(self.description)
        formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)
        groups = [group for group in self._action_groups if group is not self._optionals]
        for group in [*groups, self._optionals]:
            formatter.start_section(group.title)
            formatter.add_text(group.description)
            formatter.add_arguments(group._group_actions)
            formatter.end_section()
        formatter.add_text(self.epilog)
        return formatter.format_help()


def _boolean_argument(value: str) -> bool:
    if value.lower() in ("true", "t"):
        return True
    if value.lower() in ("false", "f"):
        return False
    raise argparse.ArgumentTypeError("expected True or False (T/F also accepted)")


def _modification_argument(value):
    names = value.split(",")
    if any(name not in ("5mC", "5hmC") for name in names) or len(names) != len(set(names)):
        raise argparse.ArgumentTypeError("choose 5mC, 5hmC, or 5mC,5hmC without duplicates; CpG context only")
    return value


def _run_selected_modifications(args):
    from .progress import print_analysis_summary, duration
    started = time.monotonic()
    modifications = args.modification.split(",")
    if len(modifications) == 1:
        manifest, files = _run_analysis(args)
        print_analysis_summary(args.output, manifest, files, time.monotonic() - started)
        return 0
    if os.path.lexists(args.output):
        raise FileExistsError(f"Output already exists: {args.output}. Choose a new run directory.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{args.output.name}.", suffix=".tmp", dir=args.output.parent))
    completed = []
    try:
        for modification in modifications:
            child = argparse.Namespace(**vars(args))
            child.modification, child.output = modification, stage / modification
            print(f"methyldiver: processing {modification} at CpGs ({len(completed) + 1}/{len(modifications)})", file=sys.stderr)
            child_started = time.monotonic()
            manifest, files = _run_analysis(child)
            if completed and manifest["inputs"] != completed[0][1]["inputs"]:
                raise ValueError("Input identities changed between modification analyses; no results were published.")
            completed.append((modification, manifest, files, time.monotonic() - child_started))
        summary = {"schema_version": "2.1", "package_version": __version__, "status": "complete",
                   "mode": "genome" if args.genome else "regions",
                   "analysis_type": "single" if args.command == "single" else "paired",
                   "modifications": modifications, "modification_context": "CpG",
                   "results": {name: {"directory": name, "analysis_summary": f"{name}/analysis_summary.json"}
                               for name in modifications}}
        (stage / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        if os.path.lexists(args.output):
            raise FileExistsError(f"Output appeared during analysis: {args.output}; it was not overwritten.")
        os.rename(stage, args.output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    for name, manifest, files, elapsed in completed:
        print_analysis_summary(args.output / name, manifest, files, elapsed)
    print(f"methyldiver: all modifications complete | elapsed {duration(time.monotonic() - started)}")
    print(f"  Output: {args.output.resolve()}")
    return 0


def _analysis_arguments(command, *, single=False):
    inputs = command.add_argument_group("Input files (required)")
    if single:
        inputs.add_argument("--input", required=True, type=Path, metavar="BAM", help="Indexed sample modBAM.")
    else:
        inputs.add_argument("--control", required=True, type=Path, metavar="BAM", help="Indexed control modBAM.")
        inputs.add_argument("--test", required=True, type=Path, metavar="BAM", help="Indexed test modBAM.")
    inputs.add_argument("--reference", required=True, type=Path, metavar="FASTA", help="Matching reference FASTA with .fai index.")

    scope = command.add_argument_group("Region selection (choose one)")
    selection = scope.add_mutually_exclusive_group(required=True)
    selection.add_argument("--genome", action="store_true", help="All BAM reference contigs; requires the [genome] extra [off].")
    selection.add_argument("--region", action="append", metavar="CHR:START-END", help="0-based, end-exclusive interval; repeat for separate regions [unset].")
    selection.add_argument("--regions-bed", type=Path, metavar="FILE", help="Nonoverlapping BED3/BED4 intervals; 0-based, end-exclusive [unset].")

    analysis = command.add_argument_group("Analysis options")
    analysis.add_argument("--modification", default="5mC", type=_modification_argument, metavar="STR",
                          help="CpG target: 5mC, 5hmC, or 5mC,5hmC for separate results [5mC].")
    analysis.add_argument("--N", dest="k", metavar="INT", default=3, type=int, help="Consecutive CpGs per unit, 1–10 [%(default)s].")
    analysis.add_argument("--call-threshold", type=float, default=0.8, metavar="FLOAT", help="Binary-call confidence, greater than 0.5 and at most 1 [%(default)s].")
    analysis.add_argument("--min-mapq", type=int, default=20, metavar="INT", help="Minimum read mapping quality, 0–254; MAPQ 255 is excluded [%(default)s].")
    analysis.add_argument("--min-base-quality", type=int, default=10, metavar="INT", help="Minimum CpG-call base quality, 0–254; missing qualities excluded [%(default)s].")
    depth_scope = "the sample" if single else "BOTH samples"
    analysis.add_argument("--min-depth", type=int, default=1, metavar="INT", help=f"Minimum complete-read depth in {depth_scope} to include a unit in the metric TSV [%(default)s].")

    output = command.add_argument_group("Output options")
    output.add_argument("--output", required=True, type=Path, metavar="DIR", help="New output directory (required; must not exist).")
    if single:
        output.add_argument("--sample-label", default="sample", metavar="STR", help="Sample label in two_metrics.tsv [%(default)s].")
    else:
        output.add_argument("--control-label", default="control", metavar="STR", help="Control label in three_metrics.tsv [%(default)s].")
        output.add_argument("--test-label", default="test", metavar="STR", help="Test label in three_metrics.tsv [%(default)s].")

    genome = command.add_argument_group("Processing options (genome or selected regions)")
    genome.add_argument("--threads", type=int, default=8, metavar="INT", help="Worker processes for genome, BED, or region chunks [%(default)s].")
    genome.add_argument("--chunk-units", type=int, metavar="INT", help="Requested units per chunk; reduced for large N; boundary units preserved [10000].")
    genome.add_argument("--chunk-manifests", type=_boolean_argument, metavar="BOOL", help="Save per-chunk metadata as .json.gz; False omits it [False].")
    genome.add_argument("--keep-patterns", type=_boolean_argument, metavar="BOOL", help="Keep complete_patterns/ and incomplete_patterns/; False removes them after TSV export [False].")

    limits = command.add_argument_group("Batch sizing and resource limits")
    limits.add_argument("--max-total-bases", type=int, metavar="INT", help="Selected-region span target per batch; larger requests split automatically; a single sparse unit may exceed it [100000].")
    limits.add_argument("--max-reads-per-region", type=int, default=100_000, metavar="INT", help="Fetched reads before filtering, per sample per region/chunk [%(default)s].")
    limits.add_argument("--max-pattern-evaluations", type=int, default=5_000_000, metavar="INT", help="Candidate read/unit patterns per sample per region/chunk [%(default)s].")


def _reader_arguments(reader):
    files = reader.add_argument_group("Input/output options")
    files.add_argument("--input", required=True, type=Path, metavar="PATH", help="complete_patterns/, incomplete_patterns/, or one .parquet file (required).")
    files.add_argument("--output", type=Path, metavar="FILE", help="New output file; .csv selects CSV [stdout].")
    files.add_argument("--format", choices=("csv", "tsv"), metavar="STR", help="csv or tsv [.csv output: csv; otherwise: tsv].")
    filters = reader.add_argument_group("Selection options")
    rows = filters.add_mutually_exclusive_group()
    rows.add_argument("--limit", type=int, default=20, metavar="INT", help="Maximum matching rows to display/export [%(default)s].")
    rows.add_argument("--all", action="store_true", help="Export all matching rows, streamed in batches [off; limit: 20].")
    filters.add_argument("--chrom", metavar="STR", help="Exact chromosome/contig name [all].")
    filters.add_argument("--sample", choices=("control", "test", "sample"), metavar="STR", help="Incomplete-pattern role: control, test, or sample [all].")
    filters.add_argument("--pattern", metavar="STR", help="Exact incomplete pattern, e.g. '1.0'; dot means missing [all].")
    filters.add_argument("--metrics-only", action="store_true", help="Plotting columns only (7 single, 12 paired); require complete metrics and the stored depth flag [off: all columns and rows].")
    labels = reader.add_argument_group("Sample labels (with --metrics-only)")
    labels.add_argument("--control-label", metavar="STR", help="Paired control label [saved label, or control].")
    labels.add_argument("--test-label", metavar="STR", help="Paired test label [saved label, or test].")
    labels.add_argument("--sample-label", metavar="STR", help="Single-sample label [saved label, or sample].")


def _discovery_arguments(command):
    files = command.add_argument_group("Input/output options")
    files.add_argument("--input", required=True, type=Path, metavar="PATH",
                       help="Paired three_metrics.tsv or its result directory; grouped by chromosome and sorted by start.")
    files.add_argument("--output", required=True, type=Path, metavar="DIR", help="New discovery output directory (required).")
    settings = command.add_argument_group("Region ranking options")
    definition = settings.add_mutually_exclusive_group()
    definition.add_argument("--window-size", type=int, default=5000, metavar="INT", help="Nonoverlapping genomic window size in bp; only fully contained CpG units are used [%(default)s].")
    definition.add_argument("--regions-bed", type=Path, metavar="BED", help="Rank nonoverlapping BED3/4 regions instead of fixed windows; only fully contained units are used [none].")
    settings.add_argument("--min-units", type=int, default=30, metavar="INT", help="Minimum saved CpG units per candidate window; not read depth [%(default)s].")
    settings.add_argument("--min-ed", type=float, default=0.5, metavar="FLOAT", help="Minimum median Epiallele distance, inclusive, 0-1 [%(default)s].")
    settings.add_argument("--top", type=int, default=100, metavar="INT", help="Maximum ranked regions exported to TSV and BED [%(default)s].")


def _run_analysis(args):
    single = args.command == "single"
    control, test = (args.input, None) if single else (args.control, args.test)
    control_label = "control" if single else args.control_label
    test_label = "test" if single else args.test_label
    sample_label = args.sample_label if single else "sample"
    from .progress import LiveProgress
    from .output import sample_labels, single_sample_labels, metric_filename
    if single:
        single_sample_labels(sample_label)
    else:
        sample_labels(control_label, test_label)
    summary_name = metric_filename(single)
    if os.path.lexists(args.output):
        raise FileExistsError(f"Output already exists: {args.output}. Choose a new run directory.")
    regions = None
    chunked_regions = False
    if not args.genome:
        from .modbam import parse_region, read_regions_bed
        from .coordinates import MAX_DENSE_COUNT_CELLS, _unit_size
        _unit_size(args.k)
        regions = read_regions_bed(args.regions_bed) if args.regions_bed else [parse_region(value) for value in args.region]
        span_limit = args.max_total_bases if args.max_total_bases is not None else 100_000
        # Conservative bound: no reference can have more than one CpG per
        # two bases. The streaming planner then counts actual reference units.
        count_bound = sum(max(0, (r.end - r.start) // 2 - args.k + 1) for r in regions)
        chunked_regions = (args.threads != 1 or args.chunk_units is not None
                           or args.chunk_manifests is not None or args.keep_patterns is not None
                           or sum(r.end - r.start for r in regions) > span_limit
                           or count_bound * (1 if single else 2) * 2 ** args.k > MAX_DENSE_COUNT_CELLS)
    if args.genome:
        if args.max_total_bases is not None:
            raise ValueError("--max-total-bases applies only to selected regions, not --genome.")
    if args.genome or chunked_regions:
        from .genome import analyze_genome
        runner = analyze_genome
        inputs = (control, test, args.reference, args.output)
        extra = {}
        if chunked_regions:
            from .genome import analyze_regions
            runner = analyze_regions
            inputs = (control, test, args.reference, regions, args.output)
            extra = {"max_total_bases": span_limit, "regions_bed": args.regions_bed}
        with LiveProgress() as reporter:
            manifest = runner(*inputs,
                threads=args.threads, chunk_units=args.chunk_units if args.chunk_units is not None else 10_000,
                chunk_manifests=args.chunk_manifests if args.chunk_manifests is not None else False,
                keep_patterns=args.keep_patterns if args.keep_patterns is not None else False,
                control_label=control_label, test_label=test_label, sample_label=sample_label,
                k=args.k, modification=args.modification, call_threshold=args.call_threshold, min_mapq=args.min_mapq,
                min_base_quality=args.min_base_quality, min_depth=args.min_depth,
                max_reads_per_region=args.max_reads_per_region,
                max_pattern_evaluations=args.max_pattern_evaluations,
                progress=reporter.analysis, status=reporter.stage, **extra)
        files = ["complete_patterns/", "incomplete_patterns/"] if manifest["keep_patterns"] else []
        if manifest["chunk_manifests"]:
            files.append("chunk_manifests/")
        files += [summary_name, "chunks.tsv", "analysis_summary.json"]
        return manifest, files
    from .modbam import analyze_modbam
    from .output import write_analysis

    result = analyze_modbam(
        control, test, args.reference, regions, k=args.k, modification=args.modification,
        call_threshold=args.call_threshold, min_mapq=args.min_mapq,
        min_base_quality=args.min_base_quality, min_depth=args.min_depth,
        max_total_bases=args.max_total_bases if args.max_total_bases is not None else 100_000,
        max_reads_per_region=args.max_reads_per_region,
        max_pattern_evaluations=args.max_pattern_evaluations,
    )
    if args.regions_bed:
        result["manifest"]["regions_bed_sha256"] = hashlib.sha256(args.regions_bed.read_bytes()).hexdigest()
    manifest = write_analysis(result, args.output, control_label=control_label, test_label=test_label, sample_label=sample_label)
    return manifest, [summary_name, "metrics.tsv", "patterns.tsv", "analysis_summary.json"]


def main(argv: list[str] | None = None) -> int:
    parser = _HelpParser(
        prog="methyldiver", usage="%(prog)s <command> [options]",
        description=f"Program: methyldiver (Methylation and epiallele diversity from long reads)\nVersion: {__version__}",
        epilog="Run 'methyldiver <command> --help' for command-specific options.\nValues in brackets show defaults.")
    parser.add_argument("--version", action="version", version=f"methyldiver {__version__}")
    subcommands = parser.add_subparsers(dest="command", required=True, title="Commands", prog=parser.prog)
    common_note = "Defaults are shown in brackets. BOOL accepts True/False or T/F.\nOnly complete read patterns contribute to metric estimates."
    compare = subcommands.add_parser("compare", allow_abbrev=False,
        usage="%(prog)s [options] --control BAM --test BAM --reference FASTA --output DIR",
        description="Compare methylation level, entropy, and Epiallele distance between two samples.",
        help="Compare methylation, entropy, and Epiallele distance between two samples.", epilog=common_note)
    _analysis_arguments(compare)
    single = subcommands.add_parser("single", allow_abbrev=False,
        usage="%(prog)s [options] --input BAM --reference FASTA --output DIR",
        description="Calculate methylation level and entropy for one sample.",
        help="Calculate methylation level and entropy for one sample.", epilog=common_note)
    _analysis_arguments(single, single=True)
    reader = subcommands.add_parser("read",
        usage="%(prog)s [options] --input PATH",
        description="View or export retained MethylDiver pattern tables.",
        help="View or export retained pattern tables as TSV or CSV.",
        epilog="Defaults are shown in brackets. Parquet input requires the [genome] extra.")
    _reader_arguments(reader)
    discover = subcommands.add_parser("discover", allow_abbrev=False,
        usage="%(prog)s --input PATH --output DIR [options]",
        description="Find high-ED genomic windows or BED regions from paired CpG-unit results.",
        help="Rank genomic windows or BED regions by median Epiallele distance.",
        epilog="Descriptive effect-size ranking, not a significance test. No additional read-depth filter is applied.\nDefaults are shown in brackets.")
    _discovery_arguments(discover)
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        parser.print_help()
        return 0
    if len(argv) == 1 and argv[0] in subcommands.choices:
        subcommands.choices[argv[0]].print_help()
        return 0
    args = parser.parse_args(argv)
    try:
        if args.command == "discover":
            from .discover import discover_regions
            from .progress import LiveProgress, duration
            started = time.monotonic()
            with LiveProgress() as reporter:
                summary = discover_regions(args.input, args.output, window_size=args.window_size,
                    min_units=args.min_units, min_ed=args.min_ed, top=args.top, regions_bed=args.regions_bed,
                    progress=lambda done, total, rows: reporter.stage("scan TSV", done, total))
            print(f"methyldiver: discovery complete | elapsed {duration(time.monotonic() - started)}")
            print(f"  Regions evaluated: {summary['windows_evaluated']:,} | passing filters: {summary['windows_passing_filters']:,}")
            print(f"  Ranked regions: {summary['n_ranked_regions']:,}")
            print("  Ranking: median ED | additional read-depth filter: none | descriptive, not a significance test")
            print(f"  Output: {args.output.resolve()}")
            print("  Files: ranked_regions.tsv, top_regions.bed, window_statistics.tsv, discovery_summary.json")
            return 0
        if args.command == "read":
            from .reader import read_results, saved_sample_labels
            if not args.metrics_only and any(value is not None for value in (args.control_label, args.test_label, args.sample_label)):
                raise ValueError("Sample label options require --metrics-only when reading results.")
            labels = saved_sample_labels(args.input) if args.metrics_only else {}
            count = read_results(args.input, output=args.output, limit=None if args.all else args.limit,
                                 chrom=args.chrom, sample=args.sample, pattern=args.pattern,
                                 output_format=args.format, metrics_only=args.metrics_only,
                                 control_label=args.control_label if args.control_label is not None else labels.get("control", "control"),
                                 test_label=args.test_label if args.test_label is not None else labels.get("test", "test"),
                                 sample_label=args.sample_label if args.sample_label is not None else labels.get("sample", "sample"))
            print(f"methyldiver: {count:,} rows written" +
                  (f" to {args.output}" if args.output else " to stdout") +
                  (" (all matching rows)" if args.all else f" (limit {args.limit}; use --all for all matching rows)"),
                  file=sys.stderr)
            return 0
        if args.command in ("compare", "single"):
            return _run_selected_modifications(args)
    except (OSError, ValueError, ImportError) as exc:
        parser.exit(2, f"methyldiver: error: {exc}\n")
    except KeyboardInterrupt:
        parser.exit(130, "methyldiver: interrupted; no completed output was published.\n")
    return 0
