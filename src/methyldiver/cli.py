"""ModBAM analysis and baseline count-table CLI, with explicit input paths."""

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import tempfile

from . import __version__
from .metrics import compare_counts, metrics


def read_count_table(path: Path) -> tuple[list[int], list[int], int, str]:
    """Read exactly pattern/control/test TSV columns for one matched unit."""
    raw = path.read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    if reader.fieldnames != ["pattern", "control", "test"]:
        raise ValueError("Expected TSV header: pattern\tcontrol\ttest. bedMethyl cannot supply epiallele patterns.")
    rows = list(reader)
    if len(rows) not in (8, 16):
        raise ValueError("Supply all 8 or 16 binary patterns, including zero-count states.")
    k = len(rows).bit_length() - 1
    expected = [format(i, f"0{k}b") for i in range(2 ** k)]
    if [row["pattern"] for row in rows] != expected:
        raise ValueError("Patterns must appear exactly once in ascending binary order.")
    control, test = [], []
    for line, row in enumerate(rows, start=2):
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Line {line}: expected exactly three tab-separated fields.")
        for column, destination in (("control", control), ("test", test)):
            text = row[column].strip()
            if not re.fullmatch(r"[0-9]+", text):
                raise ValueError(f"Line {line}: {column} count must be a nonnegative whole number.")
            destination.append(int(text))
    return control, test, k, hashlib.sha256(raw).hexdigest()


def _write_json(path: Path, text: str, overwrite: bool) -> None:
    """Publish a complete JSON file without exposing a partially written result."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        if overwrite:
            os.replace(temporary, path)
        else:
            # Same-filesystem link atomically refuses an existing destination.
            os.link(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="MethylDiver: complete-pattern 3/4-CpG metrics.")
    parser.add_argument("--version", action="version", version=f"methyldiver {__version__}")
    subcommands = parser.add_subparsers(dest="command", required=True)
    analyze = subcommands.add_parser("analyze", help="Extract 3/4-CpG patterns from control/test modBAMs and compare metrics.")
    analyze.add_argument("--control", required=True, type=Path, help="Indexed control modBAM.")
    analyze.add_argument("--test", required=True, type=Path, help="Indexed test modBAM.")
    analyze.add_argument("--reference", required=True, type=Path, help="Matching reference FASTA with .fai index.")
    selection = analyze.add_mutually_exclusive_group(required=True)
    selection.add_argument("--region", action="append", help="CHROM:START-END, 0-based and end-exclusive; repeat for separate regions.")
    selection.add_argument("--regions-bed", type=Path, help="Nonoverlapping BED3/BED4 intervals, 0-based and end-exclusive.")
    analyze.add_argument("--output", required=True, type=Path, help="NEW output directory for metrics.tsv, patterns.tsv, run.json.")
    analyze.add_argument("--k", choices=(3, 4), default=3, type=int, help="Consecutive CpGs per unit (default: 3).")
    analyze.add_argument("--call-threshold", type=float, default=0.8, help="Minimum binary call confidence (default: 0.8).")
    analyze.add_argument("--min-mapq", type=int, default=20)
    analyze.add_argument("--min-base-quality", type=int, default=10)
    analyze.add_argument("--min-depth", type=int, default=1, help="Flag units with this many complete reads per sample; all estimates retained.")
    analyze.add_argument("--max-total-bases", type=int, default=100_000, help="Maximum total requested interval span (default: 100000).")
    analyze.add_argument("--max-reads-per-region", type=int, default=100_000)
    analyze.add_argument("--max-pattern-evaluations", type=int, default=5_000_000, help="Maximum candidate read/unit patterns per sample per region.")
    command = subcommands.add_parser("metrics-counts", help="Compare one unit's control and test pattern counts.")
    command.add_argument("--input", required=True, type=Path, help="TSV with pattern, control, test columns.")
    command.add_argument("--output", type=Path, help="JSON output path; defaults to stdout only.")
    command.add_argument("--overwrite", action="store_true", help="Explicitly replace an existing output file.")
    args = parser.parse_args(argv)
    try:
        if args.command == "analyze":
            if os.path.lexists(args.output):
                raise FileExistsError(f"Output already exists: {args.output}. Choose a new run directory.")
            from .modbam import analyze_modbam, parse_region, read_regions_bed
            from .output import write_analysis

            regions = read_regions_bed(args.regions_bed) if args.regions_bed else [parse_region(value) for value in args.region]
            result = analyze_modbam(
                args.control, args.test, args.reference, regions, k=args.k,
                call_threshold=args.call_threshold, min_mapq=args.min_mapq,
                min_base_quality=args.min_base_quality, min_depth=args.min_depth,
                max_total_bases=args.max_total_bases, max_reads_per_region=args.max_reads_per_region,
                max_pattern_evaluations=args.max_pattern_evaluations,
            )
            if args.regions_bed:
                result["manifest"]["regions_bed_sha256"] = hashlib.sha256(args.regions_bed.read_bytes()).hexdigest()
            write_analysis(result, args.output)
            print(json.dumps({"output_directory": str(args.output.resolve()),
                              "n_units": result["manifest"]["n_units"],
                              "n_units_passing_min_depth": result["manifest"]["n_units_passing_min_depth"],
                              "files": ["metrics.tsv", "patterns.tsv", "run.json"]}, indent=2))
            return 0
        if args.overwrite and args.output is None:
            raise ValueError("--overwrite requires --output.")
        if args.output is not None:
            if args.input.resolve() == args.output.resolve() or (
                args.output.exists() and args.input.samefile(args.output)
            ):
                raise ValueError("Output must not replace the input count table.")
        control, test, k, digest = read_count_table(args.input)
        result = {
            "schema_version": "1.0",
            "package_version": __version__,
            "unit_size": k,
            "control": metrics(control),
            "test": metrics(test),
            "comparison": compare_counts(control, test),
            "provenance": {
                "input_sha256": digest,
                "method": "complete_read_plugin",
                "entropy_log_base": 2,
                "entropy_normalization": "divide_by_unit_size",
                "distance": "sqrt_base2_jensen_shannon_divergence",
                "direction": "test_minus_control",
                "minimum_depth_filter": None,
            },
        }
        text = json.dumps(result, indent=2, allow_nan=False) + "\n"
        if args.output is not None:
            _write_json(args.output, text, args.overwrite)
        print(text, end="")
    except (OSError, ValueError, ImportError) as exc:
        parser.exit(2, f"methyldiver: error: {exc}\n")
    return 0
