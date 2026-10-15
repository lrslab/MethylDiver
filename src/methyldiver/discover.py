"""Descriptive ED ranking in fixed windows or BED regions from paired TSVs."""

from bisect import bisect_right
import csv
from datetime import datetime, timezone
import hashlib
import heapq
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import tempfile

from . import __version__


WINDOW_COLUMNS = ["window_id", "chrom", "start", "end", "n_units", "median_ed", "mean_ed",
                  "fraction_units_ge_threshold", "mean_delta_beta", "mean_delta_entropy", "passes_filters"]


def _bed_regions(path, lengths):
    """Read nonoverlapping BED3/4 intervals without importing BAM dependencies."""
    regions, names = {}, set()
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            row = line.rstrip("\r\n").split("\t")
            if len(row) not in (3, 4):
                raise ValueError(f"BED line {number}: expected 3 or 4 tab-separated columns.")
            chrom = row[0]
            try:
                start, end = int(row[1]), int(row[2])
            except ValueError as exc:
                raise ValueError(f"BED line {number}: invalid coordinates.") from exc
            if not chrom or any(c.isspace() for c in chrom) or start < 0 or end <= start:
                raise ValueError(f"BED line {number}: invalid BED0 interval.")
            if chrom in lengths and end > lengths[chrom]:
                raise ValueError(f"BED line {number}: interval exceeds the saved contig length.")
            name = row[3] if len(row) == 4 and row[3] else f"{chrom}:{start}-{end}"
            if name in names:
                raise ValueError(f"BED line {number}: duplicate region ID {name!r}.")
            names.add(name)
            regions.setdefault(chrom, []).append((start, end, name))
    if not regions:
        raise ValueError("The regions BED contains no intervals.")
    for chrom, intervals in regions.items():
        intervals.sort()
        if any(a[1] > b[0] for a, b in zip(intervals, intervals[1:])):
            raise ValueError(f"BED regions overlap on {chrom}; merge overlapping intervals before discovery.")
    return regions


def discover_regions(source, output, *, window_size=5000, min_units=30, min_ed=0.5, top=100,
                     regions_bed=None, progress=None):
    """Rank fully contained units in fixed windows or nonoverlapping BED regions.

    This is effect-size screening, not a significance test. No additional
    read-depth filtering is possible from three_metrics.tsv. Memory holds one window's
    observations and at most `top` ranked windows, rather than the full table.
    `progress(bytes_read, file_bytes, rows)` is called periodically.
    """
    for name, value in (("window_size", window_size), ("min_units", min_units), ("top", top)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be an integer >= 1.")
    if isinstance(min_ed, bool) or not isinstance(min_ed, (int, float)) or not math.isfinite(min_ed) or not 0 <= min_ed <= 1:
        raise ValueError("min_ed must be a finite number from 0 to 1.")
    source, output = Path(source).resolve(), Path(output).resolve()
    if source.is_dir():
        source = source / "three_metrics.tsv"
    if os.path.lexists(output):
        raise FileExistsError(f"Output already exists: {output}. Choose a new directory.")
    before = source.stat()
    # Contig lengths are optional, but allow clipping the final window when
    # the original MethylDiver summary is available beside the input table.
    summary_path = source.parent / "analysis_summary.json"
    source_summary = json.loads(summary_path.read_text()) if summary_path.is_file() else {}
    lengths = {r["chrom"]: r["length"] for r in source_summary.get("contigs", [])}
    summary_hash = hashlib.sha256(summary_path.read_bytes()).hexdigest() if summary_path.is_file() else None
    bed_path = Path(regions_bed).resolve() if regions_bed is not None else None
    bed_hash = hashlib.sha256(bed_path.read_bytes()).hexdigest() if bed_path else None
    bed = _bed_regions(bed_path, lengths) if bed_path else None
    bed_starts = {chrom: [r[0] for r in rows] for chrom, rows in bed.items()} if bed else {}
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent))
    counts = {"input_rows": 0, "missing_metric_rows": 0, "boundary_crossing_rows": 0,
              "outside_regions_rows": 0, "included_units": 0,
              "windows_evaluated": 0, "windows_passing_filters": 0}
    heap = []
    started = datetime.now(timezone.utc).isoformat()
    try:
        with source.open(encoding="utf-8", newline="") as handle, (stage / "window_statistics.tsv").open("w", newline="") as table:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
            required = ("chrom", "start", "end", "epiallele_distance", "delta_beta", "delta_entropy")
            if len(set(header)) != len(header) or any(name not in header for name in required):
                raise ValueError("Input must be a paired three_metrics.tsv with coordinates, epiallele_distance, delta_beta, and delta_entropy.")
            ix = {name: header.index(name) for name in required}
            writer = csv.DictWriter(table, fieldnames=WINDOW_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            current = None
            observations = []
            previous_chrom, previous_start = None, -1
            closed_contigs = set()

            def finish_window():
                if not observations:
                    return
                chrom, start, end, region_id = current
                n = len(observations)
                values = [r[0] for r in observations]
                median = statistics.median(values)
                passed = n >= min_units and median >= min_ed
                result = {"window_id": region_id, "chrom": chrom, "start": start, "end": end,
                          "n_units": n, "median_ed": median, "mean_ed": math.fsum(values) / n,
                          "fraction_units_ge_threshold": sum(v >= min_ed for v in values) / n,
                          "mean_delta_beta": math.fsum(r[1] for r in observations) / n,
                          "mean_delta_entropy": math.fsum(r[2] for r in observations) / n,
                          "passes_filters": int(passed)}
                writer.writerow(result)
                counts["windows_evaluated"] += 1
                if passed:
                    counts["windows_passing_filters"] += 1
                    score = (median, result["mean_ed"], n, -counts["windows_evaluated"])
                    entry = (score, result)
                    if len(heap) < top:
                        heapq.heappush(heap, entry)
                    elif score > heap[0][0]:
                        heapq.heapreplace(heap, entry)

            if progress:
                progress(0, before.st_size, 0)
            for number, row in enumerate(reader, 2):
                counts["input_rows"] += 1
                if len(row) != len(header):
                    raise ValueError(f"TSV line {number}: column count does not match the header.")
                chrom = row[ix["chrom"]]
                try:
                    start, end = int(row[ix["start"]]), int(row[ix["end"]])
                except ValueError as exc:
                    raise ValueError(f"TSV line {number}: invalid coordinates.") from exc
                if not chrom or any(c.isspace() for c in chrom) or start < 0 or end <= start:
                    raise ValueError(f"TSV line {number}: invalid BED0 interval.")
                if chrom in lengths and end > lengths[chrom]:
                    raise ValueError(f"TSV line {number}: interval exceeds the saved contig length.")
                if chrom != previous_chrom:
                    if chrom in closed_contigs:
                        raise ValueError("Input must be grouped by chromosome and sorted by increasing start; a contig reappeared.")
                    if previous_chrom is not None:
                        closed_contigs.add(previous_chrom)
                    previous_chrom, previous_start = chrom, -1
                if start <= previous_start:
                    raise ValueError(f"TSV line {number}: duplicate unit start or unsorted input; provide one sorted result table.")
                previous_start = start
                texts = [row[ix[name]] for name in ("epiallele_distance", "delta_beta", "delta_entropy")]
                if any(t in ("", "NA", "NaN", "nan", "null") for t in texts):
                    counts["missing_metric_rows"] += 1
                    continue
                try:
                    values = tuple(float(t) for t in texts)
                except ValueError as exc:
                    raise ValueError(f"TSV line {number}: invalid metric value.") from exc
                if (not all(math.isfinite(v) for v in values) or not 0 <= values[0] <= 1
                        or any(not -1 <= v <= 1 for v in values[1:])):
                    raise ValueError(f"TSV line {number}: metric outside its valid range.")
                if bed is not None:
                    index = bisect_right(bed_starts.get(chrom, []), start) - 1
                    if index < 0 or start >= bed[chrom][index][1]:
                        counts["outside_regions_rows"] += 1
                        continue
                    lo, hi, name = bed[chrom][index]
                    if end > hi:
                        counts["boundary_crossing_rows"] += 1
                        continue
                    key = (chrom, lo, hi, name)
                else:
                    window_start = start // window_size * window_size
                    if end > window_start + window_size:
                        counts["boundary_crossing_rows"] += 1
                        continue
                    hi = min(window_start + window_size, lengths.get(chrom, window_start + window_size))
                    key = (chrom, window_start, hi, f"{chrom}:{window_start}-{hi}")
                if key != current:
                    finish_window()
                    observations = []
                    current = key
                observations.append(values)
                counts["included_units"] += 1
                if progress and counts["input_rows"] % 100_000 == 0:
                    progress(min(handle.buffer.tell(), before.st_size), before.st_size, counts["input_rows"])
            finish_window()
            if progress:
                progress(before.st_size, before.st_size, counts["input_rows"])
        after = source.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError("Input TSV changed during discovery.")
        if summary_hash is not None and hashlib.sha256(summary_path.read_bytes()).hexdigest() != summary_hash:
            raise ValueError("Source analysis summary changed during discovery.")
        if bed_path and hashlib.sha256(bed_path.read_bytes()).hexdigest() != bed_hash:
            raise ValueError("Regions BED changed during discovery.")
        ranked = [entry[1] for entry in sorted(heap, key=lambda e: e[0], reverse=True)]
        with (stage / "ranked_regions.tsv").open("w", newline="") as handle, (stage / "top_regions.bed").open("w") as bed_output:
            writer = csv.DictWriter(handle, fieldnames=["rank", *WINDOW_COLUMNS], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for rank, row in enumerate(ranked, 1):
                writer.writerow({"rank": rank, **row})
                bed_output.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['window_id']}\n")
        summary = {"package_version": __version__, "status": "complete", "analysis": "descriptive_ED_region_ranking" if bed else "descriptive_ED_window_ranking",
            "started_at": started, "completed_at": datetime.now(timezone.utc).isoformat(),
            "input": {"path": str(source), "bytes": before.st_size, "mtime_ns": before.st_mtime_ns},
            "source_summary_sha256": summary_hash,
            "sample_labels": source_summary.get("sample_labels"),
            "source_config": source_summary.get("config"),
            "source_metric_filter": source_summary.get("outputs", {}).get(source.name),
            "parameters": {"window_size": None if bed else window_size, "min_units": min_units, "min_ed": min_ed, "top": top,
                           "regions_bed": str(bed_path) if bed_path else None},
            "regions_bed_sha256": bed_hash,
            "regions_requested": sum(map(len, bed.values())) if bed else None,
            "regions_without_usable_units": sum(map(len, bed.values())) - counts["windows_evaluated"] if bed else None,
            "coordinate_system": "BED0_half_open", "window_policy": "supplied nonoverlapping BED intervals; full unit containment" if bed else "fixed nonoverlapping windows from coordinate zero; full unit containment",
            "terminal_windows": "exact BED boundaries retained" if bed else "clipped to saved contig lengths when available; otherwise full grid windows",
            "rank_order": "median ED descending, mean ED descending, unit count descending, input window order",
            "depth_filter": "none added; uses the input TSV's existing filtering; complete-read depths are unavailable",
            "interpretation": "effect-size screening only; overlapping CpG units are not independent replicates; no p-values or FDR",
            **counts, "n_ranked_regions": len(ranked),
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (stage / "discovery_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        if os.path.lexists(output):
            raise FileExistsError(f"Output appeared during discovery: {output}; it was not overwritten.")
        os.rename(stage, output)
        return summary
    finally:
        if stage.exists():
            shutil.rmtree(stage)
