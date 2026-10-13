"""Stream MethylDiver Parquet results as human-readable TSV or CSV."""

import csv
import json
import os
from pathlib import Path
import sys
import tempfile

from .output import (format_cell, labelled_columns, sample_labels, single_sample_labels,
                     metric_columns)


def saved_sample_labels(source):
    """Recover plotting labels from a run summary beside a Parquet dataset."""
    source = Path(source).resolve()
    directory = source if source.is_dir() else source.parent
    for folder in (directory, directory.parent):
        for name in ("analysis_summary.json", "run.json"):
            path = folder / name
            if path.is_file():
                labels = json.loads(path.read_text(encoding="utf-8")).get("sample_labels", {})
                if "sample" in labels:
                    return single_sample_labels(labels["sample"])
                return sample_labels(labels.get("control", "control"), labels.get("test", "test"))
    return sample_labels()


def read_results(source, *, output=None, limit=20, chrom=None, sample=None, pattern=None,
                 output_format=None, columns=None, drop_na_metrics=False,
                 control_label="control", test_label="test", sample_label="sample",
                 metrics_only=False, progress=None):
    """Preview/export a Parquet file or dataset; return the number of rows written.

    Output defaults to stdout. ``limit=None`` exports all matching rows in
    batches. A named output must be new and is published only on success.
    CpG anchor lists become comma-separated positions; missing values are NA.
    Format defaults to CSV for a .csv destination, otherwise TSV. ``columns``
    optionally selects fields in output order.
    ``metrics_only`` selects plotting columns for the detected sample mode and
    drops missing metrics. ``drop_na_metrics`` also drops missing metrics with
    explicit columns: two values for single-sample tables, seven for paired.
    Sample labels rename beta/entropy headings only; stored roles and values
    remain unchanged. Label discovery is handled by the CLI.
    """
    if limit is not None and (isinstance(limit, bool) or not isinstance(limit, int) or limit < 1):
        raise ValueError("limit must be an integer >= 1, or None for all rows.")
    source = Path(source).resolve()
    if not source.exists():
        raise ValueError(f"Input does not exist: {source}")
    destination = Path(output) if output is not None else None
    if output_format is None:
        output_format = "csv" if destination is not None and destination.suffix.lower() == ".csv" else "tsv"
    if output_format not in ("csv", "tsv"):
        raise ValueError("output_format must be csv or tsv.")
    if destination is not None:
        if os.path.lexists(destination):
            raise FileExistsError(f"Output already exists: {destination}. Choose a new output file.")
        if source.is_dir() and destination.resolve().is_relative_to(source):
            raise ValueError("Write the export outside the input Parquet directory.")
    try:
        import pyarrow.dataset as ds
    except ImportError as exc:
        raise ImportError("Reading Parquet requires pyarrow; install MethylDiver with the [genome] extra.") from exc
    if source.is_dir():
        paths = sorted(source.glob("*.parquet"))
        if not paths:
            raise ValueError("No Parquet files found. Supply complete_patterns/, incomplete_patterns/, or one .parquet file.")
    else:
        if source.suffix != ".parquet":
            raise ValueError("Input must be a .parquet file or a Parquet dataset directory.")
        paths = [source]
    dataset = ds.dataset([str(path) for path in paths], format="parquet")
    # Avoid accidentally reading unrelated Parquet data with misleading filters.
    if not {"region_id", "chrom", "cpg_anchors"}.issubset(dataset.schema.names):
        raise ValueError("Input is not a MethylDiver result table (missing region_id, chrom, or cpg_anchors).")
    single = "sample_beta" in dataset.schema.names
    labels = single_sample_labels(sample_label) if single else sample_labels(control_label, test_label)
    if metrics_only:
        columns = metric_columns(single)
        drop_na_metrics = True
    expression = None
    for name, value in (("chrom", chrom), ("sample", sample), ("pattern", pattern)):
        if value is None:
            continue
        if name not in dataset.schema.names:
            raise ValueError(f"--{name} requires a '{name}' column; sample/pattern filters apply to incomplete_patterns/.")
        predicate = ds.field(name) == value
        expression = predicate if expression is None else expression & predicate
    if columns is None:
        columns = dataset.schema.names
    else:
        columns = list(columns)
        missing = [name for name in columns if name not in dataset.schema.names]
        if missing:
            raise ValueError(f"Input is missing requested metric columns: {', '.join(missing)}")
    if drop_na_metrics:
        for name in metric_columns(single)[5:]:
            if name not in dataset.schema.names:
                raise ValueError(f"Input is missing required metric column: {name}")
            predicate = ds.field(name).is_valid()
            expression = predicate if expression is None else expression & predicate
    scanner = dataset.scanner(columns=columns, filter=expression, batch_size=4096, use_threads=False,
                              batch_readahead=1, fragment_readahead=1)

    def write(handle):
        writer = csv.writer(handle, delimiter="," if output_format == "csv" else "\t", lineterminator="\n")
        writer.writerow(labelled_columns(columns, labels))
        count = 0
        for batch in scanner.to_batches():
            if limit is not None:
                batch = batch.slice(0, min(batch.num_rows, limit - count))
            for row in batch.to_pylist():
                writer.writerow([format_cell(row[name]) for name in columns])
            count += batch.num_rows
            if progress is not None:
                progress(count)
            if limit is not None and count >= limit:
                break
        return count

    if destination is None:
        return write(sys.stdout)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            count = write(handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, destination)
        return count
    finally:
        os.unlink(temporary)
