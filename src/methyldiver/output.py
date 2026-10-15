"""Write complete modBAM run artifacts into a new output directory."""

import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile


METRIC_SUMMARY_COLUMNS = [
    "region_id", "chrom", "start", "end", "cpg_anchors",
    "control_beta", "test_beta", "delta_beta",
    "control_entropy", "test_entropy", "delta_entropy",
    "epiallele_distance",
]
METRIC_VALUE_COLUMNS = METRIC_SUMMARY_COLUMNS[5:]
SINGLE_METRIC_SUMMARY_COLUMNS = METRIC_SUMMARY_COLUMNS[:5] + ["sample_beta", "sample_entropy"]


def metric_columns(single=False):
    return SINGLE_METRIC_SUMMARY_COLUMNS if single else METRIC_SUMMARY_COLUMNS


def metric_filename(single=False):
    return "two_metrics.tsv" if single else "three_metrics.tsv"


def single_sample_labels(label="sample"):
    if not isinstance(label, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", label):
        raise ValueError("Sample label must start with a letter or digit and contain only letters, digits, _, ., or -.")
    return {"sample": label}


def sample_labels(control_label="control", test_label="test"):
    """Validate display labels without changing the underlying sample roles."""
    labels = {"control": control_label, "test": test_label}
    for role, label in labels.items():
        if not isinstance(label, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", label):
            raise ValueError(f"{role} label must start with a letter or digit and contain only letters, digits, _, ., or -.")
    if control_label.casefold() == test_label.casefold():
        raise ValueError("Control and test labels must be different.")
    headings = labelled_columns(METRIC_SUMMARY_COLUMNS, labels)
    if len(set(name.casefold() for name in headings)) != len(headings):
        raise ValueError("Sample labels conflict with metric column names; do not use 'delta'.")
    return labels


def labelled_columns(columns, labels):
    renamed = {f"{role}_{metric}": f"{label}_{metric}"
               for role, label in labels.items() for metric in ("beta", "entropy")}
    return [renamed.get(name, name) for name in columns]


def has_complete_metrics(row, single=False):
    """Whether all plotted metrics are available (two single, seven paired)."""
    return all(row[name] is not None for name in metric_columns(single)[5:])


DEPTH_POLICY = "filter_metric_tsv_preserve_full_counts"


def passes_metric_filter(row, single=False):
    """Compact metric exports require complete metrics and sufficient depth."""
    return has_complete_metrics(row, single) and bool(row["passes_min_depth"])


def metric_filter_summary(total, complete, exported, min_depth, single=False):
    """Report disjoint exclusions: missing metrics, then valid low-depth rows."""
    return {"rows": exported, "excluded_missing_metrics": total - complete,
            "excluded_below_min_depth": complete - exported,
            "min_depth": min_depth,
            "depth_scope": "sample" if single else "both_samples",
            "row_filter": "complete_metrics_and_passes_min_depth"}


def format_cell(value):
    """Keep missing values and CpG lists explicit in text result tables."""
    if value is None:
        return "NA"
    if isinstance(value, (list, tuple)):
        return ",".join(map(str, value))
    if isinstance(value, bool):
        return int(value)
    return value


def write_analysis(result, destination, *, control_label="control", test_label="test", sample_label="sample"):
    """Publish metric/pattern tables and a run summary after successful analysis.

    Existing directories are never overwritten. Temporary output is removed if
    writing fails. TSV missing values are NA; manifest JSON uses null.
    Labels customize plotting-table headings and are recorded in the summary;
    full metric/pattern tables retain canonical sample or control/test roles.
    """
    single = result["manifest"].get("analysis_type") == "single"
    roles = ("sample",) if single else ("control", "test")
    labels = single_sample_labels(sample_label) if single else sample_labels(control_label, test_label)
    summary_columns, summary_name = metric_columns(single), metric_filename(single)
    destination = Path(destination)
    if os.path.lexists(destination):
        raise FileExistsError(f"Output already exists: {destination}. Choose a new run directory.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent))
    try:
        keys = ["region_id", "chrom", "start", "end", "cpg_anchors", "unit_size"]
        sample_keys = ["n_candidate_reads", "n_complete_reads", "n_incomplete_reads", "callable_fraction", "beta", "entropy"]
        columns = keys + [f"{sample}_{key}" for sample in roles for key in sample_keys]
        if not single:
            columns += ["delta_beta", "delta_entropy", "epiallele_distance"]
        columns += ["passes_min_depth"]
        summary_rows = complete_rows = 0
        with (temporary / "metrics.tsv").open("w", newline="", encoding="utf-8") as handle, \
                (temporary / summary_name).open("w", newline="", encoding="utf-8") as summary:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            summary_writer = csv.writer(summary, delimiter="\t", lineterminator="\n")
            summary_writer.writerow(labelled_columns(summary_columns, labels))
            for unit in result["units"]:
                row = {key: unit[key] for key in keys}
                row["cpg_anchors"] = ",".join(map(str, unit["cpg_anchors"]))
                for sample in roles:
                    row.update({f"{sample}_{key}": unit[sample][key] for key in sample_keys})
                row.update(unit.get("comparison", {}))
                row["passes_min_depth"] = int(unit["passes_min_depth"])
                writer.writerow({key: "NA" if value is None else value for key, value in row.items()})
                complete_rows += has_complete_metrics(row, single)
                if passes_metric_filter(row, single):
                    summary_writer.writerow([format_cell(row[key]) for key in summary_columns])
                    summary_rows += 1
        with (temporary / "patterns.tsv").open("w", newline="", encoding="utf-8") as handle:
            columns = ["region_id", "chrom", "cpg_anchors", "sample", "pattern", "count", "complete"]
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for unit in result["units"]:
                for sample in roles:
                    counts = {format(i, f"0{unit['unit_size']}b"): n for i, n in enumerate(unit[sample]["counts"])}
                    counts.update(unit[sample]["incomplete_patterns"])
                    for pattern, count in sorted(counts.items()):
                        writer.writerow({"region_id": unit["region_id"], "chrom": unit["chrom"],
                                         "cpg_anchors": ",".join(map(str, unit["cpg_anchors"])),
                                         "sample": sample, "pattern": pattern, "count": count,
                                         "complete": int("." not in pattern)})
        manifest = dict(result["manifest"])
        manifest["depth_policy"] = DEPTH_POLICY
        manifest["sample_labels"] = labels
        if not single:
            manifest["labelled_direction"] = f"{test_label} minus {control_label}"
        manifest["outputs"] = {name: {"sha256": hashlib.sha256((temporary / name).read_bytes()).hexdigest()}
                               for name in ("metrics.tsv", "patterns.tsv", summary_name)}
        manifest["outputs"][summary_name].update(metric_filter_summary(
            len(result["units"]), complete_rows, summary_rows,
            manifest["config"]["min_depth"], single))
        (temporary / "analysis_summary.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        if os.path.lexists(destination):
            raise FileExistsError(f"Output appeared during analysis: {destination}; it was not overwritten.")
        os.rename(temporary, destination)
        return manifest
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
