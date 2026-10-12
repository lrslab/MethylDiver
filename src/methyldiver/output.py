"""Write complete modBAM run artifacts into a new output directory."""

import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile


def write_analysis(result, destination):
    """Publish metrics.tsv, patterns.tsv and run.json after successful analysis.

    Existing directories are never overwritten. Temporary output is removed if
    writing fails. TSV missing values are NA; manifest JSON uses null.
    """
    destination = Path(destination)
    if os.path.lexists(destination):
        raise FileExistsError(f"Output already exists: {destination}. Choose a new run directory.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent))
    try:
        keys = ["region_id", "chrom", "start", "end", "cpg_anchors", "unit_size"]
        sample_keys = ["n_candidate_reads", "n_complete_reads", "n_incomplete_reads", "callable_fraction", "beta", "entropy"]
        columns = keys + [f"{sample}_{key}" for sample in ("control", "test") for key in sample_keys]
        columns += ["delta_beta", "delta_entropy", "D_epi", "passes_min_depth"]
        with (temporary / "metrics.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for unit in result["units"]:
                row = {key: unit[key] for key in keys}
                row["cpg_anchors"] = ",".join(map(str, unit["cpg_anchors"]))
                for sample in ("control", "test"):
                    row.update({f"{sample}_{key}": unit[sample][key] for key in sample_keys})
                row.update(unit["comparison"])
                row["passes_min_depth"] = int(unit["passes_min_depth"])
                writer.writerow({key: "NA" if value is None else value for key, value in row.items()})
        with (temporary / "patterns.tsv").open("w", newline="", encoding="utf-8") as handle:
            columns = ["region_id", "chrom", "cpg_anchors", "sample", "pattern", "count", "complete"]
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for unit in result["units"]:
                for sample in ("control", "test"):
                    counts = {format(i, f"0{unit['unit_size']}b"): n for i, n in enumerate(unit[sample]["counts"])}
                    counts.update(unit[sample]["incomplete_patterns"])
                    for pattern, count in sorted(counts.items()):
                        writer.writerow({"region_id": unit["region_id"], "chrom": unit["chrom"],
                                         "cpg_anchors": ",".join(map(str, unit["cpg_anchors"])),
                                         "sample": sample, "pattern": pattern, "count": count,
                                         "complete": int("." not in pattern)})
        manifest = dict(result["manifest"])
        manifest["outputs"] = {name: {"sha256": hashlib.sha256((temporary / name).read_bytes()).hexdigest()}
                               for name in ("metrics.tsv", "patterns.tsv")}
        (temporary / "run.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        if os.path.lexists(destination):
            raise FileExistsError(f"Output appeared during analysis: {destination}; it was not overwritten.")
        os.rename(temporary, destination)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
