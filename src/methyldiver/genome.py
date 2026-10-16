"""Parallel genome and selected-region analysis using the modBAM estimators.

Parquet is optional and imported only for chunked runs. Chunking preserves
every reference-CpG unit, including units spanning sequence-scan boundaries.
"""

from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from collections import Counter
from contextlib import ExitStack
import csv
from datetime import datetime, timezone
import gzip
import json
import multiprocessing
import os
from pathlib import Path
import re
import shutil
import tempfile

import pysam

from . import __version__
from .coordinates import MAX_DENSE_COUNT_CELLS
from .modbam import (Region, _digest, _file_identity, _index_path,
                     _validate_parameters, analyze_modbam, modification_config, require_modification)
from .output import (has_complete_metrics, passes_metric_filter, metric_filter_summary, DEPTH_POLICY,
                     sample_labels, single_sample_labels,
                     metric_columns, metric_filename)
from .reader import read_results


def _arrow():
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise ImportError("Chunked genome/region analysis requires pyarrow; install MethylDiver with the [genome] extra.") from exc
    return pa, pq


def _cpgs(fasta, chrom, length, scan_bases=1_000_000):
    """Scan reference tiles with one-base look-ahead, owning each C once."""
    for start in range(0, length, scan_bases):
        stop = min(length, start + scan_bases)
        sequence = fasta.fetch(chrom, start, min(length, stop + 1)).upper()
        for match in re.finditer("CG", sequence):
            yield start + match.start()


def _plan(fasta, contigs, k, chunk_units, scan_bases=1_000_000):
    jobs, summaries = [], []
    for chrom, length in contigs:
        positions = []
        n_cpgs = n_units = 0
        for position in _cpgs(fasta, chrom, length, scan_bases):
            positions.append(position)
            n_cpgs += 1
            if len(positions) == chunk_units + k - 1:
                jobs.append({"id": f"part_{len(jobs):06d}", "chrom": chrom,
                             "start": positions[0], "end": positions[-1] + 2,
                             "first_anchor": positions[0], "last_anchor": positions[chunk_units - 1],
                             "n_units": chunk_units})
                n_units += chunk_units
                positions = positions[-(k - 1):] if k > 1 else []
        remaining = max(0, len(positions) - k + 1)
        if remaining:
            jobs.append({"id": f"part_{len(jobs):06d}", "chrom": chrom,
                         "start": positions[0], "end": positions[-1] + 2,
                         "first_anchor": positions[0], "last_anchor": positions[remaining - 1],
                         "n_units": remaining})
            n_units += remaining
        assert n_units == max(0, n_cpgs - k + 1)
        summaries.append({"chrom": chrom, "length": length, "n_cpgs": n_cpgs, "n_units": n_units})
    return jobs, summaries


def _schemas(k, single=False):
    pa, _ = _arrow()
    fields = [("region_id", pa.string()), ("chrom", pa.string()), ("start", pa.int64()),
              ("end", pa.int64()), ("cpg_anchors", pa.list_(pa.int64())), ("unit_size", pa.int32())]
    for sample in (("sample",) if single else ("control", "test")):
        fields += [(f"{sample}_{key}", pa.int64()) for key in
                   ("n_candidate_reads", "n_complete_reads", "n_incomplete_reads")]
        fields += [(f"{sample}_{key}", pa.float64()) for key in ("callable_fraction", "beta", "entropy")]
        fields += [(f"{sample}_{i:0{k}b}", pa.int64()) for i in range(2 ** k)]
    if not single:
        fields += [(key, pa.float64()) for key in ("delta_beta", "delta_entropy", "epiallele_distance")]
    fields += [("passes_min_depth", pa.bool_())]
    incomplete = pa.schema([("region_id", pa.string()), ("chrom", pa.string()),
                            ("cpg_anchors", pa.list_(pa.int64())), ("sample", pa.string()),
                            ("pattern", pa.string()), ("count", pa.int64())])
    return pa.schema(fields), incomplete


def _write_chunk(result, stage, name, chunk_manifests, keep_patterns):
    pa, pq = _arrow()
    k = result["manifest"]["config"]["k"]
    single = result["manifest"].get("analysis_type") == "single"
    metrics_schema, incomplete_schema = _schemas(k, single)
    rows, incomplete = [], []
    for unit in result["units"]:
        row = {key: unit[key] for key in
               ("region_id", "chrom", "start", "end", "cpg_anchors", "unit_size", "passes_min_depth")}
        for sample in (("sample",) if single else ("control", "test")):
            values = unit[sample]
            row.update({f"{sample}_{key}": values[key] for key in
                        ("n_candidate_reads", "n_complete_reads", "n_incomplete_reads", "callable_fraction", "beta", "entropy")})
            row.update({f"{sample}_{i:0{k}b}": count for i, count in enumerate(values["counts"])})
            for pattern, count in values["incomplete_patterns"].items():
                incomplete.append({"region_id": unit["region_id"], "chrom": unit["chrom"],
                    "cpg_anchors": unit["cpg_anchors"], "sample": sample, "pattern": pattern, "count": count})
        row.update(unit.get("comparison", {}))
        rows.append(row)
    outputs = {}
    for folder, values, schema in (("complete_patterns", rows, metrics_schema),
                                   ("incomplete_patterns", incomplete, incomplete_schema)):
        relative = f"{folder}/{name}.parquet"
        path = Path(stage) / relative
        pq.write_table(pa.Table.from_pylist(values, schema=schema), path, compression="zstd")
        outputs[relative] = {"sha256": _digest(path), "rows": len(values), "retained": keep_patterns}
    manifest = dict(result["manifest"])
    manifest["outputs"] = outputs
    manifest["n_units_with_metrics"] = sum(has_complete_metrics(row, single) for row in rows)
    manifest["n_units_exported"] = sum(passes_metric_filter(row, single) for row in rows)
    relative = manifest_hash = None
    if chunk_manifests:
        relative = f"chunk_manifests/{name}.json.gz"
        data = (json.dumps(manifest, indent=2, allow_nan=False) + "\n").encode("utf-8")
        data = gzip.compress(data, mtime=0)
        path = Path(stage) / relative
        path.write_bytes(data)
        manifest_hash = _digest(path)
    return {"id": name, "n_units": len(rows), "n_units_passing_min_depth": manifest["n_units_passing_min_depth"],
            "modification_qc": manifest["modification_qc"],
            "n_units_with_metrics": manifest["n_units_with_metrics"],
            "n_units_exported": manifest["n_units_exported"],
            "manifest": relative, "manifest_sha256": manifest_hash, "outputs": outputs}


def _run_chunk(payload):
    job, inputs, settings, stage, chunk_manifests, keep_patterns = payload
    parts = job.get("parts", [{**job, "region_id": job["id"]}])
    regions = [Region(p["chrom"], p["start"], p["end"], p["region_id"]) for p in parts]
    result = analyze_modbam(*inputs, regions, max_total_bases=sum(r.end - r.start for r in regions),
                           _allow_absent_modification=True, **settings)
    units = result["units"]
    offset = 0
    for part in parts:
        subset = units[offset:offset + part["n_units"]]
        if (len(subset) != part["n_units"] or subset[0]["start"] != part["first_anchor"]
                or subset[-1]["start"] != part["last_anchor"]
                or any(u["region_id"] != part["region_id"] for u in subset)):
            raise ValueError(f"Reference CpG units changed after planning chunk {job['id']}.")
        offset += part["n_units"]
    if len(units) != job["n_units"]:
        raise ValueError(f"Reference CpG units changed after planning chunk {job['id']}.")
    return _write_chunk(result, stage, job["id"], chunk_manifests, keep_patterns)


def _analyze_chunked(control, test, reference, output, *, threads=8, chunk_units=10_000,
                   chunk_manifests=False, keep_patterns=False,
                   control_label="control", test_label="test", sample_label="sample",
                   k=3, modification="5mC", call_threshold=0.8, min_mapq=20, min_base_quality=10, min_depth=1,
                   max_reads_per_region=100_000, max_pattern_evaluations=5_000_000,
                   progress=None, status=None, regions=None, max_total_bases=100_000, regions_bed=None):
    """Write CpG units across the genome or independent selected regions.

    ``threads`` controls worker processes (default 8). ``progress(done, total, units)`` is
    optional and runs in the parent. ``status(phase, done, total)`` reports TSV
    writing and finalization. No output is published on failure. Existing output
    directories are never reused; there is no checkpoint/resume mode.
    ``chunk_manifests=True`` saves compressed JSON (.json.gz); False omits it.
    ``keep_patterns=True`` retains both Parquet folders; False removes them
    after producing the plotting TSV. Both options default to False.
    Whole-run provenance is always saved in analysis_summary.json.
    Sample labels customize only the plotting-table headings and run summary.
    ``test=None`` runs single-sample analysis; prefer analyze_single_genome for
    that use case.
    """
    single = test is None
    selected = regions is not None
    bed_path = Path(regions_bed).resolve() if regions_bed is not None else None
    bed_hash = None
    if bed_path is not None:
        from .modbam import read_regions_bed
        bed_hash = _digest(bed_path)
        regions = read_regions_bed(bed_path)
        if bed_hash != _digest(bed_path):
            raise ValueError("Regions BED changed while reading it.")
    if selected:
        regions = list(regions)
    modification_settings = modification_config(modification)
    roles = ("sample",) if single else ("control", "test")
    labels = single_sample_labels(sample_label) if single else sample_labels(control_label, test_label)
    summary_columns, summary_name = metric_columns(single), metric_filename(single)
    output = Path(output)
    if os.path.lexists(output):
        raise FileExistsError(f"Output already exists: {output}. Choose a new run directory.")
    if not isinstance(chunk_manifests, bool):
        raise ValueError("chunk_manifests must be True or False.")
    if not isinstance(keep_patterns, bool):
        raise ValueError("keep_patterns must be True or False.")
    for name, value in (("threads", threads), ("chunk_units", chunk_units)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be an integer >= 1.")
    k, _ = _validate_parameters(k, call_threshold, min_mapq, min_base_quality, min_depth,
                                max_total_bases if selected else 1, max_reads_per_region, max_pattern_evaluations)
    requested_chunk_units = chunk_units
    chunk_units = min(chunk_units, MAX_DENSE_COUNT_CELLS // (len(roles) * (2 ** k)))
    pa, pq = _arrow()
    started = datetime.now(timezone.utc).isoformat()
    sample_paths = dict(zip(roles, [Path(control).resolve()] if single else [Path(control).resolve(), Path(test).resolve()]))
    reference = Path(reference).resolve()
    paths = (*sample_paths.values(), reference)
    inputs = (Path(control).resolve(), None if single else Path(test).resolve(), reference)
    for path in paths:
        if not path.is_file():
            raise ValueError(f"Input file does not exist: {path}")
    indexes = [*(_index_path(path) for path in sample_paths.values()), Path(str(reference) + ".fai")]
    if not indexes[-1].is_file():
        raise ValueError(f"Missing FASTA index: {indexes[-1]}.")
    identities = [_file_identity(path) for path in (*paths, *indexes)]
    index_hashes = [_digest(path) for path in indexes]
    with ExitStack() as stack:
        bams = {role: stack.enter_context(pysam.AlignmentFile(str(path), "rb", index_filename=str(index)))
                for (role, path), index in zip(sample_paths.items(), indexes)}
        fasta = stack.enter_context(pysam.FastaFile(str(reference)))
        for label, bam in bams.items():
            if not bam.is_bam or not bam.has_index():
                raise ValueError(f"{label}: input must be an indexed BAM.")
        dictionaries = [dict(zip(bam.references, bam.lengths)) for bam in bams.values()]
        ar, fr = dictionaries[0], dict(zip(fasta.references, fasta.lengths))
        if selected:
            from .regions import plan_regions, validate_regions
            validate_regions(regions, [fr, *dictionaries])
            jobs, region_reports = plan_regions(fasta, regions, k, chunk_units, max_total_bases)
            used = {r.chrom for r in regions}
            contigs = [{"chrom": chrom, "length": fr[chrom],
                        "n_units": sum(r["n_units"] for r in region_reports if r["chrom"] == chrom)}
                       for chrom in dict.fromkeys(r.chrom for r in regions)]
        else:
            if not ar or any(other != ar for other in dictionaries[1:]):
                raise ValueError("Genome mode requires matching nonempty control/test BAM reference dictionaries." if not single
                                 else "Genome mode requires a nonempty BAM reference dictionary.")
            if any(fr.get(chrom) != length for chrom, length in ar.items()):
                raise ValueError("BAM/FASTA contig names or lengths do not match.")
            jobs, contigs = _plan(fasta, list(ar.items()), k, chunk_units)
            used = set(ar)
        unused = [chrom for chrom in fasta.references if chrom not in used]
        headers = {role: bam.header.to_dict() for role, bam in bams.items()}
    settings = {"k": k, "modification": modification, "call_threshold": call_threshold, "min_mapq": min_mapq,
                "min_base_quality": min_base_quality, "min_depth": min_depth,
                "max_reads_per_region": max_reads_per_region, "max_pattern_evaluations": max_pattern_evaluations}
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent))
    try:
        for directory in ("complete_patterns", "incomplete_patterns"):
            (stage / directory).mkdir()
        if chunk_manifests:
            (stage / "chunk_manifests").mkdir()
        with (stage / "chunks.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["id"] + (["region_id"] if selected else []) + ["chrom", "start", "end", "first_anchor", "last_anchor", "n_units"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(({"id": j["id"], **p} for j in jobs for p in j["parts"]) if selected else jobs)
        summaries = []
        completed_units = 0
        def record(summary):
            nonlocal completed_units
            summaries.append(summary)
            completed_units += summary["n_units"]
            if progress is not None:
                progress(len(summaries), len(jobs), completed_units)
        if progress is not None:
            progress(0, len(jobs), 0)
        if threads == 1:
            for job in jobs:
                record(_run_chunk((job, inputs, settings, stage, chunk_manifests, keep_patterns)))
        elif jobs:
            # Spawn avoids inheriting HTSlib/Arrow state across a fork. Bound the
            # queue so a failed chunk does not leave a whole genome queued.
            with ProcessPoolExecutor(max_workers=min(threads, len(jobs)),
                                     mp_context=multiprocessing.get_context("spawn")) as pool:
                iterator = iter(jobs)
                pending = set()
                def submit_one():
                    job = next(iterator, None)
                    if job is not None:
                        pending.add(pool.submit(_run_chunk, (job, inputs, settings, stage, chunk_manifests, keep_patterns)))
                for _ in range(min(threads * 2, len(jobs))):
                    submit_one()
                try:
                    while pending:
                        done, pending = wait(pending, return_when=FIRST_COMPLETED)
                        for future in done:
                            record(future.result())
                            submit_one()
                except BaseException:
                    for future in pending:
                        future.cancel()
                    raise
        else:
            # Typed empty datasets keep the result readable when no contig
            # contains k CpGs. No artificial metric values are introduced.
            for folder, schema in zip(("complete_patterns", "incomplete_patterns"), _schemas(k, single)):
                pq.write_table(pa.Table.from_pylist([], schema=schema), stage / folder / "empty.parquet", compression="zstd")
        # Serial zero-job runs need the same typed empty datasets.
        if not jobs and not list((stage / "complete_patterns").iterdir()):
            for folder, schema in zip(("complete_patterns", "incomplete_patterns"), _schemas(k, single)):
                pq.write_table(pa.Table.from_pylist([], schema=schema), stage / folder / "empty.parquet", compression="zstd")
        if completed_units != sum(c["n_units"] for c in contigs):
            raise ValueError("Analysis unit count differs from the reference plan.")
        modification_qc = {role: dict(sum((Counter(s["modification_qc"][role]) for s in summaries), Counter()))
                           for role in roles}
        require_modification(modification_qc, modification)
        for path, before in zip((*paths, *indexes), identities):
            after = _file_identity(path)
            if any(before[key] != after[key] for key in ("bytes", "mtime_ns")):
                raise ValueError(f"Input changed during genome analysis: {path}")
        if index_hashes != [_digest(path) for path in indexes]:
            raise ValueError("An input index changed during genome analysis.")
        if bed_path is not None and _digest(bed_path) != bed_hash:
            raise ValueError("Regions BED changed during analysis.")
        complete_rows = sum(s["n_units_with_metrics"] for s in summaries)
        expected_rows = sum(s["n_units_exported"] for s in summaries)
        if status is not None:
            status("write TSV", 0, expected_rows)
        def export_progress(rows):
            if status is not None:
                status("write TSV", rows, expected_rows)
        summary_rows = read_results(stage / "complete_patterns", output=stage / summary_name,
                                limit=None, columns=summary_columns, drop_na_metrics=True, require_min_depth=True,
                                control_label=control_label, test_label=test_label, sample_label=sample_label,
                                progress=export_progress)
        if summary_rows != expected_rows:
            raise ValueError("TSV unit count differs from units passing the metric and depth filters.")
        if status is not None:
            status("finalize", None, None)
        if not keep_patterns:
            for folder in ("complete_patterns", "incomplete_patterns"):
                shutil.rmtree(stage / folder)
        manifest = {"schema_version": "2.1", "package_version": __version__, "status": "complete",
            "mode": "regions" if selected else "genome", "started_at": started, "completed_at": datetime.now(timezone.utc).isoformat(),
            "scope": ("consecutive reference-CpG units fully contained within each requested region" if selected
                      else "every consecutive reference-CpG unit on all BAM reference contigs"),
            "coordinate_system": "BED0_half_open; CpG_anchors_forward_reference_C",
            "analysis_type": "single" if single else "paired",
            "config": {**settings, **modification_settings}, "threads": threads,
            "modification_qc": modification_qc,
            "sample_labels": labels,
            **({} if single else {"direction": "test_minus_control", "labelled_direction": f"{test_label} minus {control_label}"}),
            "chunk_units": chunk_units, "chunk_manifests": chunk_manifests,
            "requested_chunk_units": requested_chunk_units,
            "keep_patterns": keep_patterns,
            "contigs": contigs, "unused_fasta_contigs": unused,
            "boundary_policy": "chunks share k-1 CpGs; every complete unit belongs to exactly one chunk",
            "depth_policy": DEPTH_POLICY,
            "output_format": f"per-unit {summary_name}; optional Parquet datasets complete_patterns and incomplete_patterns",
            "outputs": {summary_name: {"sha256": _digest(stage / summary_name),
                **metric_filter_summary(completed_units, complete_rows, summary_rows, min_depth, single)}},
            "reference_validation": "BAM/FASTA contig names and lengths; sequence identity not proven",
            "inputs": dict(zip((*roles, "reference", *(f"{role}_index" for role in roles), "reference_index"), identities)),
            "index_sha256": dict(zip((*roles, "reference"), index_hashes)),
            "headers": headers, "n_chunks": len(jobs), "n_units": completed_units,
            "n_units_passing_min_depth": sum(s["n_units_passing_min_depth"] for s in summaries),
            "chunks": sorted(summaries, key=lambda s: s["id"]),
            "pysam_version": pysam.__version__, "pyarrow_version": pa.__version__,
            "source_sha256": {name: _digest(Path(__file__).with_name(name)) for name in
                ("genome.py", "regions.py", "modbam.py", "metrics.py", "coordinates.py", "output.py", "reader.py", "progress.py", "cli.py", "__init__.py")}}
        if selected:
            manifest.update({"regions": region_reports, "n_regions": len(regions),
                             "total_region_bases": sum(r.end - r.start for r in regions),
                             "max_total_bases": max_total_bases,
                             "region_batching": "span_and_unit_count; a single sparse unit may exceed the span target",
                             "qc_counting": "read processing occurrences across chunks/regions, not unique reads"})
            if bed_path is not None:
                manifest.update({"regions_bed": str(bed_path), "regions_bed_sha256": bed_hash})
        manifest["chunks_tsv_sha256"] = _digest(stage / "chunks.tsv")
        (stage / "analysis_summary.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        if os.path.lexists(output):
            raise FileExistsError(f"Output appeared during analysis: {output}; it was not overwritten.")
        os.rename(stage, output)
        return manifest
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def analyze_genome(control, test, reference, output, **kwargs):
    """Stream all reference-CpG units across BAM contigs into a new directory."""
    return _analyze_chunked(control, test, reference, output, **kwargs)


def analyze_regions(control, test, reference, regions, output, **kwargs):
    """Stream independent selected regions, with automatic batches and workers.

    max_total_bases (default 100,000) is a per-batch span target, not a
    whole-request limit. Other options follow analyze_genome. Original region
    IDs survive chunking; units never bridge separate input regions.
    """
    if regions is None:
        raise ValueError("Supply at least one explicit Region or BED interval.")
    return _analyze_chunked(control, test, reference, output, regions=regions, **kwargs)


def analyze_single_genome(sample, reference, output, **kwargs):
    """Calculate beta and entropy from one indexed modBAM across its genome."""
    return analyze_genome(sample, None, reference, output, **kwargs)
