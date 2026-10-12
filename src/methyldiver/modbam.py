"""Region-bounded, explicit-call modBAM extraction for local epiallele metrics.

SAM tags: https://samtools.github.io/hts-specs/SAMtags.pdf
Decoder: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.modified_bases_forward

Version 1 supports binary canonical C/5mC inputs, including C+m and G-m
annotations. Implicit/skipped calls are missing, never invented probabilities.
"""

from array import array
from bisect import bisect_left
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

import pysam

from . import __version__
from .coordinates import cpg_positions, sliding_units
from .metrics import compare_counts, metrics


@dataclass(frozen=True)
class Region:
    """A region in BED0 half-open coordinates."""

    chrom: str
    start: int
    end: int
    name: str = ""

    def __post_init__(self):
        if not isinstance(self.chrom, str) or not self.chrom or any(c.isspace() for c in self.chrom):
            raise ValueError("Region chromosome must be a nonempty name without whitespace.")
        if any(isinstance(v, bool) or not isinstance(v, int) for v in (self.start, self.end)):
            raise ValueError("Region coordinates must be integers.")
        if self.start < 0 or self.end <= self.start:
            raise ValueError("Regions require 0 <= start < end (BED0, end exclusive).")
        if any(c in self.name for c in "\t\r\n"):
            raise ValueError("Region names cannot contain tabs or newlines.")

    @property
    def region_id(self):
        return self.name or f"{self.chrom}:{self.start}-{self.end}"


def parse_region(value: str) -> Region:
    """Parse CHROM:START-END using zero-based, end-exclusive coordinates."""
    match = re.fullmatch(r"(.+):([0-9]+)-([0-9]+)", value)
    if not match:
        raise ValueError("Region must be CHROM:START-END in BED0 coordinates, e.g. chr1:1000-11000.")
    chrom, start, end = match.groups()
    return Region(chrom, int(start), int(end))


def read_regions_bed(path: Path) -> list[Region]:
    regions = []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            text = line.strip()
            if not text or text.startswith(("#", "track ", "browser ")):
                continue
            columns = text.split("\t")
            if len(columns) not in (3, 4):
                raise ValueError(f"BED line {number}: expected 3 or 4 tab-separated columns.")
            try:
                regions.append(Region(columns[0], int(columns[1]), int(columns[2]), columns[3] if len(columns) == 4 else ""))
            except ValueError as exc:
                raise ValueError(f"BED line {number}: {exc}") from exc
    if not regions:
        raise ValueError("The regions BED contains no intervals.")
    return regions


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_identity(path: Path) -> dict:
    stat = path.stat()
    return {"path": str(path.resolve()), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "full_content_sha256": "not_computed"}


def _index_path(path: Path) -> Path:
    for candidate in (Path(str(path) + ".bai"), path.with_suffix(".bai"),
                      Path(str(path) + ".csi"), path.with_suffix(".csi")):
        if candidate.is_file():
            return candidate
    raise ValueError(f"Missing BAM index (.bai or .csi) for {path}.")


def _tag(read, modern, historical):
    tags = [value for name, value in read.get_tags() if name in (modern, historical)]
    if len(tags) > 1:
        raise ValueError(f"Read {read.query_name}: duplicate or conflicting {modern}/{historical} tags.")
    return tags[0] if tags else None


_GROUP = re.compile(r"([ACGTUN])([+-])([A-Za-z]+|[0-9]+)([.?]?)(?:,([0-9]+(?:,[0-9]+)*))?")


def _modification_profile(read, qc):
    """Validate tags, identify the single supported modification strand.

    Returns (decoded modifications in original read coordinates, ref strand),
    or (None, alignment strand) when modification evidence is absent.
    """
    fallback = int(read.is_reverse)
    mm = _tag(read, "MM", "Mm")
    ml = _tag(read, "ML", "Ml")
    if mm is None:
        if ml is not None:
            raise ValueError(f"Read {read.query_name}: ML is present without MM.")
        qc["reads_missing_mm"] += 1
        return None, fallback
    qc["reads_with_mm"] += 1
    if not read.query_sequence:
        raise ValueError(f"Read {read.query_name}: MM tags require a query sequence.")
    if read.has_tag("MN") and read.get_tag("MN") != read.query_length:
        raise ValueError(f"Read {read.query_name}: MN does not match query length; modification tags may be stale.")
    if not isinstance(mm, str) or (mm and not mm.endswith(";")):
        raise ValueError(f"Read {read.query_name}: malformed MM tag.")
    seen = set()
    target_strands = set()
    entries = 0
    for group in mm.split(";")[:-1]:
        match = _GROUP.fullmatch(group)
        if not match:
            raise ValueError(f"Read {read.query_name}: malformed MM group {group!r}.")
        base, sign, codes, mode, deltas = match.groups()
        modifications = [int(codes)] if codes.isdigit() else list(codes)
        entries += (len(deltas.split(",")) if deltas else 0) * len(modifications)
        for code in modifications:
            key = (base, sign, code)
            if key in seen:
                raise ValueError(f"Read {read.query_name}: repeated modification group {key}.")
            seen.add(key)
        cytosine = (base, sign) in (("C", "+"), ("G", "-"))
        if cytosine:
            if len(modifications) != 1 or modifications[0] not in ("m", 27551):
                raise ValueError(f"Read {read.query_name}: competing/unsupported cytosine modifications; only binary C/5mC is supported.")
            target_strands.add((sign == "-") ^ read.is_reverse)
            qc["mm_unknown_skip_groups" if mode == "?" else "mm_implicit_skip_groups"] += 1
        elif "m" in modifications or 27551 in modifications or base == "N":
            raise ValueError(f"Read {read.query_name}: unsupported modification base/strand profile {group.split(',')[0]}.")
    if len(target_strands) > 1:
        raise ValueError(f"Read {read.query_name}: dual-strand cytosine calls require a separate molecule policy.")
    if ml is not None:
        if not isinstance(ml, array) or ml.typecode != "B" or len(ml) != entries:
            raise ValueError(f"Read {read.query_name}: ML must be a byte array matching all MM entries.")
    decoded = read.modified_bases_forward
    if decoded is None:
        raise ValueError(f"Read {read.query_name}: HTSlib could not decode MM/ML (invalid coordinates or tag contents).")
    if not target_strands:
        qc["reads_without_5mc_profile"] += 1
        return None, fallback
    qc["reads_with_5mc_profile"] += 1
    if ml is None:
        qc["reads_missing_ml"] += 1
    return decoded, int(next(iter(target_strands)))


def decode_calls(read, anchors: set[int], threshold: float, min_base_quality: int, qc: Counter):
    """Decode explicit C/5mC probabilities and map them to canonical anchors.

    Returns (anchor -> 0/1/None, modification-reference-strand).
    Unknown/implicit calls stay absent, and failed explicit calls stay None.
    """
    decoded, ref_strand = _modification_profile(read, qc)
    if decoded is None:
        return {}, ref_strand
    reference_positions = read.get_reference_positions(full_length=True)
    calls = {}
    qualities = read.query_qualities
    for (base, original_strand, modification), observations in decoded.items():
        if modification not in ("m", 27551):
            continue
        if (base, original_strand) not in (("C", 0), ("G", 1)):
            raise ValueError(f"Read {read.query_name}: unsupported decoded 5mC profile.")
        strand = original_strand ^ int(read.is_reverse)
        for forward_query, raw_probability in observations:
            query = read.query_length - 1 - forward_query if read.is_reverse else forward_query
            if not 0 <= query < len(reference_positions):
                raise ValueError(f"Read {read.query_name}: invalid modification query position.")
            reference_position = reference_positions[query]
            if reference_position is None:
                qc["calls_unaligned"] += 1
                continue
            anchor = reference_position - strand
            if anchor not in anchors:
                qc["calls_outside_region_reference_cpgs"] += 1
                continue
            if anchor in calls:
                raise ValueError(f"Read {read.query_name}: duplicate call for CpG {anchor}.")
            calls[anchor] = None
            qc["explicit_cpg_entries"] += 1
            if raw_probability == -1:
                qc["calls_missing_probability"] += 1
                continue
            if not 0 <= raw_probability <= 255:
                raise ValueError(f"Read {read.query_name}: ML byte outside 0..255.")
            if qualities is None or qualities[query] == 255 or qualities[query] < min_base_quality:
                qc["calls_low_or_missing_base_quality"] += 1
                continue
            probability = (raw_probability + 0.5) / 256
            if probability >= threshold:
                calls[anchor] = 1
                qc["calls_methylated"] += 1
            elif probability <= 1 - threshold:
                calls[anchor] = 0
                qc["calls_unmethylated"] += 1
            else:
                qc["calls_low_modification_confidence"] += 1
    return calls, ref_strand


def _extract(bam, region, positions, units, config):
    patterns = [Counter() for _ in units]
    qc = Counter()
    if not units:
        return patterns, {"reads_fetched": 0, "status": "fewer_than_k_reference_cpgs"}
    anchors = set(positions)
    read_keys, molecule_keys = set(), set()
    for read in bam.fetch(region.chrom, region.start, region.end):
        qc["reads_fetched"] += 1
        if qc["reads_fetched"] > config["max_reads_per_region"]:
            raise ValueError(f"{region.region_id}: read limit exceeded; choose a smaller region or explicitly raise --max-reads-per-region.")
        if read.flag & 3844:
            qc["reads_excluded_flags"] += 1
            continue
        if read.mapping_quality == 255 or read.mapping_quality < config["min_mapq"]:
            qc["reads_excluded_mapq"] += 1
            continue
        if read.is_paired:
            raise ValueError("Paired-end records are unsupported; supply primary long-read molecule alignments.")
        if read.has_tag("dx") and read.get_tag("dx") != 0:
            raise ValueError("Duplex/parent records need a validated molecule-selection policy; v1 supports simplex records.")
        key = (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name)
        if not read.query_name or key in read_keys:
            raise ValueError(f"{region.region_id}: missing or duplicate primary read identity {read.query_name!r}.")
        read_keys.add(key)
        if read.has_tag("MI"):
            molecule = read.get_tag("MI")
            if molecule in molecule_keys:
                raise ValueError(f"{region.region_id}: multiple primary records for the same MI molecule.")
            molecule_keys.add(molecule)
        qc["reads_eligible"] += 1
        calls, strand = decode_calls(read, anchors, config["call_threshold"], config["min_base_quality"], qc)
        # Endpoints must span the modified positions (C on +, G on -).
        # Internal gaps remain missing; we never substitute the next called CpG.
        first = bisect_left(positions, read.reference_start - strand)
        stop = min(len(units), bisect_left(positions, read.reference_end - strand) - config["k"] + 1)
        for index in range(first, max(first, stop)):
            qc["candidate_patterns"] += 1
            if qc["candidate_patterns"] > config["max_pattern_evaluations"]:
                raise ValueError(f"{region.region_id}: pattern-evaluation limit exceeded; choose a smaller region.")
            pattern = "".join("." if calls.get(p) is None else str(calls[p]) for p in units[index])
            patterns[index][pattern] += 1
            qc["incomplete_patterns" if "." in pattern else "complete_patterns"] += 1
    if qc["reads_eligible"] and not qc["reads_with_mm"]:
        raise ValueError(f"{region.region_id}: no MM/Mm tags in eligible reads; input may not be a modBAM.")
    if qc["reads_with_mm"] and not qc["reads_with_5mc_profile"]:
        raise ValueError(f"{region.region_id}: no supported 5mC modification profile in eligible reads.")
    return patterns, dict(qc)


def _sample_result(patterns, k):
    counts = [patterns[format(i, f"0{k}b")] for i in range(2 ** k)]
    candidate = sum(patterns.values())
    result = metrics(counts)
    result.update({"counts": counts, "n_candidate_reads": candidate,
                   "n_incomplete_reads": candidate - sum(counts),
                   "callable_fraction": sum(counts) / candidate if candidate else None,
                   "incomplete_patterns": {p: n for p, n in sorted(patterns.items()) if "." in p}})
    return result


def analyze_modbam(control, test, reference, regions, *, k=3, call_threshold=0.8,
                   min_mapq=20, min_base_quality=10, min_depth=1,
                   max_total_bases=100_000, max_reads_per_region=100_000,
                   max_pattern_evaluations=5_000_000):
    """Compare indexed modBAMs over explicit nonoverlapping BED0 regions.

    Returns a run manifest and per-unit results. Depth filtering is a reporting
    flag; all counts and point estimates are retained. No biological inference
    or whole-genome enumeration is performed implicitly.
    """
    from .coordinates import _unit_size

    started_at = datetime.now(timezone.utc).isoformat()
    k = _unit_size(k)
    if isinstance(call_threshold, bool) or not isinstance(call_threshold, (int, float)) or not math.isfinite(call_threshold) or not 0.5 < call_threshold <= 1:
        raise ValueError("Call confidence threshold must be finite and in (0.5, 1].")
    integers = {"min_mapq": min_mapq, "min_base_quality": min_base_quality, "min_depth": min_depth,
                "max_total_bases": max_total_bases, "max_reads_per_region": max_reads_per_region,
                "max_pattern_evaluations": max_pattern_evaluations}
    for name, value in integers.items():
        lower = 0 if name in ("min_mapq", "min_base_quality") else 1
        if isinstance(value, bool) or not isinstance(value, int) or value < lower:
            raise ValueError(f"{name} must be an integer >= {lower}.")
    if min_mapq > 254 or min_base_quality > 254:
        raise ValueError("Mapping/base quality thresholds must be <= 254; quality 255 is missing.")
    regions = list(regions)
    if not regions or any(not isinstance(region, Region) for region in regions):
        raise ValueError("Supply at least one explicit Region or BED interval.")
    if len({r.region_id for r in regions}) != len(regions):
        raise ValueError("Region IDs must be unique.")
    ordered = sorted(regions, key=lambda r: (r.chrom, r.start, r.end))
    for a, b in zip(ordered, ordered[1:]):
        if a.chrom == b.chrom and b.start < a.end:
            raise ValueError("Regions overlap; merge overlapping intervals first to avoid counting units twice.")
    if sum(r.end - r.start for r in regions) > max_total_bases:
        raise ValueError("Selected intervals exceed --max-total-bases; choose smaller regions or explicitly raise the limit.")
    control, test, reference = map(Path, (control, test, reference))
    for path in (control, test, reference):
        if not path.is_file():
            raise ValueError(f"Input file does not exist: {path}")
    fai = Path(str(reference) + ".fai")
    if not fai.is_file():
        raise ValueError(f"Missing FASTA index: {fai}. Index the reference before running.")
    indexes = [_index_path(control), _index_path(test)]
    config = {"k": k, "call_threshold": call_threshold, **integers,
              "ml_decoding": "(byte+0.5)/256", "implicit_and_unknown_calls": "missing",
              "read_filter": "primary_mapped_nonduplicate_nonqcfail", "exclude_mapq_255": True,
              "state_space": "canonical_C_vs_5mC", "context_policy": "reference_CG_no_variant_mask",
              "unit_policy": "consecutive_reference_CpGs_full_dyads_within_each_region",
              "depth_policy": "report_flag_only_preserve_all_estimates"}
    inputs = {"control": _file_identity(control), "test": _file_identity(test), "reference": _file_identity(reference)}
    inputs["reference"]["fai_sha256"] = _digest(fai)
    results, region_reports = [], []
    with pysam.FastaFile(str(reference)) as fasta, \
            pysam.AlignmentFile(str(control), "rb", index_filename=str(indexes[0])) as control_bam, \
            pysam.AlignmentFile(str(test), "rb", index_filename=str(indexes[1])) as test_bam:
        for label, bam, index in zip(("control", "test"), (control_bam, test_bam), indexes):
            if not bam.is_bam or not bam.has_index():
                raise ValueError(f"{label}: input must be an indexed BAM (CRAM is not yet supported).")
            inputs[label]["index_path"] = str(index.resolve())
            inputs[label]["index_sha256"] = _digest(index)
            header = bam.header.to_dict()
            inputs[label]["header_sha256"] = hashlib.sha256(json.dumps(header, sort_keys=True).encode()).hexdigest()
            inputs[label]["read_groups"] = header.get("RG", [])
            inputs[label]["programs"] = header.get("PG", [])
            for region in regions:
                if region.chrom not in fasta.references or region.chrom not in bam.references:
                    raise ValueError(f"{label}: region contig {region.chrom!r} absent from BAM or FASTA.")
                if bam.get_reference_length(region.chrom) != fasta.get_reference_length(region.chrom):
                    raise ValueError(f"{label}: BAM/FASTA length mismatch for {region.chrom}.")
                if region.end > fasta.get_reference_length(region.chrom):
                    raise ValueError(f"Region {region.region_id} extends beyond its reference contig.")
        for region in regions:
            sequence = fasta.fetch(region.chrom, region.start, region.end).upper()
            positions = cpg_positions(sequence, region.start)
            units = list(sliding_units(positions, k))
            cp, cq = _extract(control_bam, region, positions, units, config)
            tp, tq = _extract(test_bam, region, positions, units, config)
            eligible = 0
            for unit, a, b in zip(units, cp, tp):
                am, bm = _sample_result(a, k), _sample_result(b, k)
                passes = min(am["n_complete_reads"], bm["n_complete_reads"]) >= min_depth
                eligible += passes
                results.append({"region_id": region.region_id, "chrom": region.chrom,
                                "start": unit[0], "end": unit[-1] + 2, "cpg_anchors": list(unit),
                                "unit_size": k, "control": am, "test": bm,
                                "comparison": compare_counts(am["counts"], bm["counts"]),
                                "passes_min_depth": passes})
            region_reports.append({**asdict(region), "region_id": region.region_id,
                                   "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                                   "n_reference_cpgs": len(positions), "n_units": len(units),
                                   "n_units_passing_min_depth": eligible, "control_qc": cq, "test_qc": tq})
    # Abort if a source changed while it was being read; this is not a full BAM hash.
    for label, path in (("control", control), ("test", test), ("reference", reference)):
        current = _file_identity(path)
        if any(current[key] != inputs[label][key] for key in ("bytes", "mtime_ns")):
            raise ValueError(f"{label} changed during analysis; output was not published.")
    manifest = {"schema_version": "1.0", "package_version": __version__, "status": "complete",
                "started_at": started_at, "completed_at": datetime.now(timezone.utc).isoformat(),
                "pysam_version": pysam.__version__, "htslib_version": pysam.__samtools_version__,
                "method": "complete_read_plugin_explicit_modbam_calls", "direction": "test_minus_control",
                "reference_validation": "selected_contig_names_and_lengths_only; sequence_identity_not_proven",
                "coordinate_system": "BED0_half_open; CpG_anchors_forward_reference_C",
                "config": config, "inputs": inputs, "regions": region_reports,
                "n_units": len(results), "n_units_passing_min_depth": sum(r["passes_min_depth"] for r in results),
                "source_sha256": {name: _digest(Path(__file__).with_name(name)) for name in
                                  ("modbam.py", "metrics.py", "coordinates.py", "output.py", "cli.py", "__init__.py")}}
    return {"manifest": manifest, "units": results}
