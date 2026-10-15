"""Region-bounded, explicit-call modBAM extraction for local epiallele metrics.

SAM tags: https://samtools.github.io/hts-specs/SAMtags.pdf
Decoder: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.modified_bases_forward

Supports 5mC and 5hmC at CpGs, using selected-modification-versus-rest states.
Implicit/skipped calls are missing, never invented probabilities.
"""

from array import array
from bisect import bisect_left
from collections import Counter
from contextlib import ExitStack
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

import pysam

from . import __version__
from .coordinates import MAX_DENSE_COUNT_CELLS, cpg_positions, sliding_units
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
MODIFICATION_CODES = {"5mC": ("m", 27551), "5hmC": ("h", 76792)}
_CODE_NAMES = {code: name for name, codes in MODIFICATION_CODES.items() for code in codes}


def modification_config(modification):
    if not isinstance(modification, str) or modification not in MODIFICATION_CODES:
        raise ValueError("modification must be 5mC or 5hmC; only CpG context is supported.")
    return {"modification": modification, "modification_context": "CpG",
            "state_space": f"{modification}_vs_not_{modification}",
            "pattern_states": {"1": modification, "0": f"not_{modification}", ".": "missing_or_uncertain"},
            "implicit_and_unknown_calls": "missing", "ml_decoding": "(byte+0.5)/256",
            "empty_mm_zero_ml_placeholder": "missing"}


def require_modification(qc_by_sample, modification):
    """Check presence across the whole requested scope, not individual chunks."""
    for sample, qc in qc_by_sample.items():
        if qc.get("reads_eligible", 0) and not qc.get("reads_with_selected_modification_profile", 0):
            raise ValueError(f"{sample}: no {modification} modification profile in eligible reads in the requested scope.")


def _modification_profile(read, qc, modification="5mC"):
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
    # Older Bonito/Remora BAMs can contain MM:Z: with ML:B:C,0.
    # There is no site to which that byte can be assigned. Treat only this
    # exact placeholder as absent evidence; all other length checks remain.
    if mm == "" and isinstance(ml, array) and ml.typecode == "B" and list(ml) == [0]:
        qc["reads_empty_mm_zero_ml_placeholder"] += 1
        qc["reads_without_selected_modification_profile"] += 1
        return None, fallback
    seen = set()
    target_strands = set()
    target_present = False
    entries = 0
    for group in mm.split(";")[:-1]:
        match = _GROUP.fullmatch(group)
        if not match:
            raise ValueError(f"Read {read.query_name}: malformed MM group {group!r}.")
        base, sign, codes, mode, deltas = match.groups()
        modifications = [int(codes)] if codes.isdigit() else list(codes)
        entries += (len(deltas.split(",")) if deltas else 0) * len(modifications)
        for code in modifications:
            key = (base, sign, _CODE_NAMES.get(code, code))
            if key in seen:
                raise ValueError(f"Read {read.query_name}: repeated modification group {key}.")
            seen.add(key)
        cytosine = (base, sign) in (("C", "+"), ("G", "-"))
        if cytosine:
            if any(code not in _CODE_NAMES for code in modifications):
                raise ValueError(f"Read {read.query_name}: unsupported cytosine modification; only 5mC and 5hmC are supported.")
            target_strands.add((sign == "-") ^ read.is_reverse)
            target_present |= any(code in MODIFICATION_CODES[modification] for code in modifications)
            qc["mm_unknown_skip_groups" if mode == "?" else "mm_implicit_skip_groups"] += 1
        elif any(code in _CODE_NAMES for code in modifications) or base == "N":
            raise ValueError(f"Read {read.query_name}: unsupported modification base/strand profile {group.split(',')[0]}.")
    if len(target_strands) > 1:
        raise ValueError(f"Read {read.query_name}: dual-strand cytosine calls require a separate molecule policy.")
    if ml is not None:
        if not isinstance(ml, array) or ml.typecode != "B" or len(ml) != entries:
            raise ValueError(f"Read {read.query_name}: ML must be a byte array matching all MM entries.")
    decoded = read.modified_bases_forward
    if decoded is None:
        raise ValueError(f"Read {read.query_name}: HTSlib could not decode MM/ML (invalid coordinates or tag contents).")
    # ML probabilities are absolute, not relative weights. Validate all
    # overlapping m/h entries before selecting a target. Allow the integer
    # boundary 256 for quantization; larger sums cannot describe alternatives.
    totals = Counter()
    for (base, strand, code), observations in decoded.items():
        if code in _CODE_NAMES:
            for position, probability in observations:
                if probability >= 0:
                    totals[base, strand, position] += probability
    if any(total > 256 for total in totals.values()):
        raise ValueError(f"Read {read.query_name}: overlapping 5mC/5hmC probabilities exceed 1 (ML byte sum >256).")
    ref_strand = int(next(iter(target_strands))) if target_strands else fallback
    if not target_present:
        qc["reads_without_selected_modification_profile"] += 1
        return None, ref_strand
    qc["reads_with_selected_modification_profile"] += 1
    if ml is None:
        qc["reads_missing_ml"] += 1
    return decoded, ref_strand


def decode_calls(read, anchors: set[int], threshold: float, min_base_quality: int, qc: Counter,
                 modification="5mC"):
    """Decode selected-modification probabilities at canonical CpG anchors.

    Returns (anchor -> 0/1/None, modification-reference-strand).
    Unknown/implicit calls stay absent, and failed explicit calls stay None.
    """
    modification_config(modification)
    decoded, ref_strand = _modification_profile(read, qc, modification)
    if decoded is None:
        return {}, ref_strand
    reference_positions = read.get_reference_positions(full_length=True)
    calls = {}
    qualities = read.query_qualities
    for (base, original_strand, code), observations in decoded.items():
        if code not in MODIFICATION_CODES[modification]:
            continue
        if (base, original_strand) not in (("C", 0), ("G", 1)):
            raise ValueError(f"Read {read.query_name}: unsupported decoded {modification} profile.")
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
                qc["calls_selected_modification"] += 1
            elif probability <= 1 - threshold:
                calls[anchor] = 0
                qc["calls_not_selected_modification"] += 1
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
        calls, strand = decode_calls(read, anchors, config["call_threshold"], config["min_base_quality"], qc,
                                     config["modification"])
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


def _validate_parameters(k, call_threshold, min_mapq, min_base_quality, min_depth,
                         max_total_bases, max_reads_per_region, max_pattern_evaluations):
    """Shared parameter validation for region and genome orchestration."""
    from .coordinates import _unit_size
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
    return k, integers


def analyze_modbam(control, test, reference, regions, *, k=3, modification="5mC", call_threshold=0.8,
                   min_mapq=20, min_base_quality=10, min_depth=1,
                   max_total_bases=100_000, max_reads_per_region=100_000,
                   max_pattern_evaluations=5_000_000, _allow_absent_modification=False):
    """Compare indexed modBAMs over explicit nonoverlapping BED0 regions.

    Returns all per-unit counts and point estimates, with a depth flag used to
    filter compact metric TSV exports. No biological inference
    or whole-genome enumeration is performed implicitly.
    ``test=None`` runs single-sample analysis; prefer analyze_single_modbam for
    that use case.
    """
    started_at = datetime.now(timezone.utc).isoformat()
    modification_settings = modification_config(modification)
    k, integers = _validate_parameters(k, call_threshold, min_mapq, min_base_quality,
        min_depth, max_total_bases, max_reads_per_region, max_pattern_evaluations)
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
    single = test is None
    sample_paths = {"sample": Path(control)} if single else {"control": Path(control), "test": Path(test)}
    reference = Path(reference)
    for path in (*sample_paths.values(), reference):
        if not path.is_file():
            raise ValueError(f"Input file does not exist: {path}")
    fai = Path(str(reference) + ".fai")
    if not fai.is_file():
        raise ValueError(f"Missing FASTA index: {fai}. Index the reference before running.")
    indexes = {role: _index_path(path) for role, path in sample_paths.items()}
    config = {"k": k, "call_threshold": call_threshold, **integers, **modification_settings,
              "read_filter": "primary_mapped_nonduplicate_nonqcfail", "exclude_mapq_255": True,
              "context_policy": "reference_CG_no_variant_mask",
              "unit_policy": "consecutive_reference_CpGs_full_dyads_within_each_region",
              "depth_policy": "filter_metric_tsv_preserve_full_counts"}
    inputs = {role: _file_identity(path) for role, path in sample_paths.items()}
    inputs["reference"] = _file_identity(reference)
    inputs["reference"]["fai_sha256"] = _digest(fai)
    results, region_reports = [], []
    with ExitStack() as stack:
        fasta = stack.enter_context(pysam.FastaFile(str(reference)))
        bams = {role: stack.enter_context(pysam.AlignmentFile(str(path), "rb", index_filename=str(indexes[role])))
                for role, path in sample_paths.items()}
        for label, bam in bams.items():
            index = indexes[label]
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
            if len(bams) * (len(results) + len(units)) * (2 ** k) > MAX_DENSE_COUNT_CELLS:
                raise ValueError("Dense-pattern count limit exceeded for this N; use smaller selected regions or --genome.")
            extracted = {role: _extract(bam, region, positions, units, config) for role, bam in bams.items()}
            eligible = 0
            for index, unit in enumerate(units):
                samples = {role: _sample_result(patterns[index], k) for role, (patterns, qc) in extracted.items()}
                passes = all(values["n_complete_reads"] >= min_depth for values in samples.values())
                eligible += passes
                result = {"region_id": region.region_id, "chrom": region.chrom,
                          "start": unit[0], "end": unit[-1] + 2, "cpg_anchors": list(unit),
                          "unit_size": k, **samples, "passes_min_depth": passes}
                if not single:
                    result["comparison"] = compare_counts(samples["control"]["counts"], samples["test"]["counts"])
                results.append(result)
            region_reports.append({**asdict(region), "region_id": region.region_id,
                                   "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                                   "n_reference_cpgs": len(positions), "n_units": len(units),
                                   "n_units_passing_min_depth": eligible,
                                   **{f"{role}_qc": qc for role, (patterns, qc) in extracted.items()}})
    modification_qc = {role: Counter() for role in sample_paths}
    for report in region_reports:
        for role in sample_paths:
            modification_qc[role].update({key: value for key, value in report[f"{role}_qc"].items()
                                          if isinstance(value, int)})
    modification_qc = {role: dict(qc) for role, qc in modification_qc.items()}
    if not _allow_absent_modification:
        require_modification(modification_qc, modification)
    # Abort if a source changed while it was being read; this is not a full BAM hash.
    for label, path in {**sample_paths, "reference": reference}.items():
        current = _file_identity(path)
        if any(current[key] != inputs[label][key] for key in ("bytes", "mtime_ns")):
            raise ValueError(f"{label} changed during analysis; output was not published.")
    manifest = {"schema_version": "2.1", "package_version": __version__, "status": "complete",
                "started_at": started_at, "completed_at": datetime.now(timezone.utc).isoformat(),
                "pysam_version": pysam.__version__, "htslib_version": pysam.__samtools_version__,
                "method": "complete_read_plugin_explicit_modbam_calls",
                "analysis_type": "single" if single else "paired",
                **({} if single else {"direction": "test_minus_control"}),
                "reference_validation": "selected_contig_names_and_lengths_only; sequence_identity_not_proven",
                "coordinate_system": "BED0_half_open; CpG_anchors_forward_reference_C",
                "config": config, "inputs": inputs, "regions": region_reports,
                "modification_qc": modification_qc,
                "n_units": len(results), "n_units_passing_min_depth": sum(r["passes_min_depth"] for r in results),
                "source_sha256": {name: _digest(Path(__file__).with_name(name)) for name in
                                  ("modbam.py", "metrics.py", "coordinates.py", "output.py", "cli.py", "__init__.py")}}
    return {"manifest": manifest, "units": results}


def analyze_single_modbam(sample, reference, regions, **kwargs):
    """Calculate beta and entropy from one indexed modBAM over explicit regions."""
    return analyze_modbam(sample, None, reference, regions, **kwargs)
