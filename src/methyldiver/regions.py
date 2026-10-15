"""Bounded planning for independent BED regions, preserving every CpG unit."""

from dataclasses import asdict
import re

from .modbam import Region


def validate_regions(regions, dictionaries):
    if not regions or any(not isinstance(r, Region) for r in regions):
        raise ValueError("Supply at least one explicit Region or BED interval.")
    if len({r.region_id for r in regions}) != len(regions):
        raise ValueError("Region IDs must be unique.")
    ordered = sorted(regions, key=lambda r: (r.chrom, r.start))
    for a, b in zip(ordered, ordered[1:]):
        if a.chrom == b.chrom and b.start < a.end:
            raise ValueError("Regions overlap; merge overlapping intervals first.")
    reference = dictionaries[0]
    for region in regions:
        if any(region.chrom not in d for d in dictionaries):
            raise ValueError(f"Region contig {region.chrom!r} absent from BAM or FASTA.")
        if any(d[region.chrom] != reference[region.chrom] for d in dictionaries[1:]):
            raise ValueError(f"BAM/FASTA length mismatch for {region.chrom}.")
        if region.end > reference[region.chrom]:
            raise ValueError(f"Region {region.region_id} extends beyond its reference contig.")


def _region_cpgs(fasta, region, scan_bases):
    for start in range(region.start, region.end, scan_bases):
        stop = min(region.end, start + scan_bases)
        sequence = fasta.fetch(region.chrom, start, min(region.end, stop + 1)).upper()
        for match in re.finditer("CG", sequence):
            yield start + match.start()


def plan_regions(fasta, regions, k, chunk_units, max_bases, scan_bases=1_000_000):
    """Pack region pieces by span and unit count; never join separate regions.

    Adjacent pieces of one region share k-1 CpGs, so boundary units occur once.
    A single sparse unit may exceed max_bases: its full span must be preserved.
    """
    jobs, reports, parts = [], [], []
    batch_span = batch_units = 0

    def flush():
        nonlocal parts, batch_span, batch_units
        if parts:
            jobs.append({"id": f"part_{len(jobs):06d}", "parts": parts, "n_units": batch_units})
            parts, batch_span, batch_units = [], 0, 0

    def emit(region, positions):
        nonlocal batch_span, batch_units
        n_units = len(positions) - k + 1
        span = positions[-1] + 2 - positions[0]
        # Repeated region IDs must not occur inside one analyze_modbam call.
        if parts and (batch_span + span > max_bases or batch_units + n_units > chunk_units
                      or parts[-1]["region_id"] == region.region_id):
            flush()
        parts.append({"region_id": region.region_id, "chrom": region.chrom,
                      "start": positions[0], "end": positions[-1] + 2,
                      "first_anchor": positions[0], "last_anchor": positions[n_units - 1],
                      "n_units": n_units})
        batch_span += span
        batch_units += n_units

    for region in regions:
        positions = []
        n_cpgs = 0
        for position in _region_cpgs(fasta, region, scan_bases):
            if len(positions) >= k and (len(positions) - k + 2 > chunk_units
                                       or position + 2 - positions[0] > max_bases):
                emit(region, positions)
                positions = positions[-(k - 1):] if k > 1 else []
            positions.append(position)
            n_cpgs += 1
        if len(positions) >= k:
            emit(region, positions)
        reports.append({**asdict(region), "region_id": region.region_id,
                        "n_cpgs": n_cpgs, "n_units": max(0, n_cpgs - k + 1)})
    flush()
    if sum(j["n_units"] for j in jobs) != sum(r["n_units"] for r in reports):
        raise ValueError("Region unit count differs from the reference plan.")
    return jobs, reports
