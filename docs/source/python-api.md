# Python API

Select `--modification 5mC` (default), `5hmC`, or `5mC,5hmC` in the CLI.
See [modification selection](modifications.md) for CpG-specific state definitions,
separate result folders, and the Python `modification` keyword.

Single-sample extraction is also available in Python:

```python
from methyldiver.modbam import Region, analyze_single_modbam
from methyldiver.output import write_analysis
from methyldiver.genome import analyze_single_genome

result = analyze_single_modbam("sample.bam", "reference.fa",
                               [Region("chr1", 1000, 2000)], k=3)
write_analysis(result, "single_region_results", sample_label="D0")
# Or process the genome, with the optional genome dependency installed:
# analyze_single_genome("sample.bam", "reference.fa", "single_genome_results",
#                       k=3, threads=24, sample_label="D0")
```

For large or parallel selected-region analysis, use the streaming API instead
of collecting all units with `analyze_modbam`:

```python
from methyldiver.modbam import read_regions_bed
from methyldiver.genome import analyze_regions
from pathlib import Path

bed = Path("promoters.bed")
analyze_regions("control.bam", "test.bam", "reference.fa",
                read_regions_bed(bed), "promoter_results",
                regions_bed=bed, threads=24, k=3, keep_patterns=True)
# For a single sample, pass None instead of the test BAM.
```

This API requires the `[genome]` extra. `max_total_bases=100_000` controls the
per-batch span target; `chunk_units=10_000` controls the requested units per
chunk, reduced automatically for large N. The original in-memory extraction
APIs keep their existing whole-request limits.

The scalar metric functions also accept complete-pattern count vectors:

```python
from methyldiver import metrics, compare_counts

control = [2, 0, 0, 0, 0, 2, 4, 12]
test = [10, 4, 3, 0, 0, 0, 0, 3]

print(metrics(control))
# {'n_complete_reads': 20, 'beta': 0.8, 'entropy': 0.5236501981515562}
print(compare_counts(control, test))
# delta_beta, delta_entropy, epiallele_distance
```

Supply `2^N` counts in ascending binary order, e.g. `000, 001, …, 111` for N=3
or `00000, …, 11111` for N=5. The count-vector length determines N.
Python functions that accept the unit size retain the keyword `k` for compatibility
(e.g. `analyze_modbam(..., k=5)`); the saved configuration also retains `config.k`.
The API accepts integer-valued floats for legacy compatibility, but rejects
booleans, strings, negative counts, fractional counts, and nonfinite values.
Scalar counts are normalized to Python integers to avoid fixed-width overflow.

If the total count is zero, β and H are `None`. If either comparison sample is
empty, all comparison metrics are `None`. Analysis tables represent missing
estimates as `NA` in TSV and null values in Parquet. Eligibility thresholds are
separate from metric calculation: scalar metric functions and in-memory
extraction retain all estimates. Compact TSV writers apply the extraction
`min_depth` flag (default threshold 1), requiring both samples for comparisons.

## Pre-called read patterns and coordinates

```python
from methyldiver import count_patterns, metrics, read_pattern
from methyldiver.coordinates import cpg_positions, sliding_units, unit_interval

positions = cpg_positions('CGACGTTTCGACG')  # [0, 3, 8, 11]
unit = next(sliding_units(positions))      # (0, 3, 8)
assert unit_interval(unit) == (0, 10)      # full-dyad BED0 interval

# Input calls are canonical-anchor → integer 0, integer 1, or None.
assert read_pattern({0: 1, 8: 0}, unit) == '1.0'
qc = count_patterns(['111', '1.0', '000', '101'])
assert qc['counts'] == [1, 0, 0, 0, 0, 1, 0, 1]
assert qc['callable_fraction'] == 0.75
print(metrics(qc['counts']))
```

Missing calls use `.` in pattern strings and never become state 0. Incomplete
patterns are retained in `incomplete_patterns` for QC and excluded from the
complete count vector. `n_candidate_reads` is the number of patterns supplied;
the caller must first enforce candidate alignment eligibility and deduplicate
molecules. Empty candidate input has a missing callable fraction.

Coordinates are zero-based forward-reference C anchors. `canonical_cpg` checks
reference CG context and maps a minus-strand modification at reference G to
the preceding C. Its strand argument is the modification's reference strand,
not simply a read's alignment orientation. `sliding_units` takes the complete
ordered reference CpG catalogue for one contig. A missing middle call does not
permit bridging to the next observed CpG. No maximum physical span is imposed.

## Optional batch calculations

```python
from methyldiver.arrays import array_metrics, compare_arrays

n, beta, entropy = array_metrics([control, test])
delta_beta, delta_entropy, distance = compare_arrays([control], [test])
```

Matrices have shape `(n_units, 2**N)` for N from 1 to 10. Every paired row must refer
to the same reference assembly, contig, and ordered CpG anchors; count arrays
alone cannot verify this. Empty matrices with an explicit valid column count
are supported. Empty sample rows produce NumPy `NaN`, while the scalar API
uses `None`. Batch row totals must fit a signed 64-bit integer. Floating
calculations use float64; totals remain integer counts.

Streaming genome and region APIs default to `threads=8`; pass `threads=1`
for serial execution. The bounded in-memory extraction APIs remain serial.
