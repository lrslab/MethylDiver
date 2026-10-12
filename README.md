# MethylDiver

MethylDiver calculates methylation level, normalized pattern entropy, and
epiallele distance directly from control and test modBAM files. Three consecutive
reference CpGs are the primary unit; four CpGs are also supported.

**Version:** 1.0.0 · **Python:** 3.10+ · **Input:** binary C/5mC modBAM

| Quantity | Definition | Interpretation |
| --- | --- | --- |
| β (`beta`) | Σ p(x) × number of methylated sites in x / k | Mean methylation across complete patterns |
| H (`entropy`) | −Σ p(x) log₂ p(x) / k | Normalized within-sample pattern diversity |
| D_epi | √JSD₂(P, Q) | Unsigned between-sample pattern-distribution distance |

All three quantities range from 0 to 1. `delta_beta` and `delta_entropy` are
**test minus control**. The two sample distributions have equal weight in JSD,
even when their sequencing depths differ. No pseudocount is added.

## Install

Clone the private repository using your authorized GitHub account:

```bash
git clone https://github.com/lxd98/MethylDiver.git
cd MethylDiver
python3 -m venv .venv
.venv/bin/python -m pip install .
.venv/bin/methyldiver --version
```

Python 3.10 or later is required. The package uses `pysam` to read modBAMs and
indexed FASTA references. The scalar metric module itself uses only Python's
standard library. No modkit executable is needed for normal analysis.
The optional NumPy batch API can be installed with:

```bash
.venv/bin/python -m pip install -e '.[arrays]'
```

The package is installed from the repository root. For development, use
`pip install -e .` inside the virtual environment. Releases are identified by Git
tags, starting with `v1.0.0`.

## Run with modBAM inputs

```bash
.venv/bin/methyldiver analyze \
  --control /path/to/control.bam \
  --test /path/to/test.bam \
  --reference /path/to/reference.fa \
  --regions-bed /path/to/regions.bed \
  --output outputs/my_analysis
```

Inputs:

- Coordinate-sorted control and test **modBAMs**, each with a `.bai` or `.csi` index.
- The matching reference FASTA and its existing `.fai` index. The command never
  creates or modifies input indexes.
- A BED3/BED4 file with nonoverlapping intervals. The optional fourth column is
  a unique region name. BED coordinates are **zero-based, end-exclusive**.

For a single interval, replace `--regions-bed` with `--region CHROM:START-END`.
This option also uses **zero-based, end-exclusive** coordinates; unlike some
other tools' region syntax, it is not one-based. `--region` can be repeated.
Use the contig names in your BAM/FASTA, which may be accession IDs rather than
`chr1`. Each interval is analyzed independently; units never bridge intervals.

The output must be a **new directory**. Existing runs are not overwritten.

| Output | Contents |
| --- | --- |
| `metrics.tsv` | One row per consecutive reference-CpG unit: anchors, BED interval, control/test β and H, Δβ, ΔH, D_epi, candidate/complete/incomplete read counts, callable fractions, and depth flag |
| `patterns.tsv` | Complete binary pattern counts (including zeros) and observed incomplete patterns, per sample and unit |
| `run.json` | Resolved settings, input/index/header identities, package/decoder versions, reference checks, region-level extraction QC, timestamps, and output hashes |

Missing estimates are `NA` in TSV and `null` in JSON. The command prints a short
summary with the output location and number of units. It emits every reference
unit, including units with no complete reads. D_epi and deltas require nonempty
complete-read counts in both samples.

Useful settings:

```text
--k 3                       Three CpGs per unit; use 4 for sensitivity analysis
--call-threshold 0.8         Minimum binary methylation-call confidence
--min-mapq 20                Exclude MAPQ 255 as unknown
--min-base-quality 10        Missing base qualities are excluded
--min-depth 1               Report whether BOTH samples pass this complete-read depth
```

`--min-depth` sets `passes_min_depth`; it does not remove rows or change the
estimators. For example, `--min-depth 15` flags units with at least 15 complete reads per
sample while retaining all counts. The default of 1 is a reporting
minimum, not a reliability guarantee.

Execution is region-bounded. By default total requested interval span is limited
to 100,000 bases, fetched reads to 100,000 per sample per region, and candidate
read/unit evaluations to 5,000,000 per sample per region. Exceeding a limit aborts
without publishing a completed output directory. The corresponding limit flags
are explicit overrides; no reads are silently downsampled. There is no implicit
whole-genome mode, parallel scheduler, or resume/cache reuse in this release.

## Try the complete workflow with synthetic BAMs

These commands create tiny indexed inputs without downloading data:

```bash
.venv/bin/python examples/make_modbam_example.py --output outputs/toy_inputs

.venv/bin/methyldiver analyze \
  --control outputs/toy_inputs/control.bam \
  --test outputs/toy_inputs/test.bam \
  --reference outputs/toy_inputs/reference.fa \
  --regions-bed outputs/toy_inputs/regions.bed \
  --output outputs/toy_analysis
```

Expected: two three-CpG units, 20 complete reads per sample per unit, β = 0.5
and H = 1/3 in both samples, Δβ = ΔH = 0, and D_epi = 1. Use new output names
when repeating these commands. With `--k 4`, there is one unit with H = 1/4.

## Modification-call policy

The direct decoder supports binary canonical C/5mC profiles: `C+m` and `G-m`
annotations (or 5mC ChEBI code `27551`), either alignment orientation, and
modern `MM`/`ML` or historical `Mm`/`Ml` tags. Explicit ML bytes are decoded as
`(byte + 0.5) / 256`. At the default threshold, state 1 requires p(5mC) ≥0.8;
state 0 requires p(5mC) ≤0.2. Intermediate calls remain missing.

Skipped `?` sites, skipped implicit `.`/omitted-mode sites, missing tags, and
missing probabilities **remain missing**. In particular, implicit low-modification
evidence is not assigned an exact numerical probability. This explicit-only
policy differs from historical modkit workflows that accept inferred canonical
calls; the selected policy is recorded in the run manifest.

Only primary mapped, nonduplicate, non-QC-failed reads passing mapping quality
are eligible. A candidate must span the unit's first and last modified positions;
internal CIGAR gaps or failed calls produce incomplete patterns. Reverse-strand
cytosines at reference G are anchored to the preceding reference C. Full CpG
dyads must lie within the requested interval. Calls are joined to reference CG
coordinates; a CpG-disrupting variant mask is not applied.

Malformed tags, stale MN lengths, duplicate read identities, repeated MI molecule
identities, competing cytosine modifications such as 5hmC, dual-strand cytosine
annotations, paired-end records, and marked ONT duplex/parent records are rejected.
No cross-platform or arbitrary multi-modification support is implied by successful
BAM reading. BAMs containing only ordinary alignments without modification tags
fail with a diagnostic when eligible reads are present.

Reference validation checks selected contig names and lengths. It records hashes
of the requested reference slices and index, but does not prove whole-reference
sequence identity from the BAM header. You must supply the actual matching
assembly. BAM/FASTA fingerprints are path/size/mtime plus index/header hashes,
not full large-file content hashes. Input paths and BAM header metadata appear
in the local manifest.

These rules follow the [SAM modification-tag specification](https://samtools.github.io/hts-specs/SAMtags.pdf)
and the [pysam modification API](https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.modified_bases_forward).

## Count-table command (also supported)

The included example contains one matched unit's counts in a tab-separated file:

```text
pattern	control	test
000	2	10
001	0	4
010	0	3
011	0	0
100	0	0
101	2	0
110	4	0
111	12	3
```

```bash
.venv/bin/methyldiver metrics-counts \
  --input examples/legacy_counts.tsv \
  --output outputs/legacy_metrics.json
```

The same command is available through `.venv/bin/python -m methyldiver`.
Omit `--output` to print JSON only. Existing outputs are refused unless
`--overwrite` is supplied. The input file cannot be overwritten, including
through a symlink or hard link.

Expected numerical results:

| Result | Value |
| --- | ---: |
| Control β | 0.8 |
| Control H | 0.5236501981515562 |
| Test β | 0.2666666666666667 |
| Test H | 0.5951584324091114 |
| Δβ | −0.5333333333333333 |
| ΔH | 0.07150823425755515 |
| D_epi | 0.7309379165994421 |

Floating-point values may differ in their last displayed digit. JSON also records
the package/schema version, input SHA-256, direction, and estimator conventions.

The CLI accepts exactly `pattern`, `control`, and `test` columns, with all eight
or sixteen binary patterns in ascending order, including zero-count patterns.
TSV counts must be nonnegative whole numbers written as decimal integers. Each
file represents **one fixed unit**, not counts pooled across different loci.

## Python API

```python
from methyldiver import metrics, compare_counts

control = [2, 0, 0, 0, 0, 2, 4, 12]
test = [10, 4, 3, 0, 0, 0, 0, 3]

print(metrics(control))
# {'n_complete_reads': 20, 'beta': 0.8, 'entropy': 0.5236501981515562}
print(compare_counts(control, test))
# delta_beta, delta_entropy, D_epi
```

Supply counts in `000, 001, …, 111` order, or `0000, …, 1111` for four CpGs.
The API accepts integer-valued floats for legacy compatibility, but rejects
booleans, strings, negative counts, fractional counts, and nonfinite values.
Scalar counts are normalized to Python integers to avoid fixed-width overflow.

If the total count is zero, β and H are `None`. If either comparison sample is
empty, all comparison metrics are `None`. The CLI serializes these as JSON
`null`. Eligibility thresholds are separate from metric calculation; this
package applies no default minimum-depth cutoff.

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

Matrices have shape `(n_units, 8)` or `(n_units, 16)`. Every paired row must refer
to the same reference assembly, contig, and ordered CpG anchors; count arrays
alone cannot verify this. Empty matrices with an explicit valid column count
are supported. Empty sample rows produce NumPy `NaN`, while the scalar/JSON API
uses `None`/`null`. Batch row totals must fit a signed 64-bit integer. Floating
calculations use float64; totals remain integer counts.

## Scientific scope

These are established complete-read plug-in estimators. This release provides
no confidence intervals,
P values, missing-pattern imputation, or biological-validation claim.

Complete-read estimates describe molecules passing the input call filters.
Pattern-dependent call failure can bias their relationship to all molecules.
Low depth depresses empirical entropy and can produce positive D_epi even
under identical underlying distributions. Choose depth thresholds appropriate for the experiment; no universal
reliability cutoff is applied.

Reads and overlapping units are not independent biological replicates.
Matching site-average bedMethyl values cannot reconstruct molecular patterns.

The modBAM command performs bounded extraction and call filtering. CRAM support,
cohort sample manifests, genome-wide chunking/resume, competing-modification
models, and resampling remain future stages.

## Tests and provenance

```bash
.venv/bin/python -m unittest discover -s tests -v
```

NumPy tests run when the optional dependency is installed; otherwise they are
explicitly skipped. Tests cover the legacy example, all pure patterns,
hidden redistribution, zero depth, invalid counts, signed direction,
depth-scaling invariance, 500 randomized comparisons against the original
scalar implementation, and 500 rows of scalar/batch agreement. Coordinate,
missingness, CLI, and output-preservation checks are included.

The additional modBAM tests cover four alignment/modification orientations,
MM/ML quantization and skipped calls, historical tag names, missing ML/MN checks,
quality filters, internal CIGAR gaps and clips, duplicate identities, interval
boundaries, unsupported states, resource limits, and complete CLI output.

The unchanged scalar baseline is in `tests/fixtures/legacy_metrics.py`.
[`tests/fixtures/baseline.json`](tests/fixtures/baseline.json) records its origin
and fingerprint. See [`NOTICE`](NOTICE) for attribution.
