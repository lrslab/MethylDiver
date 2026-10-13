# MethylDiver

MethylDiver calculates methylation level, normalized pattern entropy, and
epiallele distance directly from control and test modBAM files with `compare`.
For one sample, `single` calculates methylation level and entropy. Each unit contains
N consecutive reference CpGs, selected with `--N` (default **3**, supported **1–10**).

**All three metrics use complete read patterns only.** In a pattern, `1` means
methylated, `0` means unmethylated, and `.` means a missing or uncertain call.
Any pattern containing `.` is excluded from methylation level, entropy, and
Epiallele distance calculations for that CpG unit. Missing calls are never
converted to `0` or `1`, and the observed sites of an incomplete pattern do not
contribute to these metrics.

| Example pattern | Treatment for that CpG unit |
| --- | --- |
| `101` or `100` | Complete; included in metric calculations |
| `1.0`, `.11`, or `...` | Incomplete; excluded from metrics and retained for quality checks |

Incomplete-pattern counts are saved in `incomplete_patterns/` for genome runs
when `--keep-patterns T` is selected,
and `patterns.tsv` for selected-region runs. Exclusion applies per unit: the
same read can contribute to other units where all CpG calls are complete.
The user-adjustable `--min-depth` threshold counts **complete reads in each
sample**; it sets a reporting flag and does not remove output rows or estimates.
If a sample has no complete reads, its methylation level and entropy are missing;
Epiallele distance requires complete reads in both samples.

**Version:** 2.0.0 · **Python:** 3.10+ · **Input:** binary C/5mC long-read modBAM

See [CHANGELOG.md](CHANGELOG.md) for release notes and changes from v1.0.0,
including `compare` replacing `analyze` and `--N` replacing `--k`.

| Quantity | Definition | Interpretation |
| --- | --- | --- |
| β (`beta`) | Σ p(x) × number of methylated sites in x / N | Mean methylation across complete patterns |
| H (`entropy`) | −Σ p(x) log₂ p(x) / N | Normalized within-sample pattern diversity |
| Epiallele distance (`epiallele_distance`) | √JSD₂(P, Q) | Unsigned between-sample pattern-distribution distance |

All three quantities range from 0 to 1. `delta_beta` and `delta_entropy` are
**test minus control**. The two sample distributions have equal weight in JSD,
even when their sequencing depths differ. No pseudocount is added.

The output/API field for **Epiallele distance** is `epiallele_distance`.
Version 2.0.0 replaces the earlier `D_epi` field name (output schema 2.0);
the square-root Jensen–Shannon calculation and its values are unchanged.
When reading older output tables, rename `D_epi` to `epiallele_distance`.
The NumPy API still returns the same three arrays in the same order.

## Install

Clone the private repository using your authorized GitHub account:

```bash
git clone https://github.com/lxd98/MethylDiver.git
cd MethylDiver
python3 -m venv .venv
.venv/bin/python -m pip install .
.venv/bin/methyldiver --version
```

Activate the environment to use the command directly, like other command-line
bioinformatics programs:

```bash
source .venv/bin/activate
methyldiver --help
methyldiver compare --help
methyldiver single --help
```

Help uses a samtools-style layout: program/version, concise usage, grouped
commands/options, and defaults in brackets (for example, `--min-mapq INT ... [20]`).
Running `methyldiver` alone shows the command list; running `methyldiver compare`,
`methyldiver single`, or `methyldiver read` alone shows that command's help.

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

To run from any directory without activating `.venv`, create a user-level command
link once, from the repository directory after installation:

```bash
mkdir -p "$HOME/.local/bin"
ln -s "$PWD/.venv/bin/methyldiver" "$HOME/.local/bin/methyldiver"
export PATH="$HOME/.local/bin:$PATH"
methyldiver --help
```

Keep `~/.local/bin` on your shell's PATH (add the export to `~/.bashrc` if needed).
The command uses the installed virtual environment's Python automatically, so it
also works from Conda `base`. Keep the repository and `.venv` at their installed
paths; moving or deleting them invalidates the link. This does not install
MethylDiver's Python dependencies into `base`.

## Single-sample analysis

Use `single` for methylation level (β) and normalized pattern entropy (H) from
one indexed modBAM. Supply the matching reference FASTA with its `.fai` index.
Genome mode requires the `[genome]` extra described below.

```bash
methyldiver single \
  --input /path/to/sample.bam \
  --reference /path/to/reference.fa \
  --genome \
  --N 3 \
  --threads 24 \
  --min-depth 5 \
  --output MD_single_results
```

For selected regions, replace `--genome --threads 24` with
`--region chr1:1000-2000` or `--regions-bed regions.bed`. Region coordinates
are 0-based and end-exclusive. Both modes use sliding units of N consecutive
reference CpGs (default 3, supported 1–10). The sample is extracted once per
region or genome chunk; no comparison metrics are calculated.

The default genome output contains:

| File | Contents |
| --- | --- |
| `two_metrics.tsv` | Per-unit coordinates, methylation level, and entropy |
| `chunks.tsv` | Genome chunk coordinates and unit counts |
| `analysis_summary.json` | Input identity, parameters, sample label, and QC summary |

`two_metrics.tsv` has seven tab-separated columns:

```text
region_id  chrom  start  end  cpg_anchors  sample_beta  sample_entropy
```

CpG anchors remain comma-separated within one TSV field. `--sample-label D0`
changes the final two headings to `D0_beta` and `D0_entropy`; the default label
is `sample`. There are no delta or Epiallele distance columns.

Only complete patterns contribute to β and H; a `.` at any of the N sites
excludes that read from the unit's metrics. Rows without complete reads are
omitted from `two_metrics.tsv`. `--min-depth` defaults to 1 and flags whether
the sample reaches that complete-read depth; it does not remove valid estimates
below the chosen threshold. Filtering and confidence defaults match `compare`:
mapping quality 20, base quality 10, and modification-call confidence 0.8.

Genome options `--keep-patterns T` and `--chunk-manifests T` retain the same
optional compressed folders as paired analysis; both default to False.
Retained tables use the canonical `sample` role, e.g. `sample_beta`,
`sample_n_complete_reads`, and `sample_000`. Selected-region runs also write
`metrics.tsv` (all units, metrics, and depth/QC fields) and `patterns.tsv`
(complete and incomplete counts); they do not write `chunks.tsv`.

The reader detects single-sample metric tables automatically:

```bash
methyldiver read --input MD_single_results/complete_patterns \
  --metrics-only --all --output single_metrics_export.tsv
```

This requires `--keep-patterns T` when running `single`. Saved sample labels
are recovered automatically; `read --sample-label NAME --metrics-only` overrides
the label. To filter incomplete-pattern files by role, use `--sample sample`.

## Compare two samples genome-wide

Use `compare` for paired samples. It replaces the former `analyze` command;
the old command name is no longer accepted.

Install the genome extra once. It adds PyArrow for compressed output tables:

```bash
python -m pip install '.[genome]'

methyldiver compare \
  --control /path/to/D0.bam \
  --test /path/to/D10.bam \
  --reference /path/to/reference.fa \
  --genome \
  --threads 24 \
  --min-depth 5 \
  --output outputs/genome_run
```

Both BAMs require existing indexes and matching reference-contig dictionaries.
The FASTA requires its `.fai` index. Every consecutive reference-CpG unit on all
BAM reference contigs is analyzed, including units without complete reads.
Extra FASTA contigs absent from both BAMs are recorded but are not analyzed.

Choose sample labels for the plotting table with optional flags:

```text
--control-label DMSO --test-label DEC
```

The defaults are `control` and `test`. With the example labels, `three_metrics.tsv`
uses `DMSO_beta`, `DEC_beta`, `DMSO_entropy`, and `DEC_entropy` as column names.
The other eight columns keep their names. `delta_beta` and `delta_entropy` remain
**test minus control**, so they represent **DEC minus DMSO** in this example;
Epiallele distance remains unsigned. Labels do not change the calculations.
These options work in both genome and selected-region modes.

`analysis_summary.json` records the role-to-label mapping and labelled comparison
direction. Full metric/QC tables, pattern tables, and Python calculation results
retain the canonical `control`/`test` roles, so existing filters and joins continue
to work. Labels must be distinct (ignoring case), start with a letter or digit,
and contain only letters, digits, underscores, periods, or hyphens. `delta` is
reserved to avoid duplicate metric column names.

`--threads` selects the number of worker processes (default 1).
`--chunk-units` controls units per chunk (default 10,000). Neighboring chunks
share `N-1` CpGs so boundary-spanning units are retained exactly once.
Use `--N 4` or `--N 5`, for example, to change the unit size. Omitting `--N`
uses three CpGs. There is no fixed
physical-span cutoff.

Complete patterns have `2^N` possible states (8 for N=3; 1,024 for N=10).
Genome runs automatically reduce the requested chunk size when needed to keep
both samples' dense count vectors within 2,000,000 entries per chunk; for N=10,
the default becomes 976 units per chunk for `compare`, or 1,953 for `single`.
`analysis_summary.json` records both
`requested_chunk_units` and the effective `chunk_units`. This bounds the count
vectors, not total memory: worker count, read depth, and incomplete patterns
also affect memory use. Selected-region analysis applies the same count-entry
limit across its regions; split large requests or use `--genome` if it is exceeded.
Larger N also requires complete calls at more sites, which can reduce usable
read depth. All N sites must be called for a read to contribute to the metrics.

**`--min-depth` is user-adjustable and defaults to 1.** The example uses 5;
this only sets `passes_min_depth`. All counts and estimates are retained, so
you can apply another depth threshold afterward without extracting BAMs again.
No threshold is a universal guarantee of statistical reliability.

Genome outputs:

| Output | Contents |
| --- | --- |
| `three_metrics.tsv` | Human-readable per-unit results: coordinates, control/test β and entropy, their deltas, and Epiallele distance (12 columns) |
| `complete_patterns/` | With `--keep-patterns T`: compressed Parquet tables with every unit's metrics, QC, and complete-pattern counts such as `control_000` and `test_111` |
| `incomplete_patterns/` | With `--keep-patterns T`: compressed Parquet tables of observed patterns containing `.`, per unit and sample; excluded from metric calculations; missing rows mean zero observations |
| `chunk_manifests/` | With `--chunk-manifests T`: per-chunk provenance, extraction QC, and checksums as `.json.gz` |
| `chunks.tsv` | Chunk plan and expected unit counts; chunk intervals overlap in their boundary CpGs and must not be passed together to `--regions-bed` |
| `analysis_summary.json` | Whole-run settings, inputs, reference scope, versions, code hashes, and completion summary |

`analysis_summary.json` records how the analysis was performed and summarizes
the completed results. It does not capture live progress or errors. Older runs
used the filename `run.json`; their contents remain valid.

The default genome output contains only `three_metrics.tsv`, `chunks.tsv`, and
`analysis_summary.json`. To keep detailed pattern/QC tables, add
`--keep-patterns T`; this retains both `complete_patterns/` and
`incomplete_patterns/`. With `--keep-patterns F` (default), those folders are
removed after the plotting table is produced. To keep compressed chunk metadata,
add `--chunk-manifests T`; its default is also `F`. The two options are independent
and apply to genome mode only. Existing runs are not renamed or deleted.

`incomplete_patterns/` contains **incomplete patterns only**. Complete binary-pattern
counts and per-unit QC are in `complete_patterns/` when retained. Older runs used
`metrics/` and, for some versions, `CpG_unit_pattern/`; `methyldiver read` still
accepts those directory names.

`three_metrics.tsv` is generated automatically for both genome and selected-region
analyses. It contains **one row per CpG unit with all seven metric values
available**, not a whole-sample average. Rows with `NA` in any metric column are
omitted from this TSV. Valid zero values are retained, as are units below
`--min-depth` when their metrics are available; removing missing estimates is
not a minimum-depth filter. Columns are separated by literal tab characters.
`cpg_anchors` is an unquoted, comma-separated list of zero-based positions in
one field, for example `3041,3043,3046` between two tabs.
Coordinates use BED0 half-open intervals. The columns, in order, are:

```text
region_id	chrom	start	end	cpg_anchors	control_beta	test_beta	delta_beta	control_entropy	test_entropy	delta_entropy	epiallele_distance
```

With `--keep-patterns T`, read depths, unit size, depth flags, and other QC fields remain in `complete_patterns/`
(or `metrics.tsv` for selected-region runs); they are not included in this summary TSV.
These full result tables retain every unit, including those with missing metrics.
`analysis_summary.json` records the summary TSV row count and how many units were omitted.
Detailed pattern counts remain in the original result tables. Whole-genome tables can be large; use the
reader's chromosome filter for a smaller table to inspect.

Keep pattern tables if you want to change depth filters or inspect detailed QC
later. With the default `F`, these details must be regenerated from the BAMs.

Use `--chunk-manifests` with `--genome` to choose how per-chunk metadata is saved:

| Option | Behavior |
| --- | --- |
| `--chunk-manifests True` | Save compressed `.json.gz` files |
| `--chunk-manifests False` | Omit the `chunk_manifests/` directory (default) |

`T` and `F` are also accepted; values are case-insensitive. Plain uncompressed
chunk manifests are no longer an output option. Both choices produce the same
analysis tables. `analysis_summary.json` is always saved, including the selected boolean,
per-chunk unit counts, and Parquet checksums with a `retained` flag that identifies
whether the corresponding file was kept. With `False`, detailed per-chunk
extraction QC and reference-slice hashes are not saved. This option applies to
new runs; it does not convert or remove files from existing results.

Parquet preserves nulls for missing estimates and stores `cpg_anchors` as an
integer list. If saved with `--keep-patterns T`, read `complete_patterns/` as one dataset:

```python
import pyarrow.dataset as ds

data = ds.dataset("outputs/genome_run/complete_patterns", format="parquet")
# Read only selected columns and units with >=10 complete reads in BOTH samples.
selected = data.to_table(
    columns=["chrom", "start", "end", "delta_beta", "delta_entropy", "epiallele_distance"],
    filter=(ds.field("control_n_complete_reads") >= 10)
           & (ds.field("test_n_complete_reads") >= 10),
)
```

Live progress goes to stderr, refreshing once per second on a terminal:

```text
methyldiver: analyze 81.6% 2775/3402 | elapsed 00:12:30 | ETA ~00:02:49
```

Percentage and ETA refer to the **current phase**. Analysis progress counts
completed chunks; the TSV-writing phase counts exported rows. Planning and
finalization show elapsed time with an unknown ETA. ETA is an estimate based
on the average progress rate for that phase and can change with workload.
Elapsed time covers the whole run. Wide terminals also show analyzed unit counts.
When redirected to a log, progress is written on separate lines about every five
seconds, plus phase transitions, without terminal control characters.

After all output is saved, a plain-text completion summary goes to stdout:

```text
methyldiver: analysis complete | elapsed 00:15:40
  Units analyzed: 33,896,006 | passing depth threshold: 629
  Plotting table: 4,908,608 rows | omitted with missing metrics: 28,987,398
  Output: /path/to/MD_results
  Files: three_metrics.tsv, chunks.tsv, analysis_summary.json
```

The numbers and timings above are illustrative. Machine-readable run details
remain in `analysis_summary.json`; `compare` no longer prints a JSON object at
completion.

Results are
published only after all chunks finish. The output directory must be new;
there is no resume mode. A failed or interrupted run removes its staged output.
Read and pattern-evaluation limits apply to each chunk. Reduce `--chunk-units`
or explicitly raise the corresponding limit if needed. `--max-total-bases`
applies to selected-region mode only.

## Read or export genome results

The reader requires retained Parquet files (`--keep-patterns T` for new runs).
For default runs, load `three_metrics.tsv` directly in your plotting software.

For existing genome results, create the same per-unit TSV without rerunning BAM
analysis:

```bash
methyldiver read --input MD_results/complete_patterns \
  --metrics-only --all --output MD_results/three_metrics.tsv
```

The `.csv` extension selects CSV automatically. Use `--format csv` to select CSV
explicitly, including when printing to the terminal. Otherwise output defaults
to TSV. `--metrics-only` selects the three-metric table columns and omits rows
missing any metric value; omit it to export all columns and retain missing estimates.
Add `--chrom` to export a single chromosome/contig.

With `--metrics-only`, the reader automatically reuses labels from the adjacent
analysis summary, falling back to `control` and `test` for older runs. Override
either label when exporting existing results, without rerunning BAM analysis:

```bash
methyldiver read --input MD_results/complete_patterns --metrics-only --all \
  --control-label DMSO --test-label DEC \
  --output MD_results/three_metrics_DMSO_DEC.tsv
```

Preview the first 20 rows without writing Python code:

```bash
methyldiver read --input MD_results/incomplete_patterns
```

Filter rows using exact values (filters are combined):

```bash
methyldiver read --input MD_results/incomplete_patterns \
  --chrom chr1 --sample control --pattern '1.0' --limit 50
```

Use the chromosome/contig names present in your results. In a pattern, `1` means
methylated, `0` unmethylated, and `.` missing or uncertain. Counts are per unit
and sample; the reader does not combine counts across units.

Export all matching rows to a new TSV file:

```bash
methyldiver read --input MD_results/incomplete_patterns \
  --all --output MD_results/incomplete_patterns.tsv
```

The input may also be one `.parquet` file or the `complete_patterns/` directory.
`--sample` and `--pattern` require the corresponding columns and apply to the
incomplete-pattern table. `--chrom` works with either table.
The default limit is **20 rows**, including when using `--output`; use `--all`
for a full export. Rows are streamed in batches. TSV output includes a header,
comma-separated zero-based `cpg_anchors`, and `NA` for missing values.
No matches produce a header-only TSV. Status messages go to stderr.
Output files must be new and outside the input Parquet directory. Both CSV and
TSV follow the same row limits and filters.
Reading requires the same `[genome]` extra used for genome analysis; pandas is
not required. Existing results from earlier versions can also be read.

## Run with modBAM inputs

```bash
.venv/bin/methyldiver compare \
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
| `three_metrics.tsv` | Per-unit three-metric table with coordinates and deltas (12 columns) |
| `metrics.tsv` | One row per consecutive reference-CpG unit: anchors, BED interval, control/test β and H, Δβ, ΔH, Epiallele distance, candidate/complete/incomplete read counts, callable fractions, and depth flag |
| `patterns.tsv` | Complete binary pattern counts (including zeros) and observed incomplete patterns, per sample and unit |
| `analysis_summary.json` | Resolved settings, input/index/header identities, package/decoder versions, reference checks, region-level extraction QC, timestamps, and output hashes |

Missing estimates are `NA` in TSV. The command prints a short
summary with the output location and number of units. It emits every reference
unit, including units with no complete reads. Epiallele distance and deltas require nonempty
complete-read counts in both samples.

Useful settings:

```text
--N 3                       Consecutive CpGs per unit (1–10; default 3)
--call-threshold 0.8         Minimum binary methylation-call confidence
--min-mapq 20                Exclude MAPQ 255 as unknown
--min-base-quality 10        Missing base qualities are excluded
--min-depth 1               Report whether BOTH samples pass this complete-read depth
```

`--min-depth` sets `passes_min_depth`; it does not remove rows or change the
estimators. For example, `--min-depth 15` flags units with at least 15 complete reads per
sample while retaining all counts. The default of 1 is a reporting
minimum, not a reliability guarantee.

In selected-region mode, total requested interval span is limited by default
to 100,000 bases, fetched reads to 100,000 per sample per region, and candidate
read/unit evaluations to 5,000,000 per sample per region. Exceeding a limit aborts
without publishing a completed output directory. The corresponding limit flags
are explicit overrides; no reads are silently downsampled. Genome-wide analysis
requires the explicit `--genome` option described above.

## Try the complete workflow with synthetic BAMs

These commands create tiny indexed inputs without downloading data:

```bash
.venv/bin/python examples/make_modbam_example.py --output outputs/toy_inputs

.venv/bin/methyldiver compare \
  --control outputs/toy_inputs/control.bam \
  --test outputs/toy_inputs/test.bam \
  --reference outputs/toy_inputs/reference.fa \
  --regions-bed outputs/toy_inputs/regions.bed \
  --output outputs/toy_analysis
```

Expected: two three-CpG units, 20 complete reads per sample per unit, β = 0.5
and H = 1/3 in both samples, Δβ = ΔH = 0, and Epiallele distance = 1. Use new output names
when repeating these commands. With `--N 4`, there is one unit with H = 1/4.

## Plot paired genome results

The optional R example creates **19 separate PDF panels**, one panel per file,
plus `plotting_notes.txt`. It uses `three_metrics.tsv` from a paired genome run
and its adjacent `analysis_summary.json` for sample labels and reference lengths.
It does not require retained Parquet tables.

Install the R dependencies once, then run from the repository directory:

```bash
Rscript -e 'install.packages(c("data.table", "ggplot2", "scales", "jsonlite"), repos="https://cloud.r-project.org")'
Rscript examples/plot_three_metrics.R MD_results/three_metrics.tsv MD_results/figures
```

Use an R build with Cairo PDF support and the DejaVu Sans font. Add `--overwrite`
only when replacing existing panels. The example includes distributions,
paired methylation/entropy density plots, change-versus-distance plots,
chromosome summaries, and 1-Mb genome tracks. Metric tracks show arithmetic
means in windows containing at least 100 plotted units; this is a display
threshold, not a complete-read depth threshold. Sparse inputs can have empty
panels where these display thresholds are not met.

These are exploratory plots of all rows in the supplied table. The plotting
TSV excludes missing metrics but **does not filter on `--min-depth`** and has
no depth columns. To apply a depth filter later, run with `--keep-patterns T`
and use the retained complete-pattern metrics to select units. The example
provides no P values, significance tests, or confidence intervals.

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

## Python API

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
separate from metric calculation; this package applies no default minimum-depth cutoff.

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

## Scientific scope

These are established complete-read plug-in estimators. This release provides
no confidence intervals,
P values, missing-pattern imputation, or biological-validation claim.

Complete-read estimates describe molecules passing the input call filters.
Pattern-dependent call failure can bias their relationship to all molecules.
Low depth depresses empirical entropy and can produce positive Epiallele distance even
under identical underlying distributions. Choose depth thresholds appropriate for the experiment; no universal
reliability cutoff is applied.

Reads and overlapping units are not independent biological replicates.
Matching site-average bedMethyl values cannot reconstruct molecular patterns.

The modBAM command performs bounded extraction and call filtering. CRAM support,
cohort sample manifests, resume/checkpoint reuse, competing-modification
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

With the genome extra installed, additional tests compare serial/parallel
genome runs to unchunked region analysis, including scan/chunk boundaries,
variable CpG-unit sizes, low-depth and empty units, typed empty outputs, and failure
cleanup. Genome tests are explicitly skipped when PyArrow is unavailable.

Single-sample tests verify agreement with paired sample metrics, extraction
without duplicate reads or distance calculations, incomplete patterns, zero-depth
and empty outputs, sample labels, serial/parallel genome runs, retention options,
reader exports, and cleanup after errors. CLI tests cover `compare`, `single`,
and `read` and reject the removed `analyze` and `metrics-counts` commands.

The unchanged scalar baseline is in `tests/fixtures/legacy_metrics.py`.
[`tests/fixtures/baseline.json`](tests/fixtures/baseline.json) records its origin
and fingerprint. See [`NOTICE`](NOTICE) for attribution.
