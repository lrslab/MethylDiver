# Run with modBAM inputs

```bash
.venv/bin/methyldiver compare \
  --control /path/to/control.bam \
  --test /path/to/test.bam \
  --reference /path/to/reference.fa \
  --regions-bed /path/to/regions.bed \
  --threads 24 \
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

## Automatic batching and multiple CPUs

The full BED can be passed directly, including large promoter collections.
MethylDiver checks interval validity and reference compatibility, then plans
chunks from genomic span and reference CpG-unit counts, not BED file bytes.
It streams results to disk and uses up to `--threads` worker processes (default 8).
Both `--regions-bed` and single/repeated `--region` support parallel workers.
A small request may have fewer chunks than available workers.

`--max-total-bases` is now a per-batch span target (default 100,000 bp), not a
limit on the entire CLI request. `--chunk-units` limits each chunk to 10,000
units by default, automatically reduced for large N to keep the dense pattern
vectors within 2,000,000 entries. Separate intervals can share a batch, but
units never bridge them. A long interval is split with N-1 shared CpGs;
each complete unit occurs exactly once, retaining its original region ID.
A single unit whose CpGs span more than the span target is preserved intact.
Neither this span target nor the count-entry limit is an exact RAM allocation;
worker count, read depth, and incomplete patterns also affect memory use.

Chunked processing is selected automatically when the request exceeds the span
target or a conservative dense-count estimate. It is also selected when
`--threads` is greater than 1 or any of `--chunk-units`, `--keep-patterns`,
or `--chunk-manifests` is explicitly supplied. Chunked processing requires
the `[genome]` installation extra, including for BED inputs. With the default
eight workers, selected-region runs use this chunked path.

## Outputs

Chunked region runs use the same layout as genome runs:

| Output | Contents |
| --- | --- |
| `three_metrics.tsv` | Per-unit three-metric table; original BED region IDs are retained |
| `chunks.tsv` | One row per region piece in a chunk, linking chunk ID to original region ID |
| `analysis_summary.json` | Requested regions, planned unit counts, settings, inputs, output hashes, and aggregated extraction QC |
| `complete_patterns/`, `incomplete_patterns/` | Compressed Parquet metric/depth and pattern tables, retained only with `--keep-patterns T` |
| `chunk_manifests/` | Detailed chunk extraction QC and provenance, saved as `.json.gz` only with `--chunk-manifests T` |

Small requests with explicit `--threads 1` and no explicit chunk/output
options retain the original output layout below for compatibility. They use the in-memory extractor. Its Python `analyze_modbam`
API still has whole-request limits; use `analyze_regions` for streaming Python
analysis. For single samples, the plotting table is `two_metrics.tsv`.

| Output | Contents |
| --- | --- |
| `three_metrics.tsv` | Per-unit three-metric table with coordinates and deltas (12 columns) |
| `metrics.tsv` | One row per consecutive reference-CpG unit: anchors, BED interval, control/test β and H, Δβ, ΔH, Epiallele distance, candidate/complete/incomplete read counts, callable fractions, and depth flag |
| `patterns.tsv` | Complete binary pattern counts (including zeros) and observed incomplete patterns, per sample and unit |
| `analysis_summary.json` | Resolved settings, input/index/header identities, package/decoder versions, reference checks, region-level extraction QC, timestamps, and output hashes |

Missing estimates are `NA` in full TSV tables and null in Parquet. Full tables
retain all reference units, including units with no complete reads. Plotting
tables omit rows with missing required metrics. Epiallele distance and deltas
require nonempty complete-read counts in both samples. These are per-unit
outputs, not one aggregate estimate per promoter.

Progress shows completed chunks, percentage, elapsed time, and ETA. Read QC
counters count processing occurrences, not unique reads: the same read may
appear in multiple regions or chunks. Results are published only after the
whole run succeeds; failed runs remove their temporary output.

Useful settings:

```text
--N 3                       Consecutive CpGs per unit (1–10; default 3)
--call-threshold 0.8         Minimum binary methylation-call confidence
--min-mapq 20                Exclude MAPQ 255 as unknown
--min-base-quality 10        Missing base qualities are excluded
--min-depth 1                Minimum complete-read depth for compact metric TSV rows
```

`--min-depth` filters `three_metrics.tsv` and `two_metrics.tsv`. For example,
`--min-depth 15` retains units with at least 15 complete reads in both samples
for `compare`, or in the sample for `single`. Retained full tables and the
in-memory extractor preserve all counts, estimates and `passes_min_depth`
flags. The default is 1; no threshold guarantees statistical reliability.

Fetched reads remain limited to 100,000 per sample per region piece, and candidate
read/unit evaluations to 5,000,000 per sample per piece. Exceeding these limits
aborts without publishing output. Reduce `--chunk-units` or explicitly raise
the corresponding limit; reads are never silently downsampled. Entire-genome
analysis still requires `--genome`.
