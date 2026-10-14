# Run with modBAM inputs

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
