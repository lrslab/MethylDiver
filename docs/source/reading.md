# Read or export genome results

The reader requires retained Parquet files (`--keep-patterns T` for new runs).
For default runs, load `three_metrics.tsv` directly in your plotting software.

For existing genome results, create the same per-unit TSV without rerunning BAM
analysis:

```bash
methyldiver read --input MD_results/complete_patterns \
  --metrics-only --all --output MD_results/three_metrics_export.tsv
```

The `.csv` extension selects CSV automatically. Use `--format csv` to select CSV
explicitly, including when printing to the terminal. Otherwise output defaults
to TSV. `--metrics-only` selects the three-metric table columns and omits rows
missing any metric value or failing the stored `passes_min_depth` flag.
This also applies to retained tables from older runs, so a new export may
contain fewer rows than their original reporting-only TSV. Omit
`--metrics-only` to export all columns and retain below-threshold and missing
estimates. The reader uses the stored threshold; it does not choose a new one.
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
the selected modification, `0` not the selected modification, and `.` missing
or uncertain. A `0` does not necessarily mean unmodified cytosine. Counts are per unit
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
