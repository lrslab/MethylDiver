# Examples and plotting

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
TSV excludes missing metrics and, for new runs, units failing `--min-depth`.
It has no depth columns. The plot caption reads the saved export policy to
distinguish new filtered runs from older reporting-only runs. To change the
threshold later, retain the complete-pattern tables with `--keep-patterns T`. The example
provides no P values, significance tests, or confidence intervals.

For a two-target run, use `MD_results/5mC/three_metrics.tsv` or
`MD_results/5hmC/three_metrics.tsv` with its adjacent child summary. The plot
captions identify the target; 5hmC panels use hydroxymethylation labels.
