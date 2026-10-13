# Changelog

## 2.0.0 - 2026-10-13

This release adds genome-wide and single-sample analysis, result export, and a
simpler command-line interface. It changes command names and output fields
relative to v1.0.0.

### Commands and calculations

- `compare` replaces `analyze` for paired control/test modBAMs.
- `single` calculates methylation level and normalized pattern entropy from one
  indexed modBAM; it does not calculate differences or Epiallele distance.
- `read` previews or exports retained Parquet tables as TSV or CSV.
- The standalone `metrics-counts` command is removed; Python count-vector APIs
  remain available.
- `--N` selects 1-10 consecutive CpGs per unit, default 3. The old CLI option
  `--k` is removed; Python APIs retain the keyword `k`.
- `epiallele_distance` replaces the `D_epi` output/API field. The underlying
  square-root Jensen-Shannon calculation is unchanged.
- Genome mode processes reference CpG units in parallel chunks while retaining
  units across chunk boundaries. Large N reduces the effective chunk size to
  bound dense count-vector memory.

### Outputs and usability

- Paired plotting table: `three_metrics.tsv`; single-sample plotting table:
  `two_metrics.tsv`. Both omit rows with missing required metrics.
- Genome outputs also include `chunks.tsv` and `analysis_summary.json`.
- `--keep-patterns T` retains `complete_patterns/` and `incomplete_patterns/`
  as compressed Parquet; default False removes them after TSV export.
- `--chunk-manifests T` retains per-chunk metadata as `.json.gz`; default False
  omits it.
- Sample labels customize plotting-table headings, with canonical roles retained
  in full result tables.
- Samtools-style help displays grouped options and defaults. Terminal progress
  shows phase, percentage, elapsed time, and ETA; completion uses plain text.
- `examples/plot_three_metrics.R` produces individual PDF panels using ggplot2
  from a paired plotting table and its adjacent analysis summary.

### Depth and interpretation

`--min-depth` defaults to 1 and reports whether each unit has enough complete
reads in the sample, or in both samples for `compare`. In this release it is a
reporting flag, not an output-row filter. Valid below-threshold estimates remain
in plotting TSVs. These TSVs do not contain read-depth columns. Retain full
pattern tables when subsequent depth filtering is required.

Patterns containing `.` are excluded from all metrics for that unit. Changing
N or averaging units into genomic windows does not replace a read-depth check.

### Validation

The Python regression suite covers the original metric baseline, variable unit
sizes, modBAM decoding, paired/single-sample agreement, serial/parallel genome
runs, output retention, result readers, and CLI validation. No patient data,
patient results, or local Codex handoff documents are included in this release.

## 1.0.0

Initial release: bounded paired modBAM extraction, complete-pattern metrics,
coordinate utilities, Python APIs, synthetic examples, and baseline tests.
