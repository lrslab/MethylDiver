# Changelog

## 2.1.0 - 2026-10-14

- Merge the full user documentation into `docs/`; keep the main README concise.

- Add `--modification 5mC` (default), `5hmC`, or `5mC,5hmC` to both `compare`
  and `single`. All targets use reference CpG units; 4mC, 6mA, and non-CpG
  contexts are not analysis targets in this release.
- Decode mixed 5mC/5hmC tags in combined or separate groups, including numeric
  aliases and both supported strand orientations. Reject duplicate aliases
  and incompatible overlapping modification probabilities.
- Each target uses its own absolute probability: `1` is the selected
  modification, `0` is not that modification, and `.` is missing/uncertain.
  A `0` can represent another modification, not only unmodified cytosine.
  Missing selected calls are never inferred from other modification calls.
  Existing binary C/5mC results are unchanged.
- Selecting both targets runs separate analyses sequentially into `5mC/`
  and `5hmC/`, with a top-level summary. Publish only after both succeed.
  Selecting one target preserves the existing flat output layout.
- Record target, context, pattern-state meanings, and target-presence QC in
  run summaries. Check absent targets over the whole requested scope so a
  single chunk without the selected profile does not abort a valid genome run.
- Extraction Python APIs accept `modification="5hmC"` (one target per call).
  Rename modification-call QC counters to selected/not-selected terminology.
- Label 5hmC plots as hydroxymethylation and include the target in PDF captions.
  Detailed guidance is maintained alongside the software in `docs/`.

Metric formulas, TSV columns, N defaults, complete-pattern requirements, and
reporting-only `--min-depth` behavior are unchanged.

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
