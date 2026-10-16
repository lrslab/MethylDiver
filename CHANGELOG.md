# Changelog

## 2.3.1 - 2026-10-16

- Default to eight worker processes for `compare`, `single`, and the streaming
  genome/region Python APIs. Explicit `--threads` / `threads` values override it.
- Default selected-region CLI runs now use chunked processing and require the
  `[genome]` extra. Small explicit `--threads 1` requests retain the serial path.
- Update help, installation instructions and documentation for the public
  GitHub repository. Metrics and extraction filters are unchanged.

## 2.3.0 - 2026-10-15

- Extend `discover` with `--regions-bed` to rank nonoverlapping promoter or
  other BED3/4 regions using fully contained units and the existing ED rules.
  Retain region IDs and exact boundaries; validate overlap, IDs, and input changes.

- Make `--min-depth` filter `two_metrics.tsv` and `three_metrics.tsv` in
  single, paired, genome, BED and region analyses, including `--keep-patterns F`.
  Require the threshold in the single sample or in both comparison samples.
  The default remains 1; pass `--min-depth 5` to require five complete reads.
- Preserve all units in retained full count/QC tables and the in-memory API.
  `read --metrics-only` now also requires the stored depth flag, including
  when re-exporting older retained results.
- Record the export depth policy and separate missing-metric exclusions from
  otherwise valid units excluded for low depth. Existing results are unchanged;
  rerun analysis or re-export retained tables to produce a filtered TSV.

- Add `methyldiver discover` to stream paired metric TSVs into fixed genomic
  windows and rank high median ED candidates. Defaults: 5 kb, at least 30
  observed units, median ED >=0.5, and the top 100 candidates.
- Export all window statistics, ranked candidates, a BED4 for reanalysis,
  and a discovery summary. Exclude and count boundary-crossing units;
  reject unsorted input, duplicate unit starts, and malformed metrics.
- This is descriptive screening without p-values, FDR, or read-depth
  filtering; coverage must be checked before interpreting candidates.

## 2.2.0 - 2026-10-14

- Automatically batch large `--regions-bed` and repeated `--region` requests;
  users can supply a full promoter BED without splitting it externally.
- Enable `--threads`, `--chunk-units`, `--keep-patterns`, and
  `--chunk-manifests` for selected regions in both `compare` and `single`.
- Plan chunks from reference CpG counts and genomic span. Preserve original
  region IDs, every boundary-spanning unit, and separation between input regions.
- Treat CLI `--max-total-bases` as a per-batch span target (default 100,000).
  A single sparse unit can exceed the target to preserve its complete span.
- Stream chunk results to disk with bounded worker queues and the genome-style
  output layout. Small default serial runs keep their existing TSV layout.
- Add the streaming Python API `methyldiver.genome.analyze_regions`;
  the in-memory `analyze_modbam` API retains its existing safety limits.

## 2.1.1 - 2026-10-14

- Handle empty MM tags with a single unsigned zero ML byte, as observed in
  older Bonito/Remora COLO829 BAMs, as missing modification evidence. Count
  these reads in QC; never convert the placeholder to an unmethylated call.
- Preserve validation of all other ML length/type mismatches and stale tags.
- Add regression checks for both tag spellings, orientations, targets,
  complete-pattern exclusion, and samples containing only placeholders.

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
