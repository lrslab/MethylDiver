# Single-sample analysis

Select `--modification 5mC` (default), `5hmC`, or `5mC,5hmC` in the CLI.
See [modification selection](modifications.md) for CpG-specific state definitions,
separate result folders, and the Python `modification` keyword.

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
