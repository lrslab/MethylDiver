# Compare two samples genome-wide

Select `--modification 5mC` (default), `5hmC`, or `5mC,5hmC` in the CLI.
See [modification selection](modifications.md) for CpG-specific state definitions,
separate result folders, and the Python `modification` keyword.

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
this only sets `passes_min_depth`. Valid estimates are retained in the plotting TSV. Add `--keep-patterns T`
to retain complete-read counts and apply another depth threshold later without
extracting BAMs again.
No threshold is a universal guarantee of statistical reliability.

## Output files and retention

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
analysis tables. `analysis_summary.json` records the selected boolean,
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

## Progress and completion

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
