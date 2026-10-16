# Find regions with large Epiallele distance

`discover` ranks fixed genomic windows or supplied BED regions using an existing paired
`three_metrics.tsv`. It does not require rereading the BAMs.

```bash
methyldiver discover \
  --input MD_results \
  --window-size 5000 \
  --min-units 30 \
  --min-ed 0.5 \
  --top 100 \
  --output MD_discovery_5kb
```

To rank promoter intervals, use the same merged BED that defined the promoter
analysis, or apply it to a genome-wide metric table:

```bash
methyldiver discover \
  --input MD_results \
  --regions-bed promoters.methyldiver.bed \
  --output MD_discovery_promoters
```

`--regions-bed` accepts BED3/4 with nonoverlapping intervals and unique region
IDs. BED rows may be unsorted; overlapping intervals must be merged first.
It is mutually exclusive with `--window-size`. Only fully contained CpG units
contribute. Exact BED boundaries and IDs are retained; `window_id` in the
output contains the BED region ID. Intervals with no usable observations are
omitted from statistics and counted in the summary. Units outside the BED
and units starting inside but crossing its end are counted separately.
The median ED, minimum-unit threshold, and ranking rules are unchanged.

`--input` accepts the TSV itself or its result directory. The output directory
must be new. Only paired results with Epiallele distance and delta columns are
accepted. Sample-specific column labels such as D0/D10 or COLO829BL/COLO829
are supported. The input must be grouped by chromosome, with strictly
increasing unit start positions inside each chromosome, as in MethylDiver
genome results. Unsorted or duplicated unit starts cause an error.

## Window definition and ranking

Windows are nonoverlapping BED0 intervals anchored at coordinate zero:
`[0,5000)`, `[5000,10000)`, and so on. A CpG unit contributes only when its
entire interval is within the window. Boundary-crossing units are excluded
and counted in the summary. When the original `analysis_summary.json` is
available, the terminal window is clipped to the recorded contig length.
Without those lengths, windows retain their full grid width.

Candidates need at least `--min-units` observed units (default 30) and a
**median unit ED** of at least `--min-ed` (default 0.5). Rank by median ED,
then mean ED, then unit count, all descending. Remaining ties follow input
window order. `--top` limits the ranked TSV and BED exports (default 100).
All evaluated windows, including those failing the filters, are saved
separately with a `passes_filters` column. Windows without usable observations
are absent, not filled with zeros.

The median and mean summarize existing per-unit ED values. They are not a
new ED calculated from pooled patterns over the whole 5 kb interval.
`mean_delta_beta` and `mean_delta_entropy` are unweighted means of unit-level
test-minus-control differences; they are not pooled-read regional estimates.

## Interpretation and coverage

This is **descriptive effect-size screening**, not a differential-methylation
significance test. There are no p-values or FDR values. ED is symmetric and
does not specify whether methylation increased or decreased. Adjacent CpG
units overlap and must not be treated as independent biological replicates.

`--min-units` is a count of CpG units, not a read-depth threshold. The compact
TSV lacks complete-read depths, so `discover` applies **no additional depth
filter**. New source TSVs are already filtered by the analysis `--min-depth`;
older reporting-only TSVs are not. The source export policy is recorded in
the discovery summary when available. For older runs, use retained
complete-pattern tables or reanalyze intervals to apply the desired depth threshold. A window selected
for high ED is an illustrative candidate, not an unbiased genome-wide effect
estimate or evidence of a gene's biological function.

## Files

| File | Contents |
| --- | --- |
| `ranked_regions.tsv` | Top qualifying windows, rank, coordinates, unit count, ED statistics, and mean deltas |
| `top_regions.bed` | BED4 coordinates of the same candidates, directly usable with `compare --regions-bed` |
| `window_statistics.tsv` | Every evaluated window, including those below the thresholds |
| `discovery_summary.json` | Input identity, parameters, source configuration when available, exclusion counts, interpretation, and code hash |

The input is streamed. Memory holds one window's observations plus the top
ranked candidates (and the interval index when a BED is supplied). Progress reports the input bytes scanned. Outputs are
published only on success; existing directories are not overwritten.

Python usage:

```python
from methyldiver.discover import discover_regions

discover_regions("MD_results/three_metrics.tsv", "MD_discovery_5kb",
                 window_size=5000, min_units=30, min_ed=0.5, top=100)
```

Discovery processes the saved TSV sequentially using one worker. It has no
`--threads` option; the eight-worker default applies to `compare` and `single`.
