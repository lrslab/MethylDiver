# Select 5mC or 5hmC at CpGs

MethylDiver 2.3.0 supports **5mCG and 5hmCG**: 5mC and 5hmC at reference CpG
sites. Both `compare` and `single` accept the following selections:

| Option | Analysis |
| --- | --- |
| `--modification 5mC` | 5mC at CpGs; default when the option is omitted |
| `--modification 5hmC` | 5hmC at CpGs |
| `--modification 5mC,5hmC` | Separate 5mC and 5hmC analyses |

The reference units remain N consecutive CpGs, with `--N 3` by default.
This release does not analyze non-CpG cytosines, 4mC, or 6mA. Unrelated
adenine modification groups may coexist in the BAM and are ignored after tag
validation. Unsupported cytosine modification profiles are rejected.

## Examples

Compare both modifications between the same samples:

```bash
methyldiver compare \
  --control /path/to/control.bam \
  --test /path/to/test.bam \
  --reference /path/to/reference.fa \
  --genome \
  --threads 24 \
  --modification 5mC,5hmC \
  --output MD_results
```

Analyze hydroxymethylation in one sample:

```bash
methyldiver single \
  --input /path/to/sample.bam \
  --reference /path/to/reference.fa \
  --genome \
  --threads 24 \
  --modification 5hmC \
  --output MD_5hmC
```

## Binary pattern states

Each target is analyzed independently, using its **absolute probability** in
the ML tag. For a selected modification X and confidence threshold t:

| State | Rule | Meaning |
| --- | --- | --- |
| `1` | p(X) ≥ t | Confidently the selected modification |
| `0` | p(X) ≤ 1 − t | Confidently not the selected modification |
| `.` | Intermediate probability, absent selected call, or failed quality filter | Missing or uncertain |

The default threshold is 0.8. Explicit ML bytes are interpreted as
`(byte + 0.5) / 256`, as in earlier versions. No renormalization between 5mC
and 5hmC is performed, and their probabilities are not added together.

For example, a site with high 5hmC probability and explicitly low 5mC probability
is `1` in the 5hmC analysis and `0` in the 5mC analysis. **Here `0` does not
mean unmodified cytosine.** For a binary C/5mC input, the existing default
5mC results are unchanged.

An absent selected-modification entry is never reconstructed from the other
modification. Skipped `?`, `.`, and omitted-mode sites remain missing, even
when the other modification has a confident call. Only complete N-site
patterns enter the metrics, independently for each selected target.

Beta is the selected-modification fraction across complete patterns. Thus,
for 5hmC it measures hydroxymethylation level. Entropy and Epiallele distance
use that target's binary pattern distribution. Table column names remain
`beta`, `entropy`, and `epiallele_distance` (with their usual sample prefixes).
The target and the meaning of all three pattern states are recorded in
`analysis_summary.json` under `config`.

Because complete-call availability can differ, the two target analyses may use
different reads or units. Their beta values are not guaranteed to sum to one.
Use retained depth/QC tables when comparing the coverage of the two analyses.

## Supported tags and checks

Supported encodings include combined `C+mh` groups, separate `C+m` and `C+h`
groups, their `G-m` / `G-h` counterparts, and ChEBI codes `27551` / `76792`.
Both alignment orientations and modern `MM`/`ML` or historical `Mm`/`Ml`
tags are handled. Numeric modification codes use separate groups, as required
by the [SAM modification-tag specification](https://samtools.github.io/hts-specs/SAMtags.pdf).

Duplicate groups, including duplicate letter/numeric aliases, malformed tags,
and dual-strand cytosine annotations are rejected. At overlapping explicit
positions, the 5mC and 5hmC ML byte sum must be at most 256; higher sums are
incompatible with alternative modification probabilities. Selection does not
bypass tag validation.

One narrow compatibility exception handles an encoding observed in older
Bonito/Remora COLO829 BAMs: an empty `MM` string accompanied by exactly one
unsigned zero byte in `ML`. With no listed modification site, this byte
cannot represent an unmethylated call. MethylDiver treats the read as missing
modification evidence (`.`), counts it under
`reads_empty_mm_zero_ml_placeholder` in the modification QC, and excludes
its incomplete patterns from the metrics. Historical `Mm`/`Ml` tags receive
the same handling. Other MM/ML length or type mismatches still fail validation.

If a sample has eligible reads but no selected-modification profile anywhere
in the requested analyzable scope, the run fails with a diagnostic. Absence in
one region or chunk alone does not abort a run that has the target elsewhere.
Empty BAMs and scopes without analyzable units retain the existing empty-output
behavior. A present but empty or entirely uncertain target profile can yield no
complete metric rows; this is not evidence that the modification level is zero.

## Outputs and performance

Selecting one target preserves the usual output layout directly under
`--output`. Selecting both creates:

```text
MD_results/
├── analysis_summary.json
├── 5mC/
│   ├── three_metrics.tsv
│   ├── chunks.tsv
│   └── analysis_summary.json
└── 5hmC/
    ├── three_metrics.tsv
    ├── chunks.tsv
    └── analysis_summary.json
```

Single-sample runs use `two_metrics.tsv`. Small default serial region runs use
`metrics.tsv` and `patterns.tsv` instead of `chunks.tsv`; large or parallel
region runs use the chunked layout shown above. Retention flags apply
independently to both target folders. The top-level summary links the two
analyses; each child summary records its own parameters, inputs, metrics, and QC.

The two analyses run sequentially, rereading the inputs. `--threads` applies to
each genome analysis; selecting both does not double the worker count. This
costs more time than analyzing one target. The output is published only after
both succeed; errors remove staged results and existing output is never
overwritten. Target presence/QC counts in genome summaries count read fetches
across chunks and are not deduplicated genome-wide molecule counts.

Point the reader and plotting script at a specific target:

```bash
# Requires --keep-patterns T during analysis.
methyldiver read --input MD_results/5hmC/complete_patterns \
  --metrics-only --all --output MD_results/5hmC/export.tsv

Rscript examples/plot_three_metrics.R \
  MD_results/5hmC/three_metrics.tsv MD_results/5hmC/figures
```

The R example reads the selected target from the adjacent child summary and
labels 5hmC panels as hydroxymethylation. Older summaries without this field
are interpreted as 5mC. `--min-depth` filters each target's metric TSV
independently using its own complete-read depths.

## Python API

The extraction APIs accept **one target per call**:

```python
from methyldiver.modbam import Region, analyze_single_modbam

result = analyze_single_modbam(
    "sample.bam", "reference.fa", [Region("chr1", 1000, 2000)],
    modification="5hmC", k=3,
)
```

The same keyword works in `analyze_modbam`, `analyze_genome`, and
`analyze_single_genome`. Omitting it selects `5mC`. Scalar and NumPy metric
functions still operate on binary count vectors and need no modification flag.
