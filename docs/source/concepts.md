# CpG units and metric definitions

MethylDiver calculates methylation level, normalized pattern entropy, and
epiallele distance directly from control and test modBAM files with `compare`.
For one sample, `single` calculates methylation level and entropy. Each unit contains
N consecutive reference CpGs, selected with `--N` (default **3**, supported **1–10**).

**All three metrics use complete read patterns only.** In a pattern, `1` means
the selected modification, `0` means not the selected modification, and `.` means
a missing or uncertain call. A `0` can represent another modification; it does
not necessarily mean unmodified cytosine. See [modification selection](modifications.md).
Any pattern containing `.` is excluded from methylation level, entropy, and
Epiallele distance calculations for that CpG unit. Missing calls are never
converted to `0` or `1`, and the observed sites of an incomplete pattern do not
contribute to these metrics.

| Example pattern | Treatment for that CpG unit |
| --- | --- |
| `101` or `100` | Complete; included in metric calculations |
| `1.0`, `.11`, or `...` | Incomplete; excluded from metrics and retained for quality checks |

Incomplete-pattern counts are saved in `incomplete_patterns/` for genome runs
when `--keep-patterns T` is selected,
and `patterns.tsv` for selected-region runs. Exclusion applies per unit: the
same read can contribute to other units where all CpG calls are complete.
The user-adjustable `--min-depth` threshold counts **complete reads in each
sample**. Compact metric TSVs include only units passing this threshold
(default 1). Both samples must pass for comparisons. Retained raw tables and
the in-memory extraction API preserve all units and their depth flags.
If a sample has no complete reads, its methylation level and entropy are missing;
Epiallele distance requires complete reads in both samples.

**Version:** 2.3.0 · **Python:** 3.10+ · **Input:** 5mC/5hmC long-read modBAM

See [CHANGELOG.md](https://github.com/lxd98/MethylDiver/blob/main/CHANGELOG.md) for release notes and changes from v1.0.0,
including `compare` replacing `analyze` and `--N` replacing `--k`.

| Quantity | Definition | Interpretation |
| --- | --- | --- |
| β (`beta`) | Σ p(x) × number of selected-modification sites in x / N | Mean selected-modification level across complete patterns |
| H (`entropy`) | −Σ p(x) log₂ p(x) / N | Normalized within-sample pattern diversity |
| Epiallele distance (`epiallele_distance`) | √JSD₂(P, Q) | Unsigned between-sample pattern-distribution distance |

All three quantities range from 0 to 1. `delta_beta` and `delta_entropy` are
**test minus control**. The two sample distributions have equal weight in JSD,
even when their sequencing depths differ. No pseudocount is added.

The output/API field for **Epiallele distance** is `epiallele_distance`.
Version 2.0.0 replaces the earlier `D_epi` field name (output schema 2.0);
the square-root Jensen–Shannon calculation and its values are unchanged.
When reading older output tables, rename `D_epi` to `epiallele_distance`.
The NumPy API still returns the same three arrays in the same order.
