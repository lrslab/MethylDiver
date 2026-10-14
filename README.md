# MethylDiver

**MethylDiver** is designed for long-read methylation analysis, calculating
methylation level, normalized pattern entropy, and Epiallele distance from
molecular CpG patterns in modBAM files.

**Version 2.0.0** · **Python 3.10+** · [Documentation](https://github.com/lxd98/MethylDiver_documentation) · [Release notes](CHANGELOG.md)

| Command | Function |
| --- | --- |
| `compare` | Compare methylation, entropy, and epiallele distributions between two samples |
| `single` | Calculate methylation and entropy for one sample |
| `read` | Preview or export retained Parquet results |

## Installation

Clone the private repository with an authorized GitHub account:

```bash
git clone https://github.com/lxd98/MethylDiver.git
cd MethylDiver
python3 -m venv .venv
source .venv/bin/activate
python -m pip install '.[genome]'
methyldiver --help
```

To run without activating the environment, see the
[installation guide](https://github.com/lxd98/MethylDiver_documentation/blob/main/docs/source/installation.md).

## Quick usage

Inputs are coordinate-sorted, indexed **modBAMs containing binary C/5mC calls**
and the matching reference FASTA with its `.fai` index.

Compare two samples:

```bash
methyldiver compare \
  --control /path/to/control.bam \
  --test /path/to/test.bam \
  --reference /path/to/reference.fa \
  --genome \
  --threads 24 \
  --min-depth 5 \
  --output MD_results
```

Analyze one sample:

```bash
methyldiver single \
  --input /path/to/sample.bam \
  --reference /path/to/reference.fa \
  --genome \
  --threads 24 \
  --output MD_single_results
```

Each unit contains **3 consecutive CpGs** by default; use `--N` to choose 1–10.
For selected regions, replace `--genome --threads 24` with
`--region chr1:1000-2000` or `--regions-bed regions.bed`. Coordinates are
**zero-based, end-exclusive**. The output directory must be new.

## Output

The main result is `three_metrics.tsv` for `compare` or `two_metrics.tsv` for
`single`, with one row per CpG unit having all required metrics. Genome runs
also produce `chunks.tsv` and `analysis_summary.json`.

Only complete read patterns contribute to metrics; patterns containing `.`
are excluded. **`--min-depth` reports passing units but does not filter the
saved TSV.** Add `--keep-patterns T` to retain read depths, QC, and pattern
counts for subsequent filtering in genome runs; these details are removed by default.

## Documentation

The [full documentation](https://github.com/lxd98/MethylDiver_documentation)
covers all options and defaults, input requirements, metric definitions,
outputs, filtering, the Python API, and worked examples.

- [Synthetic example and PDF plotting](https://github.com/lxd98/MethylDiver_documentation/blob/main/docs/source/examples.md)
- [Command reference](https://github.com/lxd98/MethylDiver_documentation/blob/main/docs/source/commands.md)
- [Call filtering and interpretation](https://github.com/lxd98/MethylDiver_documentation/blob/main/docs/source/methods.md)

Use `methyldiver compare --help`, `methyldiver single --help`, or
`methyldiver read --help` for command-specific help.
