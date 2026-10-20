# MethylDiver

**MethylDiver** is designed for long-read methylation analysis, calculating
methylation level, normalized pattern entropy, and Epiallele distance from
molecular CpG patterns in modBAM files.

**Version 2.3.1** · **Python 3.10+** · [Documentation](docs/README.md) · [Release notes](CHANGELOG.md)

| Command | Function |
| --- | --- |
| `compare` | Compare methylation, entropy, and epiallele distributions between two samples |
| `single` | Calculate methylation and entropy for one sample |
| `read` | Preview or export retained Parquet results |
| `discover` | Rank high-ED genomic windows or BED regions from paired results |

## Installation

Clone the repository:

```bash
git clone https://github.com/lxd98/MethylDiver.git
cd MethylDiver
python3 -m venv .venv
source .venv/bin/activate
python -m pip install '.[genome]'
methyldiver --help
```

To run without activating the environment, see the
[installation guide](docs/source/installation.md).

## Quick usage

For a runnable real-data example, see the
[COLO829 promoter tutorial](examples/colo829_promoter/README.md). It includes
two regional modBAMs, an indexed reference, commands and verified expected outputs.

Inputs are coordinate-sorted, indexed **modBAMs containing 5mC and/or 5hmC calls**
and the matching reference FASTA with its `.fai` index.

Compare two samples:

```bash
methyldiver compare \
  --control /path/to/control.bam \
  --test /path/to/test.bam \
  --reference /path/to/reference.fa \
  --genome \
  --min-depth 5 \
  --output MD_results
```

Analyze one sample:

```bash
methyldiver single \
  --input /path/to/sample.bam \
  --reference /path/to/reference.fa \
  --genome \
  --output MD_single_results
```

Both commands use **8 worker processes** by default; override with `--threads INT`.
Fewer workers are used when there are fewer chunks.

Each unit contains **3 consecutive CpGs** by default; use `--N` to choose 1–10.
The default target is **5mC at CpGs**. Use `--modification 5hmC` for 5hmC at CpGs,
or `--modification 5mC,5hmC` to analyze both separately. Both commands support
these choices; [modification selection](docs/source/modifications.md)
explains mixed calls and binary pattern states.

For selected regions, replace `--genome` with
`--region chr1:1000-2000` or `--regions-bed regions.bed`; `--threads 24`
works for these modes too. Large BEDs are automatically processed in chunks.
Coordinates are
**zero-based, end-exclusive**. The output directory must be new.

## Output

To find high-ED windows from paired results:

```bash
methyldiver discover --input MD_results --window-size 5000 \
  --min-units 30 --min-ed 0.5 --output MD_discovery_5kb
```

This ranks median per-unit ED, without significance testing or read-depth
filtering. See [region discovery](docs/source/discover.md).

The main result is `three_metrics.tsv` for `compare` or `two_metrics.tsv` for
`single`, with one row per CpG unit passing the depth filter and having all required metrics. Chunked runs
also produce `chunks.tsv` and `analysis_summary.json`.
When both modifications are selected, each gets its own `5mC/` or `5hmC/`
subfolder, with a top-level `analysis_summary.json` linking the two analyses.

Only complete read patterns contribute to metrics; patterns containing `.`
are excluded. **`--min-depth` filters the saved metric TSV** (default 1):
`single` requires the threshold in its sample, and `compare` requires it in
**both** samples. This works with `--keep-patterns F`. Add `--keep-patterns T`
to retain all counts, depths and QC, including below-threshold units, for
subsequent re-filtering. Existing TSVs are not changed automatically.
For each target, `1` means that modification and `0` means not that modification;
`0` does not necessarily mean unmodified cytosine.

## Documentation

The [full documentation](docs/README.md)
covers all options and defaults, input requirements, metric definitions,
outputs, filtering, the Python API, and worked examples.

- [Synthetic example and PDF plotting](docs/source/examples.md)
- [Command reference](docs/source/commands.md)
- [Call filtering and interpretation](docs/source/methods.md)

Use `methyldiver compare --help`, `methyldiver single --help`,
`methyldiver read --help`, or `methyldiver discover --help` for command-specific help.
