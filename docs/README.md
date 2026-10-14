# MethylDiver documentation

User documentation for **MethylDiver 2.1.0**, a tool for methylation level,
pattern entropy, and Epiallele distance from long-read modBAMs.

[Software and quick start](https://github.com/lxd98/MethylDiver) ·
[Release notes](https://github.com/lxd98/MethylDiver/blob/main/CHANGELOG.md)

## Read the documentation

| Guide | Contents |
| --- | --- |
| [Installation](source/installation.md) | Dependencies, installation, and running without environment activation |
| [CpG units and metrics](source/concepts.md) | Complete patterns, methylation, entropy, and Epiallele distance |
| [Modification selection](source/modifications.md) | 5mCG/5hmCG, mixed calls, pattern states, and separate outputs |
| [Paired analysis](source/compare.md) | Genome analysis, sample labels, outputs, retention, and progress |
| [Single-sample analysis](source/single.md) | Methylation and entropy from one modBAM |
| [Selected regions](source/regions.md) | BED/interval input, coordinates, and region outputs |
| [Command reference](source/commands.md) | All commands, options, and default values |
| [Reading results](source/reading.md) | Parquet previews, filters, and TSV/CSV export |
| [Examples and plotting](source/examples.md) | Synthetic inputs and individual ggplot2 PDF panels |
| [Filtering and interpretation](source/methods.md) | Modification calls, read eligibility, and scientific limits |
| [Python API](source/python-api.md) | Extraction, count vectors, coordinates, and batch calculations |
| [Tests and provenance](source/development.md) | Regression checks and attribution |

## Build HTML locally

The documentation builder requires Python 3.12 or later. MethylDiver itself
supports Python 3.10 or later.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r docs/requirements.txt
.venv/bin/sphinx-build -n -W --keep-going -b html docs/source docs/_build/html
```

Open `docs/_build/html/index.html` in a browser. Run the build commands from the MethylDiver repository root. The documentation
uses Sphinx with the Read the Docs theme. The Markdown chapters
are also readable directly on GitHub.

Documentation is maintained in the software repository and uses the same
private repository access. No public documentation website is published.
