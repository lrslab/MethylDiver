# Install

Clone the private repository using your authorized GitHub account:

```bash
git clone https://github.com/lxd98/MethylDiver.git
cd MethylDiver
python3 -m venv .venv
.venv/bin/python -m pip install '.[genome]'
.venv/bin/methyldiver --version
```

Activate the environment to use the command directly, like other command-line
bioinformatics programs:

```bash
source .venv/bin/activate
methyldiver --help
methyldiver compare --help
methyldiver single --help
```

Help uses a samtools-style layout: program/version, concise usage, grouped
commands/options, and defaults in brackets (for example, `--min-mapq INT ... [20]`).
Running `methyldiver` alone shows the command list; running `methyldiver compare`,
`methyldiver single`, or `methyldiver read` alone shows that command's help.

Python 3.10 or later is required. The package uses `pysam` to read modBAMs and
indexed FASTA references. The scalar metric module itself uses only Python's
standard library. No modkit executable is needed for normal analysis.
The optional NumPy batch API can be installed with:

```bash
.venv/bin/python -m pip install -e '.[arrays]'
```

The package is installed from the repository root. For development, use
`pip install -e .` inside the virtual environment. Releases are identified by Git
tags, starting with `v1.0.0`.

To run from any directory without activating `.venv`, create a user-level command
link once, from the repository directory after installation:

```bash
mkdir -p "$HOME/.local/bin"
ln -s "$PWD/.venv/bin/methyldiver" "$HOME/.local/bin/methyldiver"
export PATH="$HOME/.local/bin:$PATH"
methyldiver --help
```

Keep `~/.local/bin` on your shell's PATH (add the export to `~/.bashrc` if needed).
The command uses the installed virtual environment's Python automatically, so it
also works from Conda `base`. Keep the repository and `.venv` at their installed
paths; moving or deleting them invalidates the link. This does not install
MethylDiver's Python dependencies into `base`.


Shell commands throughout this guide assume the **MethylDiver software repository** is the working directory. The `docs/` folder contains the guide; installation commands run from the
software repository root, one directory above `docs/`.
