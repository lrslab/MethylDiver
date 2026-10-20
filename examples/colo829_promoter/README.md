# Real-data example: COLO829 promoter

Run MethylDiver on real long-read modification calls from COLO829_1 (melanoma)
and COLO829BL_1 (matched normal cell line), covering the merged
**HOXC13;HOXC13-AS promoter**. Dataset suffixes are file labels, not biological
pairing identifiers. This is a selected demonstration region from the manuscript,
not an independent validation dataset.

## Try it

From the MethylDiver software repository directory, with MethylDiver installed:

```bash
source .venv/bin/activate
bash examples/colo829_promoter/run_example.sh MD_COLO829_example
python examples/colo829_promoter/verify_results.py MD_COLO829_example
```

The script runs `single` for each sample, `compare`, then promoter `discover`.
It produces metric tables, regional discovery tables and analysis logs; no plots
are generated, and R is not required.
The four analyses completed in a few seconds on the development workstation;
performance on other machines will vary. Use a new output directory each time.
The supplied `expected/` files are reference results, separate from your output.
Nothing is uploaded by these commands.

If MethylDiver is already available on your PATH, activation is unnecessary.
Alternatively, set `METHYLDIVER_BIN=/absolute/path/to/.venv/bin/methyldiver`
when calling the shell script. The verifier uses only Python's standard library.

## Run the commands separately

From the software repository root:

```bash
example=examples/colo829_promoter

methyldiver compare \
  --control "$example/inputs/COLO829BL_1.bam" \
  --test "$example/inputs/COLO829_1.bam" \
  --reference "$example/inputs/chr12.fa.gz" \
  --regions-bed "$example/inputs/promoter.bed" \
  --min-depth 5 \
  --output MD_promoter_compare

methyldiver single \
  --input "$example/inputs/COLO829_1.bam" \
  --reference "$example/inputs/chr12.fa.gz" \
  --regions-bed "$example/inputs/promoter.bed" \
  --min-depth 5 \
  --output MD_promoter_single

methyldiver discover \
  --input MD_promoter_compare \
  --regions-bed "$example/inputs/promoter.bed" \
  --output MD_promoter_discovery
```

For the normal sample, replace the single-sample input with `COLO829BL_1.bam`
and choose another output directory. The verifier expects the folder names
created by `run_example.sh`.

Only the nondefault depth threshold is specified. Defaults are N = 3, 5mC,
call threshold 0.8, MAPQ >=20, base quality >=10, eight workers (one chunk here),
no retained patterns and no chunk manifests. Discovery defaults to at least
30 units and median ED >=0.5. It does not apply another read-depth filter.

## Inputs and coordinates

The package is approximately **42 MB**, including both BAMs and the complete,
unmasked chr12 reference. No additional download or decompression is needed.

| File | Contents |
| --- | --- |
| `inputs/COLO829BL_1.bam` and `.bam.bai` | 87 original alignments overlapping the promoter, before analysis filters |
| `inputs/COLO829_1.bam` and `.bam.bai` | 86 original alignments overlapping the promoter, before analysis filters |
| `inputs/chr12.fa.gz`, `.fa.gz.fai`, `.fa.gz.gzi` | Complete chr12 sequence from GRCh38_no_alt.fna, BGZF compressed and indexed |
| `inputs/promoter.bed` | One merged promoter: chr12, 53936830–53941692, BED0 |
| `inputs/input_manifest.json` | Input sizes, SHA-256 checksums and extraction validation |

BED coordinates are **0-based, end-exclusive**. In a browser using 1-based
coordinates, the interval is `chr12:53,936,831-53,941,692`.
All alignments overlapping the full promoter were retained, including reads
subsequently rejected by MethylDiver's filters. Read sequences, alignments,
modification tags and original headers were preserved; reads were not trimmed,
reoriented, grouped or selected by ED. No coordinates were shifted. Alignment
counts above are not complete-read depths for every CpG unit.

Use the supplied BED or that exact `--region` interval. The BAMs contain only
regional reads, even though their original headers list the full reference;
this package is not a genome-wide input dataset. The full chr12 reference
preserves chromosome lengths and avoids synthetic padding or altered sequence.

## Expected results (MethylDiver 2.3.1)

There are 198 candidate three-CpG units in the promoter. At depth >=5:

| Analysis | Main result | Rows |
| --- | --- | ---: |
| Normal sample | `single_COLO829BL_1/two_metrics.tsv` | 192 |
| Melanoma sample | `single_COLO829_1/two_metrics.tsv` | 195 |
| Comparison | `compare/three_metrics.tsv` | 191 |
| Promoter discovery | `discover_promoter/ranked_regions.tsv` | 1 |

The comparison uses COLO829BL_1 as control and COLO829_1 as test; deltas are
test minus control. The regional median ED is **0.624061**, mean ED is
**0.623464**, mean unit delta methylation is **−0.000275**, and mean unit delta
entropy is **−0.054367**. These summarize the 191 eligible units without pooling
read patterns. Discovery is descriptive and supplies no P value or FDR.

The manuscript's example uses **115 units shared across all ten datasets**.
This two-sample tutorial uses all **191 units eligible in these two samples**,
so its regional summaries are not the manuscript's all-ten-shared summaries.

Expected metric tables and discovery outputs are in `expected/`. The verifier
checks input checksums, table columns, row counts, coordinates and values
(numeric tolerance: relative 1e-9, absolute 1e-10). Paths and elapsed times in
JSON run logs are intentionally not compared. The comparison table was also
verified byte-for-byte against a fresh run on the original full BAMs and full
reference. Single-sample values match comparison values on their common units.

## Sources and rebuilding

The BAMs come from Oxford Nanopore Technologies' public **rrms_2022.07** release:
[dataset description](https://epi2me.nanoporetech.com/rrms2022.07/).
Exact BAM, index, reference and GENCODE annotation download URLs are in
[sources.tsv](sources.tsv); checksums there refer to original local source files
(the reference and annotation after decompression), not the subsets.
No large source download is needed to run this packaged example.

The promoter was defined from GENCODE v49 gene-level TSSs, 2,000 bp upstream
and 500 bp downstream, strand-aware. Overlapping/adjacent promoter intervals
were merged irrespective of strand. The two original BED0 intervals were
`chr12:53936830-53939330` (HOXC13, +, ENSG00000123364.5) and
`chr12:53939192-53941692` (HOXC13-AS, −, ENSG00000249641.3).
The semicolon denotes two genes assigned to the merged region.

To recreate the inputs from already downloaded, indexed source files:

```bash
python examples/colo829_promoter/build_inputs.py \
  --bam-directory /path/to/original_bams \
  --reference /path/to/GRCh38_no_alt.fna \
  --output rebuilt_inputs
```

The source directory must contain `COLO829_1.bam`, `COLO829BL_1.bam` and their
indexes. The builder checks every extracted alignment against its source and
every reference base against the original chr12 sequence. Compressed-file
checksums can vary with compression-library versions even if content matches.
