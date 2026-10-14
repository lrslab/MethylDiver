# Command reference

Help captured from MethylDiver **2.1.0**. Defaults appear in brackets.

## methyldiver

```text
Program: methyldiver (Methylation and epiallele diversity from long reads)
Version: 2.1.0

Usage: methyldiver <command> [options]

Commands:
  compare                Compare methylation, entropy, and Epiallele distance between two samples.
  single                 Calculate methylation level and entropy for one sample.
  read                   View or export retained pattern tables as TSV or CSV.

General options:
  -h, --help             show this help message and exit
  --version              show program's version number and exit

Run 'methyldiver <command> --help' for command-specific options.
Values in brackets show defaults.
```

## compare

```text
Compare methylation level, entropy, and Epiallele distance between two samples.

Usage: methyldiver compare [options] --control BAM --test BAM --reference FASTA --output DIR

Input files (required):
  --control BAM                  Indexed control modBAM.
  --test BAM                     Indexed test modBAM.
  --reference FASTA              Matching reference FASTA with .fai index.

Region selection (choose one):
  --genome                       All BAM reference contigs; requires the [genome] extra [off].
  --region CHR:START-END         0-based, end-exclusive interval; repeat for separate regions
                                 [unset].
  --regions-bed FILE             Nonoverlapping BED3/BED4 intervals; 0-based, end-exclusive [unset].

Analysis options:
  --modification STR             CpG target: 5mC, 5hmC, or 5mC,5hmC for separate results [5mC].
  --N INT                        Consecutive CpGs per unit, 1–10 [3].
  --call-threshold FLOAT         Binary-call confidence, greater than 0.5 and at most 1 [0.8].
  --min-mapq INT                 Minimum read mapping quality, 0–254; MAPQ 255 is excluded [20].
  --min-base-quality INT         Minimum CpG-call base quality, 0–254; missing qualities excluded
                                 [10].
  --min-depth INT                Complete-read depth in BOTH samples for the pass flag; does not
                                 filter estimates [1].

Output options:
  --output DIR                   New output directory (required; must not exist).
  --control-label STR            Control label in three_metrics.tsv [control].
  --test-label STR               Test label in three_metrics.tsv [test].

Genome options (--genome only):
  --threads INT                  Worker processes; selected regions run serially [1].
  --chunk-units INT              Requested units per chunk; reduced for large N; boundary units
                                 preserved [10000].
  --chunk-manifests BOOL         Save per-chunk metadata as .json.gz; False omits it [False].
  --keep-patterns BOOL           Keep complete_patterns/ and incomplete_patterns/; False removes
                                 them after TSV export [False].

Resource limits (exceeding a limit stops analysis):
  --max-total-bases INT          Total interval span in selected-region mode only [100000].
  --max-reads-per-region INT     Fetched reads before filtering, per sample per region/chunk
                                 [100000].
  --max-pattern-evaluations INT  Candidate read/unit patterns per sample per region/chunk [5000000].

General options:
  -h, --help                     show this help message and exit

Defaults are shown in brackets. BOOL accepts True/False or T/F.
Only complete read patterns contribute to metric estimates.
```

## single

```text
Calculate methylation level and entropy for one sample.

Usage: methyldiver single [options] --input BAM --reference FASTA --output DIR

Input files (required):
  --input BAM                    Indexed sample modBAM.
  --reference FASTA              Matching reference FASTA with .fai index.

Region selection (choose one):
  --genome                       All BAM reference contigs; requires the [genome] extra [off].
  --region CHR:START-END         0-based, end-exclusive interval; repeat for separate regions
                                 [unset].
  --regions-bed FILE             Nonoverlapping BED3/BED4 intervals; 0-based, end-exclusive [unset].

Analysis options:
  --modification STR             CpG target: 5mC, 5hmC, or 5mC,5hmC for separate results [5mC].
  --N INT                        Consecutive CpGs per unit, 1–10 [3].
  --call-threshold FLOAT         Binary-call confidence, greater than 0.5 and at most 1 [0.8].
  --min-mapq INT                 Minimum read mapping quality, 0–254; MAPQ 255 is excluded [20].
  --min-base-quality INT         Minimum CpG-call base quality, 0–254; missing qualities excluded
                                 [10].
  --min-depth INT                Complete-read depth in the sample for the pass flag; does not
                                 filter estimates [1].

Output options:
  --output DIR                   New output directory (required; must not exist).
  --sample-label STR             Sample label in two_metrics.tsv [sample].

Genome options (--genome only):
  --threads INT                  Worker processes; selected regions run serially [1].
  --chunk-units INT              Requested units per chunk; reduced for large N; boundary units
                                 preserved [10000].
  --chunk-manifests BOOL         Save per-chunk metadata as .json.gz; False omits it [False].
  --keep-patterns BOOL           Keep complete_patterns/ and incomplete_patterns/; False removes
                                 them after TSV export [False].

Resource limits (exceeding a limit stops analysis):
  --max-total-bases INT          Total interval span in selected-region mode only [100000].
  --max-reads-per-region INT     Fetched reads before filtering, per sample per region/chunk
                                 [100000].
  --max-pattern-evaluations INT  Candidate read/unit patterns per sample per region/chunk [5000000].

General options:
  -h, --help                     show this help message and exit

Defaults are shown in brackets. BOOL accepts True/False or T/F.
Only complete read patterns contribute to metric estimates.
```

## read

```text
View or export retained MethylDiver pattern tables.

Usage: methyldiver read [options] --input PATH

Input/output options:
  --input PATH         complete_patterns/, incomplete_patterns/, or one .parquet file (required).
  --output FILE        New output file; .csv selects CSV [stdout].
  --format STR         csv or tsv [.csv output: csv; otherwise: tsv].

Selection options:
  --limit INT          Maximum matching rows to display/export [20].
  --all                Export all matching rows, streamed in batches [off; limit: 20].
  --chrom STR          Exact chromosome/contig name [all].
  --sample STR         Incomplete-pattern role: control, test, or sample [all].
  --pattern STR        Exact incomplete pattern, e.g. '1.0'; dot means missing [all].
  --metrics-only       Plotting columns only (7 single, 12 paired); omit missing metrics [off: all
                       columns and rows].

Sample labels (with --metrics-only):
  --control-label STR  Paired control label [saved label, or control].
  --test-label STR     Paired test label [saved label, or test].
  --sample-label STR   Single-sample label [saved label, or sample].

General options:
  -h, --help           show this help message and exit

Defaults are shown in brackets. Parquet input requires the [genome] extra.
```
