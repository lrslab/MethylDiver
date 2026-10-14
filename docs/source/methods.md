# Call filtering and scientific interpretation

## Modification-call policy

The decoder supports 5mC and 5hmC at reference CpGs, separately or together
in one modBAM. Use `--modification 5mC` (default), `5hmC`, or `5mC,5hmC`.
Both alignment orientations, `C+m`/`C+h` and `G-m`/`G-h` annotations,
combined `C+mh` groups, numeric ChEBI codes, and modern or historical tags
are supported. See [modification selection](modifications.md) for encodings,
probability validation, and output layouts.

Explicit ML bytes are decoded as `(byte + 0.5) / 256`. At the default threshold,
state 1 requires the selected modification's probability ≥0.8 and state 0
requires it ≤0.2. Intermediate calls remain missing. State 0 means **not the
selected modification**, which may include another cytosine modification.
No probabilities are combined or renormalized between classes.

Skipped `?` sites, skipped implicit `.`/omitted-mode sites, missing tags, and
missing probabilities **remain missing**. In particular, implicit low-modification
evidence is not assigned an exact numerical probability. This explicit-only
policy differs from historical modkit workflows that accept inferred canonical
calls; the selected policy is recorded in the run manifest.

Only primary mapped, nonduplicate, non-QC-failed reads passing mapping quality
are eligible. A candidate must span the unit's first and last modified positions;
internal CIGAR gaps or failed calls produce incomplete patterns. Reverse-strand
cytosines at reference G are anchored to the preceding reference C. Full CpG
dyads must lie within the requested interval. Calls are joined to reference CG
coordinates; a CpG-disrupting variant mask is not applied.

Malformed tags, stale MN lengths, duplicate read identities, repeated MI molecule
identities, unsupported cytosine modifications beyond 5mC/5hmC, dual-strand cytosine
annotations, paired-end records, and marked ONT duplex/parent records are rejected.
No cross-platform or arbitrary multi-modification support is implied by successful
BAM reading. BAMs containing only ordinary alignments without modification tags
fail with a diagnostic when eligible reads are present.

Reference validation checks selected contig names and lengths. It records hashes
of the requested reference slices and index, but does not prove whole-reference
sequence identity from the BAM header. You must supply the actual matching
assembly. BAM/FASTA fingerprints are path/size/mtime plus index/header hashes,
not full large-file content hashes. Input paths and BAM header metadata appear
in the local manifest.

These rules follow the [SAM modification-tag specification](https://samtools.github.io/hts-specs/SAMtags.pdf)
and the [pysam modification API](https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.modified_bases_forward).

## Scientific scope

These are established complete-read plug-in estimators. This release provides
no confidence intervals,
P values, missing-pattern imputation, or biological-validation claim.

Complete-read estimates describe molecules passing the input call filters.
Pattern-dependent call failure can bias their relationship to all molecules.
Low depth depresses empirical entropy and can produce positive Epiallele distance even
under identical underlying distributions. Choose depth thresholds appropriate for the experiment; no universal
reliability cutoff is applied.

Reads and overlapping units are not independent biological replicates.
Matching site-average bedMethyl values cannot reconstruct molecular patterns.

The modBAM command performs bounded extraction and call filtering. CRAM support,
cohort sample manifests, resume/checkpoint reuse, joint multistate
pattern models, non-CpG contexts, and resampling remain future stages.
