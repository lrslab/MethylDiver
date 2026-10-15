# Tests and provenance

```bash
.venv/bin/python -m unittest discover -s tests -v
```

NumPy tests run when the optional dependency is installed; otherwise they are
explicitly skipped. Tests cover the legacy example, all pure patterns,
hidden redistribution, zero depth, invalid counts, signed direction,
depth-scaling invariance, 500 randomized comparisons against the original
scalar implementation, and 500 rows of scalar/batch agreement. Coordinate,
missingness, CLI, and output-preservation checks are included.

The additional modBAM tests cover four alignment/modification orientations,
MM/ML quantization and skipped calls, historical tag names, missing ML/MN checks,
quality filters, internal CIGAR gaps and clips, duplicate identities, interval
boundaries, unsupported states, resource limits, and complete CLI output.

With the genome extra installed, additional tests compare serial/parallel
genome runs to unchunked region analysis, including scan/chunk boundaries,
variable CpG-unit sizes, low-depth and empty units, typed empty outputs, and failure
cleanup. Genome tests are explicitly skipped when PyArrow is unavailable.

Single-sample tests verify agreement with paired sample metrics, extraction
without duplicate reads or distance calculations, incomplete patterns, zero-depth
and empty outputs, sample labels, serial/parallel genome runs, retention options,
reader exports, and cleanup after errors. CLI tests cover `compare`, `single`,
and `read` and reject the removed `analyze` and `metrics-counts` commands.

The unchanged scalar baseline is in `tests/fixtures/legacy_metrics.py`.
[`tests/fixtures/baseline.json`](https://github.com/lxd98/MethylDiver/blob/v2.0.0/tests/fixtures/baseline.json) records its origin
and fingerprint. See [`NOTICE`](https://github.com/lxd98/MethylDiver/blob/v2.0.0/NOTICE) for attribution.

Mixed-modification regression tests cover combined/separate and numeric tags,
both alignment/modification orientations, confidence thresholds, missing target
calls, invalid probability sums, known target-specific distributions, paired/
single agreement, serial/parallel chunks, and cleanup of failed two-target runs.

Selected-region tests check automatic batching, original region IDs, boundary
units, serial/parallel equivalence, and BED validation. Depth-filter regression
tests verify the cutoff in single/paired, genome/region, and retained/discarded
pattern modes, including re-export with `read --metrics-only`.

Discovery tests cover fixed windows and supplied BED regions, full containment,
ranking and ties, empty inputs, malformed or unsorted tables, input changes,
CLI validation, and output cleanup.
