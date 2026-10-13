"""Optional NumPy batch calculations: each row is one matched CpG unit.

Install methyldiver[arrays] to use this module. Empty rows use NaN internally;
the scalar API and JSON CLI use None/null. Arrays do not carry reference keys:
the caller must match ordered units before calling compare_arrays.
"""

import numpy as np

from .metrics import _count_unit_size, checked_counts


def _checked_matrix(counts):
    raw = np.asarray(counts)
    if raw.ndim != 2:
        raise ValueError("Expected a two-dimensional count matrix with 2^N columns.")
    k = _count_unit_size(raw.shape[1])
    # Check original values before coercion so bools and strings cannot become counts.
    rows = [checked_counts(row)[0] for row in np.asarray(counts, dtype=object)]
    totals = [sum(row) for row in rows]
    if any(n > np.iinfo(np.int64).max for n in totals):
        raise ValueError("Batch complete-read totals must fit int64; use the scalar API for larger counts.")
    values = np.asarray(rows, dtype=np.float64).reshape(raw.shape)
    return values, np.asarray(totals, dtype=np.int64), k


def _calculate(values, totals, k):
    p = np.divide(values, totals[:, None], out=np.zeros_like(values), where=totals[:, None] > 0)
    logs = np.zeros_like(p)
    np.log2(p, out=logs, where=p > 0)
    entropy = -(p * logs).sum(axis=1) / k
    bits = np.array([i.bit_count() for i in range(2 ** k)], dtype=float)
    beta = (p @ bits) / k
    beta[totals == 0] = np.nan
    entropy[totals == 0] = np.nan
    return p, beta, entropy


def array_metrics(counts):
    """Return (n_complete_reads, beta, normalized_entropy) arrays."""
    values, totals, k = _checked_matrix(counts)
    _, beta, entropy = _calculate(values, totals, k)
    return totals, beta, entropy


def compare_arrays(control, test):
    """Return (delta_beta, delta_entropy, epiallele_distance), one value per matched row."""
    a, na, k = _checked_matrix(control)
    b, nb, j = _checked_matrix(test)
    if a.shape != b.shape or k != j:
        raise ValueError("Comparison matrices must have identical shapes and CpG-unit sizes.")
    p, ba, ha = _calculate(a, na, k)
    q, bb, hb = _calculate(b, nb, k)
    middle = (p + q) / 2
    terms = []
    for values in (p, q):
        ratio = np.divide(values, middle, out=np.ones_like(values), where=middle > 0)
        logs = np.zeros_like(values)
        np.log2(ratio, out=logs, where=values > 0)
        terms.append((values * logs).sum(axis=1))
    distance = np.sqrt(np.clip((terms[0] + terms[1]) / 2, 0, 1))
    distance[(na == 0) | (nb == 0)] = np.nan
    return bb - ba, hb - ha, distance
