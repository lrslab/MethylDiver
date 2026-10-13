"""Baseline definitions from the archived CK-AML epidec implementation.

See NOTICE and tests/fixtures/baseline.json for attribution and source hashes.
Counts are in binary pattern order (000 through 111 for k=3).
"""

import math
from collections.abc import Iterable
from numbers import Integral, Real
from typing import TypedDict

from .coordinates import _unit_size


def _count_unit_size(length):
    if length < 2 or length & (length - 1):
        raise ValueError("Supply exactly 2^N counts, including zero-count patterns (N from 1 to 10).")
    return _unit_size(length.bit_length() - 1)


class SampleMetrics(TypedDict):
    n_complete_reads: int
    beta: float | None
    entropy: float | None


class ComparisonMetrics(TypedDict):
    delta_beta: float | None
    delta_entropy: float | None
    epiallele_distance: float | None


def checked_counts(counts: Iterable[Real]) -> tuple[tuple[int, ...], int]:
    """Validate 2^N finite, nonnegative integer-valued counts (N from 1 to 10).

    Integral floats are accepted for compatibility with the original baseline;
    booleans, strings, missing values, and fractional counts are rejected.
    Values are normalized to Python integers to avoid fixed-width sum overflow.
    """
    try:
        values = tuple(counts)
    except TypeError as exc:
        raise ValueError("Counts must be an iterable of 2^N numbers (N from 1 to 10).") from exc
    k = _count_unit_size(len(values))
    result = []
    for i, value in enumerate(values):
        if isinstance(value, bool) or not isinstance(value, Real):
            raise ValueError(f"Count {i} must be a finite nonnegative integer, not {value!r}.")
        if isinstance(value, Integral):
            number = int(value)
        elif math.isfinite(value) and value == int(value):
            number = int(value)
        else:
            raise ValueError(f"Count {i} must be a finite nonnegative integer, not {value!r}.")
        if number < 0:
            raise ValueError(f"Count {i} must be nonnegative.")
        result.append(number)
    return tuple(result), k


def _metrics(counts: tuple[int, ...], k: int) -> SampleMetrics:
    n = sum(counts)
    if n == 0:
        return {"n_complete_reads": 0, "beta": None, "entropy": None}
    beta = sum(i.bit_count() * value for i, value in enumerate(counts)) / (k * n)
    probabilities = [value / n for value in counts]
    entropy = -math.fsum(p * math.log2(p) for p in probabilities if p > 0) / k
    return {"n_complete_reads": n, "beta": beta, "entropy": max(0.0, entropy)}


def metrics(counts: Iterable[Real]) -> SampleMetrics:
    """Return methylation level beta and Shannon pattern entropy divided by k.

    Both estimates use only the complete-read counts supplied. Empty samples
    return None, not a zero biological effect. No minimum-depth filter is applied.
    """
    values, k = checked_counts(counts)
    return _metrics(values, k)


def compare_counts(control: Iterable[Real], test: Iterable[Real]) -> ComparisonMetrics:
    """Return test-minus-control changes and square-root base-2 JSD.

    The two distributions receive equal weight regardless of read depth.
    The caller must match the reference and ordered CpG anchors before comparing.
    If either sample is empty, all comparison quantities are None.
    """
    a, k = checked_counts(control)
    b, other_k = checked_counts(test)
    if k != other_k:
        raise ValueError("Comparison requires identical CpG-unit size.")
    na, nb = sum(a), sum(b)
    if not na or not nb:
        return {"delta_beta": None, "delta_entropy": None, "epiallele_distance": None}
    am, bm = _metrics(a, k), _metrics(b, k)
    terms = []
    for ca, cb in zip(a, b):
        p, q = ca / na, cb / nb
        middle = (p + q) / 2
        if p > 0:
            terms.append(0.5 * p * math.log2(p / middle))
        if q > 0:
            terms.append(0.5 * q * math.log2(q / middle))
    jsd = math.fsum(terms)
    # Nonempty samples above guarantee numeric metrics.
    assert am["beta"] is not None and bm["beta"] is not None
    assert am["entropy"] is not None and bm["entropy"] is not None
    return {
        "delta_beta": bm["beta"] - am["beta"],
        "delta_entropy": bm["entropy"] - am["entropy"],
        "epiallele_distance": math.sqrt(max(0.0, min(1.0, jsd))),
    }
