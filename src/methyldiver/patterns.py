"""Complete binary counts and separate incomplete-pattern QC."""

from collections import Counter
from collections.abc import Iterable, Mapping
from numbers import Integral
from typing import TypedDict

from .coordinates import _positions, _unit_size


class PatternCounts(TypedDict):
    counts: list[int]
    n_candidate_reads: int
    n_complete_reads: int
    n_incomplete_reads: int
    callable_fraction: float | None
    incomplete_patterns: dict[str, int]


def read_pattern(calls: Mapping[int, int | None], unit: Iterable[int]) -> str:
    """Read states at fixed reference anchors, in genomic order.

    Calls are 0, 1, or None. An absent anchor also means missing ('.'). Do not
    remove a missing middle CpG or replace it with another observed site.
    """
    anchors = _positions(unit)
    _unit_size(len(anchors))
    pattern = []
    for position in anchors:
        value = calls.get(position)
        if value is None:
            pattern.append(".")
        elif isinstance(value, Integral) and not isinstance(value, bool) and value in (0, 1):
            pattern.append(str(int(value)))
        else:
            raise ValueError(f"Call at CpG {position} must be integer 0, 1, or None.")
    return "".join(pattern)


def count_patterns(patterns: Iterable[str], k: int = 3) -> PatternCounts:
    """Count one fixed unit's eligible candidate molecules.

    Strings contain exactly k symbols from '01.'. Incomplete patterns are
    retained for QC and never distributed over the complete count bins.
    The caller is responsible for alignment eligibility and molecule deduplication.
    """
    k = _unit_size(k)
    counts = [0] * (2 ** k)
    incomplete: Counter[str] = Counter()
    total = 0
    for pattern in patterns:
        if not isinstance(pattern, str) or len(pattern) != k or any(c not in "01." for c in pattern):
            raise ValueError(f"Expected a {k}-site pattern using only '0', '1', and '.', got {pattern!r}.")
        total += 1
        if "." in pattern:
            incomplete[pattern] += 1
        else:
            counts[int(pattern, 2)] += 1
    complete = sum(counts)
    return {
        "counts": counts,
        "n_candidate_reads": total,
        "n_complete_reads": complete,
        "n_incomplete_reads": total - complete,
        "callable_fraction": complete / total if total else None,
        "incomplete_patterns": dict(sorted(incomplete.items())),
    }
