"""Reference CpG anchors: zero-based C coordinates on the forward reference."""

from collections.abc import Iterable, Iterator
from numbers import Integral


def _coordinate(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer.")
    return int(value)


def _unit_size(k: int) -> int:
    if isinstance(k, bool) or not isinstance(k, Integral) or k not in (3, 4):
        raise ValueError("Supported CpG-unit sizes are 3 and 4.")
    return int(k)


def _positions(positions: Iterable[int]) -> tuple[int, ...]:
    values = tuple(_coordinate(p, "CpG position") for p in positions)
    if any(b <= a for a, b in zip(values, values[1:])):
        raise ValueError("CpG positions must be sorted and unique.")
    return values


def cpg_positions(sequence: str, offset: int = 0) -> list[int]:
    """Enumerate CG motifs in a reference segment, with a genomic offset."""
    offset = _coordinate(offset, "offset")
    sequence = sequence.upper()
    return [offset + i for i in range(len(sequence) - 1) if sequence[i:i + 2] == "CG"]


def canonical_cpg(sequence: str, position: int, reference_strand: str, offset: int = 0) -> int:
    """Map a modified-C reference position to a forward-reference CpG anchor.

    For a minus-strand modification at reference G, subtract one. The strand
    is the modification's reference strand, not the BAM alignment flag.
    """
    position = _coordinate(position, "position")
    offset = _coordinate(offset, "offset")
    if reference_strand not in ("+", "-"):
        raise ValueError("Modification reference strand must be '+' or '-'.")
    anchor = position - (reference_strand == "-")
    local = anchor - offset
    if local < 0 or sequence[local:local + 2].upper() != "CG":
        raise ValueError("Modification does not map to a reference CpG.")
    return anchor


def sliding_units(positions: Iterable[int], k: int = 3) -> Iterator[tuple[int, ...]]:
    """Yield consecutive reference CpGs from one contig, with stride one.

    Pass the full reference CpG catalogue for the region, not called sites from
    one read. There is no maximum physical span. Call separately per contig.
    """
    k = _unit_size(k)
    values = _positions(positions)
    for i in range(len(values) - k + 1):
        yield values[i:i + k]


def unit_interval(unit: Iterable[int]) -> tuple[int, int]:
    """Return the BED0 half-open interval containing all CpG dyads in a unit."""
    values = _positions(unit)
    _unit_size(len(values))
    return values[0], values[-1] + 2
