"""Local epiallele metrics from complete N-CpG patterns (default N=3)."""

from .metrics import checked_counts, compare_counts, metrics
from .patterns import count_patterns, read_pattern

__version__ = "2.1.0"
__all__ = [
    "__version__", "checked_counts", "metrics", "compare_counts",
    "count_patterns", "read_pattern",
]
