"""Local epiallele metrics from complete three- or four-CpG patterns."""

from .metrics import checked_counts, compare_counts, metrics
from .patterns import count_patterns, read_pattern

__version__ = "1.0.0"
__all__ = [
    "__version__", "checked_counts", "metrics", "compare_counts",
    "count_patterns", "read_pattern",
]
