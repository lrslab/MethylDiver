import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import unittest

from methyldiver import checked_counts, compare_counts, metrics


def load_baseline():
    path = Path(__file__).parent / "fixtures" / "legacy_metrics.py"
    spec = importlib.util.spec_from_file_location("legacy_metrics", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestMetrics(unittest.TestCase):
    def test_legacy_toy(self):
        a = [2, 0, 0, 0, 0, 2, 4, 12]
        b = [10, 4, 3, 0, 0, 0, 0, 3]
        self.assertAlmostEqual(metrics(a)["beta"], 0.8)
        self.assertAlmostEqual(metrics(a)["entropy"], 0.5236501981515562)
        self.assertAlmostEqual(metrics(b)["beta"], 0.2666666666666667)
        self.assertAlmostEqual(metrics(b)["entropy"], 0.5951584324091114)
        result = compare_counts(a, b)
        self.assertAlmostEqual(result["delta_beta"], -0.5333333333333333)
        self.assertAlmostEqual(result["delta_entropy"], 0.07150823425755515)
        self.assertAlmostEqual(result["D_epi"], 0.7309379165994421)

    def test_all_pure_patterns_and_uniform(self):
        for k in (3, 4):
            for pattern in range(2 ** k):
                counts = [0] * (2 ** k)
                counts[pattern] = 10
                self.assertEqual(metrics(counts)["beta"], pattern.bit_count() / k)
                self.assertEqual(metrics(counts)["entropy"], 0)
                self.assertEqual(compare_counts(counts, counts)["D_epi"], 0)
            self.assertEqual(metrics([1] * (2 ** k))["entropy"], 1)

    def test_hidden_redistribution(self):
        a = [10, 0, 0, 0, 0, 0, 0, 10]
        b = [0, 10, 0, 0, 0, 0, 10, 0]
        self.assertEqual(metrics(a)["beta"], 0.5)
        self.assertEqual(metrics(a)["entropy"], 1 / 3)
        self.assertEqual(compare_counts(a, b), {"delta_beta": 0, "delta_entropy": 0, "D_epi": 1})

    def test_empty_is_missing(self):
        for k in (3, 4):
            empty, full = [0] * (2 ** k), [1] * (2 ** k)
            self.assertEqual(metrics(empty), {"n_complete_reads": 0, "beta": None, "entropy": None})
            for a, b in ((empty, full), (full, empty), (empty, empty)):
                self.assertTrue(all(value is None for value in compare_counts(a, b).values()))

    def test_invalid_counts_and_sizes(self):
        for bad in (True, False, -1, 0.5, math.nan, math.inf, -math.inf, "1", None, 1j, [1]):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                metrics([bad] + [0] * 7)
        for counts in (None, 2, [], [1] * 7, [1] * 9, [1] * 32):
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                metrics(counts)
        with self.assertRaises(ValueError):
            compare_counts([1] * 8, [1] * 16)

    def test_integer_valued_floats_and_iterators(self):
        self.assertEqual(metrics(iter([1.0] * 8)), metrics([1] * 8))
        counts, k = checked_counts([1.0] * 8)
        self.assertEqual(k, 3)
        self.assertTrue(all(type(value) is int for value in counts))

    def test_large_integer_arithmetic(self):
        # Counts and methylated-bit sums must not wrap at int64 boundaries.
        counts = [0] * 7 + [2 ** 65]
        self.assertEqual(metrics(counts), {"n_complete_reads": 2 ** 65, "beta": 1.0, "entropy": 0.0})

    def test_randomized_baseline_and_invariants(self):
        baseline = load_baseline()
        rng = random.Random(20261012)
        for k in (3, 4):
            for _ in range(250):
                a = [rng.randrange(1000) for _ in range(2 ** k)]
                b = [rng.randrange(1000) for _ in range(2 ** k)]
                result, expected = compare_counts(a, b), baseline.compare_counts(a, b)
                reverse = compare_counts(b, a)
                scaled = compare_counts([x * 3 for x in a], [x * 7 for x in b])
                for key, value in result.items():
                    self.assertAlmostEqual(value, expected[key], delta=1e-12)
                    self.assertAlmostEqual(value, scaled[key], delta=1e-12)
                    self.assertAlmostEqual(value, reverse[key] * (1 if key == "D_epi" else -1), delta=1e-12)
                sample = metrics(a)
                for key, value in sample.items():
                    self.assertAlmostEqual(value, baseline.metrics(a)[key], delta=1e-12)
                self.assertTrue(0 <= sample["beta"] <= 1)
                self.assertTrue(0 <= sample["entropy"] <= 1)
                self.assertTrue(0 <= result["D_epi"] <= 1)
                # Relabel patterns together: H and distance are invariant, beta need not be.
                order = list(range(2 ** k))
                rng.shuffle(order)
                permuted = compare_counts([a[i] for i in order], [b[i] for i in order])
                self.assertAlmostEqual(result["D_epi"], permuted["D_epi"], delta=1e-12)
                self.assertAlmostEqual(result["delta_entropy"], permuted["delta_entropy"], delta=1e-12)

    def test_baseline_snapshot_fingerprint(self):
        root = Path(__file__).resolve().parents[1]
        provenance = json.loads((root / "tests/fixtures/baseline.json").read_text())
        source = root / "tests/fixtures/legacy_metrics.py"
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), provenance["baseline_metrics_sha256"])


if __name__ == "__main__":
    unittest.main()
