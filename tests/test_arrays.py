import unittest

from methyldiver import compare_counts, metrics

try:
    import numpy as np
except ImportError:
    np = None
else:
    from methyldiver.arrays import array_metrics, compare_arrays


@unittest.skipIf(np is None, "Optional NumPy backend not installed")
class TestArrays(unittest.TestCase):
    def test_randomized_scalar_equivalence(self):
        rng = np.random.default_rng(20261012)
        for k in (1, 2, 3, 4, 5, 10):
            a = rng.integers(0, 1000, size=(25 if k == 10 else 250, 2 ** k), dtype=np.int64)
            b = rng.integers(0, 1000, size=a.shape, dtype=np.int64)
            a[0] = 0
            b[1] = 0
            a[2] = 0
            a[2, -1] = 10
            n, beta, entropy = array_metrics(a)
            db, dh, distance = compare_arrays(a, b)
            for i, (ar, br) in enumerate(zip(a, b)):
                expected = metrics(ar)
                self.assertEqual(n[i], expected["n_complete_reads"])
                for actual, key in ((beta[i], "beta"), (entropy[i], "entropy")):
                    if expected[key] is None:
                        self.assertTrue(np.isnan(actual))
                    else:
                        self.assertAlmostEqual(actual, expected[key], delta=1e-12)
                expected = compare_counts(ar, br)
                for actual, key in ((db[i], "delta_beta"), (dh[i], "delta_entropy"), (distance[i], "epiallele_distance")):
                    if expected[key] is None:
                        self.assertTrue(np.isnan(actual))
                    else:
                        self.assertAlmostEqual(actual, expected[key], delta=1e-12)

    def test_zero_rows_and_empty_samples(self):
        for k in (1, 2, 3, 4, 5, 10):
            empty = np.empty((0, 2 ** k), dtype=int)
            for result in (*array_metrics(empty), *compare_arrays(empty, empty)):
                self.assertEqual(result.shape, (0,))
            n, beta, entropy = array_metrics(np.zeros((3, 2 ** k), dtype=int))
            self.assertTrue(np.all(n == 0))
            self.assertTrue(np.isnan(beta).all())
            self.assertTrue(np.isnan(entropy).all())

    def test_invalid_and_mismatched_matrices(self):
        invalid = [[1] * 8, [[1] * 7], [[True] * 8], [["1"] * 8]]
        invalid += [[[bad] + [0] * 7] for bad in (-1, np.nan, np.inf, 0.5, 2 ** 64)]
        invalid += [[[True, 1, 1, 1, 1, 1, 1, 1]]]
        for value in invalid:
            with self.subTest(value=value), self.assertRaises(ValueError):
                array_metrics(value)
        for a, b in (([[1] * 8], [[1] * 16]), ([[1] * 8], [[1] * 8] * 2)):
            with self.assertRaises(ValueError):
                compare_arrays(a, b)

    def test_int64_sums_and_unsigned_scalar_counts(self):
        a = np.zeros((1, 8), dtype=np.int64)
        a[0, -1] = 2 ** 62
        n, beta, entropy = array_metrics(a)
        self.assertEqual(n[0], 2 ** 62)
        self.assertEqual(beta[0], 1)
        self.assertEqual(entropy[0], 0)
        values = np.full(8, 2 ** 62, dtype=np.uint64)
        self.assertEqual(metrics(values)["n_complete_reads"], 2 ** 65)
        with self.assertRaises(ValueError):
            array_metrics([values])


if __name__ == "__main__":
    unittest.main()
