import unittest

from methyldiver import count_patterns, metrics, read_pattern
from methyldiver.coordinates import canonical_cpg, cpg_positions, sliding_units, unit_interval


class TestPatterns(unittest.TestCase):
    def test_binary_order_supported_sizes(self):
        for k in (1, 2, 3, 4, 5):
            patterns = [format(i, f"0{k}b") for i in range(2 ** k) for _ in range(i)]
            result = count_patterns(patterns, k)
            self.assertEqual(result["counts"], list(range(2 ** k)))
            self.assertEqual(result["callable_fraction"], 1)

    def test_missing_middle_is_not_bridged(self):
        calls = {10: 1, 30: 0, 40: 1}
        self.assertEqual(read_pattern(calls, (10, 20, 30)), "1.0")
        result = count_patterns(["111", "1.0", "000", "101", "1.0", "..."])
        self.assertEqual(result["counts"], [1, 0, 0, 0, 0, 1, 0, 1])
        self.assertEqual(result["n_candidate_reads"], 6)
        self.assertEqual(result["n_complete_reads"], 3)
        self.assertEqual(result["n_incomplete_reads"], 3)
        self.assertEqual(result["callable_fraction"], 0.5)
        self.assertEqual(result["incomplete_patterns"], {"...": 1, "1.0": 2})
        self.assertAlmostEqual(metrics(result["counts"])["beta"], 5 / 9)

    def test_empty_and_all_incomplete(self):
        self.assertIsNone(count_patterns([])["callable_fraction"])
        result = count_patterns(["...", ".1."])
        self.assertEqual(result["callable_fraction"], 0)
        self.assertIsNone(metrics(result["counts"])["entropy"])
        self.assertEqual(read_pattern({1: None}, (1, 2, 3)), "...")

    def test_invalid_patterns_calls_and_units(self):
        for pattern in ("", "00", "0000", "0?1", "012", None):
            with self.subTest(pattern=pattern), self.assertRaises(ValueError):
                count_patterns([pattern])
        for value in (-1, 2, True, False, 0.0, "1"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                read_pattern({1: value}, (1, 2, 3))
        for unit in ((1, 1, 3), (3, 2, 1), (-1, 2, 3), ()):
            with self.assertRaises(ValueError):
                read_pattern({}, unit)
        for k in (0, 11, True, 3.0):
            with self.assertRaises(ValueError):
                count_patterns([], k)


class TestCoordinates(unittest.TestCase):
    def test_strand_and_segment_offset(self):
        self.assertEqual(canonical_cpg("AACGTT", 102, "+", 100), 102)
        self.assertEqual(canonical_cpg("AACGTT", 103, "-", 100), 102)
        self.assertEqual(canonical_cpg("cg", 1, "-"), 0)
        for position, strand in ((103, "+"), (101, "-"), (106, "+"), (102, "?")):
            with self.assertRaises(ValueError):
                canonical_cpg("AACGTT", position, strand, 100)

    def test_catalogue_units_and_full_dyad_interval(self):
        positions = cpg_positions("cgacgtttcgacg", offset=100)
        self.assertEqual(positions, [100, 103, 108, 111])
        self.assertEqual(list(sliding_units(positions)), [(100, 103, 108), (103, 108, 111)])
        self.assertEqual(list(sliding_units(iter(positions), 4)), [(100, 103, 108, 111)])
        self.assertEqual(unit_interval((100, 103, 108)), (100, 110))
        self.assertEqual(list(sliding_units([1, 10000, 100000])), [(1, 10000, 100000)])

    def test_empty_short_and_invalid_catalogues(self):
        self.assertEqual(cpg_positions(""), [])
        self.assertEqual(cpg_positions("ATAT"), [])
        self.assertEqual(list(sliding_units([])), [])
        self.assertEqual(list(sliding_units([1, 2])), [])
        for values in ([1, 1, 2], [3, 2, 1], [-1, 2, 3], [True, 2, 3], [1.0, 2, 3]):
            with self.assertRaises(ValueError):
                list(sliding_units(values))
        with self.assertRaises(ValueError):
            cpg_positions("CG", -1)


if __name__ == "__main__":
    unittest.main()
