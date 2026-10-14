from array import array
from collections import Counter
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

from methyldiver.metrics import compare_counts
from methyldiver.modbam import Region, analyze_modbam, analyze_single_modbam, decode_calls
from test_modbam import HEADER, SEQUENCE, make_read


def mixed_read(name, methyl, hydroxy, *, reverse=False, opposite=False, separate=False):
    m = [255 if x == "1" else 0 for x in methyl]
    h = [255 if x == "1" else 0 for x in hydroxy]
    if reverse:
        m, h = m[::-1], h[::-1]
    base = "G-" if opposite else "C+"
    mm = f"{base}m?,0,0,0;{base}h?,0,0,0;" if separate else f"{base}mh?,0,0,0;"
    ml = m + h if separate else [v for pair in zip(m, h) for v in pair]
    return make_read(name, reverse=reverse, opposite=opposite, mm=mm, ml=ml)


class TestModifications(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.reference = self.root / "ref.fa"
        self.reference.write_text(">chrToy\n" + SEQUENCE + "\n")
        pysam.faidx(str(self.reference))
        self.regions = [Region("chrToy", 0, len(SEQUENCE))]

    def bam(self, name, reads, header=HEADER):
        path = self.root / name
        with pysam.AlignmentFile(str(path), "wb", header=header) as out:
            for read in reads:
                out.write(read)
        pysam.index(str(path))
        return path

    def pair(self):
        control = self.bam("control.bam", [mixed_read("c1", "100", "010"), mixed_read("c2", "000", "000")])
        test = self.bam("test.bam", [mixed_read("t1", "100", "000"), mixed_read("t2", "100", "001")])
        return control, test

    def cli(self, command, output, *extra, control=None, test=None):
        inputs = ["--input", control] if command == "single" else ["--control", control, "--test", test]
        return subprocess.run([sys.executable, "-m", "methyldiver", command, *map(str, inputs),
                               "--reference", str(self.reference), "--output", str(output), *map(str, extra)],
                              text=True, capture_output=True)

    def test_combined_separate_and_all_strand_orientations(self):
        for separate in (False, True):
            for reverse in (False, True):
                for opposite in (False, True):
                    read = mixed_read("r", "100", "010", reverse=reverse, opposite=opposite, separate=separate)
                    for target, expected in (("5mC", [1, 0, 0]), ("5hmC", [0, 1, 0])):
                        with self.subTest(separate=separate, reverse=reverse, opposite=opposite, target=target):
                            calls, strand = decode_calls(read, {2, 6, 10}, .8, 10, Counter(), target)
                            self.assertEqual(calls, dict(zip([2, 6, 10], expected)))
                            self.assertEqual(strand, int(reverse ^ opposite))

    def test_numeric_aliases_historical_tags_and_unrelated_adenine(self):
        for legacy in (False, True):
            for reverse in (False, True):
                read = mixed_read("r", "100", "010", reverse=reverse, separate=True)
                adenine_group = "T-a?,0;" if reverse else "A+a?,0;"
                read.set_tag("MM", read.get_tag("MM").replace("+m?", "+27551?").replace("+h?", "+76792?") + adenine_group)
                read.set_tag("ML", array("B", list(read.get_tag("ML")) + [255]))
                if legacy:
                    mm, ml = read.get_tag("MM"), read.get_tag("ML")
                    read.set_tag("MM", None); read.set_tag("ML", None)
                    read.set_tag("Mm", mm); read.set_tag("Ml", ml)
                for target, expected in (("5mC", [1, 0, 0]), ("5hmC", [0, 1, 0])):
                    self.assertEqual(decode_calls(read, {2, 6, 10}, .8, 10, Counter(), target)[0],
                                     dict(zip([2, 6, 10], expected)))

    def test_selected_calls_are_never_inferred_from_other_modification(self):
        for mode in ("?", ".", ""):
            read = make_read(mm=f"C+m?,0,0,0;C+h{mode},0;", ml=[0, 255, 0, 255])
            calls, _ = decode_calls(read, {2, 6, 10}, .8, 10, Counter(), "5hmC")
            self.assertEqual(calls, {2: 1})
        read = mixed_read("no_ml", "100", "010")
        read.set_tag("ML", None)
        self.assertEqual(decode_calls(read, {2, 6, 10}, .8, 10, Counter(), "5hmC")[0],
                         {2: None, 6: None, 10: None})
        self.assertEqual(decode_calls(make_read(), {2, 6, 10}, .8, 10, Counter(), "5hmC")[0], {})

    def test_hydroxymethyl_confidence_quality_and_probability_validation(self):
        for byte, state in ((0, 0), (50, 0), (51, None), (204, None), (205, 1), (255, 1)):
            read = make_read(mm="C+mh?,0,0,0;", ml=[0, byte] * 3)
            self.assertEqual(decode_calls(read, {2, 6, 10}, .8, 10, Counter(), "5hmC")[0],
                             dict.fromkeys([2, 6, 10], state))
        read = mixed_read("low", "100", "010")
        qualities = list(read.query_qualities); qualities[6] = 9
        read.query_qualities = array("B", qualities)
        self.assertEqual(decode_calls(read, {2, 6, 10}, .8, 10, Counter(), "5hmC")[0], {2: 0, 6: None, 10: 0})
        for mm, ml in (("C+mh?,0;", [255, 255]), ("C+m?,0;C+h?,0;", [205, 100]),
                       ("C+h?,0;C+76792?,0;", [0, 0]), ("C+h?,0;G-h?,0;", [0, 0]),
                       ("C-h?,0;", [255]), ("C+f?,0;", [255])):
            with self.subTest(mm=mm), self.assertRaises(ValueError):
                decode_calls(make_read(mm=mm, ml=ml), {2, 6, 10}, .8, 10, Counter(), "5hmC")
        # Quantization-compatible edge sum is permitted.
        decode_calls(make_read(mm="C+mh?,0;", ml=[255, 1]), {2}, .8, 10, Counter())

    def test_known_distributions_and_single_paired_agreement(self):
        control, test = self.pair()
        expected = {"5mC": ([1, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 2, 0, 0, 0]),
                    "5hmC": ([1, 0, 1, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0, 0])}
        for target, (c, t) in expected.items():
            result = analyze_modbam(control, test, self.reference, self.regions, modification=target)
            unit = result["units"][0]
            self.assertEqual(unit["control"]["counts"], c)
            self.assertEqual(unit["test"]["counts"], t)
            self.assertEqual(unit["comparison"], compare_counts(c, t))
            single = analyze_single_modbam(control, self.reference, self.regions, modification=target)
            self.assertEqual(single["units"][0]["sample"], unit["control"])
            self.assertEqual(result["manifest"]["config"]["pattern_states"]["0"], "not_" + target)
            self.assertNotIn("calls_unmethylated", result["manifest"]["regions"][0]["control_qc"])
        default = analyze_modbam(control, test, self.reference, self.regions)
        explicit = analyze_modbam(control, test, self.reference, self.regions, modification="5mC")
        self.assertEqual(default["units"], explicit["units"])

    def test_absent_target_and_invalid_selection(self):
        bam = self.bam("only_m.bam", [make_read()])
        with self.assertRaisesRegex(ValueError, "no 5hmC modification profile"):
            analyze_single_modbam(bam, self.reference, self.regions, modification="5hmC")
        for value in ("6mA", "4mC", "5mC,5hmC", None):
            with self.assertRaisesRegex(ValueError, "modification must"):
                analyze_single_modbam(bam, self.reference, self.regions, modification=value)
        for value in ("6mA", "all", "5mC,5mC", "5mC,", "5hmCG"):
            out = self.root / "bad"
            run = self.cli("single", out, "--region", "chrToy:0-14", "--modification", value, control=bam)
            self.assertEqual(run.returncode, 2, run.stderr)
            self.assertFalse(out.exists())

    def test_region_multiple_targets_outputs_labels_reader_scope_and_cleanup(self):
        control, test = self.pair()
        for command in ("compare", "single"):
            output = self.root / command
            labels = ["--sample-label", "D0"] if command == "single" else ["--control-label", "D0", "--test-label", "D10"]
            run = self.cli(command, output, "--region", "chrToy:0-14", "--modification", "5mC,5hmC",
                           *labels, control=control, test=test)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertEqual(set(p.name for p in output.iterdir()), {"5mC", "5hmC", "analysis_summary.json"})
            self.assertIn(str(output / "5hmC"), run.stdout)
            for target in ("5mC", "5hmC"):
                summary = json.loads((output/target/"analysis_summary.json").read_text())
                self.assertEqual(summary["config"]["modification"], target)
                name = "two_metrics.tsv" if command == "single" else "three_metrics.tsv"
                rows = list(csv.DictReader((output/target/name).read_text().splitlines(), delimiter="\t"))
                self.assertIn("D0_beta", rows[0])
        bam = self.bam("only_m.bam", [make_read()])
        failed = self.root / "failed"
        run = self.cli("single", failed, "--region", "chrToy:0-14", "--modification", "5mC,5hmC", control=bam)
        self.assertEqual(run.returncode, 2)
        self.assertIn("no 5hmC", run.stderr)
        self.assertNotIn("analysis complete", run.stdout)
        self.assertFalse(failed.exists())
        self.assertFalse(list(self.root.glob(".failed.*")))
        protected = self.root / "single/5mC/two_metrics.tsv"
        before = protected.read_bytes()
        run = self.cli("single", self.root/"single", "--region", "chrToy:0-14", "--modification", "5mC,5hmC", control=bam)
        self.assertEqual(run.returncode, 2)
        self.assertEqual(protected.read_bytes(), before)

    @unittest.skipUnless(importlib.util.find_spec("pyarrow"), "requires pyarrow")
    def test_genome_parallel_matches_region_and_keeps_target_metadata(self):
        from methyldiver.genome import analyze_genome
        control, test = self.pair()
        for target in ("5mC", "5hmC"):
            for threads in (1, 2):
                out = self.root/f"{target}_{threads}"
                manifest = analyze_genome(control, test, self.reference, out, modification=target,
                                          threads=threads, chunk_units=1, keep_patterns=True, chunk_manifests=True)
                self.assertEqual(manifest["config"]["modification"], target)
                self.assertGreater(manifest["modification_qc"]["control"]["reads_with_selected_modification_profile"], 0)
                rows = list(csv.DictReader((out/"three_metrics.tsv").read_text().splitlines(), delimiter="\t"))
                region = analyze_modbam(control, test, self.reference, self.regions, modification=target)["units"][0]
                self.assertAlmostEqual(float(rows[0]["control_beta"]), region["control"]["beta"])
                self.assertAlmostEqual(float(rows[0]["epiallele_distance"]), region["comparison"]["epiallele_distance"])
        for command in ("compare", "single"):
            out = self.root/f"multi_{command}"
            run = self.cli(command, out, "--genome", "--threads", "2", "--modification", "5mC,5hmC",
                           control=control, test=test)
            self.assertEqual(run.returncode, 0, run.stderr)
            for target in ("5mC", "5hmC"):
                self.assertTrue((out/target/"chunks.tsv").exists())
                self.assertFalse((out/target/"complete_patterns").exists())

    @unittest.skipUnless(importlib.util.find_spec("pyarrow"), "requires pyarrow")
    def test_target_absent_in_one_chunk_does_not_abort_whole_genome(self):
        from methyldiver.genome import analyze_single_genome
        sequence = SEQUENCE + "AAAA" + SEQUENCE
        self.reference = self.root/"long.fa"
        self.reference.write_text(">chrToy\n" + sequence + "\n")
        pysam.faidx(str(self.reference))
        header = {"HD": HEADER["HD"], "SQ": [{"SN": "chrToy", "LN": len(sequence)}]}
        reads = [make_read("m", header=header),
                 make_read("h", start=18, mm="C+h?,0,0,0;", ml=[255, 0, 255], header=header)]
        bam = self.bam("split.bam", reads, header)
        for threads in (1, 2):
            manifest = analyze_single_genome(bam, self.reference, self.root/f"split{threads}", modification="5hmC",
                                             threads=threads, chunk_units=1)
            self.assertEqual(manifest["n_units"], 4)
            self.assertEqual(manifest["outputs"]["two_metrics.tsv"]["rows"], 1)
        bam = self.bam("absent.bam", [make_read("m", header=header)], header)
        out = self.root/"absent"
        with self.assertRaisesRegex(ValueError, "no 5hmC"):
            analyze_single_genome(bam, self.reference, out, modification="5hmC", threads=2, chunk_units=1)
        self.assertFalse(out.exists())
        self.assertFalse(list(self.root.glob(".absent.*")))


if __name__ == "__main__":
    unittest.main()
