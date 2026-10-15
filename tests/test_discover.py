import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from methyldiver.discover import discover_regions


HEADER = ["region_id", "chrom", "start", "end", "cpg_anchors", "D0_beta", "D10_beta",
          "delta_beta", "D0_entropy", "D10_entropy", "delta_entropy", "epiallele_distance"]


class TestDiscover(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "three_metrics.tsv"

    def write(self, rows):
        with self.source.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER)
            for chrom, start, end, ed in rows:
                writer.writerow(["part_1", chrom, start, end, f"{start},{start+2},{end-2}",
                                 .2, .6, .4, .3, .1, -.2, ed])

    def read(self, output, name):
        return list(csv.DictReader((output / name).read_text().splitlines(), delimiter="\t"))

    def test_boundaries_filters_ranking_and_summary(self):
        self.write([("chr1", 0, 6, .2), ("chr1", 100, 5000, .8),
                    ("chr1", 4998, 5004, 1), ("chr1", 5000, 5006, .9),
                    ("chr1", 6000, 6006, .9), ("chr1", 10000, 10006, 1),
                    ("chr2", 0, 6, .7), ("chr2", 10, 16, "NA")])
        out = self.root / "results"
        events = []
        summary = discover_regions(self.source, out, min_units=2, min_ed=.5, top=1,
                                   progress=lambda *args: events.append(args))
        rows = self.read(out, "ranked_regions.tsv")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["window_id"], "chr1:5000-10000")
        self.assertEqual(float(rows[0]["median_ed"]), .9)
        self.assertEqual(float(rows[0]["mean_delta_beta"]), .4)
        self.assertEqual(float(rows[0]["mean_delta_entropy"]), -.2)
        all_rows = self.read(out, "window_statistics.tsv")
        self.assertEqual(len(all_rows), 4)
        self.assertEqual(float(all_rows[0]["median_ed"]), .5)
        self.assertEqual(float(all_rows[0]["fraction_units_ge_threshold"]), .5)
        self.assertEqual(summary["windows_passing_filters"], 2)
        self.assertEqual(summary["boundary_crossing_rows"], 1)
        self.assertEqual(summary["missing_metric_rows"], 1)
        self.assertEqual(summary["included_units"], 6)
        self.assertEqual((out / "top_regions.bed").read_text(), "chr1\t5000\t10000\tchr1:5000-10000\n")
        self.assertEqual(events[0][0], 0)
        self.assertEqual(events[-1][0], self.source.stat().st_size)
        self.assertTrue(summary["depth_filter"].startswith("none"))

    def test_terminal_window_clipping_and_ties(self):
        self.write([("chr1", 0, 6, 1), ("chr1", 5000, 5006, 1), ("chr2", 5000, 5100, 1)])
        (self.root / "analysis_summary.json").write_text(json.dumps({
            "contigs": [{"chrom": "chr1", "length": 8000}, {"chrom": "chr2", "length": 5100}],
            "sample_labels": {"control": "D0", "test": "D10"}, "config": {"min_depth": 5}}))
        out = self.root / "ties"
        summary = discover_regions(self.root, out, min_units=1, top=2)
        rows = self.read(out, "ranked_regions.tsv")
        self.assertEqual([r["window_id"] for r in rows], ["chr1:0-5000", "chr1:5000-8000"])
        self.assertEqual(self.read(out, "window_statistics.tsv")[-1]["end"], "5100")
        self.assertEqual(summary["sample_labels"]["test"], "D10")
        self.assertEqual(summary["source_config"]["min_depth"], 5)

    def test_rejects_unsorted_duplicates_invalid_values_and_preserves_output(self):
        for number, rows in enumerate([
            [("chr1", 10, 20, .5), ("chr1", 0, 6, .5)],
            [("chr1", 0, 6, .5), ("chr1", 0, 7, .5)],
            [("chr1", 0, 6, .5), ("chr2", 0, 6, .5), ("chr1", 10, 16, .5)],
            [("chr1", -1, 6, .5)], [("chr1", 0, 0, .5)],
            [("chr1", 0, 6, 1.1)], [("chr1", 0, 6, "inf")], [("chr1", 0, 6, "bad")],
        ]):
            self.write(rows)
            out = self.root / f"bad_{number}"
            with self.assertRaises(ValueError):
                discover_regions(self.source, out)
            self.assertFalse(out.exists())
            self.assertFalse(list(self.root.glob(f".{out.name}.*")))
        protected = self.root / "protected"
        protected.mkdir()
        (protected / "keep").write_text("original")
        with self.assertRaises(FileExistsError):
            discover_regions(self.source, protected)
        self.assertEqual((protected / "keep").read_text(), "original")

    def test_empty_input_parameter_validation_and_input_change(self):
        self.write([])
        out = self.root / "empty"
        summary = discover_regions(self.source, out)
        self.assertEqual(summary["n_ranked_regions"], 0)
        self.assertEqual(self.read(out, "ranked_regions.tsv"), [])
        self.assertEqual((out / "top_regions.bed").read_text(), "")
        for kwargs in ({"window_size": 0}, {"min_units": 0}, {"top": 0},
                       {"min_ed": float("nan")}, {"min_ed": 1.1}, {"min_ed": True}):
            with self.assertRaises(ValueError):
                discover_regions(self.source, self.root / "invalid", **kwargs)
        self.write([("chr1", 0, 6, .8)])
        def mutate(done, total, rows):
            if done:
                with self.source.open("a") as handle:
                    handle.write("changed\n")
        with self.assertRaisesRegex(ValueError, "changed"):
            discover_regions(self.source, self.root / "changed", progress=mutate)
        self.assertFalse((self.root / "changed").exists())

    def test_cli_help_success_and_single_sample_rejection(self):
        self.write([("chr1", 0, 6, .9), ("chr1", 100, 110, .8)])
        output = self.root / "cli"
        base = [sys.executable, "-m", "methyldiver", "discover"]
        help_run = subprocess.run(base + ["--help"], capture_output=True, text=True)
        self.assertIn("--window-size", help_run.stdout)
        self.assertIn("not a significance test", help_run.stdout)
        run = subprocess.run(base + ["--input", str(self.source), "--output", str(output),
                             "--min-units", "2"], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn("discovery complete", run.stdout)
        self.assertEqual(len(self.read(output, "ranked_regions.tsv")), 1)
        self.source.write_text("region_id\tchrom\tstart\tend\tsample_beta\tsample_entropy\n")
        run = subprocess.run(base + ["--input", str(self.source), "--output", str(self.root / "single")],
                             capture_output=True, text=True)
        self.assertEqual(run.returncode, 2)
        self.assertIn("paired three_metrics.tsv", run.stderr)
        self.assertFalse((self.root / "single").exists())

    def test_bed_regions_exact_boundaries_and_unobserved_regions(self):
        self.write([("chr1", 5, 11, .9), ("chr1", 10, 20, .6),
                    ("chr1", 19, 21, .9), ("chr1", 20, 25, .8),
                    ("chr1", 30, 36, .9), ("chr2", 0, 6, .8)])
        bed = self.root / "promoters.bed"
        bed.write_text("chr1\t20\t30\tP2\nchr1\t10\t20\tP1\nchr3\t0\t100\tP3\n")
        out = self.root / "bed"
        summary = discover_regions(self.source, out, regions_bed=bed, min_units=1)
        rows = self.read(out, "window_statistics.tsv")
        self.assertEqual([(r['window_id'],r['start'],r['end'],r['n_units']) for r in rows],
                         [('P1','10','20','1'),('P2','20','30','1')])
        self.assertEqual([float(r['median_ed']) for r in rows], [.6,.8])
        self.assertEqual(summary['outside_regions_rows'],3)
        self.assertEqual(summary['boundary_crossing_rows'],1)
        self.assertEqual(summary['regions_without_usable_units'],1)
        self.assertEqual(summary['regions_requested'],3)
        self.assertEqual((out/'top_regions.bed').read_text(),"chr1\t20\t30\tP2\nchr1\t10\t20\tP1\n")
        run = subprocess.run([sys.executable,'-m','methyldiver','discover','--input',str(self.source),
                              '--output',str(self.root/'cli_bed'),'--regions-bed',str(bed),'--min-units','1'],
                             capture_output=True,text=True)
        self.assertEqual(run.returncode,0,run.stderr)
        self.assertEqual(self.read(self.root/'cli_bed','window_statistics.tsv'),rows)

    def test_bed_invalid_overlap_duplicate_ids_and_change(self):
        self.write([('chr1',10,16,.9)])
        bed=self.root/'regions.bed'
        for value in ['', 'chr1\t-1\t20\n', 'chr1\t0\t20\tx\nchr1\t15\t30\ty\n',
                      'chr1\t0\t20\tx\nchr1\t20\t30\tx\n']:
            bed.write_text(value)
            with self.assertRaises(ValueError):
                discover_regions(self.source,self.root/'bad_bed',regions_bed=bed)
            self.assertFalse((self.root/'bad_bed').exists())
        bed.write_text('chr1\t0\t20\tx\n')
        def mutate(done,total,rows):
            if done:bed.write_text('chr1\t0\t21\tx\n')
        with self.assertRaisesRegex(ValueError,'BED changed'):
            discover_regions(self.source,self.root/'changed_bed',regions_bed=bed,progress=mutate)
        self.assertFalse((self.root/'changed_bed').exists())


if __name__ == "__main__":
    unittest.main()
