import csv
import importlib.util
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from methyldiver.reader import read_results


@unittest.skipUnless(importlib.util.find_spec("pyarrow"), "optional pyarrow is not installed")
class TestReader(unittest.TestCase):
    def setUp(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "incomplete_patterns"
        self.source.mkdir()
        self.rows = [
            dict(region_id="part_000000", chrom="chr1", cpg_anchors=[2, 5, 9], sample="control", pattern="1.0", count=3),
            dict(region_id="part_000000", chrom="chr1", cpg_anchors=[2, 5, 9], sample="test", pattern="1.0", count=2),
            dict(region_id="part_000001", chrom="chr2", cpg_anchors=[4, 8, 12], sample="control", pattern=".11", count=1),
        ]
        pq.write_table(pa.Table.from_pylist(self.rows[:2]), self.source / "part_000000.parquet")
        pq.write_table(pa.Table.from_pylist(self.rows[2:]), self.source / "part_000001.parquet")

    def run_cli(self, *args):
        return subprocess.run([sys.executable, "-m", "methyldiver", "read", *map(str, args)],
                              text=True, capture_output=True)

    def parse(self, text):
        return list(csv.DictReader(io.StringIO(text), delimiter="\t"))

    def test_preview_and_combined_filters(self):
        run = self.run_cli("--input", self.source, "--limit", 2)
        self.assertEqual(run.returncode, 0, run.stderr)
        rows = self.parse(run.stdout)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["cpg_anchors"], "2,5,9")
        self.assertIn("limit 2", run.stderr)
        run = self.run_cli("--input", self.source, "--chrom", "chr1", "--sample", "control", "--pattern", "1.0")
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual([r["count"] for r in self.parse(run.stdout)], ["3"])

    def test_full_export_single_file_and_empty_results(self):
        output = self.root / "export.tsv"
        run = self.run_cli("--input", self.source, "--all", "--output", output)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(run.stdout, "")
        self.assertEqual(len(self.parse(output.read_text())), 3)
        run = self.run_cli("--input", self.source / "part_000001.parquet", "--all")
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(self.parse(run.stdout)[0]["pattern"], ".11")
        # Historical directory names remain readable without conversion.
        legacy = self.root / "CpG_unit_pattern"
        legacy.symlink_to(self.source, target_is_directory=True)
        run = self.run_cli("--input", legacy, "--all")
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(len(self.parse(run.stdout)), 3)
        run = self.run_cli("--input", self.source, "--chrom", "absent")
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(len(run.stdout.splitlines()), 1)
        self.assertEqual(self.parse(run.stdout), [])

    def test_batch_boundary_default_limit_and_metrics_nulls(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        large = self.root / "metrics.parquet"
        pq.write_table(pa.Table.from_pylist([
            dict(region_id="part", chrom="chr1", cpg_anchors=[i, i+2, i+4], epiallele_distance=None)
            for i in range(4200)]), large)
        preview = self.root / "preview.tsv"
        self.assertEqual(read_results(large, output=preview), 20)
        self.assertEqual(len(self.parse(preview.read_text())), 20)
        full = self.root / "full.tsv"
        self.assertEqual(read_results(large, output=full, limit=None), 4200)
        self.assertEqual(self.parse(full.read_text())[-1]["epiallele_distance"], "NA")
        limited = self.root / "limited.tsv"
        self.assertEqual(read_results(large, output=limited, limit=4100), 4100)
        run = self.run_cli("--input", large, "--sample", "control")
        self.assertEqual(run.returncode, 2)
        self.assertIn("requires a 'sample' column", run.stderr)

    def test_invalid_inputs_and_output_preservation(self):
        output = self.root / "existing.tsv"
        output.write_text("keep")
        cases = [
            ("--input", self.source, "--output", output),
            ("--input", self.source, "--output", self.source / "export.tsv"),
            ("--input", self.root / "absent"),
            ("--input", self.source, "--limit", 0),
            ("--input", self.source, "--all", "--limit", 1),
            ("--input", self.source, "--control-label", "DMSO"),
            ("--input", self.root),
        ]
        for args in cases:
            with self.subTest(args=args):
                run = self.run_cli(*args)
                self.assertEqual(run.returncode, 2, run.stderr)
                self.assertNotIn("Traceback", run.stderr)
        self.assertEqual(output.read_text(), "keep")
        self.assertFalse((self.source / "export.tsv").exists())

    def test_failed_export_cleans_up_partial_file(self):
        output = self.root / "failed.tsv"
        with patch("methyldiver.reader.csv.writer") as factory:
            factory.return_value.writerow.side_effect = [None, OSError("write failed")]
            with self.assertRaisesRegex(OSError, "write failed"):
                read_results(self.source, output=output)
        self.assertFalse(output.exists())
        self.assertEqual(list(self.root.glob(".failed.tsv.*.tmp")), [])

    def test_csv_export_quotes_lists_and_selects_metrics(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        from methyldiver.output import METRIC_SUMMARY_COLUMNS, METRIC_VALUE_COLUMNS
        output = self.root / "patterns.csv"
        run = self.run_cli("--input", self.source, "--all", "--output", output)
        self.assertEqual(run.returncode, 0, run.stderr)
        with output.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0]["cpg_anchors"], "2,5,9")
        self.assertEqual(rows[0]["pattern"], "1.0")
        self.assertIn('"2,5,9"', output.read_text())
        metrics = self.root / "metrics.parquet"
        row = dict(region_id="part", chrom="chr1", start=2, end=11, cpg_anchors=[2, 5, 9], unit_size=3,
                   control_n_complete_reads=0, test_n_complete_reads=5, control_beta=None, test_beta=0.5,
                   delta_beta=None, control_entropy=None, test_entropy=0.7, delta_entropy=None,
                   epiallele_distance=None, passes_min_depth=False, control_000=0)
        complete = {**row, **{name: 0.0 for name in METRIC_VALUE_COLUMNS}, "control_n_complete_reads": 5, "passes_min_depth": True}
        cases = [row, {**complete, "passes_min_depth": False}] + [{**complete, name: None} for name in METRIC_VALUE_COLUMNS] + [complete]
        pq.write_table(pa.Table.from_pylist(cases), metrics)
        run = self.run_cli("--input", metrics, "--metrics-only", "--format", "csv", "--limit", 1)
        self.assertEqual(run.returncode, 0, run.stderr)
        reader = csv.DictReader(io.StringIO(run.stdout))
        self.assertEqual(reader.fieldnames, METRIC_SUMMARY_COLUMNS)
        self.assertEqual(reader.fieldnames, [
            "region_id", "chrom", "start", "end", "cpg_anchors",
            "control_beta", "test_beta", "delta_beta", "control_entropy",
            "test_entropy", "delta_entropy", "epiallele_distance"])
        exported = next(reader)
        self.assertEqual(float(exported["epiallele_distance"]), 0.0)
        self.assertEqual(list(reader), [])
        self.assertNotIn("NA", exported.values())
        self.assertNotIn("passes_min_depth", exported)
        self.assertNotIn("test_n_complete_reads", exported)
        # Unfiltered reading still exposes missing estimates in the full table.
        run = self.run_cli("--input", metrics, "--format", "csv", "--all")
        raw = list(csv.DictReader(io.StringIO(run.stdout)))
        self.assertEqual(len(raw), len(cases))
        self.assertEqual(raw[0]["epiallele_distance"], "NA")
        run = self.run_cli("--input", self.source, "--metrics-only")
        self.assertEqual(run.returncode, 2)
        self.assertIn("missing requested metric columns", run.stderr)


if __name__ == "__main__":
    unittest.main()
