import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from methyldiver import __version__


ROOT = Path(__file__).resolve().parents[1]


class TestCLI(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root / "counts.tsv"
        self.source.write_bytes((ROOT / "examples/legacy_counts.tsv").read_bytes())

    def run_cli(self, *arguments):
        return subprocess.run(
            [sys.executable, "-m", "methyldiver", *map(str, arguments)],
            capture_output=True, text=True, cwd=self.root,
        )

    def test_legacy_cli_and_output_provenance(self):
        output = self.root / "new" / "metrics.json"
        run = self.run_cli("metrics-counts", "--input", self.source, "--output", output)
        self.assertEqual(run.returncode, 0, run.stderr)
        result = json.loads(output.read_text())
        self.assertEqual(json.loads(run.stdout), result)
        self.assertEqual(result["unit_size"], 3)
        self.assertAlmostEqual(result["comparison"]["D_epi"], 0.7309379165994421)
        self.assertEqual(result["provenance"]["input_sha256"], hashlib.sha256(self.source.read_bytes()).hexdigest())
        self.assertIsNone(result["provenance"]["minimum_depth_filter"])

    def test_four_cpg_empty_group_is_json_null(self):
        self.source.write_text("pattern\tcontrol\ttest\n" + "".join(f"{i:04b}\t0\t1\n" for i in range(16)))
        run = self.run_cli("metrics-counts", "--input", self.source)
        self.assertEqual(run.returncode, 0, run.stderr)
        result = json.loads(run.stdout)
        self.assertEqual(result["unit_size"], 4)
        self.assertIsNone(result["control"]["beta"])
        self.assertIsNone(result["comparison"]["D_epi"])
        self.assertEqual(result["test"]["entropy"], 1)
        self.assertNotIn("NaN", run.stdout)

    def test_refuses_existing_output_and_allows_explicit_overwrite(self):
        output = self.root / "result.json"
        output.write_text("existing work")
        run = self.run_cli("metrics-counts", "--input", self.source, "--output", output)
        self.assertNotEqual(run.returncode, 0)
        self.assertEqual(output.read_text(), "existing work")
        run = self.run_cli("metrics-counts", "--input", self.source, "--output", output, "--overwrite")
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(json.loads(output.read_text())["control"]["beta"], 0.8)
        self.assertEqual(list(self.root.glob(".*.tmp")), [])

    def test_input_cannot_be_overwritten_through_alias(self):
        original = self.source.read_bytes()
        hard = self.root / "hard.tsv"
        os.link(self.source, hard)
        link = self.root / "link.tsv"
        link.symlink_to(self.source)
        for path in (self.source, hard, link):
            run = self.run_cli("metrics-counts", "--input", self.source, "--output", path, "--overwrite")
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("must not replace", run.stderr)
            self.assertEqual(self.source.read_bytes(), original)

    def test_bad_tables_have_clear_errors_and_no_output(self):
        original = self.source.read_text()
        rows = original.splitlines()
        cases = [
            original.replace("pattern\tcontrol\ttest", "chrom\tstart\tend"),
            "\n".join(rows[:-1]),
            "\n".join([rows[0], rows[2], rows[1], *rows[3:]]),
            original.replace("000\t2\t10", "000\t-2\t10"),
            original.replace("000\t2\t10", "000\t2.0\t10"),
            original.replace("000\t2\t10", "000\t2"),
            original.replace("000\t2\t10", "000\t2\t10\textra"),
            original.replace("000\t2\t10", "000\tNaN\t10"),
        ]
        output = self.root / "invalid.json"
        for content in cases:
            self.source.write_text(content)
            run = self.run_cli("metrics-counts", "--input", self.source, "--output", output)
            self.assertEqual(run.returncode, 2, run.stderr)
            self.assertIn("error:", run.stderr)
            self.assertNotIn("Traceback", run.stderr)
            self.assertFalse(output.exists())

    def test_help_version_and_missing_file(self):
        self.assertEqual(self.run_cli("--help").returncode, 0)
        self.assertEqual(self.run_cli("--version").stdout.strip(), f"methyldiver {__version__}")
        run = self.run_cli("metrics-counts", "--input", self.root / "absent.tsv")
        self.assertEqual(run.returncode, 2)
        self.assertNotIn("Traceback", run.stderr)


if __name__ == "__main__":
    unittest.main()
