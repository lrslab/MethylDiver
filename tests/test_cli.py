from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from methyldiver import __version__


class TestCLI(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def run_cli(self, *arguments):
        return subprocess.run(
            [sys.executable, "-m", "methyldiver", *map(str, arguments)],
            capture_output=True, text=True, cwd=self.root,
        )

    def test_help_and_version(self):
        for arguments in (("--help",), ("compare", "--help"), ("single", "--help"), ("read", "--help")):
            run = self.run_cli(*arguments)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertNotIn("metrics-counts", run.stdout)
        help_text = self.run_cli("--help").stdout
        self.assertIn("Program: methyldiver", help_text)
        self.assertIn(f"Version: {__version__}", help_text)
        self.assertIn("Commands:", help_text)
        self.assertEqual(self.run_cli().stdout, help_text)
        for command in ("compare", "single", "read"):
            run = self.run_cli(command)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertEqual(run.stdout, self.run_cli(command, "--help").stdout)
            self.assertIn(f"Usage: methyldiver {command} [options]", run.stdout)
            self.assertNotIn("<command>", run.stdout)
        self.assertEqual(self.run_cli("--version").stdout.strip(), f"methyldiver {__version__}")

    def test_removed_command_is_rejected_without_output(self):
        output = self.root / "metrics.json"
        for command in ("metrics-counts", "analyze"):
            run = self.run_cli(command, "--input", self.root / "input.bam", "--output", output)
            self.assertEqual(run.returncode, 2)
            self.assertIn("invalid choice", run.stderr)
            self.assertNotIn("Traceback", run.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
