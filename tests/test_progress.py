import io
import unittest
from unittest.mock import patch

from methyldiver.progress import LiveProgress, duration


class Terminal(io.StringIO):
    def isatty(self):
        return True


class TestProgress(unittest.TestCase):
    def test_percentage_elapsed_eta_and_phase_reset(self):
        now = [0.0]
        stream = io.StringIO()
        reporter = LiveProgress(stream=stream, clock=lambda: now[0])
        self.assertIn("ETA --:--:--", reporter._line())
        now[0] = 5
        reporter.analysis(0, 100, 0)
        now[0] = 25
        reporter.analysis(25, 100, 250000)
        line = reporter._line()
        self.assertIn("25.0% 25/100", line)
        self.assertIn("elapsed 00:00:25", line)
        self.assertIn("ETA ~00:01:00", line)
        self.assertIn("250,000 units", line)
        # Elapsed/ETA advance even with no newly completed chunks.
        now[0] = 35
        self.assertIn("elapsed 00:00:35", reporter._line())
        self.assertIn("ETA ~00:01:30", reporter._line())
        reporter.analysis(100, 100, 1000000)
        self.assertIn("100.0%", reporter._line())
        reporter.stage("write TSV", 0, 50)
        self.assertIn("ETA --:--:--", reporter._line())
        now[0] = 45
        reporter.stage("write TSV", 25, 50)
        self.assertIn("50.0% 25/50", reporter._line())
        self.assertIn("ETA ~00:00:10", reporter._line())
        reporter.stage("finalize")
        self.assertNotIn("%", reporter._line())
        self.assertNotIn("\r", stream.getvalue())
        self.assertNotIn("\x1b", stream.getvalue())

    def test_empty_work_and_long_durations(self):
        reporter = LiveProgress(stream=io.StringIO(), clock=lambda: 0)
        reporter.analysis(0, 0, 0)
        self.assertIn("100.0% 0/0", reporter._line())
        self.assertIn("ETA 00:00:00", reporter._line())
        self.assertEqual(duration(90061), "25:01:01")
        self.assertEqual(duration(-1), "00:00:00")
        reporter.analysis(99999, 100000, 999990000)
        self.assertIn("99.9%", reporter._line())

    def test_terminal_refresh_and_shutdown_on_error(self):
        stream = Terminal()
        with patch("methyldiver.progress.shutil.get_terminal_size", return_value=type("Size", (), {"columns": 80})()):
            with LiveProgress(stream=stream) as reporter:
                reporter.analysis(10, 100, 100000)
                reporter.stage("write TSV", 4908607, 4908608)
                self.assertLess(len(reporter._line()), 80)
            self.assertFalse(reporter.thread.is_alive())
            self.assertIn("\r\x1b[2K", stream.getvalue())
        failed = io.StringIO()
        with self.assertRaisesRegex(ValueError, "failed"):
            with LiveProgress(stream=failed) as reporter:
                reporter.analysis(1, 10, 10000)
                raise ValueError("failed")
        self.assertFalse(reporter.thread.is_alive())
        self.assertIn("stopped during analyze", failed.getvalue())
        self.assertNotIn("analysis complete", failed.getvalue())


if __name__ == "__main__":
    unittest.main()
