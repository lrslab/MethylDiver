"""Live terminal progress and plain-text analysis completion summaries."""

import shutil
import sys
import threading
import time


def duration(seconds):
    seconds = max(0, int(seconds))
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


class LiveProgress:
    """Refresh elapsed time even while the next chunk or I/O operation runs."""

    def __init__(self, stream=None, clock=None):
        self.stream = sys.stderr if stream is None else stream
        self.clock = time.monotonic if clock is None else clock
        self.tty = self.stream.isatty()
        self.started = self.phase_started = self.clock()
        self.phase = "planning"
        self.done = self.total = self.units = None
        self.lock = threading.Lock()
        self.stop = threading.Event()
        self.thread = None

    @property
    def elapsed(self):
        return max(0, self.clock() - self.started)

    def _line(self):
        now = self.clock()
        line = f"methyldiver: {self.phase}"
        eta = "--:--:--"
        if self.total is not None:
            fraction = self.done / self.total if self.total else 1.0
            percentage = min(99.9, 100 * fraction) if self.done < self.total else 100.0
            line += f" {percentage:5.1f}% {self.done}/{self.total}"
            if self.done == self.total:
                eta = "00:00:00"
            elif self.done and now > self.phase_started:
                eta = "~" + duration((now - self.phase_started) * (self.total - self.done) / self.done)
        line += f" | elapsed {duration(now - self.started)} | ETA {eta}"
        width = shutil.get_terminal_size((80, 24)).columns if self.tty else None
        if self.tty and len(line) >= width and self.total is not None:
            line = line.replace(f" {self.done}/{self.total}", "")
        if self.units is not None and (not self.tty or shutil.get_terminal_size((80, 24)).columns >= 120):
            line += f" | {self.units:,} units"
        if self.tty and len(line) >= width:
            line = line[:max(0, width - 4)] + "..."
        return line

    def _emit(self):
        # Carriage-return updates only on a terminal; redirected logs get lines.
        prefix, suffix = ("\r\033[2K", "") if self.tty else ("", "\n")
        self.stream.write(prefix + self._line() + suffix)
        self.stream.flush()

    def stage(self, phase, done=None, total=None):
        with self.lock:
            changed = phase != self.phase
            if changed:
                self.phase = phase
                self.phase_started = self.clock()
                self.units = None
            finished = total is not None and done == total and self.done != done
            self.done, self.total = done, total
            if changed or finished:
                self._emit()

    def analysis(self, done, total, units):
        with self.lock:
            changed = self.phase != "analyze"
            if changed:
                self.phase = "analyze"
                self.phase_started = self.clock()
            self.done, self.total, self.units = done, total, units
            if changed or done == total:
                self._emit()

    def _refresh(self):
        while not self.stop.wait(1.0 if self.tty else 5.0):
            with self.lock:
                self._emit()

    def __enter__(self):
        with self.lock:
            self._emit()
        self.thread = threading.Thread(target=self._refresh, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.stop.set()
        self.thread.join()
        with self.lock:
            if self.tty:
                self.stream.write("\r\033[2K")
            if exc_type is not None:
                self.stream.write(f"methyldiver: stopped during {self.phase} | elapsed {duration(self.elapsed)}\n")
            self.stream.flush()


def print_analysis_summary(output, manifest, files, elapsed):
    print(f"methyldiver: analysis complete | elapsed {duration(elapsed)}")
    if manifest.get("config", {}).get("modification"):
        print(f"  Modification: {manifest['config']['modification']} | context: CpG")
    print(f"  Units analyzed: {manifest['n_units']:,} | passing depth threshold: {manifest['n_units_passing_min_depth']:,}")
    outputs = manifest.get("outputs", {})
    table = outputs.get("two_metrics.tsv", outputs.get("three_metrics.tsv"))
    if table is not None:
        print(f"  Plotting table: {table['rows']:,} rows | omitted with missing metrics: {table['excluded_missing_metrics']:,}")
    print(f"  Output: {output.resolve()}")
    print(f"  Files: {', '.join(files)}")
