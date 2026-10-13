from array import array
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

from methyldiver.genome import _plan, analyze_genome
from methyldiver.modbam import Region, analyze_modbam
from methyldiver.output import METRIC_SUMMARY_COLUMNS, format_cell


@unittest.skipUnless(importlib.util.find_spec("pyarrow"), "optional genome dependency pyarrow is not installed")
class TestGenome(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        # CpG at a scan boundary, dense/sparse CpGs, and a terminal G.
        self.sequences = {"chrA": "AACGCG" + "A" * 27 + "CGTTCGCGTTTCGACG",
                          "chrEmpty": "ACGTTTCG", "chrNoReads": "CGACGACGACG"}
        self.ref = self.root / "ref.fa"
        self.ref.write_text("".join(f">{n}\n{s}\n" for n, s in self.sequences.items()))
        pysam.faidx(str(self.ref))
        self.header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                       "SQ": [{"SN": n, "LN": len(s)} for n, s in self.sequences.items()]}
        self.control = self.make_bam("control.bam", ["0000000", "1111111", "1010101"])
        self.test = self.make_bam("test.bam", ["0011001", "1100110", "1111111"])

    def make_bam(self, name, patterns, tags=True):
        path = self.root / name
        sequence = self.sequences["chrA"]
        with pysam.AlignmentFile(str(path), "wb", header=self.header) as bam:
            for i, pattern in enumerate(patterns):
                self.assertEqual(len(pattern), sequence.count("CG"))
                read = pysam.AlignedSegment(bam.header)
                read.query_name = f"r{i}"
                read.query_sequence = sequence
                read.reference_id = 0
                read.reference_start = 0
                read.mapping_quality = 60
                read.cigarstring = f"{len(sequence)}M"
                read.query_qualities = array("B", [30] * len(sequence))
                if tags:
                    read.set_tag("MM", "C+m?," + ",".join("0" for _ in pattern) + ";")
                    # Ambiguous middle call exercises sparse incomplete output.
                    values = [255 if x == "1" else 0 for x in pattern]
                    if i == 2:
                        values[3] = 128
                    read.set_tag("ML", array("B", values))
                bam.write(read)
        pysam.index(str(path))
        return path

    def run_genome(self, output, **kwargs):
        kwargs.setdefault("keep_patterns", True)
        return analyze_genome(self.control, self.test, self.ref, output, **kwargs)

    def assert_matches_monolithic(self, directory, k):
        import pyarrow.parquet as pq
        actual = pq.read_table(directory / "complete_patterns").to_pylist()
        with (directory / "three_metrics.tsv").open(newline="") as handle:
            summary = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(summary.fieldnames, METRIC_SUMMARY_COLUMNS)
            summary_rows = list(summary)
        text = (directory / "three_metrics.tsv").read_text()
        self.assertNotIn('"', text)
        self.assertEqual(text.splitlines()[0].split("\t"), METRIC_SUMMARY_COLUMNS)
        self.assertEqual(summary_rows, [
            {name: str(format_cell(row[name])) for name in METRIC_SUMMARY_COLUMNS}
            for row in actual if row["control_n_complete_reads"] > 0 and row["test_n_complete_reads"] > 0])
        incomplete = pq.read_table(directory / "incomplete_patterns").to_pylist()
        expected = analyze_modbam(self.control, self.test, self.ref,
            [Region(n, 0, len(s)) for n, s in self.sequences.items()], k=k, min_depth=5)["units"]
        self.assertEqual(len(actual), len(expected))
        got = {(r["chrom"], tuple(r["cpg_anchors"])): r for r in actual}
        self.assertEqual(len(got), len(actual))
        missing = {}
        for r in incomplete:
            key = (r["chrom"], tuple(r["cpg_anchors"]), r["sample"])
            missing.setdefault(key, {})[r["pattern"]] = r["count"]
        for e in expected:
            key = (e["chrom"], tuple(e["cpg_anchors"]))
            a = got[key]
            self.assertEqual((a["start"], a["end"]), (e["start"], e["end"]))
            self.assertEqual(a["passes_min_depth"], e["passes_min_depth"])
            for sample in ("control", "test"):
                for name in ("beta", "entropy", "n_candidate_reads", "n_complete_reads", "n_incomplete_reads", "callable_fraction"):
                    self.assertEqual(a[f"{sample}_{name}"], e[sample][name])
                self.assertEqual([a[f"{sample}_{i:0{k}b}"] for i in range(2 ** k)], e[sample]["counts"])
                self.assertEqual(missing.get((*key, sample), {}), e[sample]["incomplete_patterns"])
            for name, value in e["comparison"].items():
                self.assertEqual(a[name], value)

    def test_reference_scan_and_chunk_boundaries(self):
        for k in (1, 2, 3, 4, 5, 10):
            with pysam.FastaFile(str(self.ref)) as fa:
                jobs, contigs = _plan(fa, [(n, len(s)) for n, s in self.sequences.items()], k, 2, scan_bases=3)
                self.assertEqual(sum(j["n_units"] for j in jobs), sum(max(0, s.count("CG") - k + 1) for s in self.sequences.values()))
                actual = []
                for j in jobs:
                    from methyldiver.coordinates import cpg_positions, sliding_units
                    positions = cpg_positions(fa.fetch(j["chrom"], j["start"], j["end"]), j["start"])
                    actual.extend((j["chrom"], u) for u in sliding_units(positions, k))
                expected = [(n, u) for n, s in self.sequences.items() for u in sliding_units(cpg_positions(s), k)]
                self.assertEqual(actual, expected)

    def test_serial_and_parallel_equal_monolithic_and_preserve_low_depth(self):
        for threads, k in ((1, 1), (2, 1), (2, 2), (1, 3), (2, 4), (1, 5), (2, 5)):
            with self.subTest(threads=threads, k=k):
                directory = self.root / f"run_{threads}_{k}"
                manifest = self.run_genome(directory, threads=threads, chunk_units=2, k=k, min_depth=5)
                self.assertEqual(manifest["status"], "complete")
                self.assertEqual(manifest["n_units_passing_min_depth"], 0)
                self.assertGreater(manifest["n_units"], 0)
                self.assert_matches_monolithic(directory, k)

    def test_empty_genome_has_typed_readable_datasets(self):
        import pyarrow.parquet as pq
        # k=4 leaves zero units on a two-CpG reference.
        self.sequences = {"chrA": "ACGACG"}
        self.ref.write_text(">chrA\nACGACG\n")
        Path(str(self.ref) + ".fai").unlink()
        pysam.faidx(str(self.ref))
        self.header["SQ"] = [{"SN": "chrA", "LN": 6}]
        self.control = self.make_bam("empty_control.bam", ["01"])
        self.test = self.make_bam("empty_test.bam", ["10"])
        directory = self.root / "empty_run"
        manifest = self.run_genome(directory, k=4)
        self.assertEqual(manifest["n_units"], 0)
        self.assertEqual(pq.read_table(directory / "complete_patterns").num_rows, 0)
        self.assertEqual(pq.read_table(directory / "incomplete_patterns").num_rows, 0)
        with (directory / "three_metrics.tsv").open() as handle:
            summary = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(summary.fieldnames, METRIC_SUMMARY_COLUMNS)
            self.assertEqual(list(summary), [])

    def test_ten_cpg_extraction_and_chunk_memory_limit(self):
        from unittest.mock import patch
        self.sequences = {"chrA": "CG" * 12}
        self.ref.write_text(">chrA\n" + self.sequences["chrA"] + "\n")
        Path(str(self.ref) + ".fai").unlink()
        pysam.faidx(str(self.ref))
        self.header["SQ"] = [{"SN": "chrA", "LN": 24}]
        self.control = self.make_bam("ten_control.bam", ["0" * 12, "1" * 12])
        self.test = self.make_bam("ten_test.bam", ["01" * 6, "10" * 6])
        directory = self.root / "ten_run"
        manifest = self.run_genome(directory, k=10, min_depth=5)
        self.assertEqual(manifest["n_units"], 3)
        self.assertEqual(manifest["requested_chunk_units"], 10000)
        self.assertEqual(manifest["chunk_units"], 976)
        self.assert_matches_monolithic(directory, 10)
        # Force two units/chunk to exercise the cap across a real boundary.
        with patch("methyldiver.genome.MAX_DENSE_COUNT_CELLS", 4096), \
             patch("methyldiver.modbam.MAX_DENSE_COUNT_CELLS", 4096):
            with self.assertRaisesRegex(ValueError, "Dense-pattern count limit"):
                analyze_modbam(self.control, self.test, self.ref, [Region("chrA", 0, 24)], k=10)
            limited = self.root / "ten_limited"
            manifest = self.run_genome(limited, k=10, min_depth=5)
            self.assertEqual(manifest["chunk_units"], 2)
            self.assertEqual(manifest["n_chunks"], 2)
        self.assert_matches_monolithic(limited, 10)

    def test_cli_N_default_and_validation(self):
        base = [sys.executable, "-m", "methyldiver", "compare", "--control", str(self.control),
                "--test", str(self.test), "--reference", str(self.ref)]
        for index, (size_args, n) in enumerate((([], 3), (["--N", "1"], 1),
                                              (["--N", "5"], 5))):
            for mode in (["--genome"], ["--region", f"chrA:0-{len(self.sequences['chrA'])}"]):
                output = self.root / f"size_{index}_{mode[0][2:]}"
                run = subprocess.run(base + mode + size_args + ["--output", str(output)],
                                     capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stderr)
                manifest = json.loads((output / "analysis_summary.json").read_text())
                self.assertEqual(manifest["config"]["k"], n)
                with (output / "three_metrics.tsv").open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(len(rows), 8 - n)
                self.assertTrue(all(len(r["cpg_anchors"].split(",")) == n for r in rows))
        for value in ("0", "11", "3.5"):
            output = self.root / "bad_size"
            run = subprocess.run(base + ["--genome", "--N", value, "--output", str(output)],
                                 capture_output=True, text=True)
            self.assertEqual(run.returncode, 2, run.stderr)
            self.assertNotIn("Traceback", run.stderr)
            self.assertFalse(output.exists())
        # The removed option must not be interpreted as --keep-patterns either.
        for value in ("5", "T"):
            run = subprocess.run(base + ["--genome", "--k", value, "--output", str(self.root / "old_option")],
                                 capture_output=True, text=True)
            self.assertEqual(run.returncode, 2, run.stderr)
            self.assertIn("unrecognized arguments: --k", run.stderr)
            self.assertFalse((self.root / "old_option").exists())

    def test_custom_sample_labels_and_reexport(self):
        output = self.root / "labelled"
        base = [sys.executable, "-m", "methyldiver", "compare", "--control", str(self.control),
                "--test", str(self.test), "--reference", str(self.ref), "--output", str(output),
                "--genome", "--threads", "2", "--chunk-units", "2", "--keep-patterns", "T",
                "--control-label", "DMSO", "--test-label", "DEC"]
        run = subprocess.run(base, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        manifest = json.loads((output / "analysis_summary.json").read_text())
        self.assertEqual(manifest["sample_labels"], {"control": "DMSO", "test": "DEC"})
        self.assertEqual(manifest["direction"], "test_minus_control")
        self.assertEqual(manifest["labelled_direction"], "DEC minus DMSO")
        with (output / "three_metrics.tsv").open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames[5:10], ["DMSO_beta", "DEC_beta", "delta_beta", "DMSO_entropy", "DEC_entropy"])
            for row in reader:
                self.assertAlmostEqual(float(row["delta_beta"]), float(row["DEC_beta"]) - float(row["DMSO_beta"]))
        import pyarrow.parquet as pq
        self.assertIn("control_beta", pq.read_schema(next((output / "complete_patterns").glob("*.parquet"))).names)
        read = [sys.executable, "-m", "methyldiver", "read", "--input", str(output / "complete_patterns"), "--metrics-only", "--all"]
        run = subprocess.run(read, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(run.stdout, (output / "three_metrics.tsv").read_text())
        run = subprocess.run(read + ["--test-label", "Treatment"], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn("DMSO_beta\tTreatment_beta", run.stdout.splitlines()[0])
        self.assertEqual(run.stdout.splitlines()[1:], (output / "three_metrics.tsv").read_text().splitlines()[1:])

    def test_chunk_manifest_storage_modes_preserve_results_and_checksums(self):
        import pyarrow.parquet as pq
        baseline = None
        for mode, threads in ((True, 1), (True, 2), (False, 2)):
            with self.subTest(mode=mode):
                directory = self.root / f"{mode}_{threads}"
                self.run_genome(directory, chunk_manifests=mode, threads=threads,
                                chunk_units=2, min_depth=5)
                root = json.loads((directory / "analysis_summary.json").read_text())
                self.assertEqual(root["chunk_manifests"], mode)
                details = root["outputs"]["three_metrics.tsv"]
                self.assertEqual(details["rows"], 5)
                self.assertEqual(details["rows"] + details["excluded_missing_metrics"], root["n_units"])
                self.assertEqual(root["n_units_passing_min_depth"], 0)  # TSV retention is not a depth filter.
                self.assertEqual(details["sha256"], hashlib.sha256((directory / "three_metrics.tsv").read_bytes()).hexdigest())
                self.assert_matches_monolithic(directory, 3)
                tables = [pq.read_table(directory / folder).to_pylist()
                          for folder in ("complete_patterns", "incomplete_patterns")]
                if baseline is None:
                    baseline = tables
                self.assertEqual(tables, baseline)
                self.assertEqual((directory / "chunk_manifests").exists(), mode)
                for chunk in root["chunks"]:
                    for relative, details in chunk["outputs"].items():
                        path = directory / relative
                        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), details["sha256"])
                        self.assertEqual(pq.read_metadata(path).num_rows, details["rows"])
                    if not mode:
                        self.assertIsNone(chunk["manifest"])
                        self.assertIsNone(chunk["manifest_sha256"])
                        continue
                    path = directory / chunk["manifest"]
                    self.assertTrue(path.name.endswith(".json.gz"))
                    data = path.read_bytes()
                    self.assertEqual(hashlib.sha256(data).hexdigest(), chunk["manifest_sha256"])
                    detail = json.loads(gzip.decompress(data))
                    self.assertEqual(detail["outputs"], chunk["outputs"])
                    self.assertEqual(detail["n_units"], chunk["n_units"])
                    self.assertIn("config", detail)
                if mode:
                    self.assertEqual(len(list((directory / "chunk_manifests").iterdir())), root["n_chunks"])

    def test_retention_options_preserve_plotting_table_and_mark_discarded_files(self):
        baseline = None
        for keep, chunks in ((False, False), (False, True), (True, False), (True, True)):
            directory = self.root / f"retention_{keep}_{chunks}"
            options = {} if not keep and not chunks else {"keep_patterns": keep, "chunk_manifests": chunks}
            root = analyze_genome(self.control, self.test, self.ref, directory, threads=2, chunk_units=2, **options)
            text = (directory / "three_metrics.tsv").read_text()
            if baseline is None:
                baseline = text
            self.assertEqual(text, baseline)
            self.assertIs(root["keep_patterns"], keep)
            self.assertIs(root["chunk_manifests"], chunks)
            for folder, retained in (("complete_patterns", keep), ("incomplete_patterns", keep), ("chunk_manifests", chunks)):
                self.assertEqual((directory / folder).exists(), retained)
            for chunk in root["chunks"]:
                for relative, details in chunk["outputs"].items():
                    self.assertIs(details["retained"], keep)
                    self.assertEqual((directory / relative).exists(), keep)
                if chunks:
                    detail = json.loads(gzip.decompress((directory / chunk["manifest"]).read_bytes()))
                    self.assertEqual(detail["outputs"], chunk["outputs"])

    def test_failures_do_not_publish_and_existing_work_is_preserved(self):
        output = self.root / "failed"
        with self.assertRaisesRegex(ValueError, "read limit"):
            self.run_genome(output, threads=2, chunk_units=1, max_reads_per_region=1)
        self.assertFalse(output.exists())
        self.assertEqual(list(self.root.glob(".failed.*.tmp")), [])
        output.mkdir()
        sentinel = output / "keep.txt"
        sentinel.write_text("existing work")
        with self.assertRaises(FileExistsError):
            self.run_genome(output)
        self.assertEqual(sentinel.read_text(), "existing work")

    def test_parameter_and_reference_validation(self):
        for kwargs in ({"threads": 0}, {"chunk_units": 0}, {"call_threshold": .5},
                       {"min_depth": 0}, {"k": 0}, {"k": 11}, {"min_base_quality": 255},
                       {"chunk_manifests": "invalid"}, {"chunk_manifests": "False"},
                       {"chunk_manifests": 0}, {"chunk_manifests": 1},
                       {"control_label": "delta"}, {"test_label": "control"},
                       {"control_label": "bad\tname"}, {"test_label": ""},
                       {"keep_patterns": "False"}, {"keep_patterns": 1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.run_genome(self.root / "invalid", **kwargs)
        self.header["SQ"][0]["LN"] += 1
        self.test = self.make_bam("mismatch.bam", ["0000000"])
        with self.assertRaisesRegex(ValueError, "reference dictionaries"):
            self.run_genome(self.root / "mismatch")

    def test_command_line_genome_and_incompatible_options(self):
        output = self.root / "cli"
        base = [sys.executable, "-m", "methyldiver", "compare", "--control", str(self.control),
                "--test", str(self.test), "--reference", str(self.ref), "--output", str(output)]
        run = subprocess.run(base + ["--genome", "--threads", "2", "--chunk-units", "2", "--min-depth", "5",
                                    "--keep-patterns", "T", "--chunk-manifests", "T"], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn("methyldiver: analysis complete | elapsed", run.stdout)
        self.assertIn("Units analyzed: 7 | passing depth threshold: 0", run.stdout)
        self.assertIn("Plotting table: 5 rows", run.stdout)
        self.assertIn("chunk_manifests/", run.stdout)
        self.assertIn("incomplete_patterns/", run.stdout)
        self.assertIn("complete_patterns/", run.stdout)
        self.assertNotIn('"output_directory"', run.stdout)
        self.assertIn("analyze 100.0%", run.stderr)
        self.assertIn("write TSV 100.0%", run.stderr)
        self.assertIn("ETA", run.stderr)
        self.assertNotIn("\x1b", run.stderr)
        self.assertFalse((output / "CpG_unit_pattern").exists())
        self.assertFalse((output / "metrics").exists())
        self.assertIs(json.loads((output / "analysis_summary.json").read_text())["chunk_manifests"], True)
        self.assertTrue(list((output / "chunk_manifests").glob("*.json.gz")))
        self.assertFalse(list((output / "chunk_manifests").glob("*.json")))
        self.assert_matches_monolithic(output, 3)
        for mode in ("True", "False", "T", "F", "true", "false", "t", "f"):
            enabled = mode.lower() in ("true", "t")
            mode_output = self.root / f"cli_{mode}"
            run = subprocess.run(base[:-1] + [str(mode_output), "--genome", "--chunk-manifests", mode, "--keep-patterns", mode],
                                 capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertEqual("chunk_manifests/" in run.stdout, enabled)
            self.assertEqual((mode_output / "chunk_manifests").exists(), enabled)
            self.assertIs(json.loads((mode_output / "analysis_summary.json").read_text())["chunk_manifests"], enabled)
            self.assertEqual((mode_output / "complete_patterns").exists(), enabled)
            self.assertEqual((mode_output / "incomplete_patterns").exists(), enabled)
        default_output = self.root / "cli_default"
        run = subprocess.run(base[:-1] + [str(default_output), "--genome"], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn("Files: three_metrics.tsv, chunks.tsv, analysis_summary.json", run.stdout)
        self.assertEqual({p.name for p in default_output.iterdir()}, {"three_metrics.tsv", "chunks.tsv", "analysis_summary.json"})
        for options in (["--genome", "--region", "chrA:0-6"], ["--genome", "--max-total-bases", "100"],
                        ["--region", "chrA:0-6", "--chunk-manifests", "False"],
                        ["--region", "chrA:0-6", "--keep-patterns", "False"],
                        ["--genome", "--keep-patterns", "invalid"],
                        ["--genome", "--chunk-manifests", "invalid"],
                        ["--genome", "--chunk-manifests", "json"],
                        ["--genome", "--chunk-manifests", "gzip"],
                        ["--genome", "--chunk-manifests", "none"],
                        ["--region", "chrA:0-6", "--threads", "2"], ["--genome", "--chunk-units", "0"]):
            run = subprocess.run(base[:-1] + [str(self.root / "bad_cli")] + options, capture_output=True, text=True)
            self.assertEqual(run.returncode, 2, run.stderr)
            self.assertNotIn("Traceback", run.stderr)
            self.assertFalse((self.root / "bad_cli").exists())


if __name__ == "__main__":
    unittest.main()
