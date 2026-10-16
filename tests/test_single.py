from array import array
import csv
import gzip
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam

from methyldiver.modbam import Region, _extract, analyze_modbam, analyze_single_modbam
from methyldiver.output import SINGLE_METRIC_SUMMARY_COLUMNS, write_analysis


class TestSingle(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.sequences = {"chrA": "CG" * 12, "chrNoReads": "CG" * 4}
        self.ref = self.root / "ref.fa"
        self.ref.write_text("".join(f">{name}\n{seq}\n" for name, seq in self.sequences.items()))
        pysam.faidx(str(self.ref))
        self.header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                       "SQ": [{"SN": name, "LN": len(seq)} for name, seq in self.sequences.items()]}
        self.bam = self.make_bam("sample.bam")
        self.regions = [Region(name, 0, len(seq)) for name, seq in self.sequences.items()]

    def make_bam(self, name, all_missing=False, tags=True):
        path = self.root / name
        with pysam.AlignmentFile(str(path), "wb", header=self.header) as bam:
            for i, pattern in enumerate(("0" * 12, "1" * 12, "01" * 6)):
                read = pysam.AlignedSegment(bam.header)
                read.query_name = f"r{i}"
                read.query_sequence = self.sequences["chrA"]
                read.reference_id, read.reference_start = 0, 0
                read.mapping_quality = 60
                read.cigarstring = "24M"
                read.query_qualities = array("B", [30] * 24)
                if tags:
                    read.set_tag("MM", "C+m?," + ",".join("0" for _ in pattern) + ";")
                    values = [128 if all_missing else (255 if x == "1" else 0) for x in pattern]
                    if i == 2:
                        values[3] = 128
                    read.set_tag("ML", array("B", values))
                bam.write(read)
        pysam.index(str(path))
        return path

    def run_cli(self, command, *args):
        return subprocess.run([sys.executable, "-m", "methyldiver", command, *map(str, args)],
                              capture_output=True, text=True)

    def parse(self, path):
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_single_matches_paired_sample_without_comparison_or_duplicate_extraction(self):
        for n in (1, 3, 5, 10):
            paired = analyze_modbam(self.bam, self.bam, self.ref, self.regions, k=n, min_depth=5)
            with patch("methyldiver.modbam._extract", wraps=_extract) as extract, \
                 patch("methyldiver.modbam.compare_counts", side_effect=AssertionError("No comparison in single mode")):
                single = analyze_single_modbam(self.bam, self.ref, self.regions, k=n, min_depth=5)
            self.assertEqual(extract.call_count, len(self.regions))
            self.assertEqual(single["manifest"]["analysis_type"], "single")
            self.assertEqual(set(single["manifest"]["inputs"]), {"sample", "reference"})
            self.assertNotIn("direction", single["manifest"])
            for actual, expected in zip(single["units"], paired["units"], strict=True):
                self.assertEqual(actual["sample"], expected["control"])
                self.assertFalse(actual["passes_min_depth"])
                self.assertNotIn("comparison", actual)
                self.assertNotIn("test", actual)
            self.assertGreater(single["units"][0]["sample"]["n_complete_reads"], 0)

    def test_region_cli_and_bed_with_label_and_depth(self):
        bed = self.root / "regions.bed"
        bed.write_text("chrA\t0\t24\ttoy\nchrNoReads\t0\t8\tempty\n")
        for name, selection, label_args, label in (
            ("region", ["--region", "chrA:0-24", "--region", "chrNoReads:0-8"], [], "sample"),
            ("bed", ["--regions-bed", bed], ["--sample-label", "D0"], "D0"),
        ):
            output = self.root / name
            run = self.run_cli("single", "--input", self.bam, "--reference", self.ref,
                               *selection, *label_args, "--threads", 1, "--min-depth", 5, "--output", output)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn("Plotting table: 0 rows", run.stdout)
            self.assertEqual({p.name for p in output.iterdir()},
                             {"two_metrics.tsv", "metrics.tsv", "patterns.tsv", "analysis_summary.json"})
            rows = self.parse(output / "two_metrics.tsv")
            self.assertEqual(rows, [])
            with (output / "two_metrics.tsv").open() as handle:
                self.assertEqual(next(csv.reader(handle, delimiter="\t")),
                                 SINGLE_METRIC_SUMMARY_COLUMNS[:5] + [f"{label}_beta", f"{label}_entropy"])
            full = self.parse(output / "metrics.tsv")
            self.assertEqual(len(full), 12)
            self.assertTrue(all(row["passes_min_depth"] == "0" for row in full))
            self.assertEqual(full[-1]["sample_beta"], "NA")
            self.assertEqual({row["sample"] for row in self.parse(output / "patterns.tsv")}, {"sample"})
            manifest = json.loads((output / "analysis_summary.json").read_text())
            self.assertEqual(manifest["sample_labels"], {"sample": label})
            self.assertEqual(manifest["outputs"]["two_metrics.tsv"]["excluded_missing_metrics"], 2)
            self.assertEqual(manifest["outputs"]["two_metrics.tsv"]["excluded_below_min_depth"], 10)
            if name == "bed":
                self.assertIn("regions_bed_sha256", manifest)

    def test_missing_metrics_and_short_region(self):
        bam = self.make_bam("missing.bam", all_missing=True)
        result = analyze_single_modbam(bam, self.ref, self.regions)
        self.assertTrue(all(unit["sample"]["beta"] is None for unit in result["units"]))
        output = self.root / "missing"
        write_analysis(result, output)
        self.assertEqual(self.parse(output / "two_metrics.tsv"), [])
        self.assertGreater(len(self.parse(output / "patterns.tsv")), 0)
        short = analyze_single_modbam(self.bam, self.ref, [Region("chrA", 0, 4)])
        self.assertEqual(short["units"], [])

    def test_validation_and_failed_run_preserve_existing_work(self):
        ordinary = self.make_bam("ordinary.bam", tags=False)
        for extra in (["--sample-label", "bad name"], ["--N", "11"], ["--k", "3"],
                      ["--test", self.bam], ["--keep-patterns", "invalid"], ["--threads", "0"]):
            output = self.root / "invalid"
            run = self.run_cli("single", "--input", self.bam, "--reference", self.ref,
                               "--region", "chrA:0-24", "--output", output, *extra)
            self.assertEqual(run.returncode, 2, run.stderr)
            self.assertNotIn("Traceback", run.stderr)
            self.assertFalse(output.exists())
        run = self.run_cli("single", "--input", ordinary, "--reference", self.ref,
                           "--region", "chrA:0-24", "--output", self.root / "ordinary")
        self.assertEqual(run.returncode, 2)
        self.assertFalse((self.root / "ordinary").exists())
        existing = self.root / "existing"
        existing.mkdir()
        (existing / "keep").write_text("existing work")
        run = self.run_cli("single", "--input", self.bam, "--reference", self.ref,
                           "--region", "chrA:0-24", "--output", existing)
        self.assertEqual(run.returncode, 2)
        self.assertEqual((existing / "keep").read_text(), "existing work")

    @unittest.skipUnless(importlib.util.find_spec("pyarrow"), "optional pyarrow is not installed")
    def test_genome_serial_parallel_retention_and_reader(self):
        import pyarrow.parquet as pq
        from methyldiver.genome import analyze_single_genome
        for n, threads in ((1, 1), (3, 2), (5, 2), (10, 1)):
            output = self.root / f"genome_{n}"
            run = self.run_cli("single", "--input", self.bam, "--reference", self.ref, "--genome",
                               "--threads", threads, "--N", n, "--chunk-units", 2,
                               "--keep-patterns", "T", "--chunk-manifests", "T",
                               "--sample-label", "D0", "--min-depth", 5, "--output", output)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn("two_metrics.tsv", run.stdout)
            actual = pq.read_table(output / "complete_patterns").to_pylist()
            expected = analyze_single_modbam(self.bam, self.ref, self.regions, k=n, min_depth=5)["units"]
            self.assertEqual(len(actual), len(expected))
            for row, unit in zip(actual, expected, strict=True):
                self.assertEqual(row["cpg_anchors"], unit["cpg_anchors"])
                self.assertEqual(row["chrom"], unit["chrom"])
                for key in ("beta", "entropy", "n_complete_reads", "n_incomplete_reads"):
                    self.assertEqual(row[f"sample_{key}"], unit["sample"][key])
                self.assertEqual([row[f"sample_{i:0{n}b}"] for i in range(2 ** n)], unit["sample"]["counts"])
                self.assertNotIn("epiallele_distance", row)
            manifest = json.loads((output / "analysis_summary.json").read_text())
            self.assertEqual(set(manifest["inputs"]), {"sample", "reference", "sample_index", "reference_index"})
            self.assertEqual(set(manifest["headers"]), {"sample"})
            for path in (output / "chunk_manifests").glob("*.json.gz"):
                self.assertEqual(json.loads(gzip.decompress(path.read_bytes()))["analysis_type"], "single")
            read = self.run_cli("read", "--input", output / "complete_patterns", "--metrics-only", "--all")
            self.assertEqual(read.returncode, 0, read.stderr)
            self.assertEqual(read.stdout, (output / "two_metrics.tsv").read_text())
            relabel = self.run_cli("read", "--input", output / "complete_patterns", "--metrics-only",
                                   "--sample-label", "Tumor", "--all")
            self.assertEqual(relabel.returncode, 0, relabel.stderr)
            self.assertEqual(relabel.stdout, read.stdout.replace("D0_beta\tD0_entropy", "Tumor_beta\tTumor_entropy"))
            incomplete = self.run_cli("read", "--input", output / "incomplete_patterns", "--sample", "sample", "--all")
            self.assertEqual(incomplete.returncode, 0, incomplete.stderr)
            self.assertTrue(all("." in row["pattern"] for row in csv.DictReader(io.StringIO(incomplete.stdout), delimiter="\t")))
            default_output = self.root / f"default_{n}"
            default = analyze_single_genome(self.bam, self.ref, default_output, k=n, min_depth=5,
                                            chunk_units=2, sample_label="D0")
            self.assertEqual({p.name for p in default_output.iterdir()}, {"two_metrics.tsv", "chunks.tsv", "analysis_summary.json"})
            self.assertEqual((default_output / "two_metrics.tsv").read_bytes(), (output / "two_metrics.tsv").read_bytes())
            self.assertEqual(default["n_units_passing_min_depth"], 0)

    @unittest.skipUnless(importlib.util.find_spec("pyarrow"), "optional pyarrow is not installed")
    def test_empty_genome_and_failed_genome_cleanup(self):
        import pyarrow.parquet as pq
        from methyldiver.genome import analyze_single_genome
        empty = self.root / "empty.bam"
        with pysam.AlignmentFile(str(empty), "wb", header={"HD": self.header["HD"], "SQ": [self.header["SQ"][1]]}):
            pass
        pysam.index(str(empty))
        for threads in (1, 2):
            output = self.root / f"empty_{threads}"
            manifest = analyze_single_genome(empty, self.ref, output, k=10, threads=threads, keep_patterns=True)
            self.assertEqual(manifest["n_units"], 0)
            self.assertEqual(manifest["chunk_units"], 1953)
            table = pq.read_table(output / "complete_patterns")
            self.assertEqual(table.num_rows, 0)
            self.assertIn("sample_beta", table.schema.names)
            self.assertNotIn("control_beta", table.schema.names)
            self.assertEqual(self.parse(output / "two_metrics.tsv"), [])
        ordinary = self.make_bam("ordinary_genome.bam", tags=False)
        output = self.root / "failed_genome"
        with self.assertRaisesRegex(ValueError, "no MM/Mm tags"):
            analyze_single_genome(ordinary, self.ref, output, threads=2, chunk_units=2)
        self.assertFalse(output.exists())
        self.assertEqual(list(self.root.glob(".failed_genome.*")), [])


if __name__ == "__main__":
    unittest.main()
