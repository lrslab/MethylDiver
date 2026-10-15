import csv
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import unittest

import pysam

import test_genome as fixtures
from methyldiver.coordinates import cpg_positions, sliding_units
from methyldiver.genome import analyze_regions
from methyldiver.modbam import Region, analyze_modbam
from methyldiver.regions import plan_regions


@unittest.skipUnless(importlib.util.find_spec("pyarrow"), "optional pyarrow is not installed")
class TestRegions(unittest.TestCase):
    setUp = fixtures.TestGenome.setUp
    make_bam = fixtures.TestGenome.make_bam

    def test_planning_preserves_boundaries_and_all_units(self):
        regions = [Region("chrA", 0, 6, "first"),
                   Region("chrA", 6, len(self.sequences["chrA"]), "second"),
                   Region("chrEmpty", 0, 2, "partial_dyad"),
                   Region("chrNoReads", 0, len(self.sequences["chrNoReads"]), "no_reads")]
        with pysam.FastaFile(str(self.ref)) as fasta:
            for k in (1, 2, 3, 5, 10):
                expected = [(r.region_id, unit) for r in regions for unit in sliding_units(
                    cpg_positions(fasta.fetch(r.chrom, r.start, r.end), r.start), k)]
                for limit in (3, 10, 100):
                    jobs, reports = plan_regions(fasta, regions, k, 2, limit, scan_bases=3)
                    actual = []
                    for job in jobs:
                        self.assertLessEqual(job["n_units"], 2)
                        span = sum(p["end"] - p["start"] for p in job["parts"])
                        if span > limit:
                            self.assertEqual(job["n_units"], 1)
                        for p in job["parts"]:
                            actual.extend((p["region_id"], u) for u in sliding_units(
                                cpg_positions(fasta.fetch(p["chrom"], p["start"], p["end"]), p["start"]), k))
                    self.assertEqual(actual, expected)
                    self.assertEqual(sum(r["n_units"] for r in reports), len(expected))

    def assert_matches(self, output, expected, single=False):
        import pyarrow.parquet as pq
        rows = pq.read_table(output / "complete_patterns").to_pylist()
        incomplete = pq.read_table(output / "incomplete_patterns").to_pylist()
        got = {(r["region_id"], tuple(r["cpg_anchors"])): r for r in rows}
        self.assertEqual(len(got), len(rows))
        self.assertEqual(len(rows), len(expected))
        missing = {}
        for r in incomplete:
            missing.setdefault((r["region_id"], tuple(r["cpg_anchors"]), r["sample"]), {})[r["pattern"]] = r["count"]
        for unit in expected:
            key = (unit["region_id"], tuple(unit["cpg_anchors"]))
            row = got[key]
            self.assertEqual(row["passes_min_depth"], unit["passes_min_depth"])
            for role in (("sample",) if single else ("control", "test")):
                for field in ("beta", "entropy", "n_complete_reads", "n_candidate_reads", "n_incomplete_reads", "callable_fraction"):
                    self.assertEqual(row[f"{role}_{field}"], unit[role][field])
                k = unit["unit_size"]
                self.assertEqual([row[f"{role}_{i:0{k}b}"] for i in range(2 ** k)], unit[role]["counts"])
                self.assertEqual(missing.get((*key, role), {}), unit[role]["incomplete_patterns"])
            for name, value in unit.get("comparison", {}).items():
                self.assertEqual(row[name], value)

    def test_serial_parallel_paired_single_match_unchunked(self):
        regions = [Region("chrA", 0, len(self.sequences["chrA"]), "promoter_A"),
                   Region("chrNoReads", 0, len(self.sequences["chrNoReads"]), "promoter_B")]
        for single, k in ((False, 1), (False, 3), (True, 3), (True, 5)):
            test = None if single else self.test
            expected = analyze_modbam(self.control, test, self.ref, regions, k=k, min_depth=5)["units"]
            previous = None
            for threads in (1, 2):
                output = self.root / f"regions_{single}_{k}_{threads}"
                report = analyze_regions(self.control, test, self.ref, regions, output,
                    threads=threads, chunk_units=2, max_total_bases=10, k=k,
                    min_depth=5, keep_patterns=True, chunk_manifests=True)
                self.assertEqual(report["mode"], "regions")
                self.assertEqual(report["n_regions"], 2)
                self.assertEqual(report["n_units_passing_min_depth"], 0)
                self.assert_matches(output, expected, single)
                name = "two_metrics.tsv" if single else "three_metrics.tsv"
                data = (output / name).read_bytes()
                if previous is not None:
                    self.assertEqual(data, previous)
                previous = data
                for path in (output / "chunk_manifests").glob("*.gz"):
                    chunk = json.loads(gzip.decompress(path.read_bytes()))
                    self.assertTrue(all(r["name"].startswith("promoter_") for r in chunk["regions"]))

    def test_cli_auto_batches_full_span_and_retains_original_id(self):
        # A sparse unit crossing 100 kb must survive the automatic split.
        self.sequences["chrA"] += "A" * 100_010 + "CG"
        self.ref.write_text("".join(f">{n}\n{s}\n" for n, s in self.sequences.items()))
        Path(str(self.ref) + ".fai").unlink()
        pysam.faidx(str(self.ref))
        self.header["SQ"][0]["LN"] = len(self.sequences["chrA"])
        control = self.make_bam("long_control.bam", ["00000000", "11111111"])
        test = self.make_bam("long_test.bam", ["00110011", "11001100"])
        region = Region("chrA", 0, len(self.sequences["chrA"]), "full_promoter")
        bed = self.root / "all.bed"
        bed.write_text(f"chrA\t0\t{region.end}\tfull_promoter\n")
        output = self.root / "auto"
        run = subprocess.run([sys.executable, "-m", "methyldiver", "compare", "--control", str(control),
            "--test", str(test), "--reference", str(self.ref), "--regions-bed", str(bed),
            "--output", str(output)], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        report = json.loads((output / "analysis_summary.json").read_text())
        self.assertEqual(report["mode"], "regions")
        self.assertGreater(report["n_chunks"], 1)
        self.assertEqual(report["n_units"], 6)
        self.assertIn("regions_bed_sha256", report)
        self.assertEqual({p.name for p in output.iterdir()}, {"chunks.tsv", "three_metrics.tsv", "analysis_summary.json"})
        rows = list(csv.DictReader((output / "three_metrics.tsv").read_text().splitlines(), delimiter="\t"))
        expected = analyze_modbam(control, test, self.ref, [region], max_total_bases=region.end)["units"]
        self.assertEqual(len(rows), len(expected))
        for row, unit in zip(rows, expected):
            self.assertEqual(row["region_id"], "full_promoter")
            self.assertEqual(row["cpg_anchors"], ",".join(map(str, unit["cpg_anchors"])))
            self.assertEqual(float(row["control_beta"]), unit["control"]["beta"])

    def test_cli_parallel_region_and_bed(self):
        bed = self.root / "small.bed"
        bed.write_text(f"chrA\t0\t{len(self.sequences['chrA'])}\tpromoter\n")
        for command, scope in (("compare", ["--regions-bed", str(bed)]),
                               ("single", ["--region", f"chrA:0-{len(self.sequences['chrA'])}",
                                           "--region", "chrNoReads:0-11"])):
            output = self.root / command
            inputs = ["--input", str(self.control)] if command == "single" else ["--control", str(self.control), "--test", str(self.test)]
            run = subprocess.run([sys.executable, "-m", "methyldiver", command, *inputs,
                "--reference", str(self.ref), *scope, "--threads", "2", "--chunk-units", "2",
                "--keep-patterns", "True", "--output", str(output)], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn("100.0%", run.stderr)
            self.assertTrue((output / "complete_patterns").is_dir())

    def test_empty_invalid_and_failed_requests(self):
        regions = [Region("chrA", 0, len(self.sequences["chrA"]), "p")]
        for threads in (1, 2):
            output = self.root / f"empty_{threads}"
            report = analyze_regions(self.control, self.test, self.ref, [Region("chrEmpty", 0, 3)],
                                     output, threads=threads, keep_patterns=True)
            self.assertEqual(report["n_units"], 0)
            self.assertEqual(report["n_regions"], 1)
        for invalid in ([], [regions[0], regions[0]], [Region("chrA", 0, 10), Region("chrA", 5, 15)],
                        [Region("missing", 0, 10)], [Region("chrA", 0, 1000)]):
            with self.assertRaises(ValueError):
                analyze_regions(self.control, self.test, self.ref, invalid, self.root / "invalid", threads=2)
            self.assertFalse((self.root / "invalid").exists())
        for kwargs in ({"max_reads_per_region": 1}, {"modification": "5hmC"},
                       {"max_total_bases": 0}, {"threads": 0}):
            with self.assertRaises(ValueError):
                analyze_regions(self.control, self.test, self.ref, regions, self.root / "failed", chunk_units=1, **kwargs)
            self.assertFalse((self.root / "failed").exists())
            self.assertFalse(list(self.root.glob(".failed.*")))
        existing = self.root / "protected"
        existing.mkdir()
        (existing / "sentinel").write_text("keep")
        with self.assertRaises(FileExistsError):
            analyze_regions(self.control, self.test, self.ref, regions, existing)
        self.assertEqual((existing / "sentinel").read_text(), "keep")

    def test_dense_limit_applies_per_chunk_and_bed_change_rolls_back(self):
        from unittest.mock import patch
        region = Region("chrA", 0, len(self.sequences["chrA"]), "dense")
        with patch("methyldiver.genome.MAX_DENSE_COUNT_CELLS", 32):
            report = analyze_regions(self.control, self.test, self.ref, [region], self.root / "dense")
        self.assertEqual(report["chunk_units"], 2)
        self.assertEqual(report["n_chunks"], 3)
        bed = self.root / "mutable.bed"
        bed.write_text(f"chrA\t0\t{region.end}\tdense\n")
        def change_bed(done, total, units):
            if done:
                bed.write_text(f"chrA\t0\t{region.end}\tchanged\n")
        with self.assertRaisesRegex(ValueError, "BED changed"):
            analyze_regions(self.control, self.test, self.ref, [region], self.root / "changed",
                            regions_bed=bed, progress=change_bed)
        self.assertFalse((self.root / "changed").exists())
        self.assertFalse(list(self.root.glob(".changed.*")))


if __name__ == "__main__":
    unittest.main()
