from array import array
from collections import Counter
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

from methyldiver.modbam import Region, analyze_modbam, decode_calls, parse_region, read_regions_bed
from methyldiver.output import write_analysis


ROOT = Path(__file__).resolve().parents[1]
SEQUENCE = "AACGAACGAACGAA"
HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrToy", "LN": len(SEQUENCE)}]}


def make_read(name="read", pattern="101", reverse=False, opposite=False, sequence=SEQUENCE,
              cigar=None, start=0, mm=None, ml=None, no_tags=False, header=HEADER):
    read = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header))
    read.query_name = name
    read.query_sequence = sequence
    read.flag = 16 if reverse else 0
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigarstring = cigar or f"{len(sequence)}M"
    read.query_qualities = array("B", [30] * len(sequence))
    if not no_tags:
        group = "G-m?" if opposite else "C+m?"
        read.set_tag("MM", mm if mm is not None else group + "," + ",".join("0" for _ in pattern) + ";")
        values = pattern[::-1] if reverse else pattern
        read.set_tag("ML", array("B", ml if ml is not None else [255 if value == "1" else 0 for value in values]))
        read.set_tag("MN", len(sequence))
    return read


class TestModBAM(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.reference = self.root / "ref.fa"
        self.reference.write_text(">chrToy\n" + SEQUENCE + "\n")
        pysam.faidx(str(self.reference))

    def bam(self, name, reads, header=HEADER):
        path = self.root / name
        with pysam.AlignmentFile(str(path), "wb", header=header) as out:
            for read in sorted(reads, key=lambda r: r.reference_start):
                out.write(read)
        pysam.index(str(path))
        return path

    def analyze(self, control, test=None, regions=None, **kwargs):
        return analyze_modbam(control, test or control, self.reference,
                              regions or [Region("chrToy", 0, len(SEQUENCE))], **kwargs)

    def test_four_alignment_and_modification_orientations(self):
        reads = [make_read(f"r{reverse}{opposite}", "101", reverse, opposite)
                 for reverse in (False, True) for opposite in (False, True)]
        for read in reads:
            calls, strand = decode_calls(read, {2, 6, 10}, 0.8, 10, Counter())
            self.assertEqual(calls, {2: 1, 6: 0, 10: 1})
        bam = self.bam("orientations.bam", reads)
        result = self.analyze(bam)
        sample = result["units"][0]["control"]
        self.assertEqual(sample["counts"], [0, 0, 0, 0, 0, 4, 0, 0])
        self.assertEqual(sample["n_candidate_reads"], 4)
        self.assertEqual(sample["beta"], 2 / 3)
        self.assertEqual(result["units"][0]["comparison"]["D_epi"], 0)

    def test_known_distribution_shift_and_depth_flag(self):
        control = self.bam("control.bam", [make_read("c0", "000"), make_read("c1", "111", True)])
        test = self.bam("test.bam", [make_read("t0", "001"), make_read("t1", "110", True)])
        result = self.analyze(control, test, min_depth=3)
        unit = result["units"][0]
        self.assertEqual(unit["cpg_anchors"], [2, 6, 10])
        self.assertEqual((unit["start"], unit["end"]), (2, 12))
        self.assertEqual(unit["control"]["beta"], 0.5)
        self.assertEqual(unit["control"]["entropy"], 1 / 3)
        self.assertEqual(unit["comparison"], {"delta_beta": 0, "delta_entropy": 0, "D_epi": 1})
        self.assertFalse(unit["passes_min_depth"])
        self.assertEqual(unit["control"]["n_complete_reads"], 2)

    def test_quantization_confidence_boundaries(self):
        for byte, expected in ((0, 0), (50, 0), (51, None), (128, None), (204, None), (205, 1), (255, 1)):
            with self.subTest(byte=byte):
                read = make_read(ml=[byte] * 3)
                calls, _ = decode_calls(read, {2, 6, 10}, 0.8, 10, Counter())
                self.assertEqual(calls, dict.fromkeys((2, 6, 10), expected))

    def test_unknown_implicit_missing_tags_and_missing_ml(self):
        reads = []
        for i, mode in enumerate(("?", ".", "")):
            reads.append(make_read(f"skip{i}", mm=f"C+m{mode},0;", ml=[255]))
        absent = make_read("no_tags", no_tags=True)
        reads.append(absent)
        missing_ml = make_read("no_ml")
        missing_ml.set_tag("ML", None)
        reads.append(missing_ml)
        bam = self.bam("missing.bam", reads)
        sample = self.analyze(bam)["units"][0]["control"]
        self.assertEqual(sample["n_candidate_reads"], 5)
        self.assertEqual(sample["n_complete_reads"], 0)
        self.assertEqual(sample["incomplete_patterns"], {"...": 2, "1..": 3})
        self.assertIsNone(sample["beta"])
        ordinary = self.bam("ordinary.bam", [absent])
        with self.assertRaisesRegex(ValueError, "no MM/Mm"):
            self.analyze(ordinary)

    def test_cigar_deletion_skip_insertion_and_soft_clip(self):
        reads = [
            make_read("deletion", sequence=SEQUENCE[:6] + SEQUENCE[7:], cigar="6M1D7M", mm="C+m?,0,0;", ml=[255, 0]),
            make_read("skip", sequence=SEQUENCE[:6] + SEQUENCE[7:], cigar="6M1N7M", mm="C+m?,0,0;", ml=[255, 0]),
            make_read("insertion", sequence=SEQUENCE[:6] + "C" + SEQUENCE[6:], cigar="6M1I8M", mm="C+m?,0,0,0,0;", ml=[255, 255, 0, 255]),
            make_read("softclip", sequence="C" + SEQUENCE, cigar="1S14M", mm="C+m?,0,0,0,0;", ml=[255, 255, 0, 255]),
        ]
        sample = self.analyze(self.bam("cigars.bam", reads))["units"][0]["control"]
        self.assertEqual(sample["n_candidate_reads"], 4)
        self.assertEqual(sample["counts"][5], 2)
        self.assertEqual(sample["incomplete_patterns"], {"1.0": 2})

    def test_base_quality_and_reference_context(self):
        low = make_read("low")
        quality = list(low.query_qualities)
        quality[6] = 9
        low.query_qualities = array("B", quality)
        absent = make_read("no_qualities")
        absent.query_qualities = None
        # An extra non-reference-CpG C at position 0 cannot replace a missing site.
        variant = make_read("variant", sequence="C" + SEQUENCE[1:], mm="C+m?,0,0,0,0;", ml=[255, 255, 0, 255])
        sample = self.analyze(self.bam("quality.bam", [low, absent, variant]))["units"][0]["control"]
        self.assertEqual(sample["incomplete_patterns"], {"...": 1, "1.1": 1})
        self.assertEqual(sample["counts"][5], 1)

    def test_alignment_filters_and_partial_endpoint(self):
        reads = [make_read("good")]
        for flag in (256, 512, 1024, 2048):
            read = make_read(f"flag{flag}")
            read.flag |= flag
            reads.append(read)
        for mapq in (19, 255):
            read = make_read(f"mapq{mapq}")
            read.mapping_quality = mapq
            reads.append(read)
        partial = make_read("partial", pattern="01", start=3, sequence=SEQUENCE[3:])
        reads.append(partial)
        result = self.analyze(self.bam("filters.bam", reads))
        sample = result["units"][0]["control"]
        self.assertEqual(sample["n_candidate_reads"], 1)
        qc = result["manifest"]["regions"][0]["control_qc"]
        self.assertEqual(qc["reads_excluded_flags"], 4)
        self.assertEqual(qc["reads_excluded_mapq"], 2)

    def test_legacy_and_numeric_tag_encodings(self):
        for legacy in (False, True):
            for code in ("m", "27551"):
                read = make_read(mm=f"C+{code}?,0,0,0;", ml=[255, 0, 255])
                if legacy:
                    mm, ml = read.get_tag("MM"), read.get_tag("ML")
                    read.set_tag("MM", None)
                    read.set_tag("ML", None)
                    read.set_tag("Mm", mm)
                    read.set_tag("Ml", ml)
                calls, _ = decode_calls(read, {2, 6, 10}, 0.8, 10, Counter())
                self.assertEqual(calls, {2: 1, 6: 0, 10: 1})

    def test_unsupported_and_malformed_tags(self):
        reads = []
        for mm, ml in (("C+mh?,0;", [200, 40]), ("C+h?;", []),
                       ("C+m?,0;G-m?,0;", [255, 255]), ("C-m?,0;", [255]),
                       ("N+m?,0;", [255]), ("C+m?,0;C+m?,0;", [255, 255]),
                       ("C+m?,0;", [255, 0]), ("C+m?,0", [255]),
                       ("C+m?,1000;", [255])):
            reads.append(make_read(mm=mm, ml=ml))
        stale = make_read()
        stale.set_tag("MN", len(SEQUENCE) + 1)
        reads.append(stale)
        mixed = make_read()
        mixed.set_tag("Mm", mixed.get_tag("MM"))
        reads.append(mixed)
        orphan = make_read()
        orphan.set_tag("MM", None)
        reads.append(orphan)
        wrong_type = make_read()
        wrong_type.set_tag("ML", array("H", [255, 0, 255]))
        reads.append(wrong_type)
        for read in reads:
            with self.subTest(tags=read.get_tags()), self.assertRaises(ValueError):
                decode_calls(read, {2, 6, 10}, 0.8, 10, Counter())

    def test_duplicate_reads_molecules_and_unsupported_profiles(self):
        with self.assertRaisesRegex(ValueError, "duplicate primary"):
            self.analyze(self.bam("duplicate.bam", [make_read("same"), make_read("same")]))
        a, b = make_read("a"), make_read("b")
        a.set_tag("MI", "molecule")
        b.set_tag("MI", "molecule")
        with self.assertRaisesRegex(ValueError, "same MI"):
            self.analyze(self.bam("molecule.bam", [a, b]))
        for tag in ("paired", "duplex"):
            read = make_read()
            if tag == "paired":
                read.flag |= 1
            else:
                read.set_tag("dx", 1)
            with self.assertRaises(ValueError):
                self.analyze(self.bam(f"{tag}.bam", [read]))

    def test_empty_bam_and_region_boundary(self):
        bam = self.bam("empty.bam", [])
        unit = self.analyze(bam)["units"][0]
        self.assertIsNone(unit["control"]["callable_fraction"])
        self.assertIsNone(unit["comparison"]["D_epi"])
        result = self.analyze(bam, regions=[Region("chrToy", 0, 11)])
        # Last C is at 10, but its full CG dyad ends at 12 and is not contained.
        self.assertEqual(result["manifest"]["regions"][0]["n_reference_cpgs"], 2)
        self.assertEqual(result["units"], [])

    def test_missing_indexes_reference_mismatch_and_limits(self):
        bam = self.bam("reads.bam", [make_read("a"), make_read("b")])
        for config in ({"max_total_bases": 5}, {"max_reads_per_region": 1},
                       {"max_pattern_evaluations": 1}, {"call_threshold": 0.5},
                       {"call_threshold": float("nan")}, {"min_mapq": 255}, {"min_depth": 0}):
            with self.subTest(config=config), self.assertRaises(ValueError):
                self.analyze(bam, **config)
        for regions in ([Region("missing", 0, 1)], [Region("chrToy", 0, 15)],
                        [Region("chrToy", 0, 10, "a"), Region("chrToy", 5, 14, "b")]):
            with self.assertRaises(ValueError):
                self.analyze(bam, regions=regions)
        bad_header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chrToy", "LN": 15}]}
        mismatch = self.bam("wrong_reference.bam", [], bad_header)
        with self.assertRaisesRegex(ValueError, "length mismatch"):
            self.analyze(mismatch)
        Path(str(bam) + ".bai").unlink()
        with self.assertRaisesRegex(ValueError, "Missing BAM index"):
            self.analyze(bam)
        Path(str(self.reference) + ".fai").unlink()
        with self.assertRaisesRegex(ValueError, "Missing FASTA index"):
            self.analyze(mismatch)

    def test_region_parser_bed_and_multiple_intervals(self):
        self.assertEqual(parse_region("chrToy:0-14"), Region("chrToy", 0, 14))
        for value in ("chrToy", "chrToy:-1-5", "chrToy:10-5", "chrToy:0-0"):
            with self.assertRaises(ValueError):
                parse_region(value)
        bed = self.root / "regions.bed"
        bed.write_text("# intervals\nchrToy\t0\t7\tfirst\nchrToy\t7\t14\tsecond\n")
        regions = read_regions_bed(bed)
        bam = self.bam("regions.bam", [make_read()])
        result = self.analyze(bam, regions=regions)
        self.assertEqual(len(result["manifest"]["regions"]), 2)
        self.assertEqual(result["units"], [])  # independent interval boundaries are not bridged

    def test_output_roundtrip_and_preservation(self):
        bam = self.bam("result.bam", [make_read()])
        result = self.analyze(bam)
        output = self.root / "run"
        write_analysis(result, output)
        with (output / "metrics.tsv").open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["cpg_anchors"], "2,6,10")
        self.assertEqual(float(row["control_beta"]), 2 / 3)
        manifest = json.loads((output / "run.json").read_text())
        self.assertEqual(manifest["status"], "complete")
        for name, identity in manifest["outputs"].items():
            self.assertEqual(identity["sha256"], hashlib.sha256((output / name).read_bytes()).hexdigest())
        with (output / "patterns.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        counts = [int(row["count"]) for row in rows if row["sample"] == "control"]
        self.assertEqual(counts, [0, 0, 0, 0, 0, 1, 0, 0])
        with self.assertRaises(FileExistsError):
            write_analysis(result, output)

    def test_complete_cli_example_three_and_four_cpg(self):
        spec = importlib.util.spec_from_file_location("make_modbam_example", ROOT / "examples/make_modbam_example.py")
        example = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(example)
        inputs = example.create_example(self.root / "example")
        for k in (3, 4):
            output = self.root / f"cli_k{k}"
            command = [sys.executable, "-m", "methyldiver", "analyze", "--control", str(inputs / "control.bam"),
                       "--test", str(inputs / "test.bam"), "--reference", str(inputs / "reference.fa"),
                       "--regions-bed", str(inputs / "regions.bed"), "--output", str(output), "--k", str(k)]
            run = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            with (output / "metrics.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 5 - k)
            for row in rows:
                self.assertEqual(float(row["D_epi"]), 1)
                self.assertEqual(float(row["delta_beta"]), 0)
                self.assertEqual(float(row["delta_entropy"]), 0)
                self.assertEqual(int(row["control_n_complete_reads"]), 20)
            again = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(again.returncode, 2)
            self.assertIn("already exists", again.stderr)

    def test_cli_failure_does_not_publish_output(self):
        ordinary = self.bam("ordinary.bam", [make_read(no_tags=True)])
        output = self.root / "should_not_exist"
        run = subprocess.run([sys.executable, "-m", "methyldiver", "analyze", "--control", str(ordinary),
                              "--test", str(ordinary), "--reference", str(self.reference),
                              "--region", "chrToy:0-14", "--output", str(output)], capture_output=True, text=True)
        self.assertEqual(run.returncode, 2)
        self.assertIn("no MM/Mm", run.stderr)
        self.assertNotIn("Traceback", run.stderr)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
