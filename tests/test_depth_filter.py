"""Depth filtering must survive all CLI extraction and retention paths."""
from array import array
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam


class TestDepthFilter(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.names = ['exact', 'low_control', 'low_test', 'missing']
        self.ref = self.root / 'ref.fa'
        self.ref.write_text(''.join(f'>{name}\nCGCGCG\n' for name in self.names))
        pysam.faidx(str(self.ref))
        self.bed = self.root / 'regions.bed'
        self.bed.write_text(''.join(f'{name}\t0\t6\t{name}\n' for name in self.names))
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': n, 'LN': 6} for n in self.names]}
        for role, depths in [('control', [5, 4, 5, 0]), ('test', [5, 5, 4, 5])]:
            bam_path = self.root / f'{role}.bam'
            with pysam.AlignmentFile(str(bam_path), 'wb', header=header) as bam:
                for chrom, depth in enumerate(depths):
                    # Extra uncertain reads must not count towards --min-depth.
                    for i in range(depth + 6):
                        r = pysam.AlignedSegment(bam.header)
                        r.query_name = f'{role}_{chrom}_{i}'
                        r.query_sequence = 'CGCGCG'
                        r.reference_id = chrom
                        r.reference_start = 0
                        r.mapping_quality = 60
                        r.cigarstring = '6M'
                        r.query_qualities = array('B', [30]*6)
                        r.set_tag('MM', 'C+m?,0,0,0;')
                        r.set_tag('ML', array('B', [0, 0 if i < depth else 128, 0]))
                        bam.write(r)
            pysam.index(str(bam_path))

    def test_cli_depth_boundary_all_modes(self):
        for single in (False, True):
            for mode, threads in [('region', 1), ('region', 2), ('bed', 2), ('genome', 1), ('genome', 2)]:
                for keep in (None, 'F', 'T') if mode == 'region' and threads == 1 else ('F', 'T'):
                    with self.subTest(single=single, mode=mode, threads=threads, keep=keep):
                        out = self.root / f'out_{single}_{mode}_{threads}_{keep}'
                        cmd = [sys.executable, '-m', 'methyldiver', 'single' if single else 'compare',
                               '--reference', str(self.ref), '--output', str(out),
                               '--N', '3', '--min-depth', '5', '--threads', str(threads)]
                        if single:
                            cmd += ['--input', str(self.root/'control.bam'), '--sample-label', 'D0']
                        else:
                            cmd += ['--control', str(self.root/'control.bam'), '--test', str(self.root/'test.bam'),
                                    '--control-label', 'D0', '--test-label', 'D10']
                        if mode == 'genome':
                            cmd += ['--genome']
                        elif mode == 'bed':
                            cmd += ['--regions-bed', str(self.bed)]
                        else:
                            for name in self.names:
                                cmd += ['--region', f'{name}:0-6']
                        if keep is not None:
                            cmd += ['--keep-patterns', keep, '--chunk-units', '1']
                        run = subprocess.run(cmd, capture_output=True, text=True)
                        self.assertEqual(run.returncode, 0, run.stderr)
                        name = 'two_metrics.tsv' if single else 'three_metrics.tsv'
                        with (out/name).open() as handle:
                            rows = list(csv.DictReader(handle, delimiter='\t'))
                        self.assertEqual([r['chrom'] for r in rows], ['exact', 'low_test'] if single else ['exact'])
                        # True zero estimates at exactly the threshold must survive.
                        self.assertTrue(all(float(r['D0_beta']) == 0 for r in rows))
                        summary = json.loads((out/'analysis_summary.json').read_text())
                        details = summary['outputs'][name]
                        self.assertEqual(details['rows'], 2 if single else 1)
                        self.assertEqual(details['excluded_missing_metrics'], 1)
                        self.assertEqual(details['excluded_below_min_depth'], 1 if single else 2)
                        self.assertEqual(sum(details[k] for k in ('rows','excluded_missing_metrics','excluded_below_min_depth')), 4)
                        if keep == 'F':
                            self.assertFalse((out/'complete_patterns').exists())
                            self.assertFalse((out/'incomplete_patterns').exists())
                        if keep == 'T':
                            read = subprocess.run([sys.executable, '-m', 'methyldiver', 'read',
                                '--input', str(out/'complete_patterns'), '--metrics-only', '--all'],
                                capture_output=True, text=True)
                            self.assertEqual(read.returncode, 0, read.stderr)
                            self.assertEqual(read.stdout, (out/name).read_text())
