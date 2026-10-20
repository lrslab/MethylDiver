"""Compare a local example run with the supplied expected tables (no extra dependencies)."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def check_table(actual, expected):
    with actual.open(newline='') as a, expected.open(newline='') as e:
        ar, er = csv.DictReader(a, delimiter='\t'), csv.DictReader(e, delimiter='\t')
        assert ar.fieldnames == er.fieldnames, f'Column mismatch: {actual}'
        aa, ee = list(ar), list(er)
    assert len(aa) == len(ee), f'Row count mismatch: {actual}'
    for i, (a, e) in enumerate(zip(aa, ee), 1):
        for key in a:
            if a[key] == e[key]:
                continue
            try:
                same = math.isclose(float(a[key]), float(e[key]), rel_tol=1e-9, abs_tol=1e-10)
            except (TypeError, ValueError):
                same = False
            assert same, f'Mismatch in {actual}, row {i}, column {key}: {a[key]} != {e[key]}'
    print(f'PASS {actual.name}: {len(aa)} rows ({actual.parent.name})')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='Output directory passed to run_example.sh')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / 'inputs/input_manifest.json').read_text())
    for name, info in manifest['files'].items():
        path = root / 'inputs' / name
        h = hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                h.update(block)
        assert h.hexdigest() == info['sha256'], f'Input checksum mismatch: {name}'
    print('PASS input checksums')
    tables = [
        'single_COLO829BL_1/two_metrics.tsv',
        'single_COLO829_1/two_metrics.tsv',
        'compare/three_metrics.tsv',
        'discover_promoter/ranked_regions.tsv',
        'discover_promoter/window_statistics.tsv',
    ]
    for name in tables:
        check_table(args.output / name, root / 'expected' / name)
    bed = 'discover_promoter/top_regions.bed'
    assert (args.output / bed).read_bytes() == (root / 'expected' / bed).read_bytes(), 'Ranked BED differs'
    print('All expected result tables match.')


if __name__ == '__main__':
    main()
