"""Extract the real COLO829 promoter example from indexed public source files."""
import argparse
import hashlib
import json
from pathlib import Path

import pysam

CHROM, START, END = 'chr12', 53936830, 53941692
REGION_ID = 'HOXC13;HOXC13-AS'


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam-directory', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True,
                        help='New directory for extracted inputs; must not exist.')
    args = parser.parse_args()
    samples = ['COLO829BL_1', 'COLO829_1']
    sources = {name: args.bam_directory / f'{name}.bam' for name in samples}
    with pysam.FastaFile(str(args.reference)) as fasta:
        length = fasta.get_reference_length(CHROM)
        for path in sources.values():
            with pysam.AlignmentFile(str(path), 'rb') as bam:
                assert bam.has_index(), f'Missing BAM index: {path}'
                assert bam.get_reference_length(CHROM) == length
        args.output.mkdir(parents=True, exist_ok=False)
        reference = args.output / 'chr12.fa.gz'
        # BGZF supports indexed random access directly, without decompression.
        with pysam.BGZFile(str(reference), 'wb') as handle:
            handle.write(b'>chr12\n')
            for start in range(0, length, 60 * 10000):
                sequence = fasta.fetch(CHROM, start, min(start + 60 * 10000, length))
                handle.write(('\n'.join(sequence[i:i + 60] for i in range(0, len(sequence), 60)) + '\n').encode())
        pysam.faidx(str(reference))
        with pysam.FastaFile(str(reference)) as copied:
            for start in range(0, length, 1000000):
                stop = min(start + 1000000, length)
                assert copied.fetch(CHROM, start, stop) == fasta.fetch(CHROM, start, stop)
    counts = {}
    for sample, source in sources.items():
        target = args.output / f'{sample}.bam'
        with pysam.AlignmentFile(str(source), 'rb') as bam:
            with pysam.AlignmentFile(str(target), 'wb', header=bam.header) as out:
                count = 0
                for read in bam.fetch(CHROM, START, END):
                    out.write(read)
                    count += 1
            pysam.index(str(target))
            with pysam.AlignmentFile(str(target), 'rb') as subset:
                # Check all records, including sequence, orientation, CIGAR and tags.
                original = (r.to_string() for r in bam.fetch(CHROM, START, END))
                extracted = (r.to_string() for r in subset.fetch(CHROM, START, END))
                from itertools import zip_longest
                assert all(a == b for a, b in zip_longest(original, extracted))
        counts[sample] = count
    (args.output / 'promoter.bed').write_text(
        f'{CHROM}\t{START}\t{END}\t{REGION_ID}\n')
    metadata = {
        'region': {'chrom': CHROM, 'start': START, 'end': END, 'name': REGION_ID,
                   'coordinates': '0-based, end-exclusive'},
        'sample_alignment_counts_before_analysis_filters': counts,
        'reference': 'Complete chr12 from GRCh38_no_alt.fna; BGZF compressed; no masking or coordinate shifts.',
        'extraction': 'Every alignment overlapping the BED interval, no read subsampling or analysis filtering. Original headers and all alignment fields/tags retained.',
        'checks': ['All subset SAM records equal corresponding source SAM records.',
                   'Every chr12 reference base equal to source FASTA.'],
        'pysam_version': pysam.__version__,
        'files': {p.name: {'bytes': p.stat().st_size, 'sha256': sha256(p)}
                  for p in sorted(args.output.iterdir()) if p.is_file()},
    }
    (args.output / 'input_manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps({'alignment_counts': counts, 'input_bytes': sum(p.stat().st_size for p in args.output.iterdir())}, indent=2))


if __name__ == '__main__':
    main()
