"""Create tiny, indexed synthetic modBAM inputs with known pattern counts.

Run: python examples/make_modbam_example.py --output outputs/toy_inputs
The destination must not exist. No human data or legacy paths are used.
"""

import argparse
from array import array
from pathlib import Path

import pysam


def create_example(destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    sequence = "AACGAACGAACGAACGAA"
    reference = destination / "reference.fa"
    reference.write_text(">chrToy\n" + sequence + "\n", encoding="utf-8")
    pysam.faidx(str(reference))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chrToy", "LN": len(sequence)}]}
    for sample, patterns in (("control", ["0000", "1111"]), ("test", ["0011", "1100"])):
        path = destination / f"{sample}.bam"
        with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
            for i in range(20):
                pattern = patterns[i // 10]
                reverse = bool(i % 2)
                read = pysam.AlignedSegment(bam.header)
                read.query_name = f"{sample}_synthetic_{i}"
                read.query_sequence = sequence
                read.flag = 16 if reverse else 0
                read.reference_id = 0
                read.reference_start = 0
                read.mapping_quality = 60
                read.cigarstring = f"{len(sequence)}M"
                read.query_qualities = array("B", [30] * len(sequence))
                read.set_tag("MM", "C+m?,0,0,0,0;")
                # MM/ML follow original read orientation; reverse genomic order
                # when the alignment sequence is reverse-complemented.
                values = pattern[::-1] if reverse else pattern
                read.set_tag("ML", array("B", [255 if x == "1" else 0 for x in values]))
                read.set_tag("MN", len(sequence))
                bam.write(read)
        pysam.index(str(path))
    (destination / "regions.bed").write_text(f"chrToy\t0\t{len(sequence)}\ttoy\n", encoding="utf-8")
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(create_example(args.output).resolve())
