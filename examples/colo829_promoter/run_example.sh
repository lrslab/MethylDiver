#!/usr/bin/env bash
# Run from any directory with MethylDiver on PATH.
set -euo pipefail
example_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
output=${1:-MD_COLO829_example}
md_bin=${METHYLDIVER_BIN:-methyldiver}
if [[ -e "$output" ]]; then
    echo "Output already exists: $output. Choose a new output directory." >&2
    exit 1
fi
command -v "$md_bin" >/dev/null
mkdir -p -- "$output"
for sample in COLO829BL_1 COLO829_1; do
    "$md_bin" single \
        --input "$example_dir/inputs/$sample.bam" \
        --reference "$example_dir/inputs/chr12.fa.gz" \
        --regions-bed "$example_dir/inputs/promoter.bed" \
        --min-depth 5 \
        --output "$output/single_$sample"
done
"$md_bin" compare \
    --control "$example_dir/inputs/COLO829BL_1.bam" \
    --test "$example_dir/inputs/COLO829_1.bam" \
    --reference "$example_dir/inputs/chr12.fa.gz" \
    --regions-bed "$example_dir/inputs/promoter.bed" \
    --min-depth 5 \
    --output "$output/compare"
"$md_bin" discover \
    --input "$output/compare" \
    --regions-bed "$example_dir/inputs/promoter.bed" \
    --output "$output/discover_promoter"
printf '\nExample complete: %s\n' "$output"
