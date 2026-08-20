#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: validate_e2e_pipeline.sh OUTPUT_DIR [SAMPLE]" >&2
    exit 2
}

[[ $# -ge 1 && $# -le 2 ]] || usage
output_dir=$(realpath -s "$1")
sample=${2:-SpeedSmoke21}
sample_dir="$output_dir/$sample"
trim_dir="$sample_dir/trim"
cell_dir="$sample_dir/bam"
main_bam="$sample_dir/bam_merged/${sample}_merged_dedup.bam"
unmatched_bam="$sample_dir/bam_unmatched/${sample}_unmatched.bam"
lambda_bam="$sample_dir/bam_lambda/${sample}_Lambda.bam"
trim_stats="$trim_dir/${sample}_trim_stats.txt"
main_dup="$sample_dir/bam_merged/${sample}_merged_dedup.dupsifter.stat"
unmatched_dup="$sample_dir/bam_unmatched/${sample}_unmatched.dupsifter.stat"

for path in "$main_bam" "$unmatched_bam" "$lambda_bam" "$trim_stats" \
        "$main_dup" "$unmatched_dup"; do
    [[ -s $path ]] || usage
done

declare -A expected_fastq_md5=(
    [R1]=a4a6e04dd9dc0d0240d08b40acbf7313
    [R2]=b4963e39c59e150ee26fe425b841fccf
    [Unmatched_R1]=179a18a1b1d142d3f3142ad567ed3f2b
    [Unmatched_R2]=a42e8e046c3210601fbeeb4640f4c577
    [Lambda_R1]=6d8650ffd790dde02b38e39d37089c8b
    [Lambda_R2]=64eedf69df99127ea5bb2ed0bef1dc32
)
expected_stats_sha=d714fdef15e915440b1f0dfb2a19082593a2700ea496f003502d8be6d1a26dd5
expected_main_dup_sha=6fa2168946eec4748911fce68e31a686eee122aad26a913f1668617f380c76c9
expected_unmatched_dup_sha=e43fcbd91a9cfbbe2c15e9ba099227cd6dec3b3f0a27d246bcdba9eb5b126c58

for label in R1 R2 Unmatched_R1 Unmatched_R2 Lambda_R1 Lambda_R2; do
    fastq="$trim_dir/${sample}_${label}.fq.gz"
    observed=$(pigz -dc -p 2 "$fastq" | md5sum | cut -d' ' -f1)
    [[ $observed == "${expected_fastq_md5[$label]}" ]] || {
        echo "$label decompressed MD5 differs: $observed" >&2
        exit 1
    }
done

observed_stats_sha=$(sha256sum "$trim_stats" | cut -d' ' -f1)
[[ $observed_stats_sha == "$expected_stats_sha" ]] || {
    echo "trim stats SHA-256 differs: $observed_stats_sha" >&2
    exit 1
}
[[ $(sha256sum "$main_dup" | cut -d' ' -f1) == "$expected_main_dup_sha" ]]
[[ $(sha256sum "$unmatched_dup" | cut -d' ' -f1) == "$expected_unmatched_dup_sha" ]]
if compgen -G "$sample_dir/tmp/trim/_out_*_barcodes_all.txt" >/dev/null; then
    echo "optional per-block barcode tables were not suppressed" >&2
    exit 1
fi

samtools quickcheck -v "$main_bam" "$unmatched_bam" "$lambda_bam"
find "$cell_dir" -maxdepth 1 -type f -name '[0-9][0-9][0-9][0-9].bam' \
    -print0 | xargs -0 -r samtools quickcheck -v

cell_bams=$(find "$cell_dir" -maxdepth 1 -type f \
    -name '[0-9][0-9][0-9][0-9].bam' | wc -l)
cell_csi=$(find "$cell_dir" -maxdepth 1 -type f \
    -name '[0-9][0-9][0-9][0-9].bam.csi' | wc -l)
cell_bai=$(find "$cell_dir" -maxdepth 1 -type f \
    -name '[0-9][0-9][0-9][0-9].bam.bai' | wc -l)
align_cells=$(wc -l < "$sample_dir/align_list.txt")
[[ $cell_bams == 2500 && $cell_csi == 2500 && $cell_bai == 0 \
    && $align_cells == 2500 ]] || {
    echo "unexpected cell/index counts: bam=$cell_bams csi=$cell_csi bai=$cell_bai align=$align_cells" >&2
    exit 1
}

fofn=$(mktemp)
trap 'rm -f "$fofn"' EXIT
find "$cell_dir" -maxdepth 1 -type f -name '[0-9][0-9][0-9][0-9].bam' \
    | sort > "$fofn"

main_records=$(samtools view -c "$main_bam")
cell_records=$(samtools cat -b "$fofn" | samtools view -c)
unmatched_records=$(samtools view -c "$unmatched_bam")
unmatched_primary=$(samtools flagstat "$unmatched_bam" \
    | awk '$4 == "primary" {print $1; exit}')
lambda_records=$(samtools view -c "$lambda_bam")

[[ $main_records == 5274032 && $cell_records == "$main_records" ]] || {
    echo "matched count mismatch: merged=$main_records cells=$cell_records" >&2
    exit 1
}
[[ $unmatched_records == 3850680 && $unmatched_primary == 3728358 ]] || {
    echo "unmatched count mismatch: total=$unmatched_records primary=$unmatched_primary" >&2
    exit 1
}
(( lambda_records > 0 )) || {
    echo "lambda BAM is empty" >&2
    exit 1
}

first_cell=$(head -n 1 "$fofn")
samtools idxstats "$first_cell" >/dev/null

printf 'metric\tvalue\n'
printf 'matched_records\t%s\n' "$main_records"
printf 'summed_cell_records\t%s\n' "$cell_records"
printf 'unmatched_records\t%s\n' "$unmatched_records"
printf 'unmatched_primary\t%s\n' "$unmatched_primary"
printf 'lambda_records\t%s\n' "$lambda_records"
printf 'cell_bams\t%s\n' "$cell_bams"
printf 'cell_csi\t%s\n' "$cell_csi"
