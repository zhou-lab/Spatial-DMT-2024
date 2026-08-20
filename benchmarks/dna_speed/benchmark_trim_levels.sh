#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_trim_levels.sh R1.fq.gz R2.fq.gz WHITELIST OUTDIR

Environment:
  COMPRESSION_LEVELS  whitespace-separated gzip levels [1 3 6]
  WORKERS             parallel trim processes [allocated CPUs]
  BLOCK_SIZE          GNU Parallel block size [512M]
  CLEAN_OUTPUTS       delete per-case FASTQ shards after checksumming [1]
EOF
    exit 2
}

[[ $# -eq 4 ]] || usage
r1=$(realpath -s "$1")
r2=$(realpath -s "$2")
whitelist=$(realpath -s "$3")
outdir=$(realpath -sm "$4")
script_dir=$(realpath -s "$(dirname "$0")")
if [[ -n ${SLURM_SUBMIT_DIR:-} \
        && -x $SLURM_SUBMIT_DIR/benchmarks/dna_speed/benchmark_trim_streaming.sh ]]; then
    script_dir=$(realpath -s "$SLURM_SUBMIT_DIR/benchmarks/dna_speed")
fi
levels=${COMPRESSION_LEVELS:-"1 3 6"}
workers=${WORKERS:-${SLURM_CPUS_PER_TASK:-1}}
block_size=${BLOCK_SIZE:-512M}
clean_outputs=${CLEAN_OUTPUTS:-1}

[[ $workers =~ ^[1-9][0-9]*$ ]] || usage
[[ $clean_outputs == 0 || $clean_outputs == 1 ]] || usage
mkdir -p "$outdir"
summary="$outdir/summary.tsv"
printf 'compression_level\tworkers\tblock_size\tblocks\twall_seconds\tmax_rss_kb\tcompressed_bytes\tstats_sha256\tfastq_sha256\n' \
    > "$summary"

for level in $levels; do
    [[ $level =~ ^[0-9]$ ]] || usage
    case_dir="$outdir/level_$level"
    if [[ -e $case_dir ]]; then
        echo "Case directory already exists: $case_dir" >&2
        exit 1
    fi

    WORKERS="$workers" BLOCK_SIZE="$block_size" COMPRESSION_LEVEL="$level" \
        SKIP_BARCODE_COUNTS=1 \
        "$script_dir/benchmark_trim_streaming.sh" \
        "$r1" "$r2" "$whitelist" "$case_dir"

    awk -F '\t' -v level="$level" 'NR == 2 {
        print level "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8
    }' "$case_dir/summary.tsv" >> "$summary"

    if [[ $clean_outputs == 1 ]]; then
        find "$case_dir/output" -type f -delete
        rmdir "$case_dir/output"
    fi
done

reference_stats=$(awk -F '\t' 'NR == 2 {print $8}' "$summary")
reference_fastq=$(awk -F '\t' 'NR == 2 {print $9}' "$summary")
awk -F '\t' -v stats="$reference_stats" -v fastq="$reference_fastq" \
    'NR > 1 && ($8 != stats || $9 != fastq) {exit 1}' "$summary" || {
        echo "Compression-level outputs differ" >&2
        exit 1
    }

cat "$summary"
