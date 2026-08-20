#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_trim_blocks.sh R1.fq.gz R2.fq.gz WHITELIST OUTDIR

Environment:
  BLOCK_SIZES    whitespace-separated GNU Parallel block sizes
                 [64M 128M 256M 512M]
  WORKERS        parallel trim processes [allocated CPUs]
  CLEAN_OUTPUTS  delete per-case FASTQ shards after checksumming [1]
EOF
    exit 2
}

[[ $# -eq 4 ]] || usage
r1=$(realpath -s "$1")
r2=$(realpath -s "$2")
whitelist=$(realpath -s "$3")
outdir=$(realpath -sm "$4")
script_dir=$(realpath -s "$(dirname "$0")")
if [[ -n ${SLURM_SUBMIT_DIR:-} && -x $SLURM_SUBMIT_DIR/benchmarks/dna_speed/benchmark_trim_streaming.sh ]]; then
    script_dir=$(realpath -s "$SLURM_SUBMIT_DIR/benchmarks/dna_speed")
fi
block_sizes=${BLOCK_SIZES:-"64M 128M 256M 512M"}
workers=${WORKERS:-${SLURM_CPUS_PER_TASK:-1}}
clean_outputs=${CLEAN_OUTPUTS:-1}

[[ $workers =~ ^[1-9][0-9]*$ ]] || usage
[[ $clean_outputs == 0 || $clean_outputs == 1 ]] || usage
mkdir -p "$outdir"
summary="$outdir/summary.tsv"
printf 'block_size\tworkers\tblocks\twall_seconds\tmax_rss_kb\tcompressed_bytes\tstats_sha256\tfastq_sha256\n' > "$summary"

for block_size in $block_sizes; do
    case_dir="$outdir/$block_size"
    if [[ -e $case_dir ]]; then
        echo "Case directory already exists: $case_dir" >&2
        exit 1
    fi

    WORKERS="$workers" BLOCK_SIZE="$block_size" \
        "$script_dir/benchmark_trim_streaming.sh" \
        "$r1" "$r2" "$whitelist" "$case_dir"

    awk -F '\t' -v size="$block_size" 'NR == 2 {
        print size "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8
    }' "$case_dir/summary.tsv" >> "$summary"

    if [[ $clean_outputs == 1 ]]; then
        find "$case_dir/output" -type f -delete
        rmdir "$case_dir/output"
    fi
done

cat "$summary"
