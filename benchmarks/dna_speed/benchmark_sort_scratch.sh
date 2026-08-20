#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: benchmark_sort_scratch.sh UNSORTED_BAM OUTDIR" >&2
    exit 2
}

[[ $# -eq 2 ]] || usage
input=$(realpath -s "$1")
outdir=$(realpath -sm "$2")
[[ -r $input ]] || {
    echo "Input is not readable: $input" >&2
    exit 1
}

threads=${SORT_THREADS:-4}
memory=${SORT_MEM:-16M}
local_root=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/dna-speed-sort-${SLURM_JOB_ID:-$$}
mkdir -p "$outdir" "$outdir/nfs_tmp" "$local_root"
trap 'rm -rf "$local_root"' EXIT

summary="$outdir/summary.tsv"
printf 'temp_mode\tthreads\tmemory_per_thread\twall_seconds\tmax_rss_kb\ttemp_files\tview_md5\n' > "$summary"

run_case() {
    local mode=$1 temp_root case_dir started ended max_rss temp_files view_md5
    case_dir="$outdir/$mode"
    mkdir -p "$case_dir"
    if [[ $mode == local ]]; then
        temp_root="$local_root"
    else
        temp_root="$outdir/nfs_tmp"
    fi

    started=$(date +%s)
    /usr/bin/time -v -o "$case_dir/time.txt" \
        samtools sort -@ "$threads" -m "$memory" -T "$temp_root/sort" \
            -O bam -o "$case_dir/output.bam" "$input" \
            2> "$case_dir/stderr.txt"
    ended=$(date +%s)
    samtools quickcheck -v "$case_dir/output.bam"
    samtools view -@ "$threads" "$case_dir/output.bam" | md5sum > "$case_dir/view.md5"

    max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$case_dir/time.txt")
    temp_files=$(/usr/bin/awk '/merging from/{print $4; exit}' "$case_dir/stderr.txt")
    view_md5=$(cut -d' ' -f1 "$case_dir/view.md5")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$mode" "$threads" "$memory" \
        "$((ended - started))" "$max_rss" "${temp_files:-unknown}" "$view_md5" \
        | tee -a "$summary"
}

run_case nfs
run_case local

if [[ $(awk 'NR == 2 {print $7}' "$summary") != $(awk 'NR == 3 {print $7}' "$summary") ]]; then
    echo "NFS and local sort outputs differ" >&2
    exit 1
fi

echo "NFS and local sort records match exactly."
