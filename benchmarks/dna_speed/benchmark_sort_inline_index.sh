#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: benchmark_sort_inline_index.sh UNSORTED_BAM OUTDIR" >&2
    exit 2
}

[[ $# -eq 2 ]] || usage

input=$(realpath -s "$1")
outdir=$(realpath -sm "$2")
threads=${SLURM_CPUS_PER_TASK:-$(nproc)}
memory=${SORT_MEM:-1G}
[[ -r $input && $threads =~ ^[1-9][0-9]*$ ]] || usage
samtools quickcheck -v "$input"

mkdir -p "$outdir/results"
tmp_root=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/dna-inline-index-${SLURM_JOB_ID:-$$}
mkdir -p "$tmp_root"
trap 'rm -rf "$tmp_root"' EXIT

summary="$outdir/summary.tsv"
printf 'case\tthreads\twall_seconds\tmax_rss_kb\tbam_bytes\tindex_bytes\tview_md5\tidxstats_md5\n' > "$summary"

run_current() {
    samtools sort -@ "$threads" -m "$memory" -T "$tmp_root/current" \
        -O bam -o "$case_dir/output.bam" "$input"
    samtools index "$case_dir/output.bam" "$case_dir/output.bam.bai"
}

run_inline() {
    samtools sort -@ "$threads" -m "$memory" -T "$tmp_root/inline" \
        --write-index -O bam \
        -o "$case_dir/output.bam##idx##$case_dir/output.bam.bai" "$input"
}
export -f run_current run_inline
export input threads memory tmp_root

reference_view=""
reference_idxstats=""
for case_name in current inline; do
    case_dir="$outdir/results/$case_name"
    mkdir -p "$case_dir"
    export case_dir
    started=$(date +%s)
    /usr/bin/time -v -o "$case_dir/time.txt" \
        bash -o pipefail -c "run_$case_name" \
        2> "$case_dir/pipeline.stderr"
    ended=$(date +%s)

    samtools quickcheck -v "$case_dir/output.bam"
    view_md5=$(samtools view -@ "$threads" "$case_dir/output.bam" \
        | md5sum | cut -d' ' -f1)
    idxstats_md5=$(samtools idxstats "$case_dir/output.bam" \
        | md5sum | cut -d' ' -f1)
    if [[ -z $reference_view ]]; then
        reference_view=$view_md5
        reference_idxstats=$idxstats_md5
    elif [[ $view_md5 != "$reference_view" \
            || $idxstats_md5 != "$reference_idxstats" ]]; then
        echo "$case_name output differs from the current path" >&2
        exit 1
    fi

    max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/ {
        gsub(/^[[:space:]]+/, "", $2); print $2
    }' "$case_dir/time.txt")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$case_name" "$threads" "$((ended - started))" "$max_rss" \
        "$(stat -Lc %s "$case_dir/output.bam")" \
        "$(stat -Lc %s "$case_dir/output.bam.bai")" \
        "$view_md5" "$idxstats_md5" | tee -a "$summary"
done

echo "Inline and post-write indexing produced identical records and idxstats."
