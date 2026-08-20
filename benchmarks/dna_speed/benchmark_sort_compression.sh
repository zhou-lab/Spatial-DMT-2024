#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_sort_compression.sh UNSORTED_BAM OUTDIR

Environment:
  SORT_LEVELS  Space-separated NAME:LEVEL cases. Use LEVEL=default to omit -l.
               Default: default:default level1:1 level3:3 level6:6
  SORT_MEM     Memory per samtools sort thread [1G].
EOF
    exit 2
}

[[ $# -eq 2 ]] || usage

input=$(realpath -s "$1")
outdir=$(realpath -sm "$2")
cases=${SORT_LEVELS:-"default:default level1:1 level3:3 level6:6"}
threads=${SLURM_CPUS_PER_TASK:-$(nproc)}
memory=${SORT_MEM:-1G}

[[ -r $input && $threads =~ ^[1-9][0-9]*$ ]] || usage
samtools quickcheck -v "$input"
mkdir -p "$outdir/results"

tmp_root=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/dna-sort-compression-${SLURM_JOB_ID:-$$}
mkdir -p "$tmp_root"
trap 'rm -rf "$tmp_root"' EXIT

summary="$outdir/summary.tsv"
printf 'case\tcompression_level\tthreads\twall_seconds\tmax_rss_kb\toutput_bytes\tview_md5\n' > "$summary"

reference_md5=""
for spec in $cases; do
    IFS=: read -r name level <<< "$spec"
    [[ -n $name && ($level == default || $level =~ ^[0-9]$) ]] || usage

    case_dir="$outdir/results/$name"
    mkdir -p "$case_dir"
    output="$case_dir/output.bam"
    level_args=()
    if [[ $level != default ]]; then
        level_args=(-l "$level")
    fi

    started=$(date +%s)
    /usr/bin/time -v -o "$case_dir/time.txt" \
        samtools sort -@ "$threads" -m "$memory" "${level_args[@]}" \
            -T "$tmp_root/$name" -O bam -o "$output" "$input" \
            2> "$case_dir/pipeline.stderr"
    ended=$(date +%s)

    samtools quickcheck -v "$output"
    view_md5=$(samtools view -@ "$threads" "$output" | md5sum | cut -d' ' -f1)
    if [[ -z $reference_md5 ]]; then
        reference_md5=$view_md5
    elif [[ $view_md5 != "$reference_md5" ]]; then
        echo "$name records differ from the first case" >&2
        exit 1
    fi

    max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/ {
        gsub(/^[[:space:]]+/, "", $2); print $2
    }' "$case_dir/time.txt")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$name" "$level" "$threads" "$((ended - started))" "$max_rss" \
        "$(stat -Lc %s "$output")" "$view_md5" | tee -a "$summary"
done

echo "All sort compression cases produced identical SAM records."
