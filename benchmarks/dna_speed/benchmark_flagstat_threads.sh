#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: benchmark_flagstat_threads.sh BAM OUTDIR" >&2
    exit 2
}

[[ $# -eq 2 ]] || usage

input=$(realpath -s "$1")
outdir=$(realpath -sm "$2")
allocated=${SLURM_CPUS_PER_TASK:-$(nproc)}
[[ -r $input && $allocated =~ ^[1-9][0-9]*$ ]] || usage
samtools quickcheck -v "$input"
mkdir -p "$outdir/results"

summary="$outdir/summary.tsv"
printf 'case\tadditional_threads\twall_seconds\tmax_rss_kb\toutput_md5\n' > "$summary"

reference_md5=""
for spec in default:0 threaded:$allocated; do
    IFS=: read -r case_name extra_threads <<< "$spec"
    case_dir="$outdir/results/$case_name"
    mkdir -p "$case_dir"
    started=$(date +%s)
    /usr/bin/time -v -o "$case_dir/time.txt" \
        samtools flagstat -@ "$extra_threads" "$input" \
        > "$case_dir/flagstat.txt" 2> "$case_dir/pipeline.stderr"
    ended=$(date +%s)

    output_md5=$(md5sum "$case_dir/flagstat.txt" | cut -d' ' -f1)
    if [[ -z $reference_md5 ]]; then
        reference_md5=$output_md5
    elif [[ $output_md5 != "$reference_md5" ]]; then
        echo "Threaded flagstat output differs" >&2
        exit 1
    fi
    max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/ {
        gsub(/^[[:space:]]+/, "", $2); print $2
    }' "$case_dir/time.txt")
    printf '%s\t%s\t%s\t%s\t%s\n' "$case_name" "$extra_threads" \
        "$((ended - started))" "$max_rss" "$output_md5" | tee -a "$summary"
done

echo "Default and threaded flagstat outputs match exactly."
