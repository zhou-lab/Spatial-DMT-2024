#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_pileup_parallelism.sh BAM_LIST REF_FASTA OUTDIR

Environment:
  PILEUP_CASES  Space-separated NAME:JOBS:THREADS[:WINDOW_STEP] cases.
                WINDOW_STEP defaults to BISCUIT's 100000.
                Default: jobs2_t1:2:1 jobs1_t2:1:2

Each nonempty BAM_LIST line is a BAM path. All cases must fit within
SLURM_CPUS_PER_TASK. Non-header VCF records are checked per BAM.
EOF
    exit 2
}

[[ $# -eq 3 ]] || usage

bam_list=$(realpath -s "$1")
ref_fasta=$(realpath -s "$2")
outdir=$(realpath -sm "$3")
cases=${PILEUP_CASES:-"jobs2_t1:2:1 jobs1_t2:1:2"}
allocated_cpus=${SLURM_CPUS_PER_TASK:-$(nproc)}

[[ -r $bam_list && -r $ref_fasta ]] || usage
while read -r bam; do
    [[ -z $bam || $bam == \#* ]] && continue
    [[ -r $bam ]] || {
        echo "BAM is not readable: $bam" >&2
        exit 1
    }
    samtools quickcheck -v "$bam"
done < "$bam_list"

mkdir -p "$outdir/results"
summary="$outdir/summary.tsv"
printf 'case\tjobs\tthreads_per_cell\twindow_step\twall_seconds\tmax_rss_kb\toutput_bytes\n' > "$summary"

pileup_one() {
    local bam=$1 label
    label=$(basename "$bam" .bam)
    biscuit pileup -m 0 -a 0 -c -u -p -@ "$pileup_threads" \
        -s "$window_step" \
        "$ref_fasta" "$bam" > "$case_dir/$label.vcf"
}
export -f pileup_one
export ref_fasta

first_case=""
for spec in $cases; do
    IFS=: read -r name jobs pileup_threads window_step <<< "$spec"
    window_step=${window_step:-100000}
    [[ -n $name && $jobs =~ ^[1-9][0-9]*$ \
        && $pileup_threads =~ ^[1-9][0-9]*$ \
        && $window_step =~ ^[1-9][0-9]*$ ]] || usage
    (( jobs * pileup_threads <= allocated_cpus )) || {
        echo "$name requests $((jobs * pileup_threads)) CPUs; $allocated_cpus allocated" >&2
        exit 2
    }

    case_dir="$outdir/results/$name"
    rm -rf "$case_dir"
    mkdir -p "$case_dir"
    export case_dir pileup_threads window_step

    started=$(date +%s)
    /usr/bin/time -v -o "$case_dir/time.txt" \
        parallel --plain --halt soon,fail=1 -j "$jobs" \
            pileup_one :::: "$bam_list" \
        2> "$case_dir/pipeline.stderr"
    ended=$(date +%s)

    while read -r bam; do
        [[ -z $bam || $bam == \#* ]] && continue
        label=$(basename "$bam" .bam)
        /usr/bin/awk '!/^#/' "$case_dir/$label.vcf" \
            | LC_ALL=C sort -k1,1V -k2,2n \
            | md5sum \
            | /usr/bin/awk -v label="$label" '{print label "\t" $1}'
    done < "$bam_list" > "$case_dir/checksums.tsv"

    if [[ -z $first_case ]]; then
        first_case=$name
    else
        cmp "$outdir/results/$first_case/checksums.tsv" "$case_dir/checksums.tsv"
    fi

    max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/ {
        gsub(/^[[:space:]]+/, "", $2); print $2
    }' "$case_dir/time.txt")
    output_bytes=$(find "$case_dir" -type f -name '*.vcf' -printf '%s\n' \
        | /usr/bin/awk '{sum += $1} END {print sum + 0}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$name" "$jobs" "$pileup_threads" "$window_step" \
        "$((ended - started))" \
        "$max_rss" "$output_bytes" | tee -a "$summary"
done

echo "All pileup cases produced identical non-header VCF records."
