#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_dupsifter_boundary_fixture.sh NAME_BAM REF_FASTA OUTDIR

Environment:
  DUPSIFTER_NEW  optimized dupsifter binary
  SORT_THREADS   samtools sort threads [allocated CPUs]
  SORT_MEM       sort memory per thread [768M]
  CLEAN_INPUT    remove NAME_BAM after exact validation [1]

OUTDIR/production_sam must contain the retained BAM, time file, and stderr
from a production run that reached the fixture's expected final singleton.
The optimized cases must stop at that same read name, reject only that final
record, and produce identical complete prefixes through SAM and uncompressed
BAM transport. The production binary exits before flushing buffered output,
so its record MD5 is reported for timing context but is not an equivalence
reference for this malformed-input fixture.
EOF
    exit 2
}

[[ $# -eq 3 ]] || usage
name_bam=$(realpath -s "$1")
ref_fasta=$(realpath -s "$2")
outdir=$(realpath -sm "$3")
script_dir=$(realpath -s "$(dirname "$0")")
if [[ -n ${SLURM_SUBMIT_DIR:-} \
        && -x $SLURM_SUBMIT_DIR/tools/dupsifter_speed/dupsifter ]]; then
    script_dir=$(realpath -s "$SLURM_SUBMIT_DIR/benchmarks/dna_speed")
fi
new=${DUPSIFTER_NEW:-$(realpath -s "$script_dir/../../tools/dupsifter_speed/dupsifter")}
sort_threads=${SORT_THREADS:-${SLURM_CPUS_PER_TASK:-4}}
sort_mem=${SORT_MEM:-768M}
clean_input=${CLEAN_INPUT:-1}
local_root=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/dna-speed-boundary-${SLURM_JOB_ID:-$$}

[[ -r $name_bam && -r $ref_fasta && -x $new ]] || usage
[[ $sort_threads =~ ^[1-9][0-9]*$ ]] || usage
[[ $clean_input == 0 || $clean_input == 1 ]] || usage
mkdir -p "$outdir" "$local_root"
trap 'rm -rf "$local_root"' EXIT

elapsed_seconds() {
    awk '/Elapsed \(wall clock\)/ {
        sub(/^.*\): /, "")
        count = split($0, part, ":")
        if (count == 3) seconds = part[1] * 3600 + part[2] * 60 + part[3]
        else seconds = part[1] * 60 + part[2]
        printf "%.0f\n", seconds
    }' "$1"
}

summary="$outdir/summary.tsv"
printf 'case\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kb\toutput_bytes\tview_md5\ttermination\n' \
    > "$summary"

validate_expected_boundary() {
    local case_dir=$1
    grep -F "Can't find read 1 and/or read 2 in 1 reads with read ID:" \
        "$case_dir/pipeline.stderr" > "$case_dir/singleton_error.txt"
    samtools quickcheck -v "$case_dir/output.bam"
}

summarize_case() {
    local name=$1 case_dir="$outdir/$1" bam="$outdir/$1/output.bam"
    validate_expected_boundary "$case_dir"
    samtools flagstat -@ "$sort_threads" "$bam" > "$case_dir/flagstat.txt"
    samtools view -@ "$sort_threads" "$bam" | md5sum > "$case_dir/view.md5"

    local wall user system max_rss output_bytes view_md5
    wall=$(elapsed_seconds "$case_dir/time.txt")
    user=$(awk -F: '/User time/{sub(/^[[:space:]]+/, "", $2); print $2}' "$case_dir/time.txt")
    system=$(awk -F: '/System time/{sub(/^[[:space:]]+/, "", $2); print $2}' "$case_dir/time.txt")
    max_rss=$(awk -F: '/Maximum resident set size/{sub(/^[[:space:]]+/, "", $2); print $2}' "$case_dir/time.txt")
    output_bytes=$(stat -c %s "$bam")
    view_md5=$(cut -d' ' -f1 "$case_dir/view.md5")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\texpected_final_singleton\n' \
        "$name" "$wall" "$user" "$system" "$max_rss" "$output_bytes" \
        "$view_md5" | tee -a "$summary"
}

run_candidate() {
    local name=$1 stream_arg=$2 case_dir="$outdir/$1"
    find "$case_dir" -mindepth 1 -delete 2>/dev/null || true
    mkdir -p "$case_dir" "$local_root/$name"

    set +e
    /usr/bin/time -v -o "$case_dir/time.txt" \
        bash -o pipefail -c '
            "$1" $2 --stats-output "$3" "$4" "$5" \
                | samtools sort -@ "$6" -m "$7" -T "$8/sort" \
                    -O bam -o "$9" -
        ' _ "$new" "$stream_arg" "$case_dir/dupsifter.stat" "$ref_fasta" \
            "$name_bam" "$sort_threads" "$sort_mem" "$local_root/$name" \
            "$case_dir/output.bam" 2> "$case_dir/pipeline.stderr"
    status=$?
    set -e
    if (( status == 0 )); then
        echo "$name unexpectedly accepted the incomplete final template" >&2
        exit 1
    fi
    summarize_case "$name"
}

production_dir="$outdir/production_sam"
[[ -s $production_dir/output.bam \
    && -s $production_dir/time.txt \
    && -s $production_dir/pipeline.stderr ]] || usage
summarize_case production_sam
run_candidate fullhash_sam ""
run_candidate fullhash_ubam "--uncompressed-bam"

cmp "$outdir/production_sam/singleton_error.txt" \
    "$outdir/fullhash_sam/singleton_error.txt"
cmp "$outdir/production_sam/singleton_error.txt" \
    "$outdir/fullhash_ubam/singleton_error.txt"
cmp "$outdir/fullhash_sam/flagstat.txt" \
    "$outdir/fullhash_ubam/flagstat.txt"
cmp "$outdir/fullhash_sam/view.md5" \
    "$outdir/fullhash_ubam/view.md5"

input_records=$(samtools view -c "$name_bam")
output_records=$(awk 'NR == 1 {print $1}' \
    "$outdir/fullhash_ubam/flagstat.txt")
(( output_records == input_records - 1 )) || {
    echo "Expected $((input_records - 1)) completed records; found $output_records" >&2
    exit 1
}

find "$outdir" -mindepth 2 -maxdepth 2 -type f -name 'output.bam' -delete
if [[ $clean_input == 1 ]]; then
    find "$(dirname "$name_bam")" -maxdepth 1 -type f \
        -samefile "$name_bam" -delete
fi
cat "$summary"
echo "Graceful SAM/BAM prefixes and final-singleton errors match exactly."
