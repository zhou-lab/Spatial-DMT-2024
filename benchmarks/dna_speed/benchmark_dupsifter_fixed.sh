#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_dupsifter_fixed.sh NAME_BAM REF_FASTA OUTDIR

Environment:
  DUPSIFTER_OLD  Production dupsifter binary.
  DUPSIFTER_NEW  Candidate dupsifter binary.
  SORT_THREADS   samtools sort threads (default: 4).
  SORT_MEM       sort memory per thread (default: 768M).
  CLEAN_BAMS     delete each validated result BAM immediately [1].

The exact same name-grouped BAM is supplied to every case. Therefore SAM
record checksums and dupsifter statistics must match exactly.
EOF
    exit 2
}

[[ $# -eq 3 ]] || usage

name_bam=$(realpath -s "$1")
ref_fasta=$(realpath -s "$2")
outdir=$(realpath -sm "$3")
script_dir=$(realpath -s "$(dirname "$0")")
if [[ -n ${SLURM_SUBMIT_DIR:-} \
        && -x $SLURM_SUBMIT_DIR/benchmarks/dna_speed/benchmark_dupsifter_fixed.sh ]]; then
    script_dir=$(realpath -s "$SLURM_SUBMIT_DIR/benchmarks/dna_speed")
fi
old=${DUPSIFTER_OLD:-/mnt/isilon/zhoulab/labbin/dupsifter}
new=${DUPSIFTER_NEW:-$(realpath -s "$script_dir/../../tools/dupsifter_speed/dupsifter")}
sort_threads=${SORT_THREADS:-4}
sort_mem=${SORT_MEM:-768M}
clean_bams=${CLEAN_BAMS:-1}
local_root=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/dna-speed-dupsifter-${SLURM_JOB_ID:-$$}

[[ $clean_bams == 0 || $clean_bams == 1 ]] || usage

for input in "$name_bam" "$ref_fasta" "$old" "$new"; do
    [[ -r $input ]] || {
        echo "Required input is not readable: $input" >&2
        exit 1
    }
done

mkdir -p "$outdir" "$local_root"
trap 'rm -rf "$local_root"' EXIT

summary="$outdir/summary.tsv"
printf 'case\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kb\toutput_bytes\tview_md5\n' > "$summary"

run_case() {
    local name=$1 binary=$2 stream_arg=$3
    local case_dir="$outdir/$name"
    rm -rf "$case_dir" "$local_root/$name"
    mkdir -p "$case_dir" "$local_root/$name"

    local started ended status bam="$case_dir/output.bam"
    started=$(date +%s)
    set +e
    /usr/bin/time -v -o "$case_dir/time.txt" \
        bash -o pipefail -c '
            "$1" $2 --stats-output "$3" "$4" "$5" \
                | samtools sort -@ "$6" -m "$7" -T "$8/sort" -O bam -o "$9" -
        ' _ "$binary" "$stream_arg" "$case_dir/dupsifter.stat" "$ref_fasta" \
            "$name_bam" "$sort_threads" "$sort_mem" "$local_root/$name" "$bam" \
            2> "$case_dir/pipeline.stderr"
    status=$?
    set -e
    (( status == 0 )) || {
        echo "$name failed with status $status" >&2
        exit "$status"
    }
    ended=$(date +%s)

    samtools quickcheck -v "$bam"
    samtools flagstat -@ "$sort_threads" "$bam" > "$case_dir/flagstat.txt"
    samtools view -@ "$sort_threads" "$bam" | md5sum > "$case_dir/view.md5"

    user=$(/usr/bin/awk -F: '/User time/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$case_dir/time.txt")
    system=$(/usr/bin/awk -F: '/System time/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$case_dir/time.txt")
    max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$case_dir/time.txt")
    view_md5=$(cut -d' ' -f1 "$case_dir/view.md5")
    output_bytes=$(stat -c %s "$bam")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$((ended - started))" \
        "$user" "$system" "$max_rss" "$output_bytes" "$view_md5" \
        | tee -a "$summary"
    if [[ $clean_bams == 1 ]]; then
        rm -f "$bam"
    fi
}

run_case production_sam "$old" ""
run_case fullhash_sam "$new" ""
run_case fullhash_ubam "$new" "--uncompressed-bam"

reference_md5=$(/usr/bin/awk 'NR == 2 {print $7}' "$summary")
if ! /usr/bin/awk -F '\t' -v expected="$reference_md5" \
    'NR > 1 && $7 != expected {exit 1}' "$summary"; then
    echo "Candidate output records differ from production" >&2
    exit 1
fi

if ! cmp -s "$outdir/production_sam/dupsifter.stat" "$outdir/fullhash_sam/dupsifter.stat" \
    || ! cmp -s "$outdir/production_sam/dupsifter.stat" "$outdir/fullhash_ubam/dupsifter.stat"; then
    echo "Candidate dupsifter statistics differ from production" >&2
    exit 1
fi

echo "All fixed-input records and dupsifter statistics match exactly."
