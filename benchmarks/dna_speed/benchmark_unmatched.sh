#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_unmatched.sh R1.fq.gz R2.fq.gz BISCUIT_INDEX REF_FASTA OUTDIR

Environment:
  BENCH_CASES  Space-separated NAME:ALIGN_THREADS:SORT_THREADS:TEMP_MODE cases.
               TEMP_MODE is nfs or local.
  SORT_MEM     samtools sort memory per thread (default: 768M).

The benchmark runs cases serially so they share one node and do not contend
with each other. SLURM_CPUS_PER_TASK must cover ALIGN_THREADS + SORT_THREADS for
balanced cases. A case named baseline may intentionally reproduce the current
pipeline's oversubscription.
EOF
    exit 2
}

[[ $# -eq 5 ]] || usage

# Preserve the final path component for BISCUIT indexes: their FASTA can be a
# symlink while the index sidecars live beside the symlink itself.
r1=$(realpath -s "$1")
r2=$(realpath -s "$2")
biscuit_index=$(realpath -s "$3")
ref_fasta=$(realpath -s "$4")
outdir=$(realpath -sm "$5")

for input in "$r1" "$r2" "$biscuit_index" "$ref_fasta"; do
    [[ -r $input ]] || {
        echo "Required input is not readable: $input" >&2
        exit 1
    }
done

bench_cases=${BENCH_CASES:-\
"baseline_nfs:20:10:nfs balanced_nfs:16:4:nfs balanced_local:16:4:local sort8_local:12:8:local"}
sort_mem=${SORT_MEM:-768M}
allocated_cpus=${SLURM_CPUS_PER_TASK:-$(nproc)}
local_root=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/dna-speed-${SLURM_JOB_ID:-$$}

mkdir -p "$outdir/results" "$outdir/tmp" "$local_root"
trap 'rm -rf "$local_root"' EXIT

summary="$outdir/summary.tsv"
printf 'case\talign_threads\tsort_threads\ttemp_mode\twall_seconds\tmax_rss_kb\toutput_bytes\tview_md5\n' > "$summary"

run_case() {
    local spec=$1 name align_threads sort_threads temp_mode
    IFS=: read -r name align_threads sort_threads temp_mode <<< "$spec"
    [[ -n $name && $align_threads =~ ^[1-9][0-9]*$ && $sort_threads =~ ^[1-9][0-9]*$ ]] || {
        echo "Invalid benchmark case: $spec" >&2
        return 2
    }
    [[ $temp_mode == nfs || $temp_mode == local ]] || {
        echo "Invalid temp mode in case: $spec" >&2
        return 2
    }

    local requested_threads=$((align_threads + sort_threads))
    if (( requested_threads > allocated_cpus )) && [[ $name != baseline* ]]; then
        echo "$name requests $requested_threads runnable threads but only $allocated_cpus were allocated" >&2
        return 2
    fi

    local case_dir="$outdir/results/$name"
    local temp_dir
    if [[ $temp_mode == local ]]; then
        temp_dir="$local_root/$name"
    else
        temp_dir="$outdir/tmp/$name"
    fi
    rm -rf "$case_dir" "$temp_dir"
    mkdir -p "$case_dir" "$temp_dir"

    local bam="$case_dir/unmatched.bam"
    local dup="$case_dir/dupsifter.stat"
    local metrics="$case_dir/time.txt"
    local stderr="$case_dir/pipeline.stderr"
    local started ended wall max_rss output_bytes view_md5
    started=$(date +%s)

    set +e
    {
        /usr/bin/time -v -o "$metrics" \
            bash -o pipefail -c '
                biscuit align "$1" -b 1 -@ "$2" "$3" "$4" \
                | "$5" --stats-output "$6" "$7" - \
                | samtools sort -@ "$8" -m "$9" -T "${10}/sort" -O bam -o "${11}" -
            ' _ "$biscuit_index" "$align_threads" "$r2" "$r1" \
                "${DUPSIFTER:-/mnt/isilon/zhoulab/labbin/dupsifter}" "$dup" \
                "$ref_fasta" "$sort_threads" "$sort_mem" "$temp_dir" "$bam"
    } 2> "$stderr"
    local pipeline_status=$?
    set -e
    if (( pipeline_status != 0 )); then
        echo "$name failed with status $pipeline_status; see $stderr" >&2
        return "$pipeline_status"
    fi
    ended=$(date +%s)
    wall=$((ended - started))

    samtools quickcheck -v "$bam"
    samtools flagstat -@ "$sort_threads" "$bam" > "$case_dir/flagstat.txt"
    samtools stats -@ "$sort_threads" "$bam" > "$case_dir/stats.txt"
    view_md5=$(samtools view -@ "$sort_threads" "$bam" | md5sum | cut -d' ' -f1)
    printf '%s  alignments.sam\n' "$view_md5" > "$case_dir/view.md5"

    max_rss=$(awk -F: '/Maximum resident set size/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$metrics")
    output_bytes=$(stat -c %s "$bam")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$name" "$align_threads" "$sort_threads" "$temp_mode" "$wall" \
        "$max_rss" "$output_bytes" "$view_md5" | tee -a "$summary"
}

for spec in $bench_cases; do
    run_case "$spec"
done

first_spec=${bench_cases%% *}
reference_case=${first_spec%%:*}
reference_primary=$(awk '$1 ~ /^[0-9]+$/ && $4 == "primary" {print $1; exit}' \
    "$outdir/results/$reference_case/flagstat.txt")
for flagstat in "$outdir"/results/*/flagstat.txt; do
    primary=$(awk '$1 ~ /^[0-9]+$/ && $4 == "primary" {print $1; exit}' "$flagstat")
    if [[ $primary != "$reference_primary" ]]; then
        echo "Primary read count differs in $flagstat: expected $reference_primary, observed $primary" >&2
        exit 1
    fi
done

reference_md5=$(awk 'NR == 2 {print $8}' "$summary")
if awk -F '\t' -v expected="$reference_md5" 'NR > 1 && $8 != expected {exit 1}' "$summary"; then
    echo "All benchmark cases produced identical SAM-record checksums."
else
    # BISCUIT's multi-threaded output is not record-order deterministic, and
    # changing its thread count can make tiny differences in secondary hits.
    # Exact equivalence is tested separately against a fixed name-grouped BAM.
    echo "SAM-record checksums differ across BISCUIT thread configurations; aggregate metrics retained for review." >&2
fi
