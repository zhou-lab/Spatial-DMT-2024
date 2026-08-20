#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_chunked_unmatched.sh R1.fq.gz R2.fq.gz BISCUIT_INDEX REF_FASTA OUTDIR [REFERENCE_BAM]

Environment:
  SHARDS                    Number of paired FASTQ shards (default: 4)
  ALIGN_THREADS_PER_SHARD   BISCUIT threads per shard (default: 4)
  VIEW_THREADS_PER_SHARD    samtools view threads per shard (default: 1)
  SORT_THREADS              final samtools sort threads (default: 4)
  SORT_MEM                  final sort memory per thread (default: 768M)
  FINAL_TEMP_MODE           nfs or local (default: local)

Shard BAMs are name-grouped and concatenated in original input order. One
global dupsifter pass therefore retains duplicate detection across shards.
EOF
    exit 2
}

(( $# == 5 || $# == 6 )) || usage

r1=$(realpath -s "$1")
r2=$(realpath -s "$2")
biscuit_index=$(realpath -s "$3")
ref_fasta=$(realpath -s "$4")
outdir=$(realpath -sm "$5")
reference_bam=${6:-}
[[ -z $reference_bam ]] || reference_bam=$(realpath -s "$reference_bam")

for input in "$r1" "$r2" "$biscuit_index" "$ref_fasta"; do
    [[ -r $input ]] || {
        echo "Required input is not readable: $input" >&2
        exit 1
    }
done

shards=${SHARDS:-4}
align_threads=${ALIGN_THREADS_PER_SHARD:-4}
view_threads=${VIEW_THREADS_PER_SHARD:-1}
sort_threads=${SORT_THREADS:-4}
sort_mem=${SORT_MEM:-768M}
temp_mode=${FINAL_TEMP_MODE:-local}
allocated_cpus=${SLURM_CPUS_PER_TASK:-$(nproc)}
local_root=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/dna-speed-chunked-${SLURM_JOB_ID:-$$}

for value in "$shards" "$align_threads" "$view_threads" "$sort_threads"; do
    [[ $value =~ ^[1-9][0-9]*$ ]] || usage
done
[[ $temp_mode == nfs || $temp_mode == local ]] || usage

threads_per_worker=$((align_threads + view_threads))
parallel_workers=$((allocated_cpus / threads_per_worker))
(( parallel_workers > shards )) && parallel_workers=$shards
(( parallel_workers > 0 )) || {
    echo "Allocation cannot run one shard worker" >&2
    exit 2
}

mkdir -p "$outdir/chunks" "$outdir/shard_bams" "$outdir/results" "$outdir/tmp" "$local_root"
trap 'rm -rf "$local_root"' EXIT

# Benchmark-only fixed splitting. Production chunking reuses compressed outputs
# from dna_trim_and_tag and does not materialize uncompressed FASTQ shards.
if ! compgen -G "$outdir/chunks/r1.*" >/dev/null; then
    r1_lines=$(pigz -dc -p 4 "$r1" | wc -l)
    r2_lines=$(pigz -dc -p 4 "$r2" | wc -l)
    if (( r1_lines == 0 || r1_lines != r2_lines || r1_lines % 4 != 0 )); then
        echo "Invalid paired FASTQ line counts: R1=$r1_lines R2=$r2_lines" >&2
        exit 1
    fi
    pairs=$((r1_lines / 4))
    pairs_per_shard=$(((pairs + shards - 1) / shards))
    lines_per_shard=$((pairs_per_shard * 4))
    pigz -dc -p 4 "$r1" | split -d -a 4 -l "$lines_per_shard" - "$outdir/chunks/r1."
    pigz -dc -p 4 "$r2" | split -d -a 4 -l "$lines_per_shard" - "$outdir/chunks/r2."
fi

r1_chunks=$(find "$outdir/chunks" -maxdepth 1 -type f -name 'r1.*' | wc -l)
r2_chunks=$(find "$outdir/chunks" -maxdepth 1 -type f -name 'r2.*' | wc -l)
if (( r1_chunks != shards || r2_chunks != shards )); then
    echo "Expected $shards paired chunks, found R1=$r1_chunks R2=$r2_chunks" >&2
    exit 1
fi

align_shard() {
    set -o pipefail
    local r1_chunk=$1 cid r2_chunk bam partial
    cid=${r1_chunk##*.}
    r2_chunk="$outdir/chunks/r2.$cid"
    bam="$outdir/shard_bams/$cid.name.bam"
    partial="$bam.partial.$$"
    if [[ -s $bam ]] && samtools quickcheck "$bam"; then
        return 0
    fi
    biscuit align "$biscuit_index" -b 1 -@ "$align_threads" "$r2_chunk" "$r1_chunk" \
        2> "$outdir/shard_bams/$cid.biscuit.stderr" \
        | samtools view -@ "$view_threads" -b -o "$partial" -
    samtools quickcheck -v "$partial"
    mv "$partial" "$bam"
}
export -f align_shard
export outdir biscuit_index align_threads view_threads

metrics="$outdir/results/time.txt"
pipeline_stderr="$outdir/results/pipeline.stderr"
started=$(date +%s)
set +e
{
    /usr/bin/time -v -o "$metrics" bash -o pipefail -c '
        parallel --halt soon,fail=1 -j "$1" align_shard ::: "$2"/chunks/r1.*
        find "$2/shard_bams" -maxdepth 1 -type f -name "*.name.bam" -print \
            | sort -V > "$2/results/shards.fofn"
        samtools cat -b "$2/results/shards.fofn" \
            | "$3" --stats-output "$2/results/dupsifter.stat" "$4" - \
            | samtools sort -@ "$5" -m "$6" -T "$7/sort" -O bam \
                -o "$2/results/unmatched.bam" -
    ' _ "$parallel_workers" "$outdir" \
        "${DUPSIFTER:-/mnt/isilon/zhoulab/labbin/dupsifter}" "$ref_fasta" \
        "$sort_threads" "$sort_mem" \
        "$([[ $temp_mode == local ]] && printf %s "$local_root" || printf %s "$outdir/tmp")"
} 2> "$pipeline_stderr"
pipeline_status=$?
set -e
ended=$(date +%s)
if (( pipeline_status != 0 )); then
    echo "Chunked benchmark failed with status $pipeline_status; see $pipeline_stderr" >&2
    exit "$pipeline_status"
fi

bam="$outdir/results/unmatched.bam"
samtools quickcheck -v "$bam"
samtools flagstat -@ "$sort_threads" "$bam" > "$outdir/results/flagstat.txt"
samtools stats -@ "$sort_threads" "$bam" > "$outdir/results/stats.txt"
view_md5=$(samtools view -@ "$sort_threads" "$bam" | md5sum | cut -d' ' -f1)
printf '%s  alignments.sam\n' "$view_md5" > "$outdir/results/view.md5"

max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$metrics")
printf 'shards\talign_threads\tview_threads\tparallel_workers\tsort_threads\ttemp_mode\twall_seconds\tmax_rss_kb\toutput_bytes\tview_md5\n' \
    > "$outdir/summary.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$shards" "$align_threads" "$view_threads" "$parallel_workers" \
    "$sort_threads" "$temp_mode" "$((ended - started))" "$max_rss" \
    "$(stat -c %s "$bam")" "$view_md5" | tee -a "$outdir/summary.tsv"

if [[ -n $reference_bam ]]; then
    [[ -r $reference_bam ]] || {
        echo "Reference BAM is not readable: $reference_bam" >&2
        exit 1
    }
    reference_md5=$(samtools view -@ "$sort_threads" "$reference_bam" | md5sum | cut -d' ' -f1)
    if [[ $view_md5 != "$reference_md5" ]]; then
        echo "Chunked and monolithic SAM-record checksums differ; comparing aggregate metrics." >&2
    else
        echo "Chunked and monolithic SAM-record checksums match."
    fi

    reference_flagstat="$outdir/results/reference.flagstat.txt"
    samtools flagstat -@ "$sort_threads" "$reference_bam" > "$reference_flagstat"
    primary=$(awk '$1 ~ /^[0-9]+$/ && $4 == "primary" {print $1; exit}' \
        "$outdir/results/flagstat.txt")
    reference_primary=$(awk '$1 ~ /^[0-9]+$/ && $4 == "primary" {print $1; exit}' \
        "$reference_flagstat")
    if [[ $primary != "$reference_primary" ]]; then
        echo "Chunked primary count differs: expected $reference_primary, observed $primary" >&2
        exit 1
    fi
fi
