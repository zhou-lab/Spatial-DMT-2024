#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_trim_compression.sh R1.fq.gz R2.fq.gz WHITELIST OUTDIR

Environment:
  TRIM_CASES   Space-separated NAME:BACKEND:THREADS_PER_STREAM:WORKERS[:MATCHER]
               cases. MATCHER defaults to fuzzy.
  CHUNK_READS  Read pairs per benchmark chunk (default: 250000).
  PYTHON_BIN   Python executable used by the pipeline.
  TRIM_SCRIPT  Trim script path (needed when Slurm copies this benchmark).

Input splitting is performed once outside the timed cases. Each case processes
the same paired chunks, and decompressed FASTQ plus aggregate-stat checksums are
compared with the baseline case.
EOF
    exit 2
}

[[ $# -eq 4 ]] || usage

r1=$(realpath -s "$1")
r2=$(realpath -s "$2")
whitelist=$(realpath -s "$3")
outdir=$(realpath -sm "$4")

for input in "$r1" "$r2" "$whitelist"; do
    [[ -r $input ]] || {
        echo "Required input is not readable: $input" >&2
        exit 1
    }
done

trim_cases=${TRIM_CASES:-\
"baseline:pigz:2:20 pigz1:pigz:1:20 gzip_sync:gzip:1:20 pigz1_8workers:pigz:1:8"}
chunk_reads=${CHUNK_READS:-250000}
python_bin=${PYTHON_BIN:-/mnt/isilon/zhoulab/labsoftware/anaconda/conda_2026/bin/python}
trim_script=${TRIM_SCRIPT:-$(realpath -s "$(dirname "$0")/../../tools/spatialmeth_trimtag.py")}
allocated_cpus=${SLURM_CPUS_PER_TASK:-$(nproc)}

[[ $chunk_reads =~ ^[1-9][0-9]*$ ]] || usage
mkdir -p "$outdir/chunks" "$outdir/results"

chunk_lines=$((chunk_reads * 4))
if ! compgen -G "$outdir/chunks/r1.*" >/dev/null; then
    rm -f "$outdir/chunks/r1."* "$outdir/chunks/r2."*
    pigz -dc -p 4 "$r1" | split -d -a 4 -l "$chunk_lines" - "$outdir/chunks/r1."
    pigz -dc -p 4 "$r2" | split -d -a 4 -l "$chunk_lines" - "$outdir/chunks/r2."
fi

r1_chunks=$(find "$outdir/chunks" -maxdepth 1 -type f -name 'r1.*' | wc -l)
r2_chunks=$(find "$outdir/chunks" -maxdepth 1 -type f -name 'r2.*' | wc -l)
if (( r1_chunks == 0 || r1_chunks != r2_chunks )); then
    echo "Mismatched benchmark chunks: R1=$r1_chunks R2=$r2_chunks" >&2
    exit 1
fi

summary="$outdir/summary.tsv"
printf 'case\tbackend\tthreads_per_stream\tworkers\tmatcher\twall_seconds\tmax_rss_kb\tcompressed_bytes\tstats_sha256\tfastq_sha256\n' > "$summary"

run_case() {
    local spec=$1 name backend threads workers matcher
    IFS=: read -r name backend threads workers matcher <<< "$spec"
    matcher=${matcher:-fuzzy}
    [[ -n $name && $threads =~ ^[1-9][0-9]*$ && $workers =~ ^[1-9][0-9]*$ ]] || {
        echo "Invalid trim case: $spec" >&2
        return 2
    }
    case "$backend" in
        pigz|gzip) ;;
        *) echo "Invalid compression backend in case: $spec" >&2; return 2 ;;
    esac
    case "$matcher" in
        fuzzy|deterministic|seeded|exact-seeded) ;;
        *) echo "Invalid matcher in case: $spec" >&2; return 2 ;;
    esac
    if (( workers > allocated_cpus )); then
        echo "$name requests $workers workers but only $allocated_cpus CPUs were allocated" >&2
        return 2
    fi

    local case_dir="$outdir/results/$name"
    local output_dir="$case_dir/output"
    local metrics="$case_dir/time.txt"
    local stderr="$case_dir/process.stderr"
    rm -rf "$case_dir"
    mkdir -p "$output_dir"

    process_chunk() {
        local f=$1 cid f2
        cid=${f##*.}
        f2="$outdir/chunks/r2.$cid"
        "$python_bin" "$trim_script" "$f" "$f2" \
            -o "$output_dir/out_$cid" \
            -a CTATCTCTTATA AGATGCGAGAAGCCAACGCTTG \
               AATCATACACCAATACAAAACATCAACCAC \
               CAAATACTCTAACCTCTCAAACACATAAAT \
            --protocol spatialdmt \
            --whitelist "$whitelist" \
            --whitelist_col 4 \
            --bc-mismatch 1 \
            --matcher "$matcher" \
            --compression-backend "$backend" \
            --compression-threads "$threads"
    }
    export -f process_chunk
    export outdir output_dir python_bin trim_script whitelist backend threads matcher

    local started ended wall
    started=$(date +%s)
    set +e
    /usr/bin/time -v -o "$metrics" \
        parallel --halt soon,fail=1 -j "$workers" process_chunk ::: "$outdir"/chunks/r1.* \
        2> "$stderr"
    local process_status=$?
    set -e
    ended=$(date +%s)
    wall=$((ended - started))
    if (( process_status != 0 )); then
        echo "$name failed with status $process_status; see $stderr" >&2
        return "$process_status"
    fi

    local aggregate_stats="$case_dir/aggregate_stats.tsv"
    /usr/bin/awk '{count[$1] += $2} END {for (key in count) print key "\t" count[key]}' \
        "$output_dir"/out_*_stats.txt | LC_ALL=C sort -k1,1 > "$aggregate_stats"

    local hashes="$case_dir/fastq_checksums.tsv" suffix checksum
    : > "$hashes"
    for suffix in R1 R2 Unmatched_R1 Unmatched_R2 Lambda_R1 Lambda_R2; do
        files=()
        for f in "$outdir"/chunks/r1.*; do
            cid=${f##*.}
            files+=("$output_dir/out_${cid}_${suffix}.fq.gz")
        done
        (( ${#files[@]} == r1_chunks )) || {
            echo "$name has ${#files[@]} $suffix outputs; expected $r1_chunks" >&2
            return 1
        }
        for f in "${files[@]}"; do
            [[ -f $f ]] || {
                echo "$name is missing expected output: $f" >&2
                return 1
            }
        done
        checksum=$(pigz -dc -p 4 "${files[@]}" | md5sum | cut -d' ' -f1)
        printf '%s\t%s\n' "$suffix" "$checksum" >> "$hashes"
    done

    local max_rss compressed_bytes stats_sha fastq_sha
    max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$metrics")
    compressed_bytes=$(find "$output_dir" -type f -name '*.fq.gz' -printf '%s\n' \
        | /usr/bin/awk '{total += $1} END {print total + 0}')
    stats_sha=$(sha256sum "$aggregate_stats" | cut -d' ' -f1)
    fastq_sha=$(sha256sum "$hashes" | cut -d' ' -f1)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$name" "$backend" "$threads" "$workers" "$matcher" "$wall" "$max_rss" \
        "$compressed_bytes" "$stats_sha" "$fastq_sha" | tee -a "$summary"

    rm -rf "$output_dir"
}

for spec in $trim_cases; do
    run_case "$spec"
done

baseline_stats=$(/usr/bin/awk -F '\t' 'NR == 2 {print $9}' "$summary")
baseline_fastq=$(/usr/bin/awk -F '\t' 'NR == 2 {print $10}' "$summary")
if /usr/bin/awk -F '\t' -v stats="$baseline_stats" -v fastq="$baseline_fastq" \
        'NR > 1 && ($9 != stats || $10 != fastq) {exit 1}' "$summary"; then
    echo "All trim cases produced identical decompressed FASTQs and aggregate statistics."
else
    echo "Trim cases differ in output; inspect before adopting a change." >&2
    exit 1
fi
