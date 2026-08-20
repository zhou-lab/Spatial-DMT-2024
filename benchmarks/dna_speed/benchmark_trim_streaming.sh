#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_trim_streaming.sh R1.fq[.gz] R2.fq[.gz] WHITELIST OUTDIR

Environment:
  WORKERS       Parallel trim processes (default: allocated CPUs).
  BLOCK_SIZE    GNU Parallel input block size (default: 512M).
  PYTHON_BIN    Pipeline Python executable.
  TRIM_SCRIPT   spatialmeth_trimtag.py path.
  COMPRESSION_LEVEL
                gzip compression level [6].
  SKIP_BARCODE_COUNTS
                suppress optional observed-barcode tables [1].

R1 and R2 are decompressed concurrently, collapsed to one tab-separated line
per FASTQ record, paired into eight fields, and streamed directly to workers.
No decompressed FASTQ chunks are written to disk.
EOF
    exit 2
}

[[ $# -eq 4 ]] || usage
r1=$(realpath -s "$1")
r2=$(realpath -s "$2")
whitelist=$(realpath -s "$3")
outdir=$(realpath -sm "$4")
workers=${WORKERS:-${SLURM_CPUS_PER_TASK:-1}}
block_size=${BLOCK_SIZE:-512M}
python_bin=${PYTHON_BIN:-/mnt/isilon/zhoulab/labsoftware/anaconda/conda_2026/bin/python}
trim_script=${TRIM_SCRIPT:-$(realpath -s "$(dirname "$0")/../../tools/spatialmeth_trimtag.py")}
compression_level=${COMPRESSION_LEVEL:-6}
skip_barcode_counts=${SKIP_BARCODE_COUNTS:-1}

for input in "$r1" "$r2" "$whitelist" "$python_bin" "$trim_script"; do
    [[ -r $input ]] || {
        echo "Required input is not readable: $input" >&2
        exit 1
    }
done
[[ $workers =~ ^[1-9][0-9]*$ ]] || usage
[[ $compression_level =~ ^[0-9]$ ]] || usage
[[ $skip_barcode_counts == 0 || $skip_barcode_counts == 1 ]] || usage

output_dir="$outdir/output"
mkdir -p "$output_dir"
if find "$output_dir" -mindepth 1 -print -quit | grep -q .; then
    echo "Output directory is not empty: $output_dir" >&2
    exit 1
fi

decompress() {
    if [[ $1 == *.gz ]]; then
        pigz -dc -p 2 "$1"
    else
        cat "$1"
    fi
}

process_block() {
    local block_id=$1
    local barcode_args=()
    if [[ $skip_barcode_counts == 1 ]]; then
        barcode_args+=(--skip-barcode-counts)
    fi
    "$python_bin" "$trim_script" - - --paired-tab-stdin \
        -o "$output_dir/out_$block_id" \
        -a CTATCTCTTATA AGATGCGAGAAGCCAACGCTTG \
           AATCATACACCAATACAAAACATCAACCAC \
           CAAATACTCTAACCTCTCAAACACATAAAT \
        --protocol spatialdmt \
        --whitelist "$whitelist" \
        --whitelist_col 4 \
        --bc-mismatch 1 \
        --matcher exact-seeded \
        "${barcode_args[@]}" \
        --compression-level "$compression_level" \
        --compression-backend gzip
}

export -f decompress process_block
export python_bin trim_script output_dir whitelist compression_level \
    skip_barcode_counts PYTHONHASHSEED=0

started=$(date +%s)
/usr/bin/time -v -o "$outdir/time.txt" \
    bash -o pipefail -c '
        paste <(decompress "$1" | paste - - - -) \
              <(decompress "$2" | paste - - - -) \
        | parallel --plain --halt soon,fail=1 --pipe \
              --block "$3" --recend "\n" -j "$4" process_block {#}
    ' _ "$r1" "$r2" "$block_size" "$workers" \
    2> "$outdir/process.stderr"
ended=$(date +%s)

stats_files=("$output_dir"/out_*_stats.txt)
(( ${#stats_files[@]} > 0 )) || {
    echo "Streaming trim produced no stats files" >&2
    exit 1
}
/usr/bin/awk '{count[$1] += $2} END {for (key in count) print key "\t" count[key]}' \
    "${stats_files[@]}" | LC_ALL=C sort -k1,1 > "$outdir/aggregate_stats.tsv"

: > "$outdir/fastq_checksums.tsv"
for suffix in R1 R2 Unmatched_R1 Unmatched_R2 Lambda_R1 Lambda_R2; do
    mapfile -t files < <(
        find "$output_dir" -maxdepth 1 -type f -regextype posix-extended \
            -regex ".*/out_[0-9]+_${suffix}\\.fq\\.gz" -printf '%p\n' | sort -V
    )
    (( ${#files[@]} == ${#stats_files[@]} )) || {
        echo "$suffix has ${#files[@]} shards; expected ${#stats_files[@]}" >&2
        exit 1
    }
    checksum=$(pigz -dc -p 4 "${files[@]}" | md5sum | cut -d' ' -f1)
    printf '%s\t%s\n' "$suffix" "$checksum" >> "$outdir/fastq_checksums.tsv"
done

max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$outdir/time.txt")
stats_sha=$(sha256sum "$outdir/aggregate_stats.tsv" | cut -d' ' -f1)
fastq_sha=$(sha256sum "$outdir/fastq_checksums.tsv" | cut -d' ' -f1)
compressed_bytes=$(find "$output_dir" -type f -name '*.fq.gz' -printf '%s\n' \
    | /usr/bin/awk '{total += $1} END {print total + 0}')
printf 'workers\tblock_size\tblocks\twall_seconds\tmax_rss_kb\tcompressed_bytes\tstats_sha256\tfastq_sha256\n' \
    > "$outdir/summary.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$workers" "$block_size" "${#stats_files[@]}" "$((ended - started))" \
    "$max_rss" "$compressed_bytes" "$stats_sha" "$fastq_sha" \
    | tee -a "$outdir/summary.tsv"
