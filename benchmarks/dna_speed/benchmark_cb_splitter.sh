#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: benchmark_cb_splitter.sh INPUT.bam OUTDIR" >&2
    exit 2
}

[[ $# -eq 2 ]] || usage
input=$(realpath -s "$1")
outdir=$(realpath -sm "$2")
script_root=$(realpath -s "$(dirname "$0")/../..")
python_bin=${PYTHON_BIN:-/mnt/isilon/zhoulab/labsoftware/anaconda/conda_2026/bin/python}
python_splitter=${PYTHON_SPLITTER:-$script_root/tools/split_bam_by_cb.py}
hts_splitter=${HTS_SPLITTER:-$script_root/tools/split_bam_by_cb_hts/split_bam_by_cb_hts}
bucketed_splitter=${BUCKETED_SPLITTER:-$script_root/tools/split_bam_by_cb_bucketed.sh}
python_threads=${PYTHON_THREADS:-2}
hts_threads=${HTS_THREADS:-2}
bucket_jobs=${BUCKET_JOBS:-${SLURM_CPUS_PER_TASK:-2}}
bucket_threads=${BUCKET_THREADS:-$((bucket_jobs > 1 ? bucket_jobs - 1 : 1))}

for path in "$input" "$python_bin" "$python_splitter" "$hts_splitter" "$bucketed_splitter"; do
    [[ -r $path ]] || {
        echo "Required input is not readable: $path" >&2
        exit 1
    }
done

mkdir -p "$outdir/python" "$outdir/hts" "$outdir/hts_single" "$outdir/bucketed"

run_case() {
    local name=$1
    shift
    /usr/bin/time -v -o "$outdir/$name.time.txt" "$@" \
        2> "$outdir/$name.stderr.txt"
}

run_case python "$python_bin" "$python_splitter" "$input" "$outdir/python" \
    --threads "$python_threads"
run_case hts "$hts_splitter" -@ "$hts_threads" "$input" "$outdir/hts"
run_case hts_single "$hts_splitter" -@ 0 "$input" "$outdir/hts_single"
run_case bucketed env SPLITTER="$hts_splitter" SPLIT_JOBS="$bucket_jobs" \
    BUCKET_THREADS="$bucket_threads" "$bucketed_splitter" "$input" "$outdir/bucketed"

implementations=(python hts hts_single bucketed)
for name in "${implementations[@]}"; do
    find "$outdir/$name" -maxdepth 1 -type f -name '*.bam' -printf '%f\n' | sort \
        > "$outdir/$name.files.txt"
done
cmp "$outdir/python.files.txt" "$outdir/hts.files.txt"
cmp "$outdir/python.files.txt" "$outdir/hts_single.files.txt"
cmp "$outdir/python.files.txt" "$outdir/bucketed.files.txt"

for name in "${implementations[@]}"; do
    : > "$outdir/$name.records.md5"
done
while IFS= read -r bam; do
    for name in "${implementations[@]}"; do
        samtools view "$outdir/$name/$bam" | md5sum | sed "s|  -|  $bam|" \
            >> "$outdir/$name.records.md5"
    done
done < "$outdir/python.files.txt"
cmp "$outdir/python.records.md5" "$outdir/hts.records.md5"
cmp "$outdir/python.records.md5" "$outdir/hts_single.records.md5"
cmp "$outdir/python.records.md5" "$outdir/bucketed.records.md5"

printf 'implementation\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kb\toutput_bams\n' \
    > "$outdir/summary.tsv"
for name in "${implementations[@]}"; do
    wall=$(/usr/bin/awk '/Elapsed \(wall clock\)/{sub(/^.*\): /, ""); print}' "$outdir/$name.time.txt")
    user=$(/usr/bin/awk -F: '/User time/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$outdir/$name.time.txt")
    system=$(/usr/bin/awk -F: '/System time/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$outdir/$name.time.txt")
    max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/{gsub(/^[[:space:]]+/, "", $2); print $2}' "$outdir/$name.time.txt")
    outputs=$(wc -l < "$outdir/$name.files.txt")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$wall" "$user" "$system" \
        "$max_rss" "$outputs" | tee -a "$outdir/summary.tsv"
done

echo "Per-cell filenames and SAM records match exactly."
