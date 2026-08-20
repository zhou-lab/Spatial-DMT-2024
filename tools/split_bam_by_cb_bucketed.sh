#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: split_bam_by_cb_bucketed.sh INPUT.bam OUTPUT_DIR

Environment:
  SPLITTER          compiled split_bam_by_cb_hts binary
  SPLIT_JOBS        number of deterministic cell buckets/workers [32]
  BUCKET_THREADS    shared HTSlib threads for the bucket pass [16]
  SPLIT_TMPDIR      primary temporary root; defaults to SLURM_TMPDIR, then /tmp
  SPLIT_FALLBACK_TMPDIR
                    optional second filesystem used when primary space is low
  SPLIT_HEADROOM_GB free space retained on each temporary filesystem [20]
  SPLIT_PRIMARY_BUDGET_BYTES
                    optional cap on usable primary bytes (benchmark/testing)
  KEEP_BUCKETS      set to 1 to retain temporary bucket BAMs [0]
EOF
    exit 2
}

[[ $# -eq 2 ]] || usage
input=$(realpath -s "$1")
output=$(realpath -sm "$2")
script_dir=$(realpath -s "$(dirname "$0")")
splitter=${SPLITTER:-$script_dir/split_bam_by_cb_hts/split_bam_by_cb_hts}
jobs=${SPLIT_JOBS:-32}
bucket_threads=${BUCKET_THREADS:-16}
tmp_parent=${SPLIT_TMPDIR:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}}
fallback_parent=${SPLIT_FALLBACK_TMPDIR:-}
headroom_gb=${SPLIT_HEADROOM_GB:-20}
primary_budget_bytes=${SPLIT_PRIMARY_BUDGET_BYTES:-}
run_id="${SLURM_JOB_ID:-$$}-${SLURM_RESTART_COUNT:-0}"
tmp_dir="$tmp_parent/cb-split-$run_id"
route_dir="$tmp_dir/routes"
local_bucket_dir="$tmp_dir/buckets"
fallback_tmp_dir=""
keep_buckets=${KEEP_BUCKETS:-0}

[[ -r $input && -x $splitter ]] || usage
for value in "$jobs" "$bucket_threads" "$headroom_gb"; do
    [[ $value =~ ^[1-9][0-9]*$ ]] || usage
done
if [[ -n $primary_budget_bytes && ! $primary_budget_bytes =~ ^[0-9]+$ ]]; then
    usage
fi

mkdir -p "$output" "$route_dir" "$local_bucket_dir"
cleanup() {
    if [[ $keep_buckets != 1 ]]; then
        if [[ -n $fallback_tmp_dir ]]; then
            rm -rf "$fallback_tmp_dir"
        fi
        rm -rf "$tmp_dir"
    fi
}
trap cleanup EXIT

input_bytes=$(stat -c %s "$input")
headroom_bytes=$((headroom_gb * 1024 * 1024 * 1024))
estimated_bucket_bytes=$((input_bytes + input_bytes / 5))
output_device=$(stat -c %d "$output")
primary_device=$(stat -c %d "$tmp_parent")
primary_available=$(df --output=avail -B1 "$tmp_parent" | awk 'NR == 2 {print $1}')
primary_reserved=$headroom_bytes
if [[ $primary_device == "$output_device" ]]; then
    primary_reserved=$((primary_reserved + estimated_bucket_bytes))
fi
primary_usable=$((primary_available > primary_reserved ? primary_available - primary_reserved : 0))
if [[ -n $primary_budget_bytes ]] && (( primary_usable > primary_budget_bytes )); then
    primary_usable=$primary_budget_bytes
fi
local_jobs=$((primary_usable * jobs / estimated_bucket_bytes))
(( local_jobs > jobs )) && local_jobs=$jobs

bucket_dirs=("$local_bucket_dir")
if (( local_jobs < jobs )); then
    [[ -n $fallback_parent ]] || {
        echo "Primary scratch fits $local_jobs/$jobs buckets; SPLIT_FALLBACK_TMPDIR is required" >&2
        exit 1
    }
    mkdir -p "$fallback_parent"
    fallback_device=$(stat -c %d "$fallback_parent")
    if [[ $primary_device == "$fallback_device" ]]; then
        echo "Primary and fallback temporary roots are on the same filesystem" >&2
        exit 1
    fi
    fallback_tmp_dir="$fallback_parent/cb-split-$run_id"
    fallback_bucket_dir="$fallback_tmp_dir/buckets"
    mkdir -p "$fallback_bucket_dir"
    bucket_dirs+=("$fallback_bucket_dir")

    fallback_jobs=$((jobs - local_jobs))
    fallback_bucket_bytes=$((
        (estimated_bucket_bytes * fallback_jobs + jobs - 1) / jobs
    ))
    fallback_required=$((fallback_bucket_bytes + headroom_bytes))
    if [[ $fallback_device == "$output_device" ]]; then
        fallback_required=$((fallback_required + estimated_bucket_bytes))
    fi
    fallback_available=$(df --output=avail -B1 "$fallback_parent" | awk 'NR == 2 {print $1}')
    if (( fallback_available < fallback_required )); then
        echo "Fallback/output filesystem needs about $fallback_required bytes; $fallback_available available" >&2
        exit 1
    fi
fi

if [[ $output_device != "$primary_device" \
        && ( -z ${fallback_device:-} || $output_device != "${fallback_device:-}" ) ]]; then
    output_available=$(df --output=avail -B1 "$output" | awk 'NR == 2 {print $1}')
    output_required=$((estimated_bucket_bytes + headroom_bytes))
    if (( output_available < output_required )); then
        echo "Output filesystem needs about $output_required bytes; $output_available available" >&2
        exit 1
    fi
fi

for ((bucket = 0; bucket < jobs; bucket++)); do
    printf -v label 'bucket_%05d.bam' "$bucket"
    if (( bucket < local_jobs )); then
        destination="$local_bucket_dir/$label"
    else
        destination="$fallback_bucket_dir/$label"
    fi
    ln -s "$destination" "$route_dir/$label"
done

printf 'bucket_plan input_bytes=%s jobs=%s local_jobs=%s fallback_jobs=%s\n' \
    "$input_bytes" "$jobs" "$local_jobs" "$((jobs - local_jobs))" >&2

"$splitter" -@ "$bucket_threads" --buckets "$jobs" "$input" "$route_dir"

split_bucket() {
    "$splitter" -@ 0 "$1" "$output"
    if [[ $keep_buckets != 1 ]]; then
        rm -f "$1"
    fi
}
export -f split_bucket
export splitter output keep_buckets

find "${bucket_dirs[@]}" -maxdepth 1 -type f -name 'bucket_*.bam' -print0 \
    | sort -z \
    | parallel -0 --halt soon,fail=1 -j "$jobs" \
        split_bucket {}
