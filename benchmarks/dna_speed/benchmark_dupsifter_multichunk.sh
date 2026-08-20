#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_dupsifter_multichunk.sh REF_FASTA OUTDIR SORT_CHUNK.bam [...]

Environment:
  SORT_THREADS  threads for name and result sorts [allocated CPUs]
  SORT_MEM      memory per sort thread [768M]
  CLEAN_INPUT   delete the prepared name-grouped BAM after validation [1]

The supplied samtools-sort fragments are concatenated, globally name-sorted,
and passed to benchmark_dupsifter_fixed.sh. Result BAMs are deleted after
their SAM checksums, flagstats, and dupsifter statistics are captured.
EOF
    exit 2
}

(( $# >= 3 )) || usage
ref_fasta=$(realpath -s "$1")
outdir=$(realpath -sm "$2")
shift 2
inputs=("$@")
sort_threads=${SORT_THREADS:-${SLURM_CPUS_PER_TASK:-4}}
sort_mem=${SORT_MEM:-768M}
clean_input=${CLEAN_INPUT:-1}

script_dir=$(realpath -s "$(dirname "$0")")
if [[ -n ${SLURM_SUBMIT_DIR:-} && -x $SLURM_SUBMIT_DIR/benchmarks/dna_speed/benchmark_dupsifter_fixed.sh ]]; then
    script_dir=$(realpath -s "$SLURM_SUBMIT_DIR/benchmarks/dna_speed")
fi
repo_dir=$(realpath -s "$script_dir/../..")

[[ $sort_threads =~ ^[1-9][0-9]*$ ]] || usage
[[ $clean_input == 0 || $clean_input == 1 ]] || usage
[[ -r $ref_fasta ]] || usage
for input in "${inputs[@]}"; do
    [[ -r $input ]] || {
        echo "Input is not readable: $input" >&2
        exit 1
    }
done

mkdir -p "$outdir"
name_bam="$outdir/input.name.bam"
partial="$name_bam.partial"
local_root=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/dna-speed-name-sort-${SLURM_JOB_ID:-$$}
mkdir -p "$local_root"
trap 'rm -rf "$local_root"; rm -f "$partial"' EXIT

export SORT_THREADS="$sort_threads"
export SORT_MEM="$sort_mem"
export LOCAL_ROOT="$local_root"
export PARTIAL="$partial"
/usr/bin/time -v -o "$outdir/name_sort.time.txt" \
    bash -o pipefail -c '
        samtools cat "$@" \
            | samtools sort -n -@ "$SORT_THREADS" -m "$SORT_MEM" \
                -T "$LOCAL_ROOT/sort" -O bam -o "$PARTIAL" -
    ' _ "${inputs[@]}" \
    2> "$outdir/name_sort.stderr.txt"
mv "$partial" "$name_bam"
samtools quickcheck -v "$name_bam"

DUPSIFTER_NEW="$repo_dir/tools/dupsifter_speed/dupsifter" \
SORT_THREADS="$sort_threads" SORT_MEM="$sort_mem" CLEAN_BAMS=1 \
    "$script_dir/benchmark_dupsifter_fixed.sh" \
    "$name_bam" "$ref_fasta" "$outdir/results"

if [[ $clean_input == 1 ]]; then
    rm -f "$name_bam"
fi
