#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: prepare_fastq_subset.sh READ_PAIRS R1.fq.gz R2.fq.gz OUT_R1.fq.gz OUT_R2.fq.gz

Write the first READ_PAIRS records from a paired FASTQ input. The input files
must have the same read order, as required by the DNA pipeline.
EOF
    exit 2
}

[[ $# -eq 5 ]] || usage

read_pairs=$1
r1=$2
r2=$3
out_r1=$4
out_r2=$5

[[ $read_pairs =~ ^[1-9][0-9]*$ ]] || usage
[[ -r $r1 && -r $r2 ]] || {
    echo "Input FASTQs are not readable" >&2
    exit 1
}

mkdir -p "$(dirname "$out_r1")" "$(dirname "$out_r2")"
lines=$((read_pairs * 4))

extract_prefix() {
    local input=$1 output=$2 tmp
    tmp="${output}.tmp.$$"
    # head intentionally closes the decompressor early. Without pipefail, the
    # expected SIGPIPE from pigz does not turn a valid prefix into a failure.
    set +o pipefail
    pigz -dc -p 4 "$input" | head -n "$lines" | pigz -c -p 4 > "$tmp"
    local status=$?
    set -o pipefail
    if (( status != 0 )); then
        rm -f "$tmp"
        return "$status"
    fi
    mv "$tmp" "$output"
}

extract_prefix "$r1" "$out_r1"
extract_prefix "$r2" "$out_r2"

r1_lines=$(pigz -dc -p 4 "$out_r1" | wc -l)
r2_lines=$(pigz -dc -p 4 "$out_r2" | wc -l)
expected_lines=$((read_pairs * 4))
if (( r1_lines != expected_lines || r2_lines != expected_lines )); then
    echo "Unexpected subset size: R1=$r1_lines R2=$r2_lines expected=$expected_lines" >&2
    exit 1
fi

r1_first=$(pigz -dc "$out_r1" | sed -n '1s/[[:space:]].*$//p')
r2_first=$(pigz -dc "$out_r2" | sed -n '1s/[[:space:]].*$//p')
r1_last=$(pigz -dc "$out_r1" | tail -n 4 | sed -n '1s/[[:space:]].*$//p')
r2_last=$(pigz -dc "$out_r2" | tail -n 4 | sed -n '1s/[[:space:]].*$//p')
if [[ ${r1_first%/1} != ${r2_first%/2} || ${r1_last%/1} != ${r2_last%/2} ]]; then
    echo "R1/R2 names are not paired at subset boundaries" >&2
    exit 1
fi

printf 'read_pairs\t%s\nR1\t%s\nR2\t%s\n' "$read_pairs" "$out_r1" "$out_r2"
