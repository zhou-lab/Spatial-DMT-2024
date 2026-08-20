#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: benchmark_pileup_transport.sh VCF_LIST OUTDIR

VCF_LIST contains one uncompressed VCF path per nonempty line. The benchmark
compares one-pass streaming sort, the production two-pass header+record sort,
and direct bgzip/tabix transport. All successful cases must have identical
decompressed records.

Environment:
  TRANSPORT_CASES  Space-separated modes [stream_sorted sorted].
                   Add direct only for diagnostic ordering tests.
EOF
    exit 2
}

[[ $# -eq 2 ]] || usage

vcf_list=$(realpath -s "$1")
outdir=$(realpath -sm "$2")
jobs=${SLURM_CPUS_PER_TASK:-$(nproc)}
transports=${TRANSPORT_CASES:-"stream_sorted sorted"}
[[ -r $vcf_list && $jobs =~ ^[1-9][0-9]*$ ]] || usage

while read -r vcf; do
    [[ -z $vcf || $vcf == \#* ]] && continue
    [[ -s $vcf ]] || {
        echo "VCF is not readable or empty: $vcf" >&2
        exit 1
    }
done < "$vcf_list"

mkdir -p "$outdir/results"
summary="$outdir/summary.tsv"
printf 'case\trecord_order\twall_seconds\tmax_rss_kb\toutput_bytes\n' > "$summary"

compress_one() {
    local vcf=$1 label output
    label=$(basename "$vcf" .vcf)
    output="$case_dir/$label.vcf.gz"
    if [[ $transport == direct ]]; then
        bgzip -c "$vcf" > "$output"
    elif [[ $transport == stream_sorted ]]; then
        {
            while IFS= read -r line; do
                if [[ $line == \#* ]]; then
                    printf '%s\n' "$line"
                else
                    {
                        printf '%s\n' "$line"
                        cat
                    } | LC_ALL=C sort -k1,1V -k2,2n
                    break
                fi
            done
        } < "$vcf" | bgzip -c > "$output"
    else
        {
            /usr/bin/awk '/^#/ {print}' "$vcf"
            /usr/bin/awk '!/^#/ {print}' "$vcf" \
                | LC_ALL=C sort -k1,1V -k2,2n
        } | bgzip -c > "$output"
    fi
    tabix -p vcf "$output"
}
export -f compress_one

first_case=""
for transport in $transports; do
    [[ $transport == stream_sorted || $transport == sorted \
        || $transport == direct ]] || usage
    case_dir="$outdir/results/$transport"
    rm -rf "$case_dir"
    mkdir -p "$case_dir"
    export case_dir transport

    started=$(date +%s)
    /usr/bin/time -v -o "$case_dir/time.txt" \
        parallel --plain --halt soon,fail=1 -j "$jobs" \
            compress_one :::: "$vcf_list" \
        2> "$case_dir/pipeline.stderr"
    ended=$(date +%s)

    while read -r vcf; do
        [[ -z $vcf || $vcf == \#* ]] && continue
        label=$(basename "$vcf" .vcf)
        bgzip -dc "$case_dir/$label.vcf.gz" \
            | md5sum \
            | /usr/bin/awk -v label="$label" '{print label "\t" $1}'
    done < "$vcf_list" > "$case_dir/checksums.tsv"

    if [[ -z $first_case ]]; then
        first_case=$transport
        record_order=reference
    elif cmp "$outdir/results/$first_case/checksums.tsv" \
            "$case_dir/checksums.tsv"; then
        record_order=exact
    elif [[ $transport == direct ]]; then
        record_order=different
        echo "Direct transport indexed successfully but changed VCF record order." >&2
    else
        echo "$transport differs from $first_case" >&2
        exit 1
    fi

    max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/ {
        gsub(/^[[:space:]]+/, "", $2); print $2
    }' "$case_dir/time.txt")
    output_bytes=$(find "$case_dir" -type f -name '*.vcf.gz' -printf '%s\n' \
        | /usr/bin/awk '{sum += $1} END {print sum + 0}')
    printf '%s\t%s\t%s\t%s\t%s\n' "$transport" "$record_order" \
        "$((ended - started))" "$max_rss" "$output_bytes" \
        | tee -a "$summary"
done

echo "Required pileup transports produced identical indexed VCFs."
