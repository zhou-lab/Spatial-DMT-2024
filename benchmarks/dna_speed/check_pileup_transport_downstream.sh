#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: check_pileup_transport_downstream.sh SORTED_VCF DIRECT_VCF REF_FASTA CPG_BED YAME OUTDIR

Runs the production two-pass CpG/conversion path and a one-pass path on VCFs
in natural-sort and direct BISCUIT order. The packed .cg and conversion
summaries must match byte for byte.
EOF
    exit 2
}

[[ $# -eq 6 ]] || usage

sorted_vcf=$(realpath -s "$1")
direct_vcf=$(realpath -s "$2")
ref_fasta=$(realpath -s "$3")
cpg_bed=$(realpath -s "$4")
yame=$(realpath -s "$5")
outdir=$(realpath -sm "$6")
script_dir=$(realpath -s "$(dirname "$0")")
wanding_awk=$(realpath -s "$script_dir/../../tools/wanding.awk")
mergecg_context=$(realpath -s "$script_dir/../../tools/mergecg_context/mergecg_context")

for path in "$sorted_vcf" "$direct_vcf" "$ref_fasta" "$cpg_bed" \
        "$yame" "$wanding_awk" "$mergecg_context"; do
    [[ -r $path ]] || {
        echo "Input is not readable: $path" >&2
        exit 1
    }
done

mkdir -p "$outdir/results"
summary="$outdir/summary.tsv"
printf 'case\tvcf_order\tconversion_mode\tmerge_backend\twall_seconds\tmax_rss_kb\tcg_bytes\tconv_bytes\n' > "$summary"

process_vcf_twopass() {
    local vcf=$1 prefix=$2
    biscuit vcf2bed -k 1 -t cg "$vcf" \
        | cut -f1-5 \
        | LC_ALL=C sort -k1,1 -k2,2n \
        | biscuit mergecg "$ref_fasta" - \
        | bedtools intersect -a "$cpg_bed" -b - -sorted -loj \
        | awk -v OFS='\t' -F'\t' -f "$wanding_awk" \
            -e '{ M=round($7*$8); print M,$8-M; }' \
        | "$yame" pack -f 3 - > "$prefix.cg"

    biscuit vcf2bed -e -t c "$vcf" \
        | awk -F'\t' '
            { sum[$6] += $8; cnt[$6] += 1 }
            END {
                for (c in sum) printf "%s\t%.6f\t%d\n", c, sum[c], cnt[c]
            }' \
        | LC_ALL=C sort > "$prefix.conv.tsv"
}

merge_cpg_stream() {
    if [[ $1 == context ]]; then
        "$mergecg_context" -
    else
        biscuit mergecg "$ref_fasta" -
    fi
}

prepare_cpg_stream() {
    if [[ $1 == context ]]; then
        awk -F'\t' '$5 == "CG"'
    else
        awk -F'\t' -v OFS='\t' '$5 == "CG" {
            print $1,$2,$3,$8,$9
        }'
    fi
}

process_vcf_onepass() {
    local vcf=$1 prefix=$2 merge_backend=$3
    local context_fifo conv_pid main_status conv_status
    context_fifo="$prefix.context.fifo"
    mkfifo "$context_fifo"

    awk -F'\t' '
        { sum[$6] += $8; cnt[$6] += 1 }
        END {
            for (c in sum) printf "%s\t%.6f\t%d\n", c, sum[c], cnt[c]
        }' < "$context_fifo" \
        | LC_ALL=C sort > "$prefix.conv.tsv" &
    conv_pid=$!

    set +e
    biscuit vcf2bed -e -t c "$vcf" \
        | tee "$context_fifo" \
        | prepare_cpg_stream "$merge_backend" \
        | LC_ALL=C sort -k1,1 -k2,2n \
        | merge_cpg_stream "$merge_backend" \
        | bedtools intersect -a "$cpg_bed" -b - -sorted -loj \
        | awk -v OFS='\t' -F'\t' -f "$wanding_awk" \
            -e '{ M=round($7*$8); print M,$8-M; }' \
        | "$yame" pack -f 3 - > "$prefix.cg"
    main_status=$?
    wait "$conv_pid"
    conv_status=$?
    set -e
    rm -f "$context_fifo"
    if (( main_status != 0 || conv_status != 0 )); then
        echo "one-pass conversion failed: CpG=$main_status conversion=$conv_status" >&2
        return 1
    fi
}

run_processor() {
    "$1" "$2" "$3" "$4"
}
export -f process_vcf_twopass merge_cpg_stream prepare_cpg_stream
export -f process_vcf_onepass run_processor
export ref_fasta cpg_bed yame wanding_awk mergecg_context

cases=${DOWNSTREAM_CASES:-"sorted_twopass:sorted:twopass:biscuit sorted_onepass:sorted:onepass:biscuit direct_onepass:direct:onepass:biscuit direct_context:direct:onepass:context"}
reference_case=""
for spec in $cases; do
    IFS=: read -r case_name vcf_order conversion_mode merge_backend <<< "$spec"
    if [[ $vcf_order == sorted ]]; then
        vcf=$sorted_vcf
    else
        vcf=$direct_vcf
    fi
    if [[ $conversion_mode == twopass ]]; then
        processor=process_vcf_twopass
    else
        processor=process_vcf_onepass
    fi
    tabix -l "$vcf" >/dev/null
    case_dir="$outdir/results/$case_name"
    mkdir -p "$case_dir"
    prefix="$case_dir/cell"
    started=$(date +%s)
    /usr/bin/time -v -o "$case_dir/time.txt" \
        bash -o pipefail -c 'run_processor "$@"' _ \
            "$processor" "$vcf" "$prefix" "$merge_backend"
    ended=$(date +%s)
    max_rss=$(/usr/bin/awk -F: '/Maximum resident set size/ {
        gsub(/^[[:space:]]+/, "", $2); print $2
    }' "$case_dir/time.txt")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$case_name" "$vcf_order" "$conversion_mode" "$merge_backend" \
        "$((ended - started))" "$max_rss" "$(stat -Lc %s "$prefix.cg")" \
        "$(stat -Lc %s "$prefix.conv.tsv")" | tee -a "$summary"

    if [[ -z $reference_case ]]; then
        reference_case=$case_name
    else
        cmp "$outdir/results/$reference_case/cell.cg" "$prefix.cg"
        cmp "$outdir/results/$reference_case/cell.conv.tsv" \
            "$prefix.conv.tsv"
    fi
done

echo "All VCF transport and conversion modes produced identical downstream outputs."
