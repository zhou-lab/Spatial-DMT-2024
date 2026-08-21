#!/usr/bin/env bash

set -eo pipefail

usage() {
    cat <<'EOF'
Usage: run_biscuit_pileup_batches.sh [OPTIONS]

Required options:
  --cell-list FILE             Ordered cell names, one per line
  --bam-dir DIR                Directory containing CELL.bam files
  --tmp-dir DIR                Restart/staging directory
  --output-cg FILE             Aggregate YAME output
  --output-index FILE          Aggregate YAME index output
  --output-conversion FILE     Sample conversion-rate table
  --reference FILE             Reference FASTA
  --cpg-converter FILE         cpg_to_yame executable
  --cpg-index FILE             Shared CpG coordinate index
  --yame FILE                  yame executable
  --planner FILE               plan_pileup_batches.py
  --python FILE                Python executable
  --cpus N                     Total CPU allocation
  --batch-size N               Maximum cells per joint pileup
  --max-batch-bytes N          Target aggregate BAM bytes per batch (0=off)
  --threads-per-batch N        BISCUIT threads per concurrent batch
  --cpg-threads N              VCF reader threads per converter
  --window-step N              BISCUIT dispatch window

Optional:
  --biscuit FILE               BISCUIT executable (default: biscuit)
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

cell_list=
bam_dir=
tmp_pileup=
output_cg=
output_index=
output_conversion=
ref_fasta=
cpg_converter=
cpg_index=
yame=
planner=
python=
biscuit=biscuit
cpus=
batch_size=
max_batch_bytes=
pileup_threads=
cpg_threads=
pileup_window_step=

while (($#)); do
    case "$1" in
        --cell-list) cell_list=$2; shift 2 ;;
        --bam-dir) bam_dir=$2; shift 2 ;;
        --tmp-dir) tmp_pileup=$2; shift 2 ;;
        --output-cg) output_cg=$2; shift 2 ;;
        --output-index) output_index=$2; shift 2 ;;
        --output-conversion) output_conversion=$2; shift 2 ;;
        --reference) ref_fasta=$2; shift 2 ;;
        --cpg-converter) cpg_converter=$2; shift 2 ;;
        --cpg-index) cpg_index=$2; shift 2 ;;
        --yame) yame=$2; shift 2 ;;
        --planner) planner=$2; shift 2 ;;
        --python) python=$2; shift 2 ;;
        --biscuit) biscuit=$2; shift 2 ;;
        --cpus) cpus=$2; shift 2 ;;
        --batch-size) batch_size=$2; shift 2 ;;
        --max-batch-bytes) max_batch_bytes=$2; shift 2 ;;
        --threads-per-batch) pileup_threads=$2; shift 2 ;;
        --cpg-threads) cpg_threads=$2; shift 2 ;;
        --window-step) pileup_window_step=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown argument: $1" ;;
    esac
done

for variable in cell_list bam_dir tmp_pileup output_cg output_index \
        output_conversion ref_fasta cpg_converter cpg_index yame planner python \
        cpus batch_size max_batch_bytes pileup_threads cpg_threads \
        pileup_window_step; do
    [[ -n ${!variable} ]] || die "missing --${variable//_/-}"
done
for executable in "$cpg_converter" "$yame" "$planner" "$python"; do
    [[ -e $executable ]] || die "missing tool: $executable"
done
if [[ $biscuit == */* ]]; then
    [[ -x $biscuit ]] || die "missing BISCUIT executable: $biscuit"
else
    biscuit=$(command -v "$biscuit") \
        || die "BISCUIT executable not found in PATH"
fi
for input in "$cell_list" "$ref_fasta" "$cpg_index"; do
    [[ -s $input ]] || die "missing or empty input: $input"
done
[[ -d $bam_dir ]] || die "missing BAM directory: $bam_dir"
for value in "$cpus" "$batch_size" "$max_batch_bytes" "$pileup_threads" \
        "$cpg_threads" "$pileup_window_step"; do
    [[ $value =~ ^[0-9]+$ ]] || die "numeric options must be nonnegative integers"
done
(( cpus > 0 && batch_size > 0 && pileup_threads > 0 \
   && cpg_threads > 0 && pileup_window_step > 0 )) \
    || die "CPU, batch, thread, and window values must be positive"
(( batch_size <= 4096 )) || die "batch size cannot exceed 4096 VCF samples"

if (( pileup_threads > cpus )); then pileup_threads=$cpus; fi
if (( cpg_threads > pileup_threads )); then cpg_threads=$pileup_threads; fi
pileup_jobs=$((cpus / pileup_threads))
if (( pileup_jobs < 1 )); then pileup_jobs=1; fi
cpg_jobs=$((cpus / cpg_threads))
if (( cpg_jobs < 1 )); then cpg_jobs=1; fi

batch_dir="$tmp_pileup/batches"
mkdir -p "$tmp_pileup" "$batch_dir" "$(dirname "$output_cg")" \
    "$(dirname "$output_index")" "$(dirname "$output_conversion")"
expected_cpg_rows=$("$cpg_converter" rows "$cpg_index")
[[ $expected_cpg_rows =~ ^[0-9]+$ ]] \
    || die "invalid CpG row count: $expected_cpg_rows"

export tmp_pileup batch_dir bam_dir cpg_converter cpg_index yame ref_fasta \
    pileup_window_step pileup_threads cpg_threads expected_cpg_rows biscuit

cell_complete() {
    local bam=$1 cell cg conv actual_rows
    cell=$(basename "$bam" .bam)
    cg="$tmp_pileup/$cell.cg"
    conv="$tmp_pileup/$cell.conv.tsv"
    [[ -s $cg && -e $conv \
       && $cg -nt $bam && $conv -nt $bam \
       && $cg -nt $ref_fasta && $conv -nt $ref_fasta \
       && $cg -nt $cpg_index && $conv -nt $cpg_index \
       && $cg -nt $biscuit && $conv -nt $biscuit \
       && $cg -nt $cpg_converter && $conv -nt $cpg_converter ]] || return 1
    actual_rows=$("$yame" info "$cg" | awk -F'\t' 'NR == 2 { print $4 }')
    [[ $actual_rows == "$expected_cpg_rows" ]]
}
export -f cell_complete

# Reuse valid per-cell VCFs left by an older or interrupted pipeline. Failure
# here falls through to a fresh direct-stream batch instead of blocking it.
process_legacy_vcf() {
    set -eo pipefail
    local vcf=$1 cell bam stage actual_rows
    cell=$(basename "$vcf" .vcf.gz)
    bam="$bam_dir/$cell.bam"
    [[ -s $bam ]] || return 0
    if cell_complete "$bam"; then return 0; fi
    if [[ ! $vcf -nt $bam || ! $vcf -nt $ref_fasta \
            || ! $vcf -nt $biscuit ]]; then
        return 0
    fi

    stage="$batch_dir/legacy-$cell.stage.$$"
    mkdir "$stage"
    if ! "$cpg_converter" pack-vcf -@ "$cpg_threads" \
            -i "$cpg_index" -o "$stage/$cell.cg" \
            -c "$stage/$cell.conv.tsv" "$vcf"; then
        echo "WARNING: cannot reuse legacy VCF for $cell; scheduling pileup" >&2
        rmdir "$stage" 2>/dev/null || true
        return 0
    fi
    if ! actual_rows=$("$yame" info "$stage/$cell.cg" \
            | awk -F'\t' 'NR == 2 { print $4 }'); then
        echo "WARNING: cannot validate legacy VCF for $cell; scheduling pileup" >&2
        return 0
    fi
    if [[ $actual_rows != "$expected_cpg_rows" ]]; then
        echo "WARNING: invalid legacy CpG dimensions for $cell; scheduling pileup" >&2
        return 0
    fi
    mv "$stage/$cell.cg" "$tmp_pileup/$cell.cg"
    mv "$stage/$cell.conv.tsv" "$tmp_pileup/$cell.conv.tsv"
    rmdir "$stage"
}
export -f process_legacy_vcf

find -L "$tmp_pileup" -maxdepth 1 -name '*.vcf.gz' -type f \
    -printf '%s\t%p\n' | sort -nr | cut -f2- \
    | parallel --halt soon,fail=1 -j "$cpg_jobs" process_legacy_vcf

ref_stamp=$(stat -c '%s:%Y' "$ref_fasta")
pileup_signature="ref=$ref_fasta:$ref_stamp;args=m0,a0,c,u,p;step=$pileup_window_step;threads=$pileup_threads;batch=$batch_size;bytes=$max_batch_bytes"
manifest_list="$batch_dir/current.manifests"
"$python" "$planner" \
    --cell-list "$cell_list" --bam-dir "$bam_dir" \
    --output-dir "$batch_dir" --manifest-list "$manifest_list" \
    --batch-size "$batch_size" --max-batch-bytes "$max_batch_bytes" \
    --signature "$pileup_signature"

process_batch() {
    set -eo pipefail
    local manifest=$1 batch_name stage bam cell actual_rows complete
    local -a bams
    mapfile -t bams < "$manifest"
    ((${#bams[@]} > 0)) || return 1
    complete=1
    for bam in "${bams[@]}"; do
        if ! cell_complete "$bam"; then complete=0; fi
    done
    if (( complete )); then return 0; fi

    batch_name=$(basename "$manifest" .bams)
    stage="$batch_dir/$batch_name.stage.$$"
    mkdir "$stage"
    if ! "$biscuit" pileup -m 0 -a 0 -c -u -p \
            -@ "$pileup_threads" -s "$pileup_window_step" \
            "$ref_fasta" "${bams[@]}" \
            | "$cpg_converter" pack-vcf-batch -@ "$cpg_threads" \
                -i "$cpg_index" -o "$stage" -; then
        echo "ERROR: pileup/conversion failed for $manifest" >&2
        return 1
    fi

    for bam in "${bams[@]}"; do
        cell=$(basename "$bam" .bam)
        if [[ ! -s $stage/$cell.cg || ! -e $stage/$cell.conv.tsv ]]; then
            echo "ERROR: missing staged outputs for $cell in $manifest" >&2
            return 1
        fi
        actual_rows=$("$yame" info "$stage/$cell.cg" \
            | awk -F'\t' 'NR == 2 { print $4 }')
        if [[ $actual_rows != "$expected_cpg_rows" ]]; then
            echo "ERROR: invalid CpG dimensions for $cell: $actual_rows != $expected_cpg_rows" >&2
            return 1
        fi
    done
    for bam in "${bams[@]}"; do
        cell=$(basename "$bam" .bam)
        mv "$stage/$cell.cg" "$tmp_pileup/$cell.cg"
        mv "$stage/$cell.conv.tsv" "$tmp_pileup/$cell.conv.tsv"
    done
    rmdir "$stage"
}
export -f process_batch

parallel --halt soon,fail=1 -j "$pileup_jobs" \
    process_batch :::: "$manifest_list"

aggregate_stage="$batch_dir/aggregate.stage.$$"
mkdir "$aggregate_stage"
aggregate_cg="$aggregate_stage/$(basename "$output_cg")"
aggregate_samples="$aggregate_cg.samples"
aggregate_conversion="$aggregate_stage/$(basename "$output_conversion")"
: > "$aggregate_cg"
: > "$aggregate_samples"
while IFS= read -r cell; do
    [[ -n $cell ]] || continue
    bam="$bam_dir/$cell.bam"
    cell_complete "$bam" || die "incomplete final CpG record for $cell"
    cat "$tmp_pileup/$cell.cg" >> "$aggregate_cg"
    printf '%s\n' "$cell" >> "$aggregate_samples"
done < "$cell_list"
"$yame" index -s "$aggregate_samples" "$aggregate_cg"

{
    echo "BISCUITqc Conversion Rate by Base Average Table"
    printf 'CA\tCC\tCG\tCT\n'
    while IFS= read -r cell; do
        [[ -n $cell ]] || continue
        cat "$tmp_pileup/$cell.conv.tsv"
    done < "$cell_list" | awk -F'\t' '
        { sum[$1] += $2; cnt[$1] += $3 }
        END {
            for (c in sum) m[c] = (cnt[c] >= 20) ? sum[c]/cnt[c] : -1
            printf "%g\t%g\t%g\t%g\n",
                (("CA" in m)?m["CA"]:-1), (("CC" in m)?m["CC"]:-1),
                (("CG" in m)?m["CG"]:-1), (("CT" in m)?m["CT"]:-1)
        }'
} > "$aggregate_conversion"

mv "$aggregate_cg" "$output_cg"
mv "$aggregate_cg.idx" "$output_index"
mv "$aggregate_conversion" "$output_conversion"
mv "$aggregate_samples" "$tmp_pileup/$(basename "$output_cg").samples"
rmdir "$aggregate_stage"
