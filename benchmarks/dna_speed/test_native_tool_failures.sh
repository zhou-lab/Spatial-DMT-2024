#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: test_native_tool_failures.sh NAME_SORTED.bam REF_FASTA OUTDIR

Builds a small complete-template fixture from NAME_SORTED.bam and verifies
that dupsifter and the compiled CB splitter return nonzero for failed output
and truncated input. It also checks mergecg_context output and parse failures.
Set DUPSIFTER_NEW, CB_SPLITTER, or MERGECG_CONTEXT to override binaries.
EOF
    exit 2
}

[[ $# -eq 3 ]] || usage
name_bam=$(realpath -s "$1")
ref_fasta=$(realpath -s "$2")
outdir=$(realpath -sm "$3")
repo_dir=$(realpath -s "$(dirname "$0")/../..")
dupsifter=${DUPSIFTER_NEW:-$repo_dir/tools/dupsifter_speed/dupsifter}
splitter=${CB_SPLITTER:-$repo_dir/tools/split_bam_by_cb_hts/split_bam_by_cb_hts}
mergecg_context=${MERGECG_CONTEXT:-$repo_dir/tools/mergecg_context/mergecg_context}

[[ -r $name_bam && -r $ref_fasta && -x $dupsifter && -x $splitter \
    && -x $mergecg_context ]] || usage
mkdir -p "$outdir"
find "$outdir" -mindepth 1 -delete

fixture="$outdir/complete.name.bam"
samtools view -h "$name_bam" \
    | awk '
        /^@/ {print; next}
        $1 != previous {
            groups++
            previous = $1
        }
        groups <= 100 {print}
    ' \
    | samtools view -b -o "$fixture" -
samtools quickcheck -v "$fixture"
(( $(samtools view -c "$fixture") > 0 ))

expect_failure() {
    local label=$1
    shift
    set +e
    "$@" >"$outdir/$label.stdout" 2>"$outdir/$label.stderr"
    local status=$?
    set -e
    if (( status == 0 )); then
        echo "$label unexpectedly succeeded" >&2
        exit 1
    fi
    printf '%s\tstatus=%s\n' "$label" "$status"
}

expect_failure dupsifter_output_full \
    "$dupsifter" --uncompressed-bam -o /dev/full \
    --stats-output "$outdir/output_full.stat" "$ref_fasta" "$fixture"

fixture_bytes=$(stat -c %s "$fixture")
head -c "$((fixture_bytes / 2))" "$fixture" > "$outdir/truncated.name.bam"
expect_failure dupsifter_truncated_input \
    "$dupsifter" --uncompressed-bam -o "$outdir/truncated.output.bam" \
    --stats-output "$outdir/truncated.stat" \
    "$ref_fasta" "$outdir/truncated.name.bam"

mkdir -p "$outdir/split"
expect_failure splitter_truncated_input \
    "$splitter" -@ 2 "$outdir/truncated.name.bam" "$outdir/split"

printf '%s\n' \
    $'chr1\t9\t10\tC\tCG\tCG\tAACGT\t0.500\t2' \
    $'chr1\t10\t11\tG\tCG\tCG\tAACGT\t0.667\t3' \
    $'chr1\t20\t21\tC\tCG\tCG\tAACGT\t1.000\t1' \
    $'chr1\t30\t31\tG\tCG\tCG\tAACGT\t0.000\t2' \
    > "$outdir/context.bed"
printf '%s\n' \
    $'chr1\t9\t11\t0.600\t5\tC:0.500:2,G:0.667:3' \
    $'chr1\t20\t22\t1.000\t1\tC:1.000:1,G:.:0' \
    $'chr1\t29\t31\t0.000\t2\tC:.:0,G:0.000:2' \
    > "$outdir/context.expected"
"$mergecg_context" "$outdir/context.bed" > "$outdir/context.actual"
cmp "$outdir/context.expected" "$outdir/context.actual"

printf 'malformed\n' > "$outdir/context.malformed.bed"
expect_failure mergecg_context_malformed \
    "$mergecg_context" "$outdir/context.malformed.bed"
expect_failure mergecg_context_output_full \
    bash -o pipefail -c '"$1" "$2" > /dev/full' _ \
        "$mergecg_context" "$outdir/context.bed"

grep -Eq 'ERROR|Failed|failed|truncated' "$outdir/dupsifter_output_full.stderr"
grep -Eq 'ERROR|Failed|failed|truncated' "$outdir/dupsifter_truncated_input.stderr"
grep -Eq 'ERROR|Failed|failed|truncated' "$outdir/splitter_truncated_input.stderr"
grep -Eq 'Expected|Malformed|malformed' "$outdir/mergecg_context_malformed.stderr"
grep -Eq 'space|stdout|write' "$outdir/mergecg_context_output_full.stderr"
echo "All native-tool failure checks passed."
