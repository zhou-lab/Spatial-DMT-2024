#!/usr/bin/env bash

set -eo pipefail

root=$(cd "$(dirname "$0")" && pwd)
tool=${CPG_TO_YAME:-$root/cpg_to_yame}
yame=${YAME:-/mnt/isilon/zhoulab/labsoftware/anaconda/conda_2026/bin/yame}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/cpg_to_yame.test.XXXXXX")
tool_command=("$tool")
if [[ ${VALGRIND:-0} == 1 ]]; then
    tool_command=(
        valgrind --quiet --leak-check=full --errors-for-leak-kinds=all
        --error-exitcode=99 "$tool"
    )
fi

cleanup() {
    find "$tmp" -mindepth 1 -delete
    rmdir "$tmp"
}
trap cleanup EXIT

cat > "$tmp/cpg.bed" <<'EOF'
chr1	9	11
chr1	19	21
chr1	29	31
EOF

cat > "$tmp/joint.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
chr1	10	.	C	.	.	.	CX=CG;N5=AACAA	GT:CV:BT	0/0:2:0.500	0/0:1:0.000
chr1	11	.	G	.	.	.	CX=CG;N5=AACCA	GT:CV:BT	0/0:2:1.000	./.:0:.
chr1	21	.	G	.	.	.	CX=CG;N5=AACTA	GT:CV:BT	0/0:3:0.333	./.:0:.
EOF

cat > "$tmp/A.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A
chr1	10	.	C	.	.	.	CX=CG;N5=AACAA	GT:CV:BT	0/0:2:0.500
chr1	11	.	G	.	.	.	CX=CG;N5=AACCA	GT:CV:BT	0/0:2:1.000
chr1	21	.	G	.	.	.	CX=CG;N5=AACTA	GT:CV:BT	0/0:3:0.333
EOF

mkdir "$tmp/batch"
"${tool_command[@]}" index "$tmp/cpg.bed" "$tmp/cpg.idx"
[[ $("${tool_command[@]}" rows "$tmp/cpg.idx") == 3 ]]
"${tool_command[@]}" pack-vcf -i "$tmp/cpg.idx" -o "$tmp/A.single.cg" \
    -c "$tmp/A.single.conv.tsv" "$tmp/A.vcf"
"${tool_command[@]}" pack-vcf-batch -i "$tmp/cpg.idx" -o "$tmp/batch" \
    "$tmp/joint.vcf"

cmp "$tmp/A.single.cg" "$tmp/batch/A.cg"
cmp "$tmp/A.single.conv.tsv" "$tmp/batch/A.conv.tsv"
diff -u <(printf 'CA\t0.500000\t1\nCC\t1.000000\t1\nCT\t0.333000\t1\n') \
    "$tmp/batch/A.conv.tsv"
diff -u <(printf 'CA\t0.000000\t1\n') "$tmp/batch/B.conv.tsv"

"$yame" unpack -a -f -1 "$tmp/batch/A.cg" > "$tmp/A.unpacked"
"$yame" unpack -a -f -1 "$tmp/batch/B.cg" > "$tmp/B.unpacked"
diff -u <(printf '3\t1\n1\t2\n0\t0\n') "$tmp/A.unpacked"
diff -u <(printf '0\t1\n0\t0\n0\t0\n') "$tmp/B.unpacked"

cat > "$tmp/unsafe.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	../escape
EOF
if "${tool_command[@]}" pack-vcf-batch -i "$tmp/cpg.idx" -o "$tmp/batch" \
        "$tmp/unsafe.vcf"; then
    echo "unsafe VCF sample name was accepted" >&2
    exit 1
fi

cat > "$tmp/duplicate.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	A
EOF
if "${tool_command[@]}" pack-vcf-batch -i "$tmp/cpg.idx" -o "$tmp/batch" \
        "$tmp/duplicate.vcf"; then
    echo "duplicate VCF sample name was accepted" >&2
    exit 1
fi

printf 'cpg_to_yame tests passed\n'
