#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: validate_integrated_pileup_cell.sh BAM BASELINE_CG BASELINE_CONV OUTROOT

Runs the actual optimized dna_biscuit_pileup rule on one cell and compares
its packed CpG and conversion outputs with fixed legacy-path references.
EOF
    exit 2
}

[[ $# -eq 4 ]] || usage

bam=$(realpath -s "$1")
baseline_cg=$(realpath -s "$2")
baseline_conv=$(realpath -s "$3")
outroot=$(realpath -sm "$4")
script_dir=$(realpath -s "$(dirname "$0")")
repo_dir=$(realpath -s "$script_dir/../..")
sample=PileupRule23
workdir="$outroot/work"
sample_dir="$workdir/$sample"
bam_dir="$sample_dir/bam"
cpus=${SLURM_CPUS_PER_TASK:-2}

[[ -r $bam && -r $baseline_cg && -r $baseline_conv \
    && $cpus =~ ^[1-9][0-9]*$ ]] || usage

if [[ -r $bam.csi ]]; then
    index=$bam.csi
    index_name=0501.bam.csi
elif [[ -r $bam.bai ]]; then
    index=$bam.bai
    index_name=0501.bam.bai
else
    echo "No BAI or CSI found for $bam" >&2
    exit 1
fi

mkdir -p "$bam_dir" "$outroot/final" "$outroot/fastq"
# The production rule intentionally discovers regular per-cell BAMs only.
# Hard links avoid copying the deep fixture while preserving that file type.
ln -f -- "$bam" "$bam_dir/0501.bam"
ln -f -- "$index" "$bam_dir/$index_name"
printf '0501\n' > "$sample_dir/align_list.txt"

target="$sample_dir/pileup/$sample.cg"
snakemake \
    --snakefile "$repo_dir/main/Snakefile" \
    --profile "$repo_dir/profiles/local_HPC" \
    --cores "$cpus" \
    --jobs 2 \
    --resources cpus_per_task="$cpus" \
    --rerun-incomplete \
    --allowed-rules build_local_mergecg_context dna_biscuit_pileup \
    --config \
        IDS="$sample" \
        FASTQ_DIR="$outroot/fastq" \
        WORKDIR="$workdir" \
        OUTDIR="$outroot/final" \
        ref=mm10 \
        protocol=spatialdmt \
        DNA_PILEUP_CPUS="$cpus" \
        PILEUP_THREADS_PER_CELL=1 \
        PILEUP_WINDOW_STEP=1000000 \
    -- "$target"

per_cell_dir="$sample_dir/tmp/pileup"
cmp "$baseline_cg" "$per_cell_dir/0501.cg"
cmp "$baseline_cg" "$target"
cmp "$baseline_conv" "$per_cell_dir/0501.conv.tsv"
[[ -s $target.idx ]]
tabix -l "$per_cell_dir/0501.vcf.gz" >/dev/null
[[ -z $(find "$per_cell_dir" -maxdepth 1 -type f -name '*.raw.vcf*' \
    -print -quit) ]]

echo "Integrated pileup rule matches the fixed deep-cell references."
