#!/usr/bin/env bash
set -euo pipefail

script_dir=$(realpath -s "$(dirname "$0")")
repo_dir=$(realpath -s "$script_dir/../..")
sample=${SAMPLE:-SpeedSmoke21}
fastq_dir=${FASTQ_DIR:-/scr1/users/suc3/20260818_dna_speedopt/e2e_fastq}
output_dir=${OUTPUT_DIR:-/scr1/users/suc3/20260818_dna_speedopt/e2e_output_v1}
cpus=${CPUS:-20}
trim_mem=${TRIM_MEM:-80G}
align_mem=${ALIGN_MEM:-160G}
unmatched_mem=${UNMATCHED_MEM:-160G}

main_bam="$output_dir/$sample/bam_merged/${sample}_merged_dedup.bam"
unmatched_bam="$output_dir/$sample/bam_unmatched/${sample}_unmatched.bam"
lambda_bam="$output_dir/$sample/bam_lambda/${sample}_Lambda.bam"

mkdir -p "$output_dir"

exec snakemake \
    --snakefile "$repo_dir/main/Snakefile" \
    --profile "$repo_dir/profiles/slurm" \
    --cores "$cpus" \
    --jobs 4 \
    --resources cpus_per_task="$cpus" \
    --latency-wait 100 \
    --rerun-incomplete \
    --config \
        IDS="$sample" \
        FASTQ_DIR="$fastq_dir" \
        WORKDIR="$output_dir" \
        OUTDIR="$output_dir" \
        ref=mm10 \
        protocol=spatialdmt \
        DNA_TRIM_CPUS="$cpus" \
        DNA_TRIM_MEM="$trim_mem" \
        DNA_ALIGN_CPUS="$cpus" \
        DNA_ALIGN_MEM="$align_mem" \
        DNA_ALIGN_SORT_CPUS=8 \
        DNA_SPLIT_JOBS="$cpus" \
        DNA_SPLIT_BUCKET_THREADS=16 \
        LAMBDA_ALIGN_CPUS="$cpus" \
        LAMBDA_ALIGN_SORT_CPUS=8 \
        UNMATCHED_ALIGN_CPUS="$cpus" \
        UNMATCHED_ALIGN_MEM="$unmatched_mem" \
        UNMATCHED_SORT_CPUS=8 \
    -- "$main_bam" "$unmatched_bam" "$lambda_bam"
