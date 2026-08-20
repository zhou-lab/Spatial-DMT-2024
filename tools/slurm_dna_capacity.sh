#!/usr/bin/env bash
set -euo pipefail

partition=${1:-defq}
reserve_cpus=${2:-5}
max_rule_cpus=${3:-80}
reserve_mem_gb=${4:-40}

for value in "$reserve_cpus" "$max_rule_cpus" "$reserve_mem_gb"; do
    [[ $value =~ ^[0-9]+$ ]] || {
        echo "Usage: slurm_dna_capacity.sh [PARTITION] [RESERVE_CPUS] [MAX_RULE_CPUS] [RESERVE_MEM_GB]" >&2
        exit 2
    }
done

partition_record=$(scontrol show partition "$partition" -o)
partition_qos=$(sed -n 's/.* QoS=\([^ ]*\).*/\1/p' <<< "$partition_record")
[[ -n $partition_qos && $partition_qos != "(null)" ]] || {
    echo "Partition $partition has no partition QOS" >&2
    exit 1
}

max_tres=$(sacctmgr -n -P show qos "$partition_qos" format=MaxTRESPerUser \
    | head -1)
cpu_limit=$(tr ',' '\n' <<< "$max_tres" \
    | sed -n 's/^cpu=\([0-9][0-9]*\)$/\1/p')
mem_limit=$(tr ',' '\n' <<< "$max_tres" \
    | sed -n 's/^mem=\([0-9][0-9.]*[KMGT]\{0,1\}\)$/\1/p')
[[ -n $cpu_limit && -n $mem_limit ]] || {
    echo "Could not parse CPU and memory limits from QOS $partition_qos: $max_tres" >&2
    exit 1
}

memory_to_mb() {
    /usr/bin/awk -v value="$1" 'BEGIN {
        if (!match(value, /^([0-9]+([.][0-9]+)?)([KMGT]?)$/, part)) exit 1
        factor["K"] = 1 / 1024
        factor[""] = factor["M"] = 1
        factor["G"] = 1024
        factor["T"] = 1024 * 1024
        printf "%.0f\n", part[1] * factor[part[3]]
    }'
}

mem_limit_mb=$(memory_to_mb "$mem_limit") || {
    echo "Could not convert QOS memory limit: $mem_limit" >&2
    exit 1
}

running_cpus=$(squeue -h -u "$USER" -t RUNNING -o '%C' \
    | /usr/bin/awk '{sum += $1} END {print sum + 0}')
running_memory_mb=0
while read -r job_id; do
    [[ -n $job_id ]] || continue
    if ! job_record=$(scontrol show job "$job_id" -o 2>/dev/null); then
        continue
    fi
    alloc_tres=$(sed -n 's/.* AllocTRES=\([^ ]*\).*/\1/p' <<< "$job_record")
    job_mem=$(tr ',' '\n' <<< "$alloc_tres" \
        | sed -n 's/^mem=\([0-9][0-9.]*[KMGT]\{0,1\}\)$/\1/p')
    [[ -n $job_mem ]] || continue
    job_mem_mb=$(memory_to_mb "$job_mem") || {
        echo "Could not convert memory allocation for job $job_id: $job_mem" >&2
        exit 1
    }
    running_memory_mb=$((running_memory_mb + job_mem_mb))
done < <(squeue -h -u "$USER" -t RUNNING -o '%i')

available_cpus=$((cpu_limit - running_cpus - reserve_cpus))
(( available_cpus < 0 )) && available_cpus=0
memory_budget_mb=$((mem_limit_mb - running_memory_mb - reserve_mem_gb * 1024))
(( memory_budget_mb < 0 )) && memory_budget_mb=0
empty_memory_budget_mb=$((mem_limit_mb - reserve_mem_gb * 1024))
(( empty_memory_budget_mb < 0 )) && empty_memory_budget_mb=0

empty_rule_cpus=$(((cpu_limit - reserve_cpus) / 2))
(( empty_rule_cpus > max_rule_cpus )) && empty_rule_cpus=$max_rule_cpus
(( empty_rule_cpus < 0 )) && empty_rule_cpus=0
current_rule_cpus=$((available_cpus / 2))
(( current_rule_cpus > max_rule_cpus )) && current_rule_cpus=$max_rule_cpus
single_rule_cpus=$available_cpus
(( single_rule_cpus > max_rule_cpus )) && single_rule_cpus=$max_rule_cpus
single_rule_empty_cpus=$((cpu_limit - reserve_cpus))
(( single_rule_empty_cpus > max_rule_cpus )) && single_rule_empty_cpus=$max_rule_cpus

printf 'checked_at\t%s\n' "$(date --iso-8601=seconds)"
printf 'partition\t%s\n' "$partition"
printf 'partition_qos\t%s\n' "$partition_qos"
printf 'max_tres_per_user\t%s\n' "$max_tres"
printf 'cpu_limit\t%s\n' "$cpu_limit"
printf 'memory_limit\t%s\n' "$mem_limit"
printf 'memory_limit_mb\t%s\n' "$mem_limit_mb"
printf 'running_cpus\t%s\n' "$running_cpus"
printf 'running_memory_mb\t%s\n' "$running_memory_mb"
printf 'reserve_cpus\t%s\n' "$reserve_cpus"
printf 'reserve_mem_gb\t%s\n' "$reserve_mem_gb"
printf 'controller_memory_mb\t%s\n' "$memory_budget_mb"
printf 'controller_memory_mb_empty_account\t%s\n' "$empty_memory_budget_mb"
printf 'available_cpus_now\t%s\n' "$available_cpus"
printf 'single_rule_cpus_now\t%s\n' "$single_rule_cpus"
printf 'paired_rule_cpus_now\t%s\n' "$current_rule_cpus"
printf 'controller_cpus_now\t%s\n' "$((current_rule_cpus * 2))"
printf 'single_rule_cpus_empty_account\t%s\n' "$single_rule_empty_cpus"
printf 'paired_rule_cpus_empty_account\t%s\n' "$empty_rule_cpus"
printf 'controller_cpus_empty_account\t%s\n' "$((empty_rule_cpus * 2))"
