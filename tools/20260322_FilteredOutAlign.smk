configfile: "/mnt/isilon/zhoulab/labpipelines/references/CHOP_HPC.yaml"

ref = config[config["ref"]]
IDS = config["ID"].split(",")

rule all:
  input:
    expand("tmp/filteredout_bam/{id}/{id}_filteredout.bam", id=IDS),
    expand("tmp/filteredout_bam/{id}/{id}_filteredout.bam.bai", id=IDS),
    expand("tmp/filteredout_bam/{id}/{id}_filteredout.bam.flagstat", id=IDS),
    expand("tmp/filteredout/{id}/{id}_filteredout_stats.txt", id=IDS)


rule extract_filteredout_all:
  input:
    expand("tmp/filteredout/{id}/{id}_filteredout_R1.fq.gz", id=IDS),
    expand("tmp/filteredout/{id}/{id}_filteredout_R2.fq.gz", id=IDS)


rule extract_filteredout_reads:
  input:
    r1 = "fastq/{id}_1.fq.gz",
    r2 = "fastq/{id}_2.fq.gz",
  output:
    r1_out = "tmp/filteredout/{id}/{id}_filteredout_R1.fq.gz",
    r2_out = "tmp/filteredout/{id}/{id}_filteredout_R2.fq.gz",
    stats = "tmp/filteredout/{id}/{id}_filteredout_stats.txt",
  resources:
    mem = "80G",
    time = "48:00:00",
    cpus_per_task = 12,
  shell:
    """
    ulimit -n 65535
    set -xeo pipefail
    mkdir -p tmp/filteredout/{wildcards.id}/
    pigz -dc {input.r1} | split -d -a 4 -l 40000000 - tmp/filteredout/{wildcards.id}/_chunk_1
    pigz -dc {input.r2} | split -d -a 4 -l 40000000 - tmp/filteredout/{wildcards.id}/_chunk_2

    source activate_2022 Bio

    parallel --halt soon,fail=1 -j {resources.cpus_per_task} '
        f={{}}; b=$(basename "$f"); cid=${{b#_chunk_1}}
        d=$(dirname "$f"); f2="$d/_chunk_2$cid"
        python3 /home/fuh1/qc_check/spatialmeth_extract_linkerfail.py \
            "$f" "$f2" \
            -o tmp/filteredout/{wildcards.id}/_out_"$cid"
    ' ::: tmp/filteredout/{wildcards.id}/_chunk_1*

    :> {output.r1_out}
    :> {output.r2_out}
    for f in $(ls tmp/filteredout/{wildcards.id}/_out_*_filteredout_R1.fq.gz | sort -V); do
        cat "$f" >> {output.r1_out}
        cat "${{f%_filteredout_R1.fq.gz}}_filteredout_R2.fq.gz" >> {output.r2_out}
    done

    cat tmp/filteredout/{wildcards.id}/_out_*_filteredout_stats.txt \
      | sort -k1,1 \
      | bedtools groupby -g 1 -c 2 -o sum > {output.stats}

    rm -f tmp/filteredout/{wildcards.id}/_chunk_*
    rm -f tmp/filteredout/{wildcards.id}/_out_*_filteredout_R1.fq.gz tmp/filteredout/{wildcards.id}/_out_*_filteredout_R2.fq.gz
    rm -f tmp/filteredout/{wildcards.id}/_out_*_filteredout_stats.txt
    """


rule align_filteredout_all:
  input:
    expand("tmp/filteredout_bam/{id}/{id}_filteredout.bam.flagstat", id=IDS)


rule biscuit_align_filteredout_pe:
  input:
    r1 = "tmp/filteredout/{id}/{id}_filteredout_R1.fq.gz",
    r2 = "tmp/filteredout/{id}/{id}_filteredout_R2.fq.gz",
  output:
    bam = "tmp/filteredout_bam/{id}/{id}_filteredout.bam",
    bai = "tmp/filteredout_bam/{id}/{id}_filteredout.bam.bai",
    flagstat = "tmp/filteredout_bam/{id}/{id}_filteredout.bam.flagstat",
    dup = "tmp/filteredout_bam/{id}/{id}_filteredout.dupsifter.stat",
  resources:
    mem = "120G",
    time = "48:00:00",
    cpus_per_task = 24,
  params:
    biscuit_index = ref["REF_BASE"] + ref["BISCUIT_INDEX"],
    ref_fasta = ref["REF_BASE"] + ref["REF_FASTA"],
  shell:
    """
    set -xe
    mkdir -p tmp/filteredout_bam/{wildcards.id}/
    biscuit align {params.biscuit_index} -b 1 -@ {resources.cpus_per_task} {input.r2} {input.r1} | \
      dupsifter --stats-output {output.dup} {params.ref_fasta} | \
      samtools sort -T {output.bam}_tmp -O bam -o {output.bam}
    samtools index {output.bam}
    samtools flagstat {output.bam} > {output.flagstat}
    """

rule biscuit_pileup_filteredout_all:
  input:
    expand("tmp/pileup_filteredout/{id}/{id}_filteredout.vcf.gz", id=IDS)


rule biscuit_pileup_filteredout:
  input:
    bam = "tmp/filteredout_bam/{id}/{id}_filteredout.bam",
    bai = "tmp/filteredout_bam/{id}/{id}_filteredout.bam.bai",
  output:
    vcf = "tmp/pileup_filteredout/{id}/{id}_filteredout.vcf.gz",
    tbi = "tmp/pileup_filteredout/{id}/{id}_filteredout.vcf.gz.tbi",
  resources:
    mem = "120G",
    time = "24:00:00",
    cpus_per_task = 8,
  params:
    ref_fasta = ref["REF_BASE"] + ref["REF_FASTA"],
  shell:
    """
    set -xe
    mkdir -p tmp/pileup_filteredout/{wildcards.id}/
    biscuit pileup -m 0 -a 0 -c -u -p -@ {resources.cpus_per_task} {params.ref_fasta} {input.bam} | \
      bgzip -c > {output.vcf}
    tabix -p vcf {output.vcf}
    """


rule biscuit_qc_filteredout_all:
  input:
    expand("tmp/BISCUITqc/{id}/{id}_filteredout_mapq_table.txt", id=IDS)


rule biscuit_qc_filteredout:
  input:
    vcf = "tmp/pileup_filteredout/{id}/{id}_filteredout.vcf.gz",
    bam = "tmp/filteredout_bam/{id}/{id}_filteredout.bam",
  output:
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_all_base_botgc_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_all_base_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_all_base_topgc_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_all_cpg_botgc_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_all_cpg_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_all_cpg_topgc_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_q40_base_botgc_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_q40_base_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_q40_base_topgc_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_q40_cpg_botgc_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_q40_cpg_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_covdist_q40_cpg_topgc_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_cv_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_mapq_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_strand_table.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_dup_report.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_totalBaseConversionRate.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_totalReadConversionRate.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_CpHRetentionByReadPos.txt',
    'tmp/BISCUITqc/{id}/{id}_filteredout_CpGRetentionByReadPos.txt',
  params:
    ref_fasta = ref["REF_BASE"] + ref["REF_FASTA"],
    assets = ref["REF_BASE"] + ref['BISCUIT_QCASSET'],
  resources:
    mem = "10G",
    time = "24:00:00",
    cpus_per_task = 2,
  shell:
    """
    set -xe
    mkdir -p tmp/BISCUITqc/{wildcards.id}
    QC.sh -o tmp/BISCUITqc/{wildcards.id} --vcf {input.vcf} {params.assets} {params.ref_fasta} {wildcards.id}_filteredout {input.bam}
    """
