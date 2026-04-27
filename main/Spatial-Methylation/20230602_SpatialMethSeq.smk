configfile: "/mnt/isilon/zhoulab/labpipelines/references/CHOP_HPC.yaml"
ref = config[config["ref"]]
IDS = config["ID"].split(",")
IDS_AUTO, = glob_wildcards("fastq/{id,[^/]+}_1.fq.gz")
FEATS = ["ChromHMM.20220414", "Win100k.20220228", "Win1m.20230709"]
FEATS_CH = ["TrinucCHxChromHMM.20220414", "GeneCAC.20220322", "Dinuc.20231030", "TrinucCH.20220321"]

rule trim_all:
	input: expand("tmp/trim/{id}/{id}_barcodes.txt", id=IDS)

rule trim_adapters:
	input:
		r1 = "fastq/{id}_1.fq.gz",
		r2 = "fastq/{id}_2.fq.gz",
	output:
		"tmp/trim/{id}/{id}_R1.fq.gz",
		"tmp/trim/{id}/{id}_R2.fq.gz",
		"tmp/trim/{id}/{id}_stats.txt",
		"tmp/trim/{id}/{id}_barcodes.txt",
	resources:
		mem = "120G",
		time = "96:00:00",
		cpus_per_task = 24,
	shell:
		"""
		set -xe
		mkdir -p tmp/trim/{wildcards.id}/
                zcat fastq/{wildcards.id}_1.fq.gz | split -l 40000000 - tmp/trim/{wildcards.id}/_chunk_1
                zcat fastq/{wildcards.id}_2.fq.gz | split -l 40000000 - tmp/trim/{wildcards.id}/_chunk_2
                source activate_2022 Bio
                parallel -j 24 'f={{}};b=$(basename $f);chunkID=${{b#_chunk_1}};d=$(dirname $f);f2=$d/_chunk_2$chunkID; /mnt/isilon/zhoulab/labtools/pyutils/spatialmeth_trimadapters.py $f $f2 -o tmp/trim/{wildcards.id}/_out_$chunkID -a CTATCTCTTATA AGATGCGAGAAGCCAACGCTTG AATCATACACCAATACAAAACATCAACCAC CAAATACTCTAACCTCTCAAACACATAAAT' ::: tmp/trim/{wildcards.id}/_chunk_1*

                :>tmp/trim/{wildcards.id}/{wildcards.id}_R1.fq.gz
                :>tmp/trim/{wildcards.id}/{wildcards.id}_R2.fq.gz
                for f in tmp/trim/{wildcards.id}/_out*R1.fq.gz; do
                	cat $f >>tmp/trim/{wildcards.id}/{wildcards.id}_R1.fq.gz
                	cat ${{f%_R1.fq.gz}}_R2.fq.gz >>tmp/trim/{wildcards.id}/{wildcards.id}_R2.fq.gz
                done
                cat tmp/trim/{wildcards.id}/_out_*_barcodes.txt | sort -k1,1n | \
                	bedtools groupby -g 1 -c 2 -o sum | sort -k2,2nr >tmp/trim/{wildcards.id}/{wildcards.id}_barcodes.txt
                cat tmp/trim/{wildcards.id}/_out_*_stats.txt | sort -k1,1 | \
                	bedtools groupby -g 1 -c 2 -o sum >tmp/trim/{wildcards.id}/{wildcards.id}_stats.txt
                rm -f tmp/trim/{wildcards.id}/_out*
                rm -f tmp/trim/{wildcards.id}/_chunk*
		"""

rule demultiplex_all:
	input: expand("tmp/dmux/{id}/{id}_R1_dmux.txt", id=IDS)
                
rule demultiplex:
	input:
		r1 = "tmp/trim/{id}/{id}_R1.fq.gz",
		r2 = "tmp/trim/{id}/{id}_R2.fq.gz",
                bc = "barcodes/spatial_barcodes.txt",
		# bc = "tmp/trim/{id}/{id}_barcodes.txt",
	output: "tmp/dmux/{id}/{id}_R1_dmux.txt", "tmp/dmux/{id}/{id}_R2_dmux.txt"
	resources:
		mem = "60G",
		time = "168:00:00",
		cpus_per_task = 12,
	shell:
		"""
		set -xe
		mkdir -p tmp/dmux/{wildcards.id}/
                dmux() {{
			awk -v OFS="\\t" -v ind=$1 -F"\\t" '
                	NR==FNR{{barcodes[$4]=0;}}
                	NR!=FNR{{barcode=substr($1,2,16);
                		if (barcode in barcodes) {{
	        			print $1"\\n"$2"\\n"$3"\\n"$4 > \
                			"tmp/dmux/{wildcards.id}/"barcode"_"ind".fq";
                			 barcodes[barcode] += 1;
                		}} }}
                	END {{
                		for(barcode in barcodes) {{ print barcode,barcodes[barcode]; }}
                	}}' {input.bc} \
                <(zcat tmp/trim/{wildcards.id}/{wildcards.id}_$1.fq.gz | paste - - - -) > \
                tmp/dmux/{wildcards.id}/{wildcards.id}_$1_dmux.txt;
                }}
                export -f dmux
                parallel -j 2 'dmux {{}}' ::: R1 R2
		parallel -j 12 'gzip -f {{}}' ::: tmp/dmux/{wildcards.id}/*.fq
		"""
		
# IDS2 = ["DMETLA1_CKDL230017286-1A_H5WY3DSX7_L2"]
# IDS3, BARCODES, = glob_wildcards("tmp/dmux/{id}/{barcode,[^/]+}_R1.fq.gz")

rule biscuit_align_all:
	input: expand("tmp/bam/{id}/align.txt", id=IDS)

rule biscuit_align_mm10_pe:
	input:
		r1 = "tmp/dmux/{id}/{id}_R1_dmux.txt",
		r2 = "tmp/dmux/{id}/{id}_R2_dmux.txt",
	output: "tmp/bam/{id}/align.txt",
	resources:
		mem = "150G",
		time = "96:00:00",
		cpus_per_task = 48,
	params:
		biscuit_index = ref["REF_BASE"]+ref["BISCUIT_INDEX"],
		ref_fasta = ref["REF_BASE"]+ref["REF_FASTA"],
	shell: # spatial library is directionary, read 2 is C-less and read 1 is G-less
		"""
		set -xe
                mkdir -p tmp/bam/{wildcards.id}/
                :>{output}
                for fq1 in tmp/dmux/{wildcards.id}/*_R1.fq.gz; do
			barcode=$(basename $fq1 _R1.fq.gz);
			out_dup=tmp/bam/{wildcards.id}/${{barcode}}.dupsifter.stat
			out_bam=tmp/bam/{wildcards.id}/${{barcode}}.bam
			out_flg=tmp/bam/{wildcards.id}/${{barcode}}.bam.flagstat
			biscuit align {params.biscuit_index} -b 1 -@ {resources.cpus_per_task} tmp/dmux/{wildcards.id}/${{barcode}}_R2.fq.gz tmp/dmux/{wildcards.id}/${{barcode}}_R1.fq.gz | dupsifter --stats-output $out_dup {params.ref_fasta} | samtools sort -T ${{out_bam}}_tmp -O bam -o $out_bam
			samtools index $out_bam
			samtools flagstat $out_bam >$out_flg
			echo $barcode >>{output}
		done
		"""

rule biscuit_pileup_all:
	input: expand("pileup/{id}/{id}.cg", id=IDS)

rule biscuit_pileup:
	input: "tmp/bam/{id}/align.txt"
	output:
		cg = "pileup/{id}/{id}.cg",
                idx = "pileup/{id}/{id}.cg.idx",
	resources:
		mem = "120G",
		time = "48:00:00",
		cpus_per_task = 48,
	params:
		ref_fasta = ref["REF_BASE"]+ref["REF_FASTA"],
		cpg_bed = ref["REF_BASE"]+ref["CPG_BED"],
	shell:
		"""
		set -xe
		mkdir -p tmp/pileup/{wildcards.id}/
		for bam in tmp/bam/{wildcards.id}/*.bam; do
			barcode=$(basename $bam .bam)
                	out_vcf=tmp/pileup/{wildcards.id}/$barcode.vcf.gz
                	biscuit pileup -m 0 -a 0 -c -u -p -@ {resources.cpus_per_task} {params.ref_fasta} \
                		tmp/bam/{wildcards.id}/$barcode.bam | bgzip -c >$out_vcf
			tabix -p vcf $out_vcf
                done

                vcf2cg() {{
                	biscuit vcf2bed -k 1 -t cg $1 | cut -f1-5 | \
                		LC_ALL=C sort -k1,1 -k2,2n -T tmp/pileup/{wildcards.id} | \
                		biscuit mergecg {params.ref_fasta} - | \
                		bedtools intersect -a {params.cpg_bed} -b - -sorted -loj | \
                		awk -v OFS="\\t" -F"\\t" -f wanding.awk -e '{{
                			M=round($7*$8); print M,$8-M;}}' | \
                		yame pack -f 3 - ;
                }}
                export -f vcf2cg
                parallel -j {resources.cpus_per_task} 'b={{/}}; vcf2cg {{}} \
                	>tmp/pileup/{wildcards.id}/${{b%.vcf.gz}}.cg' ::: \
                	tmp/pileup/{wildcards.id}/*.vcf.gz
                
                :>{output.cg}
                for cg in tmp/pileup/{wildcards.id}/*.cg; do
			barcode=$(basename $cg .cg)
			cat $cg >>{output.cg}
			yame index -1 $barcode {output.cg}
                done
		"""

BARCODES = ["AATGTGATAATGTGAT"]
rule biscuit_qc_all:
	input: expand("tmp/BISCUITqc/{id}/{id}_{barcode}_mapq_table.txt", id=IDS, barcode=BARCODES)
	output: "multiqc_report.html"
	resources:
		mem = "5G",
		time = "1:00:00",
		cpus_per_task = 2,
	shell:
		"""
		source activate_2022 multiqc
                rm -rf multiqc_*
		multiqc tmp/BISCUITqc/
		"""

rule biscuit_qc:
	input:
		cg = "pileup/{id}/{id}.cg",
		bam = "tmp/bam/{id}/{barcode}.bam",
	output:
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_all_base_botgc_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_all_base_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_all_base_topgc_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_all_cpg_botgc_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_all_cpg_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_all_cpg_topgc_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_q40_base_botgc_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_q40_base_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_q40_base_topgc_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_q40_cpg_botgc_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_q40_cpg_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_covdist_q40_cpg_topgc_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_cv_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_mapq_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_strand_table.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_dup_report.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_totalBaseConversionRate.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_totalReadConversionRate.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_CpHRetentionByReadPos.txt',
		'tmp/BISCUITqc/{id}/{id}_{barcode}_CpGRetentionByReadPos.txt',
	params:
		ref_fasta = ref["REF_BASE"]+ref["REF_FASTA"],
		assets = ref["REF_BASE"]+ref['BISCUIT_QCASSET'],
	resources:
		mem = "10G",
		time = "24:00:00",
		cpus_per_task = 2,
	shell:
		"""
		set -xe
		mkdir -p tmp/BISCUITqc/{wildcards.id}
		QC.sh -o tmp/BISCUITqc/{wildcards.id} --vcf tmp/pileup/{wildcards.id}/{wildcards.barcode}.vcf.gz {params.assets} {params.ref_fasta} {wildcards.id}_{wildcards.barcode} {input.bam}
		"""

rule feature_mean_all:
	input: expand("features/{id}/{feat}.txt.gz", id=IDS, feat=FEATS)

rule feature_mean:
	input:
		cg = "pileup/{id}/{id}.cg",
		feat = "/mnt/isilon/zhou_lab/projects/20191221_references/mm10/features/{feat}.bed.gz",
	output:
		txt = "features/{id}/{feat}.txt.gz",
	resources:
		mem = "120G",
		time = "6:00:00",
		cpus_per_task = 24,
	shell:
		"""
		set -xe
		mkdir -p features/{wildcards.id}/
                tmp=tmp/features/{wildcards.id}_{wildcards.feat}/
                rm -f $tmp/*
                mkdir -p $tmp
                cut -f1 {input.cg}.idx | parallel -j {resources.cpus_per_task} \
                	'yame subset {input.cg} {{}} | \
                		yame summary -q {{}} -m ~/references/mm10/KYCGKB_mm10/{wildcards.feat}.cm - | \
                		gzip -c >tmp/features/{wildcards.id}_{wildcards.feat}/{{}}.txt.gz'
                zcat $tmp/*.txt.gz | awk 'NR==1||$1!="QFile"' | gzip -c >features/{wildcards.id}/{wildcards.feat}.txt.gz
		"""

rule biscuit_pileup_allc_all:
	input: expand("pileup/{id}/{id}.allc", id=IDS)

rule biscuit_pileup_allc:
	input: "pileup/{id}"
	output:
            allc = "pileup/{id}/{id}.allc",
	resources:
		mem = "250G",
		time = "168:00:00",
		cpus_per_task = 15,
	shell:
		"""
		set -xe
		mkdir -p tmp/pileup/{wildcards.id}

                BetaCov2MU () {{ awk -F"\t" -v OFS="\t" -f wanding.awk -e \
                	'$4=="."{{print 0,0;next;}}{{M=round($4*$5); print M,$5-M;}}' $1; }};
                export -f BetaCov2MU

                ## the following might get into a memory outflow
                find tmp/pileup/{wildcards.id} -name '*.vcf.gz' | parallel -j 15 '
                barcode=$(basename {{}} .vcf.gz);
                biscuit vcf2bed -k 1 -t c {{}} | \
                	LC_ALL=C sort -k1,1 -k2,2n -T tmp/pileup/{wildcards.id} >tmp/pileup/{wildcards.id}/allc_$barcode.bed;
                bedtools intersect -a ~/references/mm10/annotation/allc/allc.bed.gz \
                	-b tmp/pileup/{wildcards.id}/allc_$barcode.bed -sorted -loj | \
                	cut -f1,2,3,10,11 | BetaCov2MU | yame pack -f 3 - >tmp/pileup/{wildcards.id}/$barcode.allc;'

                :>{output.allc}
                for f in tmp/pileup/{wildcards.id}/*.allc; do
                	barcode=$(basename $f .allc)
                	cat $f >>{output.allc}
                	yame index -1 $barcode {output.allc}
                done
		"""

## The feature means of the allc files shouldn't be done through full expansion
rule feature_mean_allc_all:
	input: expand("features_allc/{id}/{feat}.txt.gz", id=IDS, feat=FEATS_CH)

rule feature_mean_allc:
	input:
            	align = "tmp/bam/{id}/align.txt",
		feat = "/mnt/isilon/zhou_lab/projects/20191221_references/mm10/features/{feat}.bed.gz",
	output: "features_allc/{id}/{feat}.txt.gz",
	resources:
		mem = "150G",
		time = "72:00:00",
		cpus_per_task = 24,
	params:
		allc_cr = ref["REF_BASE"]+ref["ALLC_CR"],
	shell:
		"""
                set -xe
		mkdir -p features_allc/{wildcards.id}/
                mkdir -p tmp/features_allc/{wildcards.id}/{wildcards.feat}/

                tallyGroup() {{ awk -v barcode=$2 -v OFS="\\t" -F"\\t" '{{ cnt[$2]+=1; sum[$2]+=$1; }}
                END{{ for(k in cnt){{print k,sum[k]/cnt[k],cnt[k],"{wildcards.id}",barcode; }} }}' $1; }}
                export -f tallyGroup

                find tmp/pileup/{wildcards.id} -name '*.vcf.gz' | parallel -j {resources.cpus_per_task} '
                barcode=$(basename {{}} .vcf.gz);
                biscuit vcf2bed -k 1 -t c {{}} | bedtools intersect -a - -b {input.feat} -sorted -wo | \
                	cut -f4,9 | tallyGroup - "$barcode" | sort -k1,1 \
                	>tmp/features_allc/{wildcards.id}/{wildcards.feat}/$barcode.tmp'
                
                cat tmp/features_allc/{wildcards.id}/{wildcards.feat}/*.tmp | gzip -c \
                	>features_allc/{wildcards.id}/{wildcards.feat}.txt.gz
		"""
