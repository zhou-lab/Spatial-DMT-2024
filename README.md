# Spatial joint profiling of DNA methylome and transcriptome
## Introduction
This repository aims to share the raw data processing and downstream data analysis & visualization codes used in the Spatial DNA Methylation and RNA Transcriptomic Sequencing (Spatial-DMT-seq) project.Next Generation Sequencing (NGS) NGS was conducted on an Illumina NovaSeq 6000 sequencer (pair-end, 150-base-pair mode)

### Spatial_Methylation-seq/Spatial_RNA_seq analysis pipeline
![]( https://github.com/zhou-lab/Spatial-MT-seq-2024/blob/ea6589f2ce5ea9fc3a55fe506de7e67b97ec7cd0/workflow/Analysis_pipeline.jpg)

## Data analysis instructions
### Dependiencies
* [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html). 
* [Biopython](https://biopython.org/docs/1.75/api/index.html).
* [Python3]( https://docs.python.org/3/using/unix.html).
* [Biscuit](https://huishenlab.github.io/biscuit/#download-and-install).
* [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html).
* [awk](https://manpages.ubuntu.com/manpages/trusty/man1/awk.1posix.html).
* [STAR](https://github.com/alexdobin/STAR).

### 1. Data preprocessing
Next Generation Sequencing (NGS) NGS was conducted on an Illumina NovaSeq 6000 sequencer (pair-end, 150-base-pair mode) Read 1 contains the genome sequences, and Read 2 contains the spatial Barcode A & B.
 
The preprocessing pipeline for preprocessing both spatial transcriptomics and DNA methylation data was built upon the Snakemake workflow management. To run them,please update the correct directory for the snakemake file.

For RNA: Change all the directories in the Snakefile To obtain RNA count matrices: `sbatch Snakemake.sh `

Descriptions:

(1)Setup of Directories and Files. This automates the generation of directories to hold each sample's raw and processed data. 

(2)`filter_primer`: Uses bbduk to eliminate sequences containing primers from the reads.

(3)`filter_L1` & `filter_L2`: Applies additional filters to remove specific linker sequences from the reads.

(4)`fq_process`: Extracts and reformats data using fastq_process.py to obtain segments that could represent barcodes and UMI. 

(5)`star_solo`: Aligns reads to a reference genome (mm10) using STAR.

For methylation: Change the config ID to the data ID number

To obtain CG levels: `runSnakemake --config ID=SpMETSLE17DM ref=mm10 --snakefile /mnt/isilon/zhoulab/labpipelines/Snakefiles/20230602_SpatialMethSeq.smk feature_mean_all `

To obtain Biscuit QC results: `runSnakemake --config ID=SpMETSLE17DM ref=mm10 --snakefile /mnt/isilon/zhoulab/labpipelines/Snakefiles/20230602_SpatialMethSeq.smk biscuit_qc_all `

To obtain CH levels: `runSnakemake --config ID=SpMETSLE17DM ref=mm10 --snakefile /mnt/isilon/zhoulab/labpipelines/Snakefiles/20230602_SpatialMethSeq.smk feature_mean_allc_all `

Descriptions:

(1)`trim_all`: Trim the fastq files using spatialmeth_trimadapters.py 

(2)`demultiplex_all`: Split all the reads based on the barcodes, obtain 2500 fastq files

(3)`biscuit_align_all `: Aligns reads to a reference genome (mm10) using BISCUIT

(4)`biscuit_pileup_all`: Identify all the CG and call the methylation at those sites.

(5)`biscuit_qc_all`: Qualification check of alignment and methylation calling. 

(6)`feature_mean_all`: Run enrichment testing using YAME

(7)`biscuit_pileup_allc_all`: Identify all the CH and call the methylation at those sites.

For image: Identifies the location of pixels on tissue from the brightfield image using tissue_positions_list.py.

### 2. Downstream data analysis and visualization 

The data analysis and visualization were completed with R language (4.4.0).

QC_check.Rmd: Contains all the code to generate QC plots

WNN_DNAmRNA_clustering: Contains all the code to integrate two modalities into Seurat objects, identify differential marker genes and VMRs, and visualize them on spatial maps

Integration_E11E13.Rmd: Contains all the code to integrate between day 11 and day 13 mouse embryos’ data.
 
Utility_function.R: Contains all functions used including mapping VMR to overlapping genes and iterative PCA. 


