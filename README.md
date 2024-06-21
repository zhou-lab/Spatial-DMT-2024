# Spatial joint profiling of DNA methylome and transcriptome (Spatial-DMT)
## Introduction
This repository is for raw data processing and downstream data analysis & visualization code in the **"Spatial joint profiling of DNA methylome and transcriptome in mammalian tissues"**  manuscript.

### Computational workflow for spatial-DMT data processing and analysis
![](https://github.com/dyxmvp/Spatial-DMT-2024/blob/main/workflow/Analysis_pipeline.jpg)

## Data analysis instructions
### Dependiencies
* [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
* [Biopython](https://biopython.org/docs/1.75/api/index.html)
* [Python3]( https://docs.python.org/3/using/unix.html)
* [BISCUIT](https://huishenlab.github.io/biscuit/#download-and-install)
* [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
* [awk](https://manpages.ubuntu.com/manpages/trusty/man1/awk.1posix.html)
* [STAR](https://github.com/alexdobin/STAR)
* [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)

### 1. Data preprocessing
Next Generation Sequencing (NGS) was performed using an Illumina NovaSeq 6000 sequencer (150bp paired-end mode). Read 1 contains the genome sequences, and Read 2 contains the spatial Barcode A & B and UMIs (mRNA).
 
The preprocessing pipeline for both spatial DNA methylation and spatial RNA data was built upon the Snakemake workflow management.

#### - For RNA data: 
Change all the directories in the Snakefile to obtain RNA count matrices:

    sbatch Snakemake.sh
    
#### Descriptions:

(1) Setup of directories and files: Automate the generation of directories to store each sample's raw and processed data. 

(2)`filter_primer`: Use bbduk to filter sequences containing primers from the reads.

(3)`filter_L1` & `filter_L2`: Apply additional filters to select reads with specific linker sequences.

(4)`fq_process`: Extract spatial barcodes and UMIs and reformats data. 

(5)`star_solo`: Align reads to a reference genome (e.g. mm10) using STAR.

#### - For DNA methylation data: 

Change the config ID to the data ID number.
To obtain BISCUIT QC results: 

    runSnakemake --config ID=SpMETSLE17DM ref=mm10 --snakefile /mnt/isilon/zhoulab/labpipelines/Snakefiles/20230602_SpatialMethSeq.smk biscuit_qc_all
    
To obtain CG levels:

    runSnakemake --config ID=SpMETSLE17DM ref=mm10 --snakefile /mnt/isilon/zhoulab/labpipelines/Snakefiles/20230602_SpatialMethSeq.smk feature_mean_all`

To obtain CH levels:

    runSnakemake --config ID=SpMETSLE17DM ref=mm10 --snakefile /mnt/isilon/zhoulab/labpipelines/Snakefiles/20230602_SpatialMethSeq.smk feature_mean_allc_all

#### Descriptions:

(1)`trim_all`: Trim the fastq files using spatialmeth_trimadapters.py. 

(2)`demultiplex_all`: Split all the reads based on the barcodes, obtain 2500 fastq files.

(3)`biscuit_align_all `: Align reads to a reference genome (e.g. mm10) using BISCUIT.

(4)`biscuit_pileup_all`: Identify all the CG and call the methylation at those sites.

(5)`biscuit_qc_all`: Quality check for alignment and methylation calling. 

(6)`feature_mean_all`: Obtain average methylation over selected windows.

(7)`biscuit_pileup_allc_all`: Identify all the CH and call the methylation at those sites.

#### - For Image processing: 
Identify the location of pixels on tissue from the brightfield image using tissue_positions_list.py.

### 2. Downstream data analysis and visualization 

The data analysis and visualization were performed using R(4.4.0).

`QC_check.Rmd`: Contains all the code to generate QC plots

`WNN_DNAmRNA_clustering`: Contains all the code to integrate two modalities into Seurat objects, identify differential marker genes and VMRs, and visualize them on spatial maps

`Integration_E11E13.Rmd`: Contains all the code to integrate between day 11 and day 13 mouse embryos’ data.
 
`Utility_function.R`: Contains all functions used including mapping VMR to overlapping genes and iterative PCA. 


