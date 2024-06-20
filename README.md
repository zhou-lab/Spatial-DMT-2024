# Spatial joint profiling of DNA methylome and transcriptome
### Introduction
This repository aims to share the raw data processing and downstream data analysis & visualization codes used in the Spatial-MT-seq project.

Next Generation Sequencing (NGS) NGS was conducted on an Illumina NovaSeq 6000 sequencer (pair-end, 150-base-pair mode)
![]( https://github.com/zhou-lab/Spatial-MT-seq-2024/blob/ea6589f2ce5ea9fc3a55fe506de7e67b97ec7cd0/workflow/Experiment_pipeline.jpg)

### Spatial_Methylation-seq/Spatial_RNA_seq analysis pipeline
![]( https://github.com/zhou-lab/Spatial-MT-seq-2024/blob/ea6589f2ce5ea9fc3a55fe506de7e67b97ec7cd0/workflow/Analysis_pipeline.jpg)

### Running Spatial-MT snakemake files
## Dependiencies
* [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html). 
* [Biopython](https://biopython.org/docs/1.75/api/index.html).
* [Python3]( https://docs.python.org/3/using/unix.html).
* [Biscuit](https://huishenlab.github.io/biscuit/#download-and-install).
* [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html).
* [awk](https://manpages.ubuntu.com/manpages/trusty/man1/awk.1posix.html).
* [STAR](https://github.com/alexdobin/STAR).

## Run the pipeline
For RNA: Change all the directories in the Snakefile
To obtain RNA count matrices: sbatch Snakemake.sh
For methylation: Change the config ID to the data ID number
To obtain CG levels: runSnakemake --config ID=SpMETSLE17DM ref=mm10 --snakefile /mnt/isilon/zhoulab/labpipelines/Snakefiles/20230602_SpatialMethSeq.smk feature_mean_all
To obtain Biscuit QC results: runSnakemake --config ID=SpMETSLE17DM ref=mm10 --snakefile /mnt/isilon/zhoulab/labpipelines/Snakefiles/20230602_SpatialMethSeq.smk biscuit_qc_all
To obtain CH levels: runSnakemake --config ID=SpMETSLE17DM ref=mm10 --snakefile /mnt/isilon/zhoulab/labpipelines/Snakefiles/20230602_SpatialMethSeq.smk feature_mean_allc_all

### Data downstream analysis and visualization
The data analysis and visualization were completed with R language (4.4.0). 

