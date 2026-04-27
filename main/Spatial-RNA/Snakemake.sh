#!/bin/bash
### need to run in the st environment

cores=45
snakemake -j $cores -s Snakefile --rerun-incomplete --resources --cluster 'sbatch -t 60 --mem=30g -c 45'


