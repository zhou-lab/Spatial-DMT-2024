# Spatial joint profiling of DNA methylome and transcriptome (Spatial-DMT)
## Introduction
This repository is for raw data processing and downstream data analysis & visualization code in the **"Spatial joint profiling of DNA methylome and transcriptome in mammalian tissues"**  manuscript.

### Computational workflow for spatial-DMT data processing and analysis
![](https://github.com/dyxmvp/Spatial-DMT-2024/blob/main/workflow/Analysis_pipeline.jpg)

## Data analysis instructions
### Dependencies

All dependencies can be installed via [conda](https://docs.conda.io/en/latest/) (bioconda/conda-forge channels) except BBTools, which is available via [Homebrew](https://brew.sh/) on macOS or direct download on Linux.

**Install via conda:**
```bash
conda install -c conda-forge -c bioconda \
    snakemake=9.19.0 \
    star=2.7.11b \
    biscuit=1.7.1 \
    samtools=1.22.1 \
    htslib=1.22.1 \
    bedtools=2.31.1 \
    dupsifter=1.3.0 \
    yame=1.9 \
    pigz=2.8 \
    parallel=20260422
```

**Install BBTools (provides `bbduk.sh`) separately:**
```bash
# macOS (Homebrew)
brew install bbtools   # tested: 39.81b

# Linux
# Download from https://jgi.doe.gov/data-and-tools/software-tools/bbtools/
```

| Tool | Version tested | Source |
|------|---------------|--------|
| [Snakemake](https://snakemake.readthedocs.io/) | 9.19.0 | conda-forge |
| [STAR](https://github.com/alexdobin/STAR) / STARsolo | 2.7.11b | bioconda |
| [BISCUIT](https://huishenlab.github.io/biscuit/) | 1.7.1 | bioconda |
| [samtools](https://www.htslib.org/) | 1.22.1 | bioconda |
| bgzip / tabix (htslib) | 1.22.1 | bioconda |
| [bedtools](https://bedtools.readthedocs.io/) | 2.31.1 | bioconda |
| [dupsifter](https://github.com/huishenlab/dupsifter) | 1.3.0 | bioconda |
| [yame](https://github.com/zhou-lab/yame) | 1.9 | bioconda |
| [pigz](https://zlib.net/pigz/) | 2.8 | conda-forge |
| [GNU parallel](https://www.gnu.org/software/parallel/) | 20260422 | conda-forge |
| [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (bbduk) | 39.81b | Homebrew / direct |
| [Python](https://www.python.org/) | 3.x | conda-forge |
| [Biopython](https://biopython.org/) | — | conda-forge |

### 1. Data preprocessing

Next Generation Sequencing (NGS) was performed using an Illumina NovaSeq 6000 sequencer (150bp paired-end mode). Read 1 contains the genome sequences, and Read 2 contains the spatial Barcode A & B and UMIs (mRNA).

The preprocessing pipeline processes both spatial DNA methylation and spatial RNA data in a single Snakemake workflow (`Data_preprocess/Snakefile`).

#### Input data structure

```
fastq/
  {sample}/
    DNA_1.fq.gz   # DNA methylation read 1 (genomic sequence)
    DNA_2.fq.gz   # DNA methylation read 2 (spatial barcode + adapter)
    RNA_1.fq.gz   # RNA read 1 (cDNA sequence)
    RNA_2.fq.gz   # RNA read 2 (spatial barcode + UMI)
barcodes/
  spatial_barcodes.txt       # full spatial barcode whitelist (multi-column)
  RNA_whitelist_singlecol.txt  # single-column barcode list for STARsolo
```

Samples are auto-detected from `fastq/*/DNA_1.fq*`, or specified explicitly with `--config IDS=sample1,sample2`.

#### Running the pipeline

**Local / HPC (no scheduler):**
```bash
snakemake --profile profiles/local_HPC \
    --snakefile Data_preprocess/Snakefile \
    --config ref=hg38
```

**SLURM cluster:**
```bash
snakemake --profile profiles/slurm \
    --snakefile Data_preprocess/Snakefile \
    --config ref=hg38
```

Supported values for `ref`: `hg38`, `mm10` (must match a key in `Data_preprocess/genomes.yaml`).

Key profile config options (set in `profiles/local_HPC/config.yaml` or overridden with `--config`):

| Key | Default | Description |
|-----|---------|-------------|
| `ref` | — | Genome reference key (`hg38` or `mm10`) |
| `REF_DIR` | `/mnt/isilon/zhou_lab/projects/20191221_references/` | Root directory of reference files |
| `BBDUK` | `bbduk.sh` | Path to bbduk.sh (defaults to PATH) |
| `PYTHON` | `python` | Python interpreter with biopython, fuzzysearch, matplotlib, and pandas installed |

#### Pipeline steps

**DNA methylation branch:**

1. `dna_trim_and_demux` — Trim adapters and demultiplex reads by spatial barcode using `spatialmeth_trimadapters.py`; outputs per-barcode FASTQ pairs and merged trimmed FASTQs
2. `dna_biscuit_align` — Bisulfite-aware alignment of each barcode's reads to the genome with BISCUIT, duplicate marking with dupsifter, BAM sort and index
3. `dna_biscuit_pileup` — CpG methylation pileup per barcode BAM, VCF→BED conversion, strand merging, CpG BED intersection, and packing into a yame `.cg` file indexed by barcode

**RNA branch (runs in parallel with DNA):**

4. `rna_filter_primer` — bbduk: retain R2 reads containing the spatial primer sequence
5. `rna_filter_L1` — bbduk: retain reads containing linker 1 sequence
6. `rna_filter_L2` — bbduk: retain reads containing linker 2 sequence
7. `rna_fq_process` — Extract barcode (BC2+BC1, 16 bp) and UMI (10 bp) from fixed positions in R2 and reformat for STARsolo
8. `rna_star_solo` — STARsolo alignment and per-barcode gene expression quantification

#### Outputs

```
pipeline_output/{sample}/
  trim/          # merged trimmed DNA FASTQs (non-lambda and lambda)
  dmux/          # per-barcode demuxed DNA FASTQ pairs (.fq.gz)
  bam/           # per-barcode BAMs with index and dupsifter stats
  pileup/        # {sample}.cg and {sample}.cg.idx (yame packed methylation)
  rna_processed/
    tmp/         # intermediate filtered FASTQs
    qc/          # bbduk filter stats per step
    align/       # Aligned.out.sam, STAR logs, Solo.out/ (count matrix)
```

### 2. Quality control

Run the QC pipeline after preprocessing completes. It produces a MultiQC report aggregating STAR alignment and bbduk filter statistics, spatial heatmaps of per-barcode metrics (DNA read depth, mapping rate, duplication rate), and a self-contained HTML report embedding all plots.

```bash
snakemake --profile profiles/local_HPC \
    --snakefile Data_preprocess/Snakefile \
    --config ref=hg38 \
    -- qc
```

#### QC outputs

```
qc_output/{sample}/
  table/
    {sample}_barcode_counts.tsv          # per-barcode DNA reads, mapping rate, dup rate
    {sample}_spatial_barcode_counts.tsv  # counts mapped to (x, y) spatial grid
    {sample}_barcodes_bool_index.tsv     # all observed barcodes with spatial flag
  plots/
    spatial_dna_reads.png                # read depth heatmap
    spatial_log_dna_reads.png            # log10 read depth heatmap
    spatial_mapping_rate.png             # per-barcode mapping rate
    spatial_dup_rate.png                 # per-barcode duplication rate
    barcode_rank.png                     # barcode rank vs count
  multiqc_report.html                    # MultiQC report (STAR + bbduk stats)
  qc_report.html                         # self-contained HTML with embedded spatial maps
```

#### Cleaning temporary files

Intermediate VCF and per-barcode `.cg` files from the pileup step are stored under `pipeline_output/{sample}/tmp/`. The spatial plot PNGs are stored under `qc_output/{sample}/plots/` (the HTML report embeds them, so the PNGs can be removed). Run after both preprocessing and QC are complete:

```bash
snakemake --profile profiles/local_HPC \
    --snakefile Data_preprocess/Snakefile \
    --config ref=hg38 \
    -- clean
```

#### - For Image processing:
This python script processes a phase contrast tissue image to generate spatial metadata similar to 10x Visium outputs. It first converts the image to grayscale and applies adaptive thresholding to create a binary mask that distinguishes tissue from background. The image is then divided into a fixed 50×50 grid, and each grid cell is evaluated to determine whether it contains tissue based on a pixel intensity threshold. Detected tissue spots are matched with predefined spatial barcodes, and the results are compiled into output files, including a tissue positions CSV, a scale factors JSON, and a visualization image with highlighted tissue regions.

### 2. Downstream data analysis and visualization 

The data analysis and visualization were performed using R(4.4.0).

`QC_check.Rmd`: Contains all the code to generate QC plots

`WNN_DNAmRNA_clustering`: Contains all the code to integrate two modalities into Seurat objects, identify differential marker genes and VMRs, and visualize them on spatial maps

`Integration_E11E13.Rmd`: Contains all the code to integrate between day 11 and day 13 mouse embryos’ data.
 
`Utility_function.R`: Contains all functions used including mapping VMR to overlapping genes and iterative PCA. 


