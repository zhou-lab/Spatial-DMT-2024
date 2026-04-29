# Spatial joint profiling of DNA methylome and transcriptome (Spatial-DMT)

This repository contains the data preprocessing and QC pipeline, plus downstream analysis and visualization code, for [Lee, Fu et al., *Nature* 2025](https://www.nature.com/articles/s41586-025-09478-x). The Snakemake workflow (`main/Snakefile`) processes raw sequencing data through four targets:

- **`all`** (default) — adapter trimming, spatial demultiplexing, bisulfite alignment, CpG methylation pileup, lambda spike-in alignment, and RNA alignment with STARsolo
- **`qc`** — BISCUIT QC tables, spatial heatmaps, MultiQC report, self-contained HTML report, and per-feature mean methylation summaries
- **`allc`** — all-cytosine pileup and non-CpG feature methylation summaries (run before `clean`)
- **`clean`** — remove intermediate VCF/cg files and spatial plot PNGs

![](https://github.com/dyxmvp/Spatial-DMT-2024/blob/main/workflow/Analysis_pipeline.jpg)

## Pipeline instructions
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
    parallel=20260422 \
    multiqc \
    matplotlib \
    pandas
```

**For SLURM execution only**, install the cluster-generic executor plugin into your Snakemake environment (not required for local/HPC runs):
```bash
pip install snakemake-executor-plugin-cluster-generic
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
| [snakemake-executor-plugin-cluster-generic](https://github.com/snakemake/snakemake-executor-plugin-cluster-generic) | 1.0.9 | pip (SLURM only) |
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
| [MultiQC](https://multiqc.info/) | 1.34 | conda-forge |
| [Python](https://www.python.org/) | 3.x | conda-forge |
| [Biopython](https://biopython.org/) | — | conda-forge |
| matplotlib | — | conda-forge |
| pandas | — | conda-forge |

### 1. Data preprocessing

Next Generation Sequencing (NGS) was performed using an Illumina NovaSeq 6000 sequencer (150bp paired-end mode). Read 1 contains the genome sequences, and Read 2 contains the spatial Barcode A & B and UMIs (mRNA).

The preprocessing pipeline processes both spatial DNA methylation and spatial RNA data in a single Snakemake workflow (`main/Snakefile`).

#### Input data structure

The `barcodes/` directory is provided with this repository. Users only need to supply the raw FASTQ files:

```
fastq/
  {sample}/
    DNA_1.fq.gz   # DNA methylation read 1 (genomic sequence)
    DNA_2.fq.gz   # DNA methylation read 2 (spatial barcode + adapter)
    RNA_1.fq.gz   # RNA read 1 (cDNA sequence)
    RNA_2.fq.gz   # RNA read 2 (spatial barcode + UMI)
```

Samples are auto-detected from `fastq/*/DNA_1.fq*`, or specified explicitly with `--config IDS=sample1,sample2`.

#### Running the pipeline

**Local / HPC (no scheduler):**
```bash
snakemake --profile profiles/local_HPC \
    --snakefile main/Snakefile \
    --config ref=hg38
```

**SLURM cluster:**
```bash
snakemake --profile profiles/slurm \
    --snakefile main/Snakefile \
    --config ref=hg38
```

Supported values for `ref`: `hg38`, `mm10` (must match a key defined in your `profiles/*/site.yaml`).

Required profile config keys (set in `profiles/local_HPC/config.yaml` or `profiles/slurm/config.yaml`):

| Key | Description |
|-----|-------------|
| `ref` | Genome reference key (`hg38` or `mm10`) |
| `REF_DIR` | Root directory of reference files |
| `BBDUK` | Path to `bbduk.sh` executable |
| `PYTHON` | Python interpreter with biopython, fuzzysearch, matplotlib, and pandas installed |

#### Pipeline steps

**DNA methylation branch:**

1. `dna_trim_and_demux` — Trim adapters and demultiplex reads by spatial barcode; outputs per-barcode FASTQ pairs, merged trimmed FASTQs, and lambda spike-in FASTQs
2. `dna_biscuit_align` — Bisulfite-aware alignment of each barcode's reads to the genome with BISCUIT, duplicate marking with dupsifter, BAM sort and index
3. `dna_biscuit_pileup` — CpG methylation pileup per barcode BAM, VCF→BED conversion, strand merging, CpG BED intersection, and packing into a yame `.cg` file indexed by barcode
4. `dna_biscuit_align_lambda` — Align lambda spike-in reads to the lambda genome for bisulfite conversion efficiency estimation
5. `dna_biscuit_pileup_lambda` — Pileup on lambda BAM; produces allc.bed for conversion efficiency calculation

**RNA branch (runs in parallel with DNA):**

6. `rna_filter_primer` — bbduk: retain R2 reads containing the spatial primer sequence
7. `rna_filter_L1` — bbduk: retain reads containing linker 1 sequence
8. `rna_filter_L2` — bbduk: retain reads containing linker 2 sequence
9. `rna_fq_process` — Extract barcode (BC2+BC1, 16 bp) and UMI (10 bp) from fixed positions in R2 and reformat for STARsolo
10. `rna_star_solo` — STARsolo alignment and per-barcode gene expression quantification

#### Outputs

```
pipeline_output/{sample}/
  trim/              # merged trimmed DNA FASTQs (non-lambda and lambda)
  dmux/              # per-barcode demuxed DNA FASTQ pairs (.fq.gz)
  bam/               # per-barcode BAMs with index and dupsifter stats
  bam_lambda/        # lambda BAM with flagstat and dupsifter stats
  pileup/            # {sample}.cg and {sample}.cg.idx (yame packed CpG methylation)
  pileup_lambda/     # lambda VCF, allc.bed, and cg.bed
  tmp/               # all intermediates removed by clean: pileup VCFs, per-barcode .cg files, filtered FASTQs
  rna_bbduk/         # bbduk filter stats per step
  rna_STAR/          # Aligned.out.sam, STAR logs, Solo.out/ (count matrix)
```

### 2. Quality control

Run after preprocessing. Produces BISCUIT QC tables for every barcode, a MultiQC report aggregating BISCUIT, STAR, and bbduk statistics, spatial heatmaps of per-barcode metrics, a self-contained HTML report, and per-feature mean methylation summaries.

```bash
snakemake --profile profiles/local_HPC \
    --snakefile main/Snakefile \
    --config ref=hg38 \
    -- qc
```

**QC steps:**

- `dna_biscuit_qc` — Run `QC.sh` (BISCUIT QC pipeline) on every per-barcode BAM; tables aggregated by MultiQC
- `qc_barcode_counts` — Per-barcode DNA read counts, mapping rate, and duplication rate
- `qc_spatial_match` — Map barcode counts to (x, y) spatial grid coordinates
- `qc_spatial_plots` — Spatial heatmaps of read depth, mapping rate, and duplication rate
- `qc_multiqc` — MultiQC report aggregating BISCUIT QC tables, STAR logs, and bbduk stats
- `qc_html_report` — Self-contained HTML report with embedded spatial plots
- `feature_mean` — Per-barcode mean methylation summarized over genomic features (ChromHMM, windows, chromosomes)

#### QC outputs

```
pipeline_output/{sample}/qc/
  biscuit_qc/        # per-barcode BISCUIT QC tables (input to MultiQC)
  features/          # per-feature mean methylation tables (.txt.gz per feature set)
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
  multiqc_report.html                    # MultiQC report (BISCUIT + STAR + bbduk)
  qc_report.html                         # self-contained HTML with embedded spatial maps
```

### 3. AllC (non-CpG methylation)

Optional. Run after preprocessing and **before** clean (requires tmp/ pileup VCFs). Produces all-cytosine pileup and per-feature non-CpG methylation summaries.

```bash
snakemake --profile profiles/local_HPC \
    --snakefile main/Snakefile \
    --config ref=hg38 \
    -- allc
```

#### AllC outputs

```
pipeline_output/{sample}/
  pileup/{sample}.allc        # yame-packed all-cytosine methylation indexed by barcode
pipeline_output/{sample}/qc/
  features_allc/              # per-feature non-CpG mean methylation tables (.txt.gz)
```

### 4. Cleaning temporary files

Removes `pipeline_output/{sample}/tmp/` and `pipeline_output/{sample}/qc/tmp/` (pileup VCFs, per-barcode `.cg` files, intermediate filtered FASTQs, per-feature intermediate tables) and spatial plot PNGs from `pipeline_output/{sample}/qc/plots/`. Run after QC (and allc if needed) are complete:

```bash
snakemake --profile profiles/local_HPC \
    --snakefile main/Snakefile \
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


