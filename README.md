# Spatial joint profiling of DNA methylome and transcriptome (Spatial-DMT)

This repository contains the data preprocessing and QC pipeline, plus downstream analysis and visualization code, for [Lee, Fu et al., *Nature* 2025](https://www.nature.com/articles/s41586-025-09478-x). The Snakemake workflow (`main/Snakefile`) processes raw sequencing data through these main targets:

To access the processed `.cg` files: https://huggingface.co/datasets/HBBWS/Spatial-DMT

- **`final_output`** (default) — adapter+structural-prefix trim with CB-tagged read names, single-pass bisulfite alignment with cell-aware dedup (dupsifter `-B`), per-cell BAM split, CpG methylation pileup, lambda spike-in alignment, RNA alignment and gene matrix generation, QC reports, and collection of final deliverables under `final_output/`
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
    biscuit=1.8 \
    samtools=1.22.1 \
    htslib=1.22.1 \
    bedtools=2.31.1 \
    dupsifter=1.3.0 \
    yame=1.9 \
    pigz=2.8 \
    parallel=20260422 \
    multiqc \
    r-base \
    r-seurat \
    biopython \
    fuzzysearch \
    regex \
    pysam \
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
| [BISCUIT](https://huishenlab.github.io/biscuit/) | 1.8 (minimum) | bioconda |
| [samtools](https://www.htslib.org/) | 1.22.1 | bioconda |
| bgzip / tabix (htslib) | 1.22.1 | bioconda |
| [bedtools](https://bedtools.readthedocs.io/) | 2.31.1 | bioconda |
| [dupsifter](https://github.com/huishenlab/dupsifter) | 1.3.0 | bioconda |
| [yame](https://github.com/zhou-lab/yame) | 1.9 | bioconda |
| [pigz](https://zlib.net/pigz/) | 2.8 | conda-forge |
| [GNU parallel](https://www.gnu.org/software/parallel/) | 20260422 | conda-forge |
| [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (bbduk) | 39.81b | Homebrew / direct |
| [MultiQC](https://multiqc.info/) | 1.34 | conda-forge |
| [R](https://www.r-project.org/) | 4.x | conda-forge |
| [Seurat](https://satijalab.org/seurat/) | 5.x | conda-forge |
| [Python](https://www.python.org/) | 3.x | conda-forge |
| [Biopython](https://biopython.org/) | — | conda-forge |
| [fuzzysearch](https://github.com/taleinat/fuzzysearch) | — | conda-forge (used by `spatialmeth_trimtag.py` for R1 adapter scan) |
| [regex](https://pypi.org/project/regex/) | — | conda-forge (smcseq fuzzy linker patterns) |
| [pysam](https://github.com/pysam-developers/pysam) | — | conda-forge (used by `split_bam_by_cb.py` and `barcode_summary.py`) |
| matplotlib | — | conda-forge |
| pandas | — | conda-forge |

### 1. Data preprocessing

Next Generation Sequencing (NGS) was performed using an Illumina NovaSeq 6000 sequencer (150bp paired-end mode). Read 1 contains the genome sequences, and Read 2 contains the spatial Barcode A & B and UMIs (mRNA).

The preprocessing pipeline processes spatial DNA methylation and (optionally) spatial RNA data in a single Snakemake workflow (`main/Snakefile`).

#### Input data structure

The `barcodes/` directory is provided with this repository. Users only need to supply the raw FASTQ files:

```
fastq/
  {sample}/
    DNA_1.fq.gz   # DNA methylation read 1 (genomic sequence)
    DNA_2.fq.gz   # DNA methylation read 2 (spatial barcode + adapter)
    RNA_1.fq.gz   # RNA read 1 (cDNA sequence)  [optional]
    RNA_2.fq.gz   # RNA read 2 (spatial barcode + UMI)  [optional]
```

Samples are auto-detected from `fastq/*/DNA_1.fq*`, or specified explicitly with `--config IDS=sample1,sample2`. RNA processing is automatically skipped for samples where `RNA_1.fq*` is absent; the QC report will contain DNA-only metrics for those samples.

#### Running the pipeline

All invocations use **absolute** paths for `--profile` and `--snakefile` so the command works from any cwd (your scratch workdir, your home, anywhere). `site.yaml` is loaded automatically via each profile's `configfile:` directive (scalar form → profile-relative path resolution).

**Local / HPC (no scheduler):**
```bash
REPO=/absolute/path/to/Spatial-DMT-2024
snakemake --profile "$REPO/profiles/local_HPC" \
    --snakefile "$REPO/main/Snakefile" \
    --config ref=hg38
```

**SLURM cluster:**
```bash
REPO=/absolute/path/to/Spatial-DMT-2024
snakemake --profile "$REPO/profiles/slurm" \
    --snakefile "$REPO/main/Snakefile" \
    --config ref=hg38
```

The default target is `final_output`. To run it explicitly, append the target after `--`:

```bash
REPO=/absolute/path/to/Spatial-DMT-2024
snakemake --profile "$REPO/profiles/slurm" \
    --snakefile "$REPO/main/Snakefile" \
    --config ref=hg38 \
    -- final_output
```

The examples above run the default `spatialdmt` protocol on `hg38`. For the **`smcseq`** protocol (96×96 mouse embryo, mm10), add `protocol=smcseq` — that single flag auto-selects the 22-mer whitelist (`barcodes/smcseq_barcodes_22mer.txt`) and derives the 96×96 grid; nothing else changes:

```bash
REPO=/absolute/path/to/Spatial-DMT-2024
snakemake --profile "$REPO/profiles/slurm" \
    --snakefile "$REPO/main/Snakefile" \
    --config ref=mm10 protocol=smcseq \
        FASTQ_DIR=/scr1/users/you/run/fastq \
        WORKDIR=/scr1/users/you/run/pipeline_output \
        OUTDIR=/results/run
```

Supported values for `ref`: `hg38`, `mm10` (must match a key defined in your `profiles/*/site.yaml`). The published `smcseq` mouse-embryo dataset uses `ref=mm10`.

Required profile config keys (set in `profiles/local_HPC/config.yaml` or `profiles/slurm/config.yaml`):

| Key | Description |
|-----|-------------|
| `ref` | Genome reference key (`hg38` or `mm10`) |
| `REF_DIR` | Root directory of reference files |
| `BBDUK` | Path to `bbduk.sh` executable |
| `PYTHON` | Python interpreter with biopython, fuzzysearch, matplotlib, and pandas installed |

Optional path config keys (override via `--config KEY=VALUE` or in profile config.yaml). Defaults are cwd-relative, so the historical "run snakemake from the workdir" pattern still works unchanged.

| Key | Default | Description |
|-----|---------|-------------|
| `FASTQ_DIR` | `fastq` | Root containing `{sample}/{DNA,RNA}_{1,2}.fq.gz`. Accepts absolute paths. |
| `WORKDIR` | `pipeline_output` | Intermediate outputs (merged tagged FASTQs, merged BAM, per-cell BAMs, VCFs, qc/). Accepts absolute paths. |
| `OUTDIR` | `final_output` | Per-sample deliverable bundle (`.cg`, qc reports, gene matrix). Accepts absolute paths. |
| `protocol` | `spatialdmt` | R2 chemistry: `spatialdmt` (50×50 grid, 8+8 bp barcodes) or `smcseq` (96×96 grid, 11+11 bp barcodes). See [R2 chemistry](#r2-chemistry-protocol-switch) for grammar details. |
| `DUPSIFTER` | `~/software/dupsifter/v1.3.0/bin/dupsifter` | Path to a dupsifter ≥ 1.3.0 binary (the `-B` barcode-aware dedup flag was added in 1.3.0). |

Example with absolute paths (keeps intermediates and deliverables on separate disks):
```bash
snakemake --profile "$REPO/profiles/slurm" \
    --snakefile "$REPO/main/Snakefile" \
    --config ref=mm10 \
        FASTQ_DIR=/scratch/run42/fastq \
        WORKDIR=/scratch/run42/pipeline_output \
        OUTDIR=/results/run42
```

#### R2 chemistry (protocol switch)

R2 carries the spatial barcode flanked by linker sequences. Two grammars are supported via the `protocol` config key.

**`spatialdmt`** (default, 50×50 grid):

```
R2 = [8 BC1][30 linker1][8 BC2][30 linker2][19 AGATGTGTATAAGAGATAG][genomic]
       │       │           │       │           │
       │       │           │       │           └── Tn5 mosaic-end + flank        (fuzzy ≤1)
       │       │           │       └── ATTTATGTGTTTGAGAGGTTAGAGTATTTG            (fuzzy ≤2)
       │       │           └── 8 bp immediately before linker2
       │       └── GTGGTTGATGTTTTGTATTGGTGTATGATT                                 (fuzzy ≤2)
       └── 8 bp immediately before linker1; concat with BC2 = 16-mer whitelist key

Barcode key     : 16 bp (8 + 8)
Whitelist       : barcodes/spatial_barcodes.txt  (2500 rows = 50 × 50, col 4)
Total R2 prefix : 95 bp before genomic
```

**`smcseq`** (96×96 grid):

```
R2 = [19 GGTGTAGTGGGTTTGGAGG][11 BC1][30 linker1][11 BC2][30 linker2][8 AGATGTGT][1 UMI][genomic]
       │                       │       │           │       │           │           │
       │                       │       │           │       │           │           └── discarded
       │                       │       │           │       │           └── Tn5 ME anchor   (fuzzy ≤2)
       │                       │       │           │       └── linker2 (W): CCC.....C..CC....C...C...CC...
       │                       │       │           │                  (C): TTT.....T..TT....T...T...TT...   (fuzzy ≤4)
       │                       │       │           └── 11 bp immediately before linker2
       │                       │       └── linker1 (W): ..CC.C...C.....C.C.C..C...C...
       │                       │                  (C): ..TT.T...T.....T.T.T..T...T...   (fuzzy ≤4)
       │                       └── 11 bp immediately before linker1; concat with BC2 = 22-mer key
       └── primer anchored at pos 0                                                       (fuzzy ≤3)

Barcode key     : 22 bp (11 + 11)
Whitelist       : barcodes/smcseq_barcodes_22mer.txt  (9216 rows = 96 × 96, col 1)
Total R2 prefix : 110 bp before genomic
Strand handling : both linker variants (W = C-preserved, C = C→T BS-converted) are tried
                  per read; BISCUIT auto-handles strand at alignment, no Watson/Crick split
```

Switch with `--config protocol=smcseq` — see the [smcseq run example](#running-the-pipeline) above for a full invocation.

#### Pipeline steps

**Architecture note**: demux-after-align with a separate unmatched track.

Matched reads (BC1+BC2 hit the whitelist, exact or Hamming-1) get an 8-char ACGT CB tag encoding their (X, Y) coord (`tools/cb_codec.py`; the ACGT encoding is required because `dupsifter -B` silently collapses any non-{ACGTN} CB tag to a single bucket), and flow through one merged FASTQ → one biscuit alignment → cell-aware `dupsifter -B` → per-cell BAM split (`tools/split_bam_by_cb.py`) → per-cell pileup.

Unmatched reads — structure-fail (`_SF` suffix in read name) and whitelist-miss (`_WM`) — are split off at the trim step into a parallel `Unmatched_R{1,2}.fq.gz` pair, aligned independently into `bam_unmatched/`, and summarised by `mqc_unmatched_diagnostics`. They never enter the merged BAM, so per-cell stats stay clean.

This replaces the original per-cell-FASTQ-then-per-cell-BISCUIT path which spawned ~9216 BISCUIT invocations per sample.

**DNA methylation branch (matched track):**

1. `dna_trim_and_tag` — Adapter trim (R1) + structural-prefix parse (R2: primer + BC1 + linker + BC2 + linker + Tn5 ME [+ UMI for smcseq]). Look up BC1+BC2 in the whitelist (Hamming-1 tolerated), encode the matched (X, Y) coord as an 8-char ACGT string, rewrite each read name as `<origid>_<CB>_<UMI>`. Emits one merged R1/R2 tagged FASTQ pair for matched reads, plus a separate `Unmatched_R{1,2}.fq.gz` pair for structure-fail + whitelist-miss reads (`tools/spatialmeth_trimtag.py`; chunked-parallel via GNU parallel, pigz output).
2. `dna_biscuit_align` — One `biscuit align -9` over the matched-track FASTQ pair → `dupsifter -B` (cell-aware, BS-correct dedup via the CB tag) → `samtools sort` → per-cell BAM split keyed on CB. Outputs `bam_merged/{sample}_merged_dedup.bam` and `bam/<XXYY>.bam` per cell — every CB is a real whitelist coord, no sentinel cell.
3. `dna_biscuit_pileup` — CpG methylation pileup per per-cell BAM (parallel via GNU parallel), VCF→BED conversion, strand merging, CpG BED intersection, packing into a yame `.cg` file indexed by cell coord. Also writes per-cell `tmp/pileup/{XXYY}.conv.tsv` (CA/CC/CG/CT β-sum + count) and pools them into the sample-level `qc/biscuit_qc/{sample}_totalBaseConversionRate.txt` — one pileup per cell, no second sample-level pileup anywhere.
4. `dna_biscuit_align_lambda` — Align reads to the lambda genome for bisulfite conversion efficiency estimation. **Inputs differ by protocol**: `spatialdmt` uses the dedicated `trim/{sample}_Lambda_R{1,2}.fq.gz` files emitted by the trim step (barcode-85 reads only); `smcseq` aligns the FULL matched-track FASTQ pair to lambda since there's no dedicated lambda barcode in that protocol.
5. `dna_biscuit_pileup_lambda` — Pileup on lambda BAM; produces allc.bed for conversion efficiency calculation, plus `{sample}_Lambda.cg`.

**DNA methylation branch (unmatched track):**

6. `dna_biscuit_align_unmatched` — `biscuit align` + `dupsifter` (no `-B`, no `-9`) on the Unmatched FASTQ pair → `samtools sort`. Single sample-level BAM at `bam_unmatched/{sample}_unmatched.bam`.
7. `mqc_unmatched_biscuit` — Runs `tools/spatialmeth_qc.sh` on the unmatched BAM as ONE pseudo-sample, surfacing it as its own MultiQC entry parallel to (not mixed with) the matched-track sample-level entry.
8. `mqc_unmatched_diagnostics` — Bucket unmatched reads by `_SF`/`_WM` suffix; for SF reads, re-run the R2 anchor regexes to identify which structural anchor (primer, linker1, linker2, Tn5 ME, ...) failed first. Emits a per-sample diagnostic TSV and a single MultiQC custom-content `_mqc.tsv` (`tools/unmatched_diagnostics.py`).

**RNA branch (runs in parallel with DNA; skipped if `RNA_1.fq*` is absent):**

9.  `rna_filter_primer` — bbduk: retain R2 reads containing the spatial primer sequence
10. `rna_filter_L1` — bbduk: retain reads containing linker 1 sequence
11. `rna_filter_L2` — bbduk: retain reads containing linker 2 sequence
12. `rna_fq_process` — Extract barcode (BC2+BC1, 16 bp) and UMI (10 bp) from fixed positions in R2 and reformat for STARsolo
13. `rna_star_solo` — STARsolo alignment and per-barcode gene expression quantification; outputs raw and filtered count matrices under `Solo.out/GeneFull/`

#### Outputs

```
pipeline_output/{sample}/
  trim/              # matched tagged DNA FASTQs ({sample}_R{1,2}.fq.gz) +
                     #   Unmatched FASTQs ({sample}_Unmatched_R{1,2}.fq.gz, _SF/_WM tagged) +
                     #   Lambda_R{1,2} only for spatialdmt
  bam_merged/        # one big {sample}_merged_dedup.bam (matched track, post-align + dupsifter -B)
  bam/               # per-cell BAMs named <XXYY>.bam (from split_bam_by_cb.py); every CB is a
                     #   real whitelist coord, no sentinel bucket
  bam_unmatched/     # {sample}_unmatched.bam — sample-level alignment of the Unmatched fastqs
                     #   (no -9, no -B; position-only dedup since there's no real cell here)
  bam_lambda/        # lambda BAM with flagstat and dupsifter stats
  pileup/            # {sample}.cg + {sample}.cg.idx (yame packed CpG methylation indexed by cell)
                     #   and {sample}_Lambda.cg
  align_list.txt     # list of per-cell coord tags (XXYY) that have BAMs; consumed by pileup + QC
  tmp/               # all intermediates removed by clean: per-chunk trim outputs, pileup VCFs,
                     #   per-cell .cg + .conv.tsv shards, filtered RNA FASTQs, lambda pileup intermediates
  rna_bbduk/         # bbduk filter stats per step (RNA samples only)
  rna_STAR/          # Aligned.out.sam, STAR logs, Solo.out/ (count matrix; RNA samples only)
```

The `final_output` target copies the primary deliverables into one sample-organized directory:

```
final_output/{sample}/
  {sample}.cg              # genomic CpG methylation matrix
  {sample}.cg.idx          # yame index for {sample}.cg
  {sample}_Lambda.cg       # lambda spike-in CpG methylation matrix
  qc_report.html           # self-contained QC report
  multiqc_report.html      # MultiQC report
  multiqc_report_data/     # MultiQC parsed data and assets
  Solo.out/                # STARsolo count matrices (raw/ and filtered/), RNA samples only
```

### 2. Quality control

Run after preprocessing. Produces BISCUIT QC tables for every barcode, a MultiQC report aggregating BISCUIT and (if available) STAR and bbduk statistics, spatial heatmaps of per-barcode metrics, and a self-contained HTML report.

```bash
REPO=/absolute/path/to/Spatial-DMT-2024
snakemake --profile "$REPO/profiles/local_HPC" \
    --snakefile "$REPO/main/Snakefile" \
    --config ref=hg38 \
    -- qc
```

QC rules are named by their **consumer**: `mqc_*` produce inputs to MultiQC (either standard BISCUIT tables or custom-content `_mqc.tsv`); `spatial_*` produce per-cell tables and heatmaps that feed the spatial-map HTML; `qc_report` is the final aggregator.

**MultiQC inputs (matched track):**

- `mqc_biscuit` — One `tools/spatialmeth_qc.sh` invocation on the merged dedup BAM (sample-level — not per-cell). Emits `biscuit qc` tables (mapq, strand, dup, CpH/CpG retention by read pos). Does NOT pileup — the conversion-rate table is produced by `dna_biscuit_pileup` from the per-cell VCFs (single pileup per cell, never twice). covdist/cv intentionally NOT produced — see `spatialmeth_qc.sh` header note on contig mismatch.
- `mqc_curated` — Sample-level curated MultiQC custom-content table (13 cols: Sample, reads/mapped/mapping rate from flagstat, dup rate from dupsifter `-B` stats, lambda mapped + retention, CpG coverage/depth from yame summary on `{sample}.cg`, and CA/CC/CG/CT retention from the pooled per-cell pileup VCFs — CG is the genuine CpG methylation signal, CA/CC/CT are residual non-CpG retention ≈ 0 for complete bisulfite conversion). See `tools/write_curated_metrics_mqc.py`. Mirrors the Ultima WGBS `final_report` curated-metrics pattern.

**MultiQC inputs (unmatched track):**

- `mqc_unmatched_biscuit` — `tools/spatialmeth_qc.sh` on the unmatched BAM as ONE pseudo-sample → its own MultiQC entry.
- `mqc_unmatched_diagnostics` — SF/WM bucket breakdown + R2 anchor-fail analysis; emits MultiQC custom-content `_mqc.tsv` (`tools/unmatched_diagnostics.py`).

**Spatial map (per-cell):**

- `spatial_summary` — Single pysam pass over `bam_merged/{sample}_merged_dedup.bam` to aggregate per-CB (total, mapped, dup) counts, decode CB → (X, Y) via `cb_codec`, emit per-cell TSV (8 cols incl. `mapping_rate` and `dup_rate` — NA when undefined) and four spatial heatmaps + barcode-rank plot. Grid auto-derived from whitelist.
- `spatial_features` — Per-cell mean methylation summarized over genomic features (ChromHMM, windows, chromosomes).

**Aggregation:**

- `qc_report` — MultiQC report aggregating BISCUIT-QC tables (matched + unmatched), STAR logs, bbduk stats, and the two custom-content `_mqc.tsv` (curated + unmatched_summary), plus a self-contained HTML report with embedded spatial plots.
- `qc` — umbrella target; transitively pulls all rules above via `qc_report` + `spatial_features` deps.

#### QC outputs

```
pipeline_output/{sample}/qc/
  biscuit_qc/                # sample-level BISCUIT-QC tables on merged BAM (mapq, strand, dup,
                             #   CpH/CpG retention by read pos) +
                             #   {sample}_totalBaseConversionRate.txt (CA/CC/CG/CT, from dna_biscuit_pileup)
  biscuit_qc_merged/         # MultiQC custom-content TSVs:
                             #   {sample}_curated_metrics_mqc.tsv
                             #   {sample}_unmatched_summary_mqc.tsv
  biscuit_qc_unmatched/      # BISCUIT-QC tables on the unmatched BAM (pseudo-sample)
  unmatched_diagnostics.tsv  # long-format SF/WM breakdown + anchor-fail counts
  features/                  # per-feature mean methylation tables (.txt.gz per feature set)
  table/
    {sample}_spatial_barcode_counts.tsv  # per-cell row (cb, x, y, dna_reads, mapped_reads,
                                         #   dup_reads, mapping_rate, dup_rate) — NA when undefined
  plots/
    spatial_dna_reads.png                # read-count heatmap (linear)
    spatial_log_dna_reads.png            # read-count heatmap (log10)
    spatial_mapping_rate.png             # per-cell mapping rate
    spatial_dup_rate.png                 # per-cell duplication rate
    barcode_rank.png                     # barcode rank vs read count
  html/
    multiqc_report.html                  # MultiQC report (BISCUIT + STAR + bbduk + custom-content)
    qc_report.html                       # self-contained HTML with embedded spatial maps
```

### 3. AllC (non-CpG methylation)

Optional. Run after preprocessing and **before** clean (requires tmp/ pileup VCFs). Produces all-cytosine pileup and per-feature non-CpG methylation summaries.

```bash
REPO=/absolute/path/to/Spatial-DMT-2024
snakemake --profile "$REPO/profiles/local_HPC" \
    --snakefile "$REPO/main/Snakefile" \
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
REPO=/absolute/path/to/Spatial-DMT-2024
snakemake --profile "$REPO/profiles/local_HPC" \
    --snakefile "$REPO/main/Snakefile" \
    --config ref=hg38 \
    -- clean
```

#### - For Image processing:
This python script processes a phase contrast tissue image to generate spatial metadata similar to 10x Visium outputs. It first converts the image to grayscale and applies adaptive thresholding to create a binary mask that distinguishes tissue from background. The image is then divided into a grid matching the protocol (50×50 for `spatialdmt`, 96×96 for `smcseq`), and each grid cell is evaluated to determine whether it contains tissue based on a pixel intensity threshold. Detected tissue spots are matched with predefined spatial barcodes, and the results are compiled into output files, including a tissue positions CSV, a scale factors JSON, and a visualization image with highlighted tissue regions.

### 2. Downstream data analysis and visualization 

The data analysis and visualization were performed using R(4.4.0).

`QC_check.Rmd`: Contains all the code to generate QC plots

`WNN_DNAmRNA_clustering`: Contains all the code to integrate two modalities into Seurat objects, identify differential marker genes and VMRs, and visualize them on spatial maps

`Integration_E11E13.Rmd`: Contains all the code to integrate between day 11 and day 13 mouse embryos’ data.
 
`Utility_function.R`: Contains all functions used including mapping VMR to overlapping genes and iterative PCA. 
