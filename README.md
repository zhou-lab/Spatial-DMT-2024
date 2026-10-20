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

Create the tested environment from the checked-in specification. Micromamba,
Mamba, and Conda all accept this file; Micromamba is shown because it does not
require a base Conda installation.

```bash
micromamba create -f environment.yml
micromamba activate spatial-dmt
export PYTHONNOUSERSITE=1
```

The environment includes BBTools (`bbduk.sh`), the Slurm executor plugin, a C
compiler, and HTSlib development metadata. The workflow builds its modified
vendored dupsifter, BAM splitter, and CpG converter automatically against the
HTSlib in the active environment. Do not install the Bioconda `dupsifter`
package; its current dependency constraints conflict with BISCUIT 1.8 and the
pipeline does not use it.

`PYTHONNOUSERSITE=1` prevents packages in `~/.local` from silently overriding
the versions in the release environment. The Snakefile also exports this for
all rule subprocesses, but it should be set before launching Snakemake so the
controller itself is isolated.

The checked-in Slurm profile submits rule jobs with `--export=ALL`. This is
required when a controller was itself submitted with a restrictive Slurm
export policy; without it, compute nodes can lose the active environment's
`PATH` and fail to locate specialized tools such as YAME.

| Tool | Version tested | Source |
|------|---------------|--------|
| [Snakemake](https://snakemake.readthedocs.io/) | 9.19.0 | conda-forge |
| [snakemake-executor-plugin-cluster-generic](https://github.com/snakemake/snakemake-executor-plugin-cluster-generic) | 1.0.9 | bioconda (SLURM only) |
| [STAR](https://github.com/alexdobin/STAR) / STARsolo | 2.7.11b | bioconda |
| [BISCUIT](https://huishenlab.github.io/biscuit/) | 1.8 (minimum) | bioconda |
| [samtools](https://www.htslib.org/) | 1.23.1 | bioconda |
| bgzip / tabix (htslib) | 1.23.1 | bioconda |
| [bedtools](https://bedtools.readthedocs.io/) | 2.31.1 | bioconda |
| [dupsifter](https://github.com/huishenlab/dupsifter) | 1.3.x-derived vendored build | built by the workflow against HTSlib |
| [yame](https://github.com/zhou-lab/yame) | 1.10 | bioconda |
| [pigz](https://zlib.net/pigz/) | 2.8 | conda-forge |
| [GNU parallel](https://www.gnu.org/software/parallel/) | 20260422 | conda-forge |
| [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (bbduk) | 39.81 | bioconda |
| [MultiQC](https://multiqc.info/) | 1.35 | bioconda |
| [Python](https://www.python.org/) | 3.12 | conda-forge |
| [Biopython](https://biopython.org/) | — | conda-forge |
| [fuzzysearch](https://github.com/taleinat/fuzzysearch) | — | conda-forge (used by `spatialmeth_trimtag.py` for R1 adapter scan) |
| [regex](https://pypi.org/project/regex/) | — | conda-forge (smcseq fuzzy linker patterns) |
| [pysam](https://github.com/pysam-developers/pysam) | — | conda-forge (used by `barcode_summary.py`; the Python BAM splitter remains as a benchmark reference) |
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
python -s -m snakemake --profile "$REPO/profiles/local_HPC" \
    --snakefile "$REPO/main/Snakefile" \
    --config ref=hg38
```

**SLURM cluster:**
```bash
REPO=/absolute/path/to/Spatial-DMT-2024
capacity=$("$REPO/tools/slurm_dna_capacity.sh")
capacity_value() {
    awk -F '\t' -v key="$1" '$1 == key { print $2 }' <<< "$capacity"
}
rule_cpus=$(capacity_value paired_rule_cpus_now)
controller_cpus=$(capacity_value controller_cpus_now)
controller_mem=$(capacity_value controller_memory_mb)
python -s -m snakemake --profile "$REPO/profiles/slurm" \
    --snakefile "$REPO/main/Snakefile" \
    --cores "$controller_cpus" --jobs 4 \
    --resources cpus_per_task="$controller_cpus" mem_mb="$controller_mem" \
    --config ref=hg38 \
        DNA_TRIM_CPUS="$rule_cpus" \
        DNA_ALIGN_CPUS="$rule_cpus" \
        DNA_PILEUP_CPUS="$rule_cpus"
```

With an otherwise empty account under a 144-CPU cap, the helper returns a
138-CPU controller and 69 CPUs per heavy rule, leaving capacity for independent
DNA and RNA branches or multiple samples to progress together. Heavy DNA rules mirror their allocation in
both native Snakemake `threads` and the profile's `cpus_per_task` resource;
set both controller limits to the same budget. The native thread ledger stops
over-submission, while `cpus_per_task` supplies the Slurm allocation.
On the lab cluster, run `tools/slurm_dna_capacity.sh` immediately before
submission because the `defq` per-user limit changes dynamically. The defaults
fit two heavy rules plus six spare CPUs under the current 144-CPU limit. Set
the two controller values and `DNA_{TRIM,ALIGN,PILEUP}_CPUS` to the helper's
reported budgets if the live limit or existing allocation differs. When
optional BAM-backed unmatched diagnostics are enabled, also set
`UNMATCHED_ALIGN_CPUS` to that rule budget. Set the global `mem_mb` resource to its reported
`controller_memory_mb` value so concurrent samples also respect the QOS memory
cap after current job allocations; the helper reserves 40 GiB by default. For
an unmatched-only restart, use `single_rule_cpus_now` instead of the paired
budget.

The default target is `final_output`. To run it explicitly, append the target after `--`:

```bash
REPO=/absolute/path/to/Spatial-DMT-2024
python -s -m snakemake --profile "$REPO/profiles/slurm" \
    --snakefile "$REPO/main/Snakefile" \
    --config ref=hg38 \
    -- final_output
```

The examples above run the default `spatialdmt` protocol on `hg38`. For the **`smcseq`** protocol (96×96 mouse embryo, mm10), add `protocol=smcseq`. This selects the 22-mer DNA whitelist and the separate 16-mer spatial-RNA whitelist, and derives a 96×96 grid for both modalities:

```bash
REPO=/absolute/path/to/Spatial-DMT-2024
python -s -m snakemake --profile "$REPO/profiles/slurm" \
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
| `REF_DIR` | Root directory of reference files; see [reference preparation](docs/reference_setup.md) |

Optional path config keys (override via `--config KEY=VALUE` or in profile config.yaml). Defaults are cwd-relative, so the historical "run snakemake from the workdir" pattern still works unchanged.

| Key | Default | Description |
|-----|---------|-------------|
| `FASTQ_DIR` | `fastq` | Root containing `{sample}/{DNA,RNA}_{1,2}.fq.gz`. Accepts absolute paths. |
| `WORKDIR` | `pipeline_output` | Intermediate outputs (merged tagged FASTQs, merged BAM, per-cell BAMs, VCFs, qc/). Accepts absolute paths. |
| `OUTDIR` | `final_output` | Per-sample deliverable bundle (`.cg`, qc reports, gene matrix). Accepts absolute paths. |
| `protocol` | `spatialdmt` | R2 chemistry: `spatialdmt` (50×50 grid, 8+8 bp barcodes) or `smcseq` (96×96 grid, 11+11 bp barcodes). See [R2 chemistry](#r2-chemistry-protocol-switch) for grammar details. |
| `PYTHON` | `python` | Python executable from the active environment. |
| `BBDUK` | `bbduk.sh` | BBTools adapter/filter executable. |
| `YAME` | `yame` | YAME executable. |
| `STAR_BIN` | `STAR` | STAR/STARsolo executable. No environment module is loaded implicitly. |
| `DUPSIFTER` | `tools/dupsifter_speed/dupsifter` | Path to dupsifter. The default optimized vendored build is created automatically and retains the 1.3.x `-B` behavior. |
| `CB_SPLITTER` | `tools/split_bam_by_cb_hts/split_bam_by_cb_hts` | Compiled CB-tag BAM splitter; built automatically. |
| `MERGECG_CONTEXT` | `tools/mergecg_context/mergecg_context` | Context-aware single-cell CpG strand merger; built automatically and avoids repeated reference loads. |
| `DNA_TRIM_CPUS` | `69` | CPUs allocated to streaming trim/barcode assignment. |
| `DNA_TRIM_OVERHEAD_CPUS` | `4` | Allocated trim CPUs reserved for paired decompression and stream dispatch. |
| `DNA_TRIM_WORKERS` | `DNA_TRIM_CPUS - DNA_TRIM_OVERHEAD_CPUS` | Python trim workers; override only after measuring total job CPU use. |
| `DNA_TRIM_BLOCK_SIZE` | `auto` | Newline-aligned stream block size. Auto selects 64, 256, or 512 MiB from input size. |
| `DNA_TRIM_COMPRESSION_LEVEL` | `1` | Gzip level for trim FASTQs. This time-first default was 5% faster than level 3 with 14% larger output in the fixed deep-concurrency benchmark. |
| `DNA_TRIM_COMPRESSION_BACKEND` | `gzip` | Per-worker compressor (`gzip` or `pigz`). The default avoids spawning several compressor processes per trim worker. |
| `DNA_TRIM_MEM` | `160G` | Memory for deep streaming trim at the default worker count. |
| `DNA_ALIGN_CPUS` | `69` | CPUs allocated to matched BISCUIT alignment. |
| `DNA_ALIGN_MEM` | `400G` | Memory for deep matched alignment, deduplication, and sort. |
| `DNA_ALIGN_SORT_CPUS` | `16` | `samtools sort` threads inside the matched allocation. |
| `DNA_PILEUP_CPUS` | `DNA_ALIGN_CPUS` | CPUs allocated across independent per-cell BISCUIT pileups. |
| `DNA_PILEUP_MEM` | `80G` | Memory for per-cell pileup and conversion processing. |
| `PILEUP_THREADS_PER_CELL` | `1` | BISCUIT threads per cell; one-thread cells maximize aggregate throughput and load balancing. |
| `PILEUP_WINDOW_STEP` | `1000000` | Bases per BISCUIT pileup dispatch window; 1 Mb was fastest without the memory growth of larger windows. |
| `LAMBDA_ALIGN_CPUS` | `48` | CPUs allocated to lambda BISCUIT alignment. |
| `LAMBDA_ALIGN_SORT_CPUS` | `24` | `samtools sort` threads inside the lambda allocation. |
| `LAMBDA_PILEUP_MAX_MAPPED` | `2000000` | Deterministic cap on mapped lambda records supplied to pileup, preventing numerical overflow at extreme spike-in depth; full-BAM alignment metrics are retained and sampling is recorded in `qc/table/`. |
| `UNMATCHED_ALIGN_CPUS` | `69` | CPUs allocated to unmatched BISCUIT alignment. |
| `UNMATCHED_ALIGN_MEM` | `400G` | Memory for deep unmatched alignment, deduplication, and sort. |
| `UNMATCHED_SORT_CPUS` | `16` | `samtools sort` threads inside the unmatched allocation. |
| `UNMATCHED_DIAGNOSTICS_WITH_BAM` | `false` | Add unmatched-read genomic mapping and duplication fields by aligning the unmatched FASTQs. FASTQ-level failure decomposition is always generated. |
| `SORT_LOCAL_MIN_GB` | `250` | Minimum free local scratch before sort can use it. Input-aware checks may require more. |
| `SORT_LOCAL_INPUT_FACTOR` | `3` | Estimated BAM/sort-scratch bytes per compressed alignment-input byte. |
| `SORT_LOCAL_HEADROOM_GB` | `20` | Free space retained on the sort-scratch filesystem. |
| `STORAGE_HEADROOM_GB` | `20` | Free space retained by trim, alignment, and bucket-split preflight checks. |

Example with absolute paths (keeps intermediates and deliverables on separate disks):
```bash
python -s -m snakemake --profile "$REPO/profiles/slurm" \
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

Matched reads (BC1+BC2 hit the whitelist, exact or Hamming-1) get an 8-char ACGT CB tag encoding their (X, Y) coord (`tools/cb_codec.py`; the ACGT encoding is required because `dupsifter -B` silently collapses any non-{ACGTN} CB tag to a single bucket), and flow through one merged FASTQ → one biscuit alignment → cell-aware `dupsifter -B` → bounded bucketed per-cell BAM split (`tools/split_bam_by_cb_bucketed.sh`) → per-cell pileup.

Unmatched reads — structure-fail (`_SF` suffix in read name) and whitelist-miss (`_WM`) — are split off at the trim step into a parallel `Unmatched_R{1,2}.fq.gz` pair and summarized directly by `mqc_unmatched_diagnostics`. They never enter the merged BAM, so per-cell stats stay clean. Set `UNMATCHED_DIAGNOSTICS_WITH_BAM=true` only when genomic mapping and duplication fields are needed for this diagnostic bucket.

This replaces the original per-cell-FASTQ-then-per-cell-BISCUIT path which spawned ~9216 BISCUIT invocations per sample.

**DNA methylation branch (matched track):**

1. `dna_trim_and_tag` — Adapter trim (R1) + structural-prefix parse (R2: primer + BC1 + linker + BC2 + linker + Tn5 ME [+ UMI for smcseq]). Look up BC1+BC2 in the whitelist (Hamming-1 tolerated), encode the matched (X, Y) coord as an 8-char ACGT string, rewrite each read name as `<origid>_<CB>_<UMI>`. Emits one merged R1/R2 tagged FASTQ pair for matched reads, plus a separate `Unmatched_R{1,2}.fq.gz` pair for structure-fail + whitelist-miss reads (`tools/spatialmeth_trimtag.py`; paired blocks are dispatched by the bounded `tools/stream_parallel.py` coordinator, without raw temporary FASTQs).
2. `dna_biscuit_align` — One `biscuit align -9` over the matched-track FASTQ pair → `dupsifter -B` (cell-aware, BS-correct dedup via the CB tag) → `samtools sort` → per-cell BAM split keyed on CB. Outputs `bam_merged/{sample}_merged_dedup.bam` and `bam/<XXYY>.bam` per cell — every CB is a real whitelist coord, no sentinel cell.
3. `dna_biscuit_pileup` — CpG methylation pileup per per-cell BAM (parallel via GNU parallel), VCF→BED conversion, strand merging, CpG BED intersection, packing into a yame `.cg` file indexed by cell coord. Also writes per-cell `tmp/pileup/{XXYY}.conv.tsv` (CA/CC/CG/CT β-sum + count) and pools them into the sample-level `qc/biscuit_qc/{sample}_totalBaseConversionRate.txt` — one pileup per cell, no second sample-level pileup anywhere.
4. `dna_biscuit_align_lambda` — Align reads to the lambda genome for bisulfite conversion efficiency estimation. **Inputs differ by protocol**: `spatialdmt` uses the dedicated `trim/{sample}_Lambda_R{1,2}.fq.gz` files emitted by the trim step (barcode-85 reads only); `smcseq` aligns the FULL matched-track FASTQ pair to lambda since there's no dedicated lambda barcode in that protocol.
5. `dna_biscuit_pileup_lambda` — Pileup on lambda BAM; produces allc.bed for conversion efficiency calculation, plus `{sample}_Lambda.cg`.

**DNA methylation branch (unmatched track):**

6. `mqc_unmatched_diagnostics` — Bucket unmatched reads by `_SF`/`_WM` suffix; for SF reads, re-run the R2 anchor regexes to identify which structural anchor (primer, linker1, linker2, Tn5 ME, ...) failed first. Emits a per-sample diagnostic TSV and a single MultiQC custom-content `_mqc.tsv` (`tools/unmatched_diagnostics.py`). This FASTQ-level analysis is part of the default `final_output` target.
7. `dna_biscuit_align_unmatched` — Optional when `UNMATCHED_DIAGNOSTICS_WITH_BAM=true`: `biscuit align` + `dupsifter` (no `-B`, no `-9`) on the unmatched FASTQ pair → `samtools sort`, producing `bam_unmatched/{sample}_unmatched.bam` and adding genomic mapping and duplication fields to the diagnostic table.
8. `mqc_unmatched_biscuit` — Optional explicit target that runs `tools/spatialmeth_qc.sh` on the unmatched BAM as one pseudo-sample. It is not required by `final_output`.

**RNA branch (runs in parallel with DNA; skipped if `RNA_1.fq*` is absent):**

9.  `rna_fq_process` — Protocol-specific RNA structure parsing and STARsolo read preparation. `spatialdmt` retains the existing sequential bbduk filters and fixed-position extraction. `smcseq` follows LiuLab's released spatial-RNA regex: `[1 bp][14 bp motif][8 bp BC1][30 bp motif][8 bp BC2][49 bp motif][10 bp UMI]`, and uses the standard ordered 96-entry SPLiT-seq barcode list to generate a 96×96 RNA whitelist.
10. `rna_star_solo` — STARsolo alignment and per-barcode gene expression quantification; outputs raw and filtered count matrices under `Solo.out/GeneFull/`

The `smcseq` RNA grammar is ported from LiuLab's
[`spatial_RNA_analysis/Extract_barcode.sh`](https://github.com/LiuLab888/spatial-5mC-project/blob/main/spatial_RNA_analysis/Extract_barcode.sh)
and [`St_pipeline.sh`](https://github.com/LiuLab888/spatial-5mC-project/blob/main/spatial_RNA_analysis/St_pipeline.sh).
The ordered 96-entry component list is the
[SPLiT-seq Round1/Round2 barcode set](https://github.com/Teichlab/scg_lib_structs/tree/main/data/SPLiT-seq).

#### Outputs

```
pipeline_output/{sample}/
  trim/              # matched tagged DNA FASTQs ({sample}_R{1,2}.fq.gz) +
                     #   Unmatched FASTQs ({sample}_Unmatched_R{1,2}.fq.gz, _SF/_WM tagged) +
                     #   Lambda_R{1,2} only for spatialdmt
  bam_merged/        # one big {sample}_merged_dedup.bam (matched track, post-align + dupsifter -B)
  bam/               # per-cell BAMs named <XXYY>.bam plus compact CSI indexes (from the bucketed
                     #   C splitter); every CB is a real whitelist coord, no sentinel bucket
  bam_unmatched/     # optional {sample}_unmatched.bam when UNMATCHED_DIAGNOSTICS_WITH_BAM=true
                     #   (no -9, no -B; position-only dedup since there is no real cell here)
  bam_lambda/        # lambda BAM with flagstat and dupsifter stats
  pileup/            # {sample}.cg + {sample}.cg.idx (yame packed CpG methylation indexed by cell)
                     #   and {sample}_Lambda.cg
  align_list.txt     # list of per-cell coord tags (XXYY) that have BAMs; consumed by pileup + QC
  tmp/               # all intermediates removed by clean: per-block trim outputs, pileup VCFs,
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
  {sample}_lambda_pileup_sampling.tsv # lambda pileup depth and sampling provenance
  pipeline_version.tsv     # release, protocol, reference, and assignment label
  qc_report.html           # self-contained QC report
  multiqc_report.html      # MultiQC report
  multiqc_report_data/     # MultiQC parsed data and assets
  Solo.out/                # STARsolo count matrices (raw/ and filtered/), RNA samples only
```

The repository-level `VERSION` file is the canonical pipeline release label.
Every `final_output/{sample}/pipeline_version.tsv` records that label together
with the sample, protocol, and reference. Files created by Snakemake use the
`workflow_generated` assignment; compatible historical outputs may be labeled
`backfilled_compatible` after their workflow and output schema are audited.
Release tags and `VERSION` must match, and a released tag must not be moved.

### 2. Quality control

Run after preprocessing. Produces BISCUIT QC tables for every barcode, a MultiQC report aggregating BISCUIT and (if available) STAR and bbduk statistics, spatial heatmaps of per-barcode metrics, and a self-contained HTML report.

```bash
REPO=/absolute/path/to/Spatial-DMT-2024
python -s -m snakemake --profile "$REPO/profiles/local_HPC" \
    --snakefile "$REPO/main/Snakefile" \
    --config ref=hg38 \
    -- qc
```

QC rules are named by their **consumer**: `mqc_*` produce inputs to MultiQC (either standard BISCUIT tables or custom-content `_mqc.tsv`); `spatial_*` produce per-cell tables and heatmaps that feed the spatial-map HTML; `qc_report` is the final aggregator.

**MultiQC inputs (matched track):**

- `mqc_biscuit` — One `tools/spatialmeth_qc.sh` invocation on the merged dedup BAM (sample-level — not per-cell). Emits `biscuit qc` tables (mapq, strand, dup, CpH/CpG retention by read pos). Does NOT pileup — the conversion-rate table is produced by `dna_biscuit_pileup` from the per-cell VCFs (single pileup per cell, never twice). covdist/cv intentionally NOT produced — see `spatialmeth_qc.sh` header note on contig mismatch.
- `mqc_curated` — Sample-level curated MultiQC custom-content table (Sample, pre-filter reads/mapped/mapping rate from the merged BAM flagstat, QC-BAM reads and bsconv-retained fraction when applicable, dup rate from dupsifter `-B` stats, lambda mapped + retention, CpG coverage/depth from yame summary on `{sample}.cg`, and CA/CC/CG/CT retention from the pooled per-cell pileup VCFs — CG is the genuine CpG methylation signal, CA/CC/CT are residual non-CpG retention ≈ 0 for complete bisulfite conversion). See `tools/write_curated_metrics_mqc.py`. Mirrors the Ultima WGBS `final_report` curated-metrics pattern.

**MultiQC inputs (unmatched track):**

- `mqc_unmatched_diagnostics` — Default SF/WM bucket breakdown + R2 anchor-fail analysis from unmatched FASTQ; emits MultiQC custom-content `_mqc.tsv` (`tools/unmatched_diagnostics.py`). Mapping and duplication fields are included when `UNMATCHED_DIAGNOSTICS_WITH_BAM=true`.
- `mqc_unmatched_biscuit` — Optional explicit `tools/spatialmeth_qc.sh` target on the unmatched BAM as one pseudo-sample.

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
python -s -m snakemake --profile "$REPO/profiles/local_HPC" \
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
python -s -m snakemake --profile "$REPO/profiles/local_HPC" \
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
