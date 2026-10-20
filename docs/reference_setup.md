# Reference preparation

The pipeline does not download genomes automatically. Prepare one reference
root per site, point `REF_DIR` to it, and copy one of the checked-in
`profiles/*/site.yaml` files before editing paths. The checked-in Zhou Lab
profile is also a working example of the expected layout.

## Core DNA resources

For each genome, the default `final_output` target needs a FASTA, its FASTA
index, a BISCUIT index, and a coordinate-sorted CpG BED. The same resources are
needed for the Lambda control.

```bash
REF_DIR=/path/to/references
GENOME=mm10                         # or hg38
FASTA="$REF_DIR/$GENOME/$GENOME.fa"

samtools faidx "$FASTA"
mkdir -p "$REF_DIR/$GENOME/biscuit"
biscuit index -p "$REF_DIR/$GENOME/biscuit/$GENOME.fa" "$FASTA"

mkdir -p "$REF_DIR/$GENOME/annotation/cpg"
python /path/to/Spatial-DMT-2024/tools/build_cpg_bed.py "$FASTA" \
  | bgzip -c > "$REF_DIR/$GENOME/annotation/cpg/cpg.bed.gz"
tabix -f -p bed "$REF_DIR/$GENOME/annotation/cpg/cpg.bed.gz"
```

Run the same commands for the Lambda FASTA. Set `REF_FASTA`, `BISCUIT_INDEX`,
and `CPG_BED` under the corresponding genome key in `site.yaml`; paths are
relative to `REF_DIR/REF_BASE`.

## BISCUIT QC assets

The `final_output` target also runs BISCUIT alignment QC and requires the three
genome-matched files below in the directory selected by `BISCUIT_QCASSET`:

```text
cpg.bed.gz
windows100bp.gc_content.bot10p.bed.gz
windows100bp.gc_content.top10p.bed.gz
```

BISCUIT release assets are available for hg19, hg38, and mm10. For another
reference build, generate them from the exact FASTA used for alignment with
the `build_biscuit_QC_assets.pl` script installed with BISCUIT:

```bash
QC_DIR="$REF_DIR/$GENOME/biscuit/QC"
mkdir -p "$QC_DIR"
build_biscuit_QC_assets.pl --outdir "$QC_DIR" --ref "$FASTA"
```

See the [BISCUIT QC documentation](https://huishenlab.github.io/biscuit/docs/alignment/QC.html)
for release downloads and the asset-generation procedure. Do not reuse assets
from a different genome build or contig set.

## RNA resources

Samples with `RNA_1` and `RNA_2` inputs also require a STAR index. Build it
with the same FASTA naming convention used by the DNA branch and a matching
GTF annotation:

```bash
mkdir -p "$REF_DIR/$GENOME/STAR"
STAR --runMode genomeGenerate \
  --runThreadN 24 \
  --genomeDir "$REF_DIR/$GENOME/STAR" \
  --genomeFastaFiles "$FASTA" \
  --sjdbGTFfile /path/to/annotation.gtf \
  --sjdbOverhang 149
```

`--sjdbOverhang 149` matches the published 150-bp libraries. Use read length
minus one for a different library design.

## ChromHMM and feature masks

The HTML report includes TssA/Tss and Tx methylation and therefore requires
the `CHROMHMM_CM` YAME mask named in `site.yaml`. The explicit `qc` target also
requires the BED and YAME masks under `features/` and `KYCGKB_<genome>/` for
the feature sets listed in `main/Snakefile`. These are annotation-derived
resources, not executables or hidden software dependencies. Build them from
the same CpG universe as `CPG_BED`, or use the matching resources distributed
with the analysis site. A mismatched CpG universe must not be substituted.

Before a large run, use a dry run to identify every configured resource:

```bash
python -s -m snakemake --profile /path/to/profile \
  --snakefile /path/to/Spatial-DMT-2024/main/Snakefile \
  --config ref=mm10 IDS=sample_name --dry-run
```
