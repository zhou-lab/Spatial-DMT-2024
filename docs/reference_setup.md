# Reference preparation

This guide constructs every reference asset required by the standard
`final_output` target for `mm10` and `hg38`. Commands use pinned public files
and the two compact ChromHMM masks distributed in this repository. No Zhou Lab
filesystem is required.

The optional extended `qc` target needs additional user-selected feature
collections; those are described separately at the end.

## Prerequisites and paths

Activate the release environment, then set the repository and reference roots:

```bash
set -euo pipefail
REPO=/absolute/path/to/Spatial-DMT-2024
REF_DIR=/absolute/path/to/references
mkdir -p "$REF_DIR"
```

The commands require `curl`, `gzip`, `unzip`, `awk`, `md5sum`, `sha256sum`,
`samtools`, `biscuit`, `bgzip`, `tabix`, `bedtools`, `yame`, and `STAR`.

## 1. Genome FASTA

The validated profiles use UCSC hg38 and mm10 assemblies. Download and verify
the requested genome:

```bash
GENOME=mm10  # change to hg38 for the human reference

case "$GENOME" in
  mm10)
    FASTA_URL=https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
    FASTA_MD5=db005b65828db31735f384e4c5787be5
    ;;
  hg38)
    FASTA_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    FASTA_MD5=1c9dcaddfa41027f17cd8f7a82c7293b
    ;;
  *) echo "Unsupported standard genome: $GENOME" >&2; exit 2 ;;
esac

GENOME_DIR="$REF_DIR/$GENOME"
mkdir -p "$GENOME_DIR"
curl --fail --location "$FASTA_URL" -o "$GENOME_DIR/$GENOME.fa.gz"
printf '%s  %s\n' "$FASTA_MD5" "$GENOME_DIR/$GENOME.fa.gz" | md5sum --check -
gzip --decompress --stdout "$GENOME_DIR/$GENOME.fa.gz" > "$GENOME_DIR/$GENOME.fa"
samtools faidx "$GENOME_DIR/$GENOME.fa"
```

Run this section once with `GENOME=mm10` and once with `GENOME=hg38` when both
profiles are needed.

## 2. Lambda spike-in FASTA

The control reference is the NCBI RefSeq sequence `J02459.1`, assembly
`GCA_000840245.1`:

```bash
LAMBDA_DIR="$REF_DIR/Lambda"
LAMBDA_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/840/245/GCA_000840245.1_ViralProj14204/GCA_000840245.1_ViralProj14204_genomic.fna.gz
LAMBDA_MD5=5a7c51ae916edc3ebd2e32967b639bb3

mkdir -p "$LAMBDA_DIR"
curl --fail --location "$LAMBDA_URL" -o "$LAMBDA_DIR/J02459.fa.gz"
printf '%s  %s\n' "$LAMBDA_MD5" "$LAMBDA_DIR/J02459.fa.gz" | md5sum --check -
gzip --decompress --stdout "$LAMBDA_DIR/J02459.fa.gz" > "$LAMBDA_DIR/J02459.fa"
samtools faidx "$LAMBDA_DIR/J02459.fa"
```

## 3. BISCUIT indexes

Build the genomic and Lambda BISCUIT indexes:

```bash
for GENOME in mm10 hg38; do
  FASTA="$REF_DIR/$GENOME/$GENOME.fa"
  [[ -s "$FASTA" ]] || continue
  mkdir -p "$REF_DIR/$GENOME/biscuit"
  biscuit index -p "$REF_DIR/$GENOME/biscuit/$GENOME.fa" "$FASTA"
done

mkdir -p "$REF_DIR/Lambda/biscuit"
biscuit index -p "$REF_DIR/Lambda/biscuit/J02459.fa" \
  "$REF_DIR/Lambda/J02459.fa"
```

These output prefixes match `BISCUIT_INDEX` in the checked-in profiles.

## 4. BISCUIT QC assets and pipeline CpG BED

BISCUIT publishes the three QC assets required by `spatialmeth_qc.sh`. Release
`v1.10.2.20260818` is pinned here; its assets are compatible with BISCUIT 1.8
and later.

```bash
BISCUIT_QC_RELEASE=v1.10.2.20260818

for GENOME in mm10 hg38; do
  case "$GENOME" in
    mm10) QC_SHA256=54977edc5ef59103e7affe4b17aae24441cab35b8492a3457090705a43bbaa79 ;;
    hg38) QC_SHA256=14653c6f2fd630e84b434e02956bf72c6f67314d89b20cf0354a6083e54f4672 ;;
  esac

  ARCHIVE="$REF_DIR/$GENOME/biscuit/${GENOME}_biscuit_qc_assets.zip"
  QC_DIR="$REF_DIR/$GENOME/biscuit/QC"
  URL="https://github.com/huishenlab/biscuit/releases/download/${BISCUIT_QC_RELEASE}/${GENOME}_biscuit_qc_assets.zip"

  mkdir -p "$QC_DIR"
  curl --fail --location "$URL" -o "$ARCHIVE"
  printf '%s  %s\n' "$QC_SHA256" "$ARCHIVE" | sha256sum --check -
  unzip -j -o "$ARCHIVE" "$GENOME/*" -d "$QC_DIR"

  for NAME in cpg.bed.gz \
      windows100bp.gc_content.bot10p.bed.gz \
      windows100bp.gc_content.top10p.bed.gz; do
    [[ -s "$QC_DIR/$NAME" ]] || { echo "Missing QC asset: $QC_DIR/$NAME" >&2; exit 1; }
  done
done
```

The BISCUIT QC CpG files include non-primary contigs. The pipeline methylation
matrix and bundled ChromHMM masks use primary chromosomes only. Derive the
pipeline CpG BED from the downloaded QC asset without changing its coordinate
order:

```bash
for GENOME in mm10 hg38; do
  QC_CPG="$REF_DIR/$GENOME/biscuit/QC/cpg.bed.gz"
  CPG_DIR="$REF_DIR/$GENOME/annotation/cpg"
  mkdir -p "$CPG_DIR"

  if [[ "$GENOME" == hg38 ]]; then
    CPG_BED="$CPG_DIR/cpg_nocontig.bed.gz"
  else
    CPG_BED="$CPG_DIR/cpg.bed.gz"
  fi

  gzip --decompress --stdout "$QC_CPG" \
    | awk '$1 ~ /^chr([0-9]+|X|Y|M)$/' \
    | bgzip --stdout > "$CPG_BED"
  tabix --force --preset bed "$CPG_BED"
done
```

Expected primary-chromosome row counts are 21,867,837 for mm10 and 29,401,795
for hg38. For another assembly, generate the QC directory with the
`build_biscuit_QC_assets.pl` program installed by the release environment:

```bash
QC_DIR="$REF_DIR/$GENOME/biscuit/QC"
build_biscuit_QC_assets.pl --outdir "$QC_DIR" --ref "$REF_DIR/$GENOME/$GENOME.fa"
```

See the [BISCUIT QC documentation](https://huishenlab.github.io/biscuit/docs/alignment/QC.html)
for the asset definitions.

## 5. ChromHMM YAME masks

The exact compact masks used by the validated profiles are included under
`reference_assets/`. Copy them into the paths already named by
`profiles/*/site.yaml`:

```bash
mkdir -p "$REF_DIR/mm10/KYCGKB_mm10" "$REF_DIR/hg38/KYCGKB_hg38"
cp "$REPO/reference_assets/mm10/ChromHMM.20220414.cm" \
  "$REF_DIR/mm10/KYCGKB_mm10/ChromHMM.20220414.cm"
cp "$REPO/reference_assets/hg38/ChromHMM.20220303.cm" \
  "$REF_DIR/hg38/KYCGKB_hg38/ChromHMM.20220303.cm"

printf '%s  %s\n' \
  49a6a884892d92145b6f1e2ef5120e8338f1228a9ecdc9dd5edf9fbdb0d66bd8 \
  "$REF_DIR/mm10/KYCGKB_mm10/ChromHMM.20220414.cm" \
  | sha256sum --check -
printf '%s  %s\n' \
  c38ff4081fad37d4678e97c72264e1b80fd2cfb2a5377a5392c07e37f4f3e160 \
  "$REF_DIR/hg38/KYCGKB_hg38/ChromHMM.20220303.cm" \
  | sha256sum --check -
```

Confirm that each mask and CpG BED have the same row count:

```bash
for GENOME in mm10 hg38; do
  if [[ "$GENOME" == mm10 ]]; then
    CPG_BED="$REF_DIR/mm10/annotation/cpg/cpg.bed.gz"
    MASK="$REF_DIR/mm10/KYCGKB_mm10/ChromHMM.20220414.cm"
  else
    CPG_BED="$REF_DIR/hg38/annotation/cpg/cpg_nocontig.bed.gz"
    MASK="$REF_DIR/hg38/KYCGKB_hg38/ChromHMM.20220303.cm"
  fi
  CPG_ROWS=$(gzip --decompress --stdout "$CPG_BED" | wc -l)
  MASK_ROWS=$(yame info "$MASK" | awk 'NR == 2 {print $4}')
  [[ "$CPG_ROWS" == "$MASK_ROWS" ]] \
    || { echo "$GENOME row mismatch: CpG=$CPG_ROWS mask=$MASK_ROWS" >&2; exit 1; }
done
```

For a different assembly or ChromHMM segmentation, start from a sorted,
non-overlapping BED whose fourth column is the state label. Map one state to
every CpG, then pack the state vector with YAME:

```bash
CHROMHMM_BED=/path/to/chromhmm_primary_contigs.bed.gz
CPG_BED=/path/to/primary_contig_cpg.bed.gz
STATE_TXT=/path/to/ChromHMM.states.txt
CHROMHMM_CM=/path/to/ChromHMM.cm

bedtools map -a <(gzip --decompress --stdout "$CPG_BED") \
  -b <(gzip --decompress --stdout "$CHROMHMM_BED") \
  -c 4 -o first -null NA \
  | cut -f4 | sed 's/^ChromHMM;//' > "$STATE_TXT"
yame pack -f s "$STATE_TXT" "$CHROMHMM_CM"
```

## 6. RNA annotation and STAR indexes

RNA samples require a STAR index. The validated mm10 run used Ensembl release
82, converted from Ensembl contig names to UCSC names to match `mm10.fa`:

```bash
GTF_DIR="$REF_DIR/mm10/annotation/gtf"
GTF_GZ="$GTF_DIR/Mus_musculus.GRCm38.82.gtf.gz"
GTF="$GTF_DIR/Mus_musculus.GRCm38.82.gtf.UCSCnaming"
mkdir -p "$GTF_DIR"
curl --fail --location \
  https://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/Mus_musculus.GRCm38.82.gtf.gz \
  -o "$GTF_GZ"
printf '%s  %s\n' 4d1020ffe486d3911fc1d9ea80817121 "$GTF_GZ" | md5sum --check -
gzip --decompress --stdout "$GTF_GZ" \
  | sed -E '/^#/! s/^([0-9]+|X|Y)\t/chr\1\t/; /^MT\t/ s/^MT\t/chrM\t/' \
  > "$GTF"
```

For hg38, this guide pins the UCSC-compatible GENCODE v44 comprehensive gene
annotation:

```bash
GTF_DIR="$REF_DIR/hg38/annotation/gtf"
GTF_GZ="$GTF_DIR/gencode.v44.annotation.gtf.gz"
GTF="$GTF_DIR/gencode.v44.annotation.gtf"
mkdir -p "$GTF_DIR"
curl --fail --location \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz \
  -o "$GTF_GZ"
printf '%s  %s\n' ee330cfe6d0654ba9b9cf434d5c1bfb1 "$GTF_GZ" | md5sum --check -
gzip --decompress --stdout "$GTF_GZ" > "$GTF"
```

Build STAR with the annotation corresponding to each genome:

```bash
GENOME=mm10  # or hg38
if [[ "$GENOME" == mm10 ]]; then
  GTF="$REF_DIR/mm10/annotation/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming"
else
  GTF="$REF_DIR/hg38/annotation/gtf/gencode.v44.annotation.gtf"
fi

mkdir -p "$REF_DIR/$GENOME/STAR"
STAR --runMode genomeGenerate \
  --runThreadN 24 \
  --genomeDir "$REF_DIR/$GENOME/STAR" \
  --genomeFastaFiles "$REF_DIR/$GENOME/$GENOME.fa" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 149
```

`--sjdbOverhang 149` matches the validated 150-bp libraries. Use read length
minus one for a different library design; that intentionally defines a
different STAR reference.

## 7. Profile and preflight

Copy a checked-in profile and edit only `REF_DIR` when using the directory
layout above. The existing `REF_FASTA`, `BISCUIT_INDEX`, `CPG_BED`,
`BISCUIT_QCASSET`, and `CHROMHMM_CM` values then resolve without changes.

Before production, resolve the complete DAG:

```bash
python -s -m snakemake --profile /path/to/copied/profile \
  --snakefile "$REPO/main/Snakefile" \
  --config ref=mm10 IDS=sample_name \
  --dry-run -- final_output
```

## Optional extended feature QC

The standard `final_output` target does not require the broad annotation
library under `features/` and `KYCGKB_<genome>/`. Those files are consumed only
by the explicit extended `qc` target. They represent user-selected genomic
feature sets rather than protocol references, so they are not silently
downloaded or substituted. To use that target, provide each BED and matching
YAME mask named in `FEATS` in `main/Snakefile`, built against the same ordered
primary-chromosome CpG BED as the sample `.cg` matrix.
