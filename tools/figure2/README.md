# Figure 2

This directory builds the data-driven Figure 2 for the Spatial-DMT QC paper.
The plotting code reads revised pipeline tables first and records every plotted
metric's source, definition, and fallback status in tab-separated metadata.

## Run

The existing `qc_env` micromamba environment contains the plotting stack. The
external `yame` and `samtools` executables must be available on `PATH`.

```bash
micromamba run -n qc_env python tools/figure2/build_figure2.py \
  --repo-root /home/suc3/Lab/Spatial-DMT-2024 \
  --output-dir /home/suc3/Lab/Spatial-DMT-2024/results/figure2
```

To create an independent environment:

```bash
micromamba create -f tools/figure2/environment.yml
```

An independent MultiQC-style comparison of the six revised samples with
published SmC-seq QC and its direct scNMT/bulk validation benchmarks is built
with:

```bash
micromamba run -n qc_env python tools/figure2/build_panel_a_external.py \
  --repo-root /home/suc3/Lab/Spatial-DMT-2024 \
  --output-dir /home/suc3/Lab/Spatial-DMT-2024/results/figure2
```

This comparison harmonizes current input pairs, SmC assigned reads, and bulk
input reads into one input column; it likewise combines spatial mapping rates
and CpG-yield columns while retaining each cell's original unit in the output
TSV. Current-workflow global mCG is the median per-pixel YAME Beta. The revised
16D `.cg` was unavailable after cleanup, so only that cell is an explicitly
marked legacy fallback. RNA metrics are omitted and missing or incomparable
cells remain blank. The validation block records the published SmC-seq Pearson
correlations at 1 Mb, 100 kb, and 10 kb resolution and the E7.5 1-kb
consistency fractions.

A separate Panel A containing only final outputs from the protocol-aware
unified pipeline is built with:

```bash
micromamba run -n qc_env python tools/figure2/build_panel_a_unified.py \
  --repo-root /home/suc3/Lab/Spatial-DMT-2024 \
  --output-dir /home/suc3/Lab/Spatial-DMT-2024/results/figure2
```

This version excludes literature-only bulk/single-cell rows and methylation
validation coefficients. It compares the Spatial-DMT and SmC-seq protocol
branches using the shared curated metrics definitions, retains global mCG, and
writes `panel_a_unified_pipeline_audit.tsv` for omitted metrics or incomplete
runs. No published-paper values or legacy-pipeline outputs are mixed into it.

The independent unified-pipeline retention comparison is built with:

```bash
micromamba run -n qc_env python tools/figure2/build_panel_b_published_band.py \
  --repo-root /home/suc3/Lab/Spatial-DMT-2024 \
  --smcseq-trim-run-root /scr1/users/suc3/20260820_smcseq_retention_speedopt \
  --output-dir /home/suc3/Lab/Spatial-DMT-2024/results/figure2
```

The unified Panel A `QC BAM / pre-filter %` value is always computed from
primary-record flagstats (`num_reads_qc_bam / num_reads_merged`). A bsconv log
retention fraction is retained in the curated table but is not substituted for
that flagstat ratio, because the log can count non-primary records.

It loads completed Spatial-DMT and SmC-seq samples exclusively from outputs
produced by the shared Snakefile. The current Spatial-DMT and SmC-seq
trajectories are normalized to raw DNA reads (`2 x raw read pairs`). Spatial-DMT includes a measured
linker-positive checkpoint from `Reads_total - Unmatched_SF`; both genomic
mates are aligned, so paired-end primary BAM records are already read counts.
The optimized SmC-seq trim-only reruns provide a measured single-R1
checkpoint. Their raw totals are checked against the retained final-output
denominator, and their passed-pair totals are checked against the
pre-alignment record counts.
Only genomic R1 is aligned for SmC-seq. Spatial-DMT has no
bsconv branch, so its mapped and QC-retained checkpoints coincide; SmC-seq uses
the post-bsconv QC BAM at the final checkpoint. Cohort medians have bootstrap
95% confidence intervals. A single E6.5 reproduced trajectory is overlaid
without a confidence interval and retains its reported designed checkpoint.
The E6.5 reproduced source table reports structure-pass, mapped, and
CpH-passed counts for one selected mate per input pair. Its checkpoints
therefore use the same `2 x raw read pairs` denominator as live SmC-seq
(38.02% designed, 7.81% mapped, and 7.22% quality-filtered); 50% is the
single-mate upper bound.
Cohort sizes are omitted from the visible figure. `SpMETSLE21DM` and
`SpMETSLE23DM` are loaded from their completed final reports; their trim,
alignment, pileup, unmatched-diagnostics, and report-generation branches are
complete.

The builder writes the canonical
`Figure2_experimental_panel_b_spmet_band.png` plus the metrics and audit TSVs.
Use `--output-name Figure2_experimental_panel_b_spmet_band_v7` for a versioned
revision. For Spatial-DMT, the linker-positive checkpoint is defined as
`(Reads_total - Unmatched_SF) / Reads_total`; `Unmatched_SF` is the structure
failure count, and the denominator is the raw `Reads_total`, independent of
lambda alignment. This includes structure-positive reads that are subsequently
removed by barcode or read-length checks. The plot labels this checkpoint
`Linker-positive reads`. It does not create legacy filename aliases or
PDF/SVG copies.

Use --output-name Figure2_experimental_panel_b_spmet_band_v7 for the
annotation-free version. Its visible checkpoints are Raw DNA reads, Trim
passed reads, Mapped reads, and Quality filtered reads. For Spatial-DMT, the
second checkpoint is the requested linker-positive quantity even though the
visible label remains Trim passed reads. The plotted fractions use a raw-read
denominator (2 x raw read pairs); SmC-seq trim, mapped, and quality-filtered
values are single genomic-R1 records over that two-read denominator. The audit
TSV also records literal `Reads_passed` and its pair fraction so the two
Spatial-DMT checkpoints are not conflated.

For the SmC-seq branch, Reads_total is counted once per input R1/R2 pair, while
each passing pair emits one genomic R1 record. Therefore the Trim passed reads
checkpoint uses Reads_passed / (2 x Reads_total), applying the single-read
penalty. Mapped and quality-filtered checkpoints use the same raw-read
denominator.

The standalone spatial panels C and D are built with:

```bash
micromamba run -n qc_env python tools/figure2/build_spatial_panels_cd.py \
  --repo-root /home/suc3/Lab/Spatial-DMT-2024 \
  --output-dir /home/suc3/Lab/Spatial-DMT-2024/results/figure2
```

Panel C compares observed low-read scar bands in revised Spatial-DMT
`SpMETSLHSK09D` and SmC-seq `E65_10um_run1` outputs. The per-pixel TSVs for
the cleaned SmC-seq run are no longer available, so the script reproducibly
extracts the named `DNA reads per barcode (log10)` PNG from the final
self-contained `qc_report.html`; the Spatial-DMT example is loaded the same
way for consistent rendering. Red boxes mark selected observed bands and do
not imply automatic detection. Panel D reproduces the 14D sound and 16D
incorrect spatial maps from the retained barcode-count tables with one shared
logarithmic normalization. Exact report names, tables, fallback status, and
highlighted intervals are written to
`figure2_spatial_panels_CD_provenance.tsv`.
For visual orientation matching in Panel D, the displayed 16D matrix is rotated
90 degrees counterclockwise; the exact source-to-display coordinate transform
is recorded in the provenance table.

The spatial CpG-methylation panel is built with:

**2026-09-18 audit:** `Figure2_spatial_methylation_v1` is not suitable for
paper-equivalent spatial comparisons. Its estimator differs from the published
call-weighted methylation; several samples require different physical barcode
coordinates; and three legacy indexes scramble record names. See
`results/figure2/smcseq_paper_comparison/README.md` and the comparison artifacts.
Four samples show strong agreement after harmonization; the remaining four
are unresolved or lack a direct reference. Do not infer validation from v1.

```bash
micromamba run -n qc_env python tools/figure2/build_spatial_methylation_panel.py \
  --repo-root /home/suc3/Lab/Spatial-DMT-2024 \
  --run-root /scr1/users/suc3/20260528_runNewSpatialPipe \
  --output-dir /home/suc3/Lab/Spatial-DMT-2024/results/figure2 \
  --reuse-tsv
```

The retained `qc_report.html` files do not contain spatial methylation maps;
they contain aggregate methylation summaries plus spatial read, mapping, and
duplication maps. This panel summarizes the indexed YAME `.cg` records at each
spatial pixel for the five current Spatial-DMT samples and eight current SmC-seq
samples. Pixels with fewer than 100 covered CpGs are gray, and all maps share a
0–1 CpG-methylation scale. PNG, editable SVG, PDF, and a provenance TSV are
written as `Figure2_spatial_methylation_v1.*`.

For the Spatial-DMT-only read-count/methylation comparison, run:

```bash
micromamba run -n qc_env python tools/figure2/build_spmet_count_methylation_panel.py
```

This writes `Figure2_spmet_count_methylation_v1.png` with one sample per row,
DNA read counts on the left, and CpG methylation on the right. Both columns use
shared pooled 2nd–98th percentile display ranges; the narrowed methylation scale
makes within-cohort spatial differences visible while preserving one common
normalization across all five samples.

The corrected SmC-seq comparison is built with:

```bash
micromamba run -n qc_env python tools/figure2/build_smcseq_count_methylation_panel.py
```

It writes `Figure2_smcseq_count_methylation_corrected_v2.png` from the audited
barcode joins in `results/figure2/smcseq_paper_comparison/`. Barcode sequence
is joined to the authors' sample-specific physical A/B coordinates, both axes
use the published reversed orientation, and current methylation is calculated
as call-weighted `sum(M) / sum(M + U)`. Only the four samples with strong
matched-pixel validation are shown: both E5.5 replicates, E7.5 replicate 1,
and E8.5 20 µm. The left column uses published raw counts because the current
QC count tables were cleaned and the retained report maps use the corrupted
local geometry. The two E6.5 runs and E7.5 replicate 2 have unresolved geometry
or identity, and placenta lacks a direct published counterpart.

`--skip-expensive` reuses the existing deep-sample CpG cache but does not
silently substitute values when that cache is absent. `--force-recompute`
replaces it. Output formats can be selected with, for example,
`--format pdf,svg,png,tiff`.

## Definitions and provenance

- Revised spatial BAM-fragment tables define mapping, duplication, covered
  pixels, and unique-fragment output.
- Aggregate unique CpGs use revised curated metrics where available. The deep
  sample's missing aggregate is recomputed in one pass from its revised `.cg`
  file using `yame rowop -o musum` followed by `yame summary`. Per-pixel unique
  fragments are used instead of an impractical serial scan of every `.cg`
  index; the figure does not fabricate missing per-pixel CpG values.
- RNA complexity uses STARsolo raw spatial counts for pixels with at least one
  UMI; RNA mapping uses the STARsolo `GeneFull/Summary.csv` definition.
- Panel A uses a MultiQC general-statistics layout with in-cell bars. Count-bar
  lengths use log10 scaling so the deeply sequenced 14DM sample does not make
  shallower count bars disappear.
- Revised shallow-run values whose temporary scratch sources were cleaned are
  restored from `revised_shallow_metrics_snapshot.tsv`. The snapshot preserves
  the audited numeric value, metric definition, and original upstream path for
  each cell. It replaces only missing or legacy-fallback cells; a live revised
  source always takes precedence.
- Panel B compares the Spatial-DMT and SmC-seq branches as fractions of raw DNA
  reads. Spatial-DMT raw pair counts come from live trim statistics and are
  multiplied by two; its linker-positive checkpoint is
  `Reads_total - Unmatched_SF`, while mapped and final read counts come from
  the unified curated table. SmC-seq raw pairs are cross-checked against final
  MultiQC general statistics and are also multiplied by two; its measured trim
  checkpoint is one emitted genomic R1 per passing pair, while mapped and final
  counts use the unified curated table directly because that branch aligns only
  genomic R1. Faint lines
  show individual samples; thick lines and ribbons show cohort medians and
  sample-level bootstrap 95% confidence intervals. The E6.5 reproduced run is
  shown as one line without a confidence interval; its source table uses one
  selected-mate record per input pair over the two-read denominator, matching
  the live SmC-seq read-unit convention.
- Panel C compares revised 14D and 16D raw DNA read counts per pixel with one
  shared normalization and no color bar.
- Panel D shows the revised `SpMETSLHSK09D` spatial coverage map with a manual
  red annotation around the low-read scar. It makes no automatic-detection
  claim.

Missing values remain `NaN`/`NA`; they are never replaced by zero. Legacy
fallbacks, if any, are listed in `figure2_provenance.tsv`.

## Automatic scar detection and methylation comparison

Run the paired-line scar analysis with:

```bash
/home/suc3/micromamba/envs/qc_env/bin/python \
  tools/figure2/analyze_scar_methylation.py \
  --run-root /scr1/users/suc3/20260528_runNewSpatialPipe \
  --output-dir results/figure2/scar_methylation
```

The detector compares each pixel on a candidate row or column with the
geometric mean of its two immediate neighbors across that line. It requires
broad paired support and center occupancy, so fully blank structural lines in
the corrected SmC-seq grid are excluded rather than reported as scars. The
analysis selects the strongest supported line in one Spatial-DMT and one
validated corrected-geometry SmC-seq sample, then compares its call-weighted
CpG methylation with the two adjacent lines at matched positions. Outputs are
the comparison PNG, a metrics TSV, and the ranked candidate TSV under
`results/figure2/scar_methylation/`; detailed source and interpretation notes
are in that directory's `README.md`.

The direction-swapped follow-up uses additional Spatial-DMT column scars and
corrected-geometry SmC-seq row scars:

```bash
/home/suc3/micromamba/envs/qc_env/bin/python \
  tools/figure2/analyze_scar_methylation_orthogonal.py \
  --run-root /scr1/users/suc3/20260528_runNewSpatialPipe \
  --comparison-dir results/figure2/smcseq_paper_comparison \
  --output-dir results/figure2/scar_methylation_orthogonal
```

It recomputes call-weighted Spatial-DMT methylation from the indexed `.cg`
files and writes the expanded comparison PNG, metrics, candidate rankings, and
source notes under `results/figure2/scar_methylation_orthogonal/`.

## Outputs

The build writes `Figure2.pdf`, `Figure2.svg`, `Figure2.png`, `Figure2.tiff`,
`figure2_metrics.tsv`, `figure2_provenance.tsv`, and
`figure2_selected_samples.tsv` under `results/figure2/`. Expensive intermediate
tables are retained under `results/figure2/cache/` for exact reruns.

## Caption draft

**Figure 2 | Spatially aware quality-control outputs across samples, cohorts,
modalities, and failure modes.** **A**, MultiQC-style cohort summary of revised
read-level, methylation, spatial-fragment, and RNA metrics; count-bar lengths use
log10 scaling. **B**, Retention through directly comparable unified-pipeline
checkpoints for the Spatial-DMT and SmC-seq branches, with a single reproduced
E6.5 processing trajectory. Lines show cohort medians, ribbons show sample-level
bootstrap 95% confidence intervals, and faint trajectories show individual
samples; the reproduced E6.5 trajectory is shown without a confidence interval.
**C**, Revised 14D and 16D raw DNA read counts per pixel, shown with identical
hidden color normalization. **D**, Revised spatial coverage map with a manually
annotated low-read horizontal scar in `SpMETSLHSK09D`.
