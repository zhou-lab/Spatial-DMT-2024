# CpG-Coverage Occupancy Pilot

## Question

Can the number of covered CpGs per 10 um spatial pixel identify discrete
singlet, doublet, and higher-occupancy states?

The analysis calls fitted groups `coverage components`, not cell counts. A
component is only evidence for a cell-count state if its center follows integer
scaling, is well separated, remains stable with sequencing depth, and is not
explained by spatial barcode effects.

## Inputs

| Sample | Role | DNA fragments in QC BAM |
| --- | --- | ---: |
| E75_10um_rep2 | Shallow | 13.617 M |
| E65_10um_run2 | Deeper | 338.675 M |
| E65_10um_run1 | Even deeper, same-sample replicate | 535.794 M |

The E75 table was regenerated directly with `yame summary` from the retained
indexed `.cg`. Existing per-pixel CpG tables were reused for the E65 runs.

## Method

1. Populate the full 96 by 96 grid and set missing indexes to zero CpGs.
2. Use Otsu's threshold on positive `log10(covered CpGs + 1)` values to separate
   low-coverage background from foreground pixels.
3. Fit one-dimensional Gaussian mixtures to foreground log counts.
4. Select the component count by BIC and sort components by count.
5. Compare component centers with consecutive integer multiples.
6. Compare matched E65 pixels across the two sequencing runs.
7. Quantify row and column effects and inspect spatial component maps.
8. Repeat model selection with a maximum of eight rather than four components.

## Results

| Sample | Foreground pixels | Best K, cap 4 | Component centers | Ratios to C1 | Row/column adjusted R2 |
| --- | ---: | ---: | --- | --- | ---: |
| E75_10um_rep2 | 1,879 | 3 | 263, 1,557, 4,405 | 1.00, 5.92, 16.73 | 0.237 |
| E65_10um_run2 | 7,158 | 4 | 3,388, 7,246, 13,684, 23,247 | 1.00, 2.14, 4.04, 6.86 | 0.302 |
| E65_10um_run1 | 7,014 | 4 | 4,832, 9,793, 17,733, 29,560 | 1.00, 2.03, 3.67, 6.12 | 0.288 |

Matched E65 pixels were highly reproducible:

- Spearman correlation across all pixels: 0.9991
- Spearman correlation across joint foreground pixels: 0.9992
- Foreground/background agreement: 98.4%
- Four-component adjusted Rand index: 0.789

This reproducibility is not sufficient evidence for cell counts. The maps show
axis-aligned stripes spanning tissue regions: column effects are especially
visible in E75, while row bands recur in both E65 runs. These patterns are
consistent with spatial-barcode or capture-efficiency effects.

When the component cap was increased from four to eight, both E65 runs selected
all eight components. Their cross-run component ARI decreased to 0.617. This is
evidence that the GMM is partitioning a continuous coverage distribution rather
than discovering a fixed number of biological occupancy states.

An idealized synthetic control containing log-normal 1x, 2x, and 3x modes was
recovered as three components centered at 1,000, 1,999, and 2,996 CpGs. The
method can therefore recover integer modes when they are present under its
model assumptions.

## Conclusion

Per-pixel covered-CpG count works well as a tissue/background and relative
coverage metric, but these samples do not support interpreting its fitted modes
as singlets and doublets. The simple method is confounded by barcode-dependent
capture efficiency and continuous biological/technical variation.

A defensible cell-occupancy model would need at least per-barcode original DNA
molecule counts, correction for row/column capture effects, and calibration
against pixels with independently measured nuclear occupancy. CpG coverage can
remain one input feature, but should not be the sole cell-count measurement.

## Spatial-DMT comparison

The same pilot was run on retained 50 by 50 Spatial-DMT data:

| Sample | Foreground pixels | Best K, cap 8 | Component centers | High/low ratio | Row/column adjusted R2 |
| --- | ---: | ---: | --- | ---: | ---: |
| SpMETSLE14D | 1,841 | 2 | 5,048, 8,565 CpGs | 1.697 | 0.682 |
| SpMETSLE16D | 2,120 | 2 | 5,253, 9,168 CpGs | 1.745 | 0.529 |
| SpMETSLE14DM | 1,836 | 2 | 217,444, 366,771 CpGs | 1.687 | 0.702 |
| SpMETSLE14DM nonduplicate mapped reads | 1,860 | 2 | 318,864, 528,236 reads | 1.657 | 0.708 |

Unlike SmC-seq, all Spatial-DMT inputs selected two components even when the
model was allowed up to eight. The shallow SpMETSLE14D and deep SpMETSLE14DM
profiles were nearly identical across pixels: Spearman correlation was 0.9970,
the median deep/shallow count ratio was 42.6, and component assignments agreed
for 96.3% of shared foreground pixels.

This is a reproducible low/high DNA-recovery signal, but not evidence of
singlets and doublets. The centers differ by about 1.7-fold rather than
two-fold, and additive row/column effects explain 53% to 70% of foreground
log-count variance. After removing those effects, SpMETSLE16D selects one
component. SpMETSLE14D and SpMETSLE14DM retain a smaller low-count group, but
its high/low separation contracts to about 1.30-fold and 1.37-fold,
respectively.

These three Spatial-DMT comparators are 50 um, where a tissue pixel can contain
many cells. The fitted classes should therefore be treated as low/high tissue
or DNA capture states. External nuclear occupancy measurements are still needed
for a calibrated singlet/doublet test.

## Recent production deep samples

The pilot was rerun on two non-experimental production samples that completed
on August 5, 2026. SpMETSLE17DM is the E11#2 50 um sample and SpMETSLB14DM is
the P21 20 um sample. Both analyses allowed up to eight components.

| Sample | Raw best K | Raw centers or center | Raw ratios | Grid adjusted R2 |
| --- | ---: | --- | --- | ---: |
| SpMETSLE17DM | 4 | 184,693 / 272,532 / 399,457 / 598,149 | 1.000 / 1.476 / 2.163 / 3.239 | 0.773 |
| SpMETSLB14DM | 1 | 189,067 | 1.000 | 0.688 |

Neither raw result supports 1x/2x/4x occupancy. SpMETSLE17DM partitions a
continuous distribution into four noninteger modes, while SpMETSLB14DM is
unimodal. The maps show strong column-aligned coverage structure.

Robust median-polish normalization was therefore applied to foreground
log-counts while retaining the original foreground mask. It reduced adjusted
row/column R2 from 0.773 to 0.002 for SpMETSLE17DM and from 0.688 to -0.003 for
SpMETSLB14DM, demonstrating that the targeted grid artifact was removed.

| Sample | Normalized best K | Normalized centers or center | Normalized ratios |
| --- | ---: | --- | --- |
| SpMETSLE17DM | 3 | 207,967 / 297,332 / 320,495 | 1.000 / 1.430 / 1.541 |
| SpMETSLB14DM | 2 | 180,469 / 232,626 | 1.000 / 1.289 |

Normalization did not expose integer occupancy. The residual components are
close noninteger subdivisions of skewed distributions and their spatial maps
follow tissue contours. The negative result is therefore not explained by the
row/column artifact alone.

As a positive control, 1x/2x/4x states were simulated on the real
SpMETSLE17DM foreground mask using the normalized lowest center of 207,967 CpGs
and a conservative within-mode CV of 0.18. BIC selected K=3 at centers 205,745,
405,033, and 815,377 CpGs, with ARI 0.921 and 97.2% assignment accuracy. Thus,
integer occupancy would remain detectable at this production sample's depth
and residual variability if it were present.

The sample-level packed `.cg.idx` files produced by the previous Snakefile
were found to have incorrect barcode offsets: repeated `yame index -1` calls
do not correctly index a growing multi-sample file. All production results in
this section were therefore regenerated directly from the 2,500 retained
per-cell `.cg` files. Their covered-CpG counts correlate with mapped reads at
Spearman 0.9972 for SpMETSLE17DM and 0.9989 for SpMETSLB14DM. Invalid
packed-index exploratory outputs are isolated under
`production_deep/invalid_packed_index_results/`. The Snakefile now concatenates
all blocks and calls `yame index -s` once with the ordered barcode list.

### Library-complexity alternative

The same production samples were tested using estimated total DNA library
complexity instead of covered-CpG count. No upstream processing was rerun. The
analysis used retained per-barcode BAM intermediates and their fragment-level
mapped and duplicate counts. Direct BAM scans at coordinates 0101, 2525, and
5050 exactly reproduced the table counts for both samples.

For each pixel, total library molecules `L` were estimated by solving the
standard random-sampling model
`unique = L * (1 - exp(-mapped / L))`. This is a duplicate-based extrapolation,
not a direct molecule count, and assumes molecules are sampled uniformly and
duplicate marking reflects repeated observations of the same original
molecule.

| Sample | Median foreground complexity | Raw K and ratios | Grid adjusted R2 | Normalized K and ratios | Normalized R2 |
| --- | ---: | --- | ---: | --- | ---: |
| SpMETSLE17DM | 1,558,403 | K=2; 1.000 / 1.738 | 0.811 | K=2; 1.000 / 1.118 | -0.009 |
| SpMETSLB14DM | 339,178 | K=1; 1.000 | 0.710 | K=2; 1.000 / 1.293 | -0.004 |

The apparently near-twofold raw split in SpMETSLE17DM is spatially organized
and collapses after row/column normalization. SpMETSLB14DM is unimodal before
normalization and has only a weak 1.29-fold residual subdivision afterward.
Neither result supports 1x/2x/4x occupancy.

Within foreground pixels, duplication rates vary over only 15.5% to 16.8% in
SpMETSLE17DM and 50.3% to 52.6% in SpMETSLB14DM from the 10th to 90th
percentiles. Consequently, estimated complexity remains almost collinear with
unique fragments, mapped fragments, and CpG coverage (Spearman 0.990 to 1.000).
It does not add enough independent information to overcome the same spatial
capture artifact seen with CpG counts.

### Published SmC-seq 10 um library complexity

The library-complexity test was repeated on all five published SmC-seq 10 um
libraries using per-pixel `Mapped_reads` and `Unique_reads` from Supplementary
Table 2. This avoids rerunning cleaned BAM intermediates and directly uses the
published duplicate information. The analysis covered 4,975 reported pixels;
foreground selection retained 4,766 pixels.

| Sample | Foreground pixels | Median library observed | Raw K and ratios | Raw grid R2 | Normalized K and ratios | Normalized grid R2 |
| --- | ---: | ---: | --- | ---: | --- | ---: |
| E5.5 rep. 1 | 264 | 99.1% | K=1; 1.000 | 0.700 | K=2; 1.000 / 1.123 | -0.047 |
| E5.5 rep. 2 | 267 | 99.2% | K=1; 1.000 | 0.719 | K=1; 1.000 | -0.061 |
| E6.5 | 1,553 | 89.4% | K=4; 1.000 / 1.430 / 2.304 / 3.191 | 0.775 | K=4; 1.000 / 1.071 / 1.212 / 20.348* | -0.009 |
| E7.5 rep. 1 | 1,384 | 99.9% | K=2; 1.000 / 1.700 | 0.570 | K=4; 1.000 / 1.593 / 1.905 / 41.135* | 0.019 |
| E7.5 rep. 2 | 1,298 | 28.4% | K=3; 1.000 / 1.437 / 2.194 | 0.753 | K=2; 1.000 / 1.082 | -0.058 |

`*` The highest normalized E6.5 and E7.5-rep.-1 components each contain one
pixel and are outliers, not occupancy classes.

Raw library complexity is strongly organized by spatial rows and columns in
every sample. Median polish reduces adjusted grid R2 from 0.57-0.78 to
approximately zero. After that correction, four of five libraries are
unimodal or have meaningful component centers separated by only 1.07-1.21x.
E7.5 replicate 1 retains a candidate 1.9x subgroup, but also an intermediate
1.59x group; neither subgroup is reproduced in E7.5 replicate 2. Therefore,
10 um capture produces a more uniform per-pixel library than larger pixels,
but read-derived complexity alone does not provide a reproducible
singlet/doublet assignment. The E7.5-rep.-1 split is a candidate for validation
against independent nuclei or cell-segmentation labels, not a validated cell
count.

The combined summary is
`smcseq_10um_library_complexity/smcseq_10um_library_complexity_summary.png`.
Per-sample distributions, spatial maps, assignments, and exact model summaries
are under `smcseq_10um_library_complexity/results/`.

`SpMETSLE23DM` was an experimental sample and is excluded from the evidence
set. Its previous exploratory outputs are retained only for provenance under
`analysis/cpg_occupancy/excluded_experimental/SpMETSLE23DM/`.

## Simulation-based accessibility

A finite-genome simulation tested whether CpG sparsity alone prevents recovery
of known 1x, 2x, and 4x occupancy states. It used 21,867,837 possible CpGs,
the observed low-component depths and foreground pixel counts, occupancy
weights of 0.5/0.3/0.2, and 20 replicates per condition. A mean-one log-normal
pixel-level capture multiplier was added over a sweep of capture coefficients
of variation (CV). This shared multiplier is conservative relative to
independent cell-to-cell yield variation because it does not average down in
multi-cell pixels.

Finite-genome saturation was negligible at the observed depths:

| Baseline | Expected 1x/2x/4x CpGs | Loss from linearity at 4x |
| --- | --- | ---: |
| SpMETSLE14D, 5,048 CpGs | 5,047 / 10,093 / 20,181 | 0.035% |
| SpMETSLE16D, 5,253 CpGs | 5,252 / 10,503 / 21,001 | 0.036% |
| SpMETSLE14DM, 217,444 CpGs | 216,367 / 430,593 / 852,707 | 1.47% |

At zero capture variation, 100 covered CpGs per 1x pixel were already enough
for perfect classification in these simulations. At CV 0.24, accuracy plateaued
near 91% once depth reached a few hundred CpGs, and all 60 observed-depth
replicates selected the true three-component model by BIC. At CV 0.28,
known-center accuracy was 86% to 87% and BIC selected three components in 75%
to 95% of replicates. At CV 0.35, BIC selected three components in only 0% to
5% of replicates despite additional sequencing depth.

Thus, the observed CpG sparsity is not the theoretical limit. If exact integer
occupancy states existed and multiplicative variation stayed near the observed
within-mode CV of 0.22 to 0.28, count-based inference should detect them with
useful but imperfect accuracy. The practical limit is pixel-level capture and
biological variability: deeper sequencing cannot rescue the model once that
CV approaches roughly 0.3 to 0.35. This simulation establishes accessibility
under its assumptions; it does not establish that the components in the real
50 um data represent cell counts.

## Reproduction

The analysis script is `tools/cpg_coverage_occupancy.py`. Its isolated
environment is `/scr1/users/suc3/envs/cpg_occupancy`.

Primary outputs are under `analysis/cpg_occupancy/results/`; the eight-component
sensitivity analysis is under `analysis/cpg_occupancy/sensitivity_k8/`.
Spatial-DMT results are under `analysis/cpg_occupancy/spatialdmt/`.
Recent production deep-sample results are under
`analysis/cpg_occupancy/production_deep/`; grid normalization is implemented
by `tools/cpg_grid_normalized_occupancy.py`.
Library-complexity estimation and occupancy analysis are implemented by
`tools/library_complexity_occupancy.py`; their outputs are under
`analysis/cpg_occupancy/production_deep/library_complexity_results/`.
Published SmC-seq pixel counts are prepared by
`tools/prepare_smcseq_complexity.py`, and the 10 um summary figure is generated
by `tools/plot_smcseq_10um_complexity.py`.
Simulation results are under `analysis/cpg_occupancy/simulation/`; the driver
is `tools/simulate_cpg_occupancy_accessibility.py`.
