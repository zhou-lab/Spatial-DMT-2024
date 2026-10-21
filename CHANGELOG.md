# Changelog

## 2.0.3 - 2026-10-21

- Added pinned URLs and checksums for UCSC hg38/mm10 FASTAs, the NCBI Lambda
  reference, BISCUIT QC archives, and RNA annotations.
- Documented complete BISCUIT, CpG BED, Lambda, and STAR reference construction
  using the directory layout expected by the checked-in profiles.
- Bundled the exact compact hg38 and mm10 ChromHMM YAME masks used by the
  validated reports, with row-count and checksum verification instructions.
- Distinguished assets required by the standard `final_output` target from
  user-selected feature collections required only by the extended `qc` target.
- Made no workflow-logic or metric-definition changes.

## 2.0.2 - 2026-10-20

- Added a fully solved Micromamba environment for the workflow, including
  protocol-specific tools, the Slurm executor plugin, BBTools, and native-build
  dependencies.
- Removed hidden site-software defaults for Python, YAME, STAR, BBTools, and
  HTSlib; vendored native programs now build against HTSlib from the active
  environment.
- Isolated Snakemake and rule subprocesses from user-site Python packages and
  documented reference preparation and clean-environment execution.
- Made Slurm rule submissions explicitly inherit the controller environment so
  active-environment tools remain available under restrictive parent export
  policies.
- Consolidated sample-level and ChromHMM methylation summaries into one YAME
  aggregation pass and corrected the inline MultiQC configuration syntax.
- Replaced GNU Parallel's trim byte-stream launcher with a bounded,
  newline-aligned coordinator that propagates worker failures reliably.
- Added deterministic lambda-pileup depth capping with an explicit sampling
  manifest, avoiding BISCUIT numerical overflow in deeply sequenced controls
  while retaining full-BAM alignment and duplication metrics; the manifest is
  copied into each sample's final deliverables.
- Made the batched pileup runner resolve documented executable command names
  through `PATH` as well as accepting explicit paths.
- Clarified that unmatched-read failure decomposition is FASTQ-backed by
  default, with genomic alignment and BAM-derived fields available as an
  opt-in diagnostic.
- Added tests for reference CpG coordinates and precomputed YAME summaries.

## 2.0.1 - 2026-10-09

- Added protocol-aware Spatial-DMT and SmC-seq DNA processing with optional RNA co-profiling.
- Added matched and unmatched read tracks, barcode-aware duplicate removal, per-pixel methylation pileup, and detailed failure-read diagnostics.
- Added curated MultiQC and self-contained HTML reports with explicit fragment-, record-, pre-filter-, and post-filter metric definitions.
- Corrected sample-level YAME aggregation, designed-barcode counts, linker/raw denominators, duplicate accounting, lambda totals, and SmC-seq mapping metrics.
- Integrated validated streaming trim, native barcode splitting, resumable size-balanced pileup, and configurable Slurm resource allocation for deep libraries.
- Added a release-based `pipeline_version.tsv` artifact and an audited compatibility manifest for existing revised-pipeline outputs.
