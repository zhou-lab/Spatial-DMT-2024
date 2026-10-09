# Changelog

## 2.0.1 - 2026-10-09

- Added protocol-aware Spatial-DMT and SmC-seq DNA processing with optional RNA co-profiling.
- Added matched and unmatched read tracks, barcode-aware duplicate removal, per-pixel methylation pileup, and detailed failure-read diagnostics.
- Added curated MultiQC and self-contained HTML reports with explicit fragment-, record-, pre-filter-, and post-filter metric definitions.
- Corrected sample-level YAME aggregation, designed-barcode counts, linker/raw denominators, duplicate accounting, lambda totals, and SmC-seq mapping metrics.
- Integrated validated streaming trim, native barcode splitting, resumable size-balanced pileup, and configurable Slurm resource allocation for deep libraries.
- Added a release-based `pipeline_version.tsv` artifact and an audited compatibility manifest for existing revised-pipeline outputs.
