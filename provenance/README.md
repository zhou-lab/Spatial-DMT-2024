# Pipeline release provenance

`pipeline_2.0.1_samples.tsv` records completed revised-pipeline outputs that
are compatible with release `2.0.1`. Historical outputs predate the version
artifact, so they use the explicit `backfilled_compatible` assignment rather
than claiming that the release tag existed when the computation ran.

The backfill is justified for runs that use the revised protocol-aware output
schema and have the complete final DNA deliverable set. Later changes included
performance-equivalent implementations and report/QC corrections; those runs
remain part of the `2.0.1` result cohort. Incomplete final directories are
listed but are not labeled. Technical and experimental runs are labeled for
provenance but remain explicitly classified so they can be excluded from
biological summaries.

Future runs receive the same `pipeline_version.tsv` schema directly from the
Snakemake `pipeline_version` rule with assignment `workflow_generated`.
