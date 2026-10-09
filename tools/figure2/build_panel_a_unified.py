#!/usr/bin/env python3
"""Build a MultiQC-style Panel A from unified-pipeline final outputs."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from plot_style import apply_style


RUN_ROOT = Path("/scr1/users/suc3/20260528_runNewSpatialPipe")


@dataclass(frozen=True)
class SampleSpec:
    sample: str
    display_name: str
    protocol: str
    group: str
    color: str


SAMPLES = [
    SampleSpec("SpMETSLE14D", "E11#1, 50 um (shallow)", "spatialdmt", "Spatial-DMT branch", "#0072B2"),
    SampleSpec("SpMETSLE14DM", "E11#1, 50 um (deep)", "spatialdmt", "Spatial-DMT branch", "#0072B2"),
    SampleSpec("SpMETSLB14DM", "P21 brain, 20 um", "spatialdmt", "Spatial-DMT branch", "#009E73"),
    SampleSpec("SpMETSLE17DM", "E13, 50 um", "spatialdmt", "Spatial-DMT branch", "#0072B2"),
    SampleSpec("SpMETSLE16D", "E13.5, 50 um", "spatialdmt", "Spatial-DMT branch", "#0072B2"),
    SampleSpec("E55_10um_rep1", "E5.5, 10 um rep. 1", "smcseq", "SmC-seq branch", "#D55E00"),
    SampleSpec("E55_10um_rep2", "E5.5, 10 um rep. 2", "smcseq", "SmC-seq branch", "#D55E00"),
    SampleSpec("E65_10um_run1", "E6.5, 10 um run 1", "smcseq", "SmC-seq branch", "#D55E00"),
    SampleSpec("E65_10um_run2", "E6.5, 10 um run 2", "smcseq", "SmC-seq branch", "#D55E00"),
    SampleSpec("E75_10um_rep1", "E7.5, 10 um rep. 1", "smcseq", "SmC-seq branch", "#D55E00"),
    SampleSpec("E75_10um_rep2", "E7.5, 10 um rep. 2", "smcseq", "SmC-seq branch", "#D55E00"),
    SampleSpec("E85_20um", "E8.5 embryo, 20 um", "smcseq", "SmC-seq branch", "#D55E00"),
    SampleSpec("E85_placenta_10um", "E8.5 placenta, 10 um", "smcseq", "SmC-seq branch", "#D55E00"),
]


# metric, display label, scaling, color, unit, definition
METRICS = [
    (
        "prefilter_primary_reads",
        "Pre-filter\nprimary records",
        "count",
        "#2C7FB8",
        "primary BAM records",
        "primary records in pre-filter merged deduplicated BAM",
    ),
    (
        "mapping_rate",
        "Mapping\n%",
        "percent",
        "#2CA25F",
        "% pre-filter primary records",
        "pre-filter primary mapped divided by pre-filter primary records",
    ),
    (
        "qc_retained_fraction",
        "QC BAM /\npre-filter %",
        "percent",
        "#66C2A4",
        "% pre-filter primary records",
        "QC-BAM primary records divided by pre-filter primary records",
    ),
    (
        "conversion_efficiency",
        "Non-CpG\nconv. %",
        "percent",
        "#31A354",
        "%",
        "100 * (1 - mean(CA, CC, CT retention))",
    ),
    (
        "indexed_pixels",
        "Indexed\npixels",
        "count",
        "#969696",
        "pixels",
        "number of entries in final sample .cg.idx",
    ),
    (
        "covered_cpgs",
        "Covered\nCpGs",
        "count",
        "#756BB1",
        "CpG sites",
        "aggregate CpGs with M+U > 0 after summing indexed pixels",
    ),
    (
        "cpg_depth",
        "CpG\ndepth",
        "count",
        "#9E9AC8",
        "x",
        "aggregate CpG depth reported by yame summary",
    ),
    (
        "global_mcg",
        "Global\nmCG %",
        "percent",
        "#E6AB02",
        "%",
        "record-pooled mean CG beta from per-pixel pileup VCFs",
    ),
]

GROUP_HEADERS = [
    ("Input", 0, 0),
    ("Processing", 1, 3),
    ("Yield", 4, 6),
    ("Methylation", 7, 7),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--run-root", type=Path, default=RUN_ROOT)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def read_curated(path: Path, sample: str) -> dict[str, float]:
    frame = pd.read_csv(path, sep="\t", comment="#")
    row = frame.loc[frame["Sample"].eq(sample)]
    if len(row) != 1:
        raise ValueError(f"Expected one curated row for {sample} in {path}, found {len(row)}")
    result: dict[str, float] = {}
    for key, value in row.iloc[0].items():
        if key == "Sample":
            continue
        result[key] = float(pd.to_numeric(value, errors="coerce"))
    return result


def metric_row(
    spec: SampleSpec,
    order: int,
    metric: str,
    value: float,
    source_file: Path,
    source_type: str,
    definition: str,
) -> dict[str, object]:
    unit = next(item[4] for item in METRICS if item[0] == metric)
    return {
        "row_order": order,
        "sample": spec.sample,
        "display_name": spec.display_name,
        "protocol": spec.protocol,
        "group": spec.group,
        "metric": metric,
        "value": value,
        "unit": unit,
        "source_file": str(source_file),
        "source_type": source_type,
        "definition": definition,
    }


def load_metrics(repo_root: Path, run_root: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    snakefile = repo_root / "main" / "Snakefile"
    definitions = {item[0]: item[5] for item in METRICS}

    for order, spec in enumerate(SAMPLES, start=1):
        final_dir = run_root / "final_output" / spec.sample
        curated = final_dir / "multiqc_report_data" / "multiqc_curated_metrics.txt"
        cg_idx = final_dir / f"{spec.sample}.cg.idx"
        if not curated.is_file():
            raise FileNotFoundError(f"Missing curated metrics for selected sample: {curated}")
        if not cg_idx.is_file():
            raise FileNotFoundError(f"Missing final CpG index for selected sample: {cg_idx}")

        data = read_curated(curated, spec.sample)
        non_cpg = np.array(
            [
                data.get("CA_retention_biscuitqc", np.nan),
                data.get("CC_retention_biscuitqc", np.nan),
                data.get("CT_retention_biscuitqc", np.nan),
            ],
            dtype=float,
        )
        if not np.isfinite(non_cpg).all():
            raise ValueError(f"Incomplete non-CpG retention metrics for {spec.sample}")

        values = {
            "prefilter_primary_reads": data.get("num_reads_merged", np.nan) * 1e6,
            "mapping_rate": data.get("mapping_rate_merged", np.nan) * 100,
            "conversion_efficiency": 100 * (1 - float(non_cpg.mean())),
            "indexed_pixels": float(sum(1 for line in cg_idx.open() if line.strip())),
            "covered_cpgs": data.get("cpg_coverage_yame", np.nan) * 1e6,
            "cpg_depth": data.get("cpg_depth_yame", np.nan),
            "global_mcg": data.get("CG_retention_biscuitqc", np.nan) * 100,
        }

        # This panel's label is explicitly a primary-record retention ratio.
        # Do not substitute the bsconv log's processed/remaining ratio: that
        # log can count non-primary records and therefore need not equal the
        # flagstat primary-record ratio shown here.
        prefilter_primary = data.get("num_reads_merged", np.nan)
        qc_primary = data.get("num_reads_qc_bam", np.nan)
        if np.isfinite(prefilter_primary) and np.isfinite(qc_primary) and prefilter_primary > 0:
            values["qc_retained_fraction"] = 100 * qc_primary / prefilter_primary
            retained_source = curated
            retained_type = "unified_pipeline_final_flagstat_primary"
        elif spec.protocol == "spatialdmt":
            # Older Spatial-DMT final tables predate num_reads_qc_bam. Its
            # QC BAM is the same file as the pre-filter BAM by definition.
            values["qc_retained_fraction"] = 100.0
            retained_source = snakefile
            retained_type = "derived_pipeline_definition"
        else:
            raise ValueError(f"Incomplete QC-BAM retention counts for {spec.sample}")

        for metric, value in values.items():
            source = cg_idx if metric == "indexed_pixels" else curated
            source_type = "unified_pipeline_final"
            definition = definitions[metric]
            if metric == "qc_retained_fraction":
                source = retained_source
                source_type = retained_type
                if retained_type == "derived_pipeline_definition":
                    definition += "; Spatial-DMT QC BAM equals pre-filter BAM (Snakefile lines 180-185)"
            rows.append(
                metric_row(spec, order, metric, value, source, source_type, definition)
            )

    metrics = pd.DataFrame(rows)
    expected = len(SAMPLES) * len(METRICS)
    if len(metrics) != expected:
        raise AssertionError(f"Expected {expected} metric rows, found {len(metrics)}")
    if metrics.duplicated(["sample", "metric"]).any():
        raise ValueError("Duplicate sample/metric rows in unified Panel A")
    if not np.isfinite(metrics["value"]).all() or (metrics["value"] < 0).any():
        bad = metrics.loc[~np.isfinite(metrics["value"]) | (metrics["value"] < 0)]
        raise ValueError(f"Invalid plotted values:\n{bad}")
    return metrics


def build_audit(repo_root: Path, run_root: Path) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "item": "duplication_rate",
                "status": "not_plotted",
                "reason": (
                    "Archived Spatial-DMT curated rows undercount paired-end duplicates because the "
                    "parser omitted the dupsifter 'both reads marked as duplicates' field. The parser "
                    "is corrected for future reports, but cleaned shallow runs lack raw dupstat files."
                ),
                "snakefile_location": (
                    "mqc_curated consumes merged_dupstat: main/Snakefile:926-961; parser: "
                    "tools/write_curated_metrics_mqc.py:173-198"
                ),
                "reporting_gap": (
                    "copy_final_sample_outputs does not copy *_merged_dedup.dupsifter.stat: "
                    "main/Snakefile:1133-1157"
                ),
            },
            {
                "item": "raw_read_pairs_and_structure_pass",
                "status": "not_plotted",
                "reason": (
                    "Cleaned final archives do not retain standalone trim-stat files for all SmC-seq runs; "
                    "older MultiQC general-stat tables can fall back to post-bsconv spatial counts."
                ),
                "snakefile_location": "qc_report inputs/write_stats_mqc: main/Snakefile:1022-1068",
                "reporting_gap": (
                    "copy_final_sample_outputs copies CG, HTML, and MultiQC data but not the trim-stat "
                    "or stats TSV as explicit final deliverables: main/Snakefile:1133-1157"
                ),
            },
            {
                "item": "median_usable_reads_per_pixel",
                "status": "not_plotted",
                "reason": "Not present in the canonical curated table for all archived samples.",
                "snakefile_location": "spatial_summary creates the source table: main/Snakefile:975-994",
                "reporting_gap": (
                    "copy_final_sample_outputs does not copy qc/table/*_spatial_barcode_counts.tsv: "
                    "main/Snakefile:1133-1157"
                ),
            },
            {
                "item": "E65_20um",
                "status": "excluded",
                "reason": (
                    "No unified final curated metrics exist; final_output contains only the separate "
                    "original-reproduction result."
                ),
                "snakefile_location": "mqc_curated: main/Snakefile:926-961",
                "reporting_gap": str(
                    run_root / "final_output" / "E65_20um" / "multiqc_report_data" / "multiqc_curated_metrics.txt"
                ),
            },
            {
                "item": "E75_20um",
                "status": "excluded",
                "reason": (
                    "Archived curated row lacks bsconv_retained_fraction and reports 100% pre-filter mapping; "
                    "it is not a complete output of the current SmC-seq branch."
                ),
                "snakefile_location": (
                    "dna_bsconv_filter: main/Snakefile:357-393; mqc_curated bsconv input: "
                    "main/Snakefile:943-961"
                ),
                "reporting_gap": str(
                    run_root / "final_output" / "E75_20um" / "multiqc_report_data" / "multiqc_curated_metrics.txt"
                ),
            },
        ]
    )


def format_value(metric: str, value: float) -> str:
    if metric in {
        "mapping_rate",
        "qc_retained_fraction",
        "conversion_efficiency",
        "global_mcg",
    }:
        return f"{value:.1f}"
    if metric == "cpg_depth":
        return f"{value:.1f}x"
    if value >= 1e9:
        return f"{value / 1e9:.2f}B"
    if value >= 1e6:
        return f"{value / 1e6:.1f}M"
    if value >= 1e3:
        return f"{value / 1e3:.1f}k"
    return f"{value:.0f}"


def bar_fraction(values: np.ndarray, kind: str) -> np.ndarray:
    result = np.zeros(values.shape, dtype=float)
    if kind == "percent":
        return np.clip(values / 100.0, 0, 1)
    maximum = float(np.max(values))
    if maximum > 0:
        result = np.log10(values + 1) / np.log10(maximum + 1)
    return result


def build_plot(metrics: pd.DataFrame) -> plt.Figure:
    metadata = (
        metrics[["row_order", "sample", "display_name", "protocol", "group"]]
        .drop_duplicates("sample")
        .sort_values("row_order")
        .reset_index(drop=True)
    )
    metric_names = [item[0] for item in METRICS]
    matrix = (
        metrics.pivot(index="sample", columns="metric", values="value")
        .reindex(index=metadata["sample"], columns=metric_names)
        .to_numpy(dtype=float)
    )
    fractions = np.column_stack(
        [bar_fraction(matrix[:, j], METRICS[j][2]) for j in range(len(METRICS))]
    )

    apply_style()
    fig, ax = plt.subplots(figsize=(13.2, 8.1))
    nrows, ncols = matrix.shape
    label_left = -3.45
    ax.set_xlim(label_left, ncols - 0.44)
    ax.set_ylim(nrows - 0.05, -2.15)
    ax.axis("off")

    ax.text(label_left, -1.92, "A", fontsize=12, weight="bold", ha="left", va="center")
    ax.text(
        label_left + 0.35,
        -1.92,
        "Unified-pipeline QC across Spatial-DMT and SmC-seq",
        fontsize=10.5,
        weight="bold",
        ha="left",
        va="center",
    )

    for title, start, end in GROUP_HEADERS:
        center = (start + end) / 2
        ax.text(center, -1.46, title, ha="center", va="center", fontsize=7.3, color="#555555")
        ax.plot([start - 0.43, end + 0.43], [-1.27, -1.27], color="#AFAFAF", lw=0.7)
    for j, (_, label, _, _, _, _) in enumerate(METRICS):
        ax.text(j, -0.91, label, ha="center", va="center", fontsize=7.1, weight="bold")

    previous_group = None
    for i, row in metadata.iterrows():
        if row["group"] != previous_group:
            if i > 0:
                ax.plot(
                    [label_left + 1.65, ncols - 0.5],
                    [i - 0.5, i - 0.5],
                    color="#8F8F8F",
                    lw=0.8,
                )
            ax.text(
                label_left,
                i - 0.5,
                row["group"],
                ha="left",
                va="center",
                fontsize=6.7,
                weight="bold",
                color="#555555",
            )
            previous_group = row["group"]

        spec = next(item for item in SAMPLES if item.sample == row["sample"])
        ax.add_patch(
            Rectangle(
                (label_left, i - 0.29),
                0.13,
                0.58,
                facecolor=spec.color,
                edgecolor="none",
            )
        )
        ax.text(
            label_left + 0.23,
            i,
            row["display_name"],
            ha="left",
            va="center",
            fontsize=7.0,
            weight="bold" if row["protocol"] == "spatialdmt" else "normal",
        )
        ax.text(
            -0.57,
            i,
            row["protocol"],
            ha="right",
            va="center",
            fontsize=5.8,
            color="#777777",
        )

        for j, (metric, _, _, color, _, _) in enumerate(METRICS):
            ax.add_patch(
                Rectangle(
                    (j - 0.44, i - 0.29),
                    0.88,
                    0.58,
                    facecolor="white",
                    edgecolor="#D9D9D9",
                    linewidth=0.5,
                    zorder=-2,
                )
            )
            ax.add_patch(
                Rectangle(
                    (j - 0.44, i - 0.29),
                    0.88 * fractions[i, j],
                    0.58,
                    facecolor=color,
                    alpha=0.27,
                    edgecolor="none",
                    zorder=-1,
                )
            )
            ax.text(
                j,
                i,
                format_value(metric, matrix[i, j]),
                ha="center",
                va="center",
                fontsize=6.7,
                color="#202020",
            )

    ax.text(
        label_left,
        nrows - 0.02,
        "All cells use final outputs from the same protocol-aware Snakefile; no literature values or legacy files are used.",
        ha="left",
        va="top",
        fontsize=6.1,
        color="#555555",
    )
    ax.text(
        ncols - 0.5,
        nrows - 0.02,
        "Percent bars use 0-100%; count and depth bars use log10 scaling within each column.",
        ha="right",
        va="top",
        fontsize=6.1,
        color="#555555",
    )
    return fig


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    run_root = args.run_root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else repo_root / "results" / "figure2"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics = load_metrics(repo_root, run_root)
    audit = build_audit(repo_root, run_root)
    metrics_path = output_dir / "panel_a_unified_pipeline_metrics.tsv"
    audit_path = output_dir / "panel_a_unified_pipeline_audit.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    audit.to_csv(audit_path, sep="\t", index=False)

    fig = build_plot(metrics)
    stem = output_dir / "panel_a_unified_pipeline"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    n_spatial = sum(spec.protocol == "spatialdmt" for spec in SAMPLES)
    n_smc = sum(spec.protocol == "smcseq" for spec in SAMPLES)
    print(f"Panel A: {n_spatial} Spatial-DMT + {n_smc} SmC-seq unified-pipeline runs")
    print(f"Metrics: {len(METRICS)} columns, {len(metrics)} sourced values")
    print(f"Excluded/incomplete items: {len(audit)}; see {audit_path}")
    print(f"Wrote: {stem}.pdf/.svg/.png")


if __name__ == "__main__":
    main()
