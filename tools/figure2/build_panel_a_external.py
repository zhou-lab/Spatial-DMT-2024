#!/usr/bin/env python3
"""Build an independent MultiQC-style comparison with published benchmarks."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from plot_style import CONTEXT_COLORS, apply_style  # noqa: E402


CURRENT_ORDER = [
    "SpMETSLE14D",
    "SpMETSLE16D",
    "SpMETSLE14DM",
    "SpMETSLB16D",
    "SpMETSLHSK03D",
    "SpMETSLHSK05D",
]

CURRENT_METRIC_MAP = {
    "raw_read_pairs": "input_pairs",
    "pass_fraction": "current_pass_fraction",
    "mapping_rate": "spatial_mapping_rate",
    "duplication_rate": "duplication_rate",
    "conversion_efficiency": "conversion_efficiency",
    "unique_cpgs": "cpg_yield",
    "median_unique_fragments": "usable_reads_per_unit",
}

EXTERNAL_METRIC_MAP = {
    "smc_assigned_reads": "input_pairs",
    "bulk_input_reads": "input_pairs",
    "smc_mapping_rate": "spatial_mapping_rate",
    "cpgs_per_unit": "cpg_yield",
    "total_cpg_sites": "cpg_yield",
}

CURRENT_CG_SOURCES = {
    "SpMETSLE14D": (None, False),
    "SpMETSLE16D": (
        Path("/mnt/isilon/zhou_lab/projects/20230531_Spatial_methseq/pileup")
        / "SpMETSLE16D" / "SpMETSLE16D.cg",
        True,
    ),
    "SpMETSLE14DM": (
        Path("/scr1/users/suc3/20260528_runNewSpatialPipe/pipeline_output")
        / "SpMETSLE14DM" / "pileup" / "SpMETSLE14DM.cg",
        False,
    ),
    "SpMETSLB16D": (
        Path("/scr1/users/suc3/20260704_runSpMET_latest_mm10/pipeline_output")
        / "SpMETSLB16D" / "pileup" / "SpMETSLB16D.cg",
        False,
    ),
    "SpMETSLHSK03D": (
        Path("/scr1/users/suc3/20260624_runNewSpatialPipe_hg38/pipeline_output")
        / "SpMETSLHSK03D" / "pileup" / "SpMETSLHSK03D.cg",
        False,
    ),
    "SpMETSLHSK05D": (
        Path("/scr1/users/suc3/20260706_runSpMETSLHSK_more_hg38/pipeline_output")
        / "SpMETSLHSK05D" / "pileup" / "SpMETSLHSK05D.cg",
        False,
    ),
}

METRICS = [
    ("input_pairs", "Input reads\nor pairs", "count", "#2C7FB8"),
    ("current_pass_fraction", "Current\npass %", "percent", "#41B6C4"),
    ("spatial_mapping_rate", "Spatial\nmap %", "percent", "#2CA25F"),
    ("duplication_rate", "Dup.\n%", "percent", "#DE2D26"),
    ("conversion_efficiency", "Conv. /\ndeam. %", "percent", "#31A354"),
    ("units_profiled", "Pixels, cells\nor samples", "count", "#BDBDBD"),
    ("usable_reads_per_unit", "Usable reads\n/pixel or cell", "count", "#636363"),
    ("cpg_yield", "CpG sites\n(total or /unit)", "count", "#756BB1"),
    ("global_mcg", "Global\nmCG %", "percent", "#E6AB02"),
    ("correlation_1mb", "Pearson r\n1 Mb", "correlation", "#1B9E77"),
    ("correlation_100kb", "Pearson r\n100 kb", "correlation", "#66A61E"),
    ("correlation_10kb", "Pearson r\n10 kb", "correlation", "#A6D854"),
    ("consistent_1kb_fraction", "Consistent\n1-kb bins %", "percent", "#7570B3"),
]

GROUP_HEADERS = [
    ("Input", 0, 0),
    ("Processing", 1, 4),
    ("Yield", 5, 7),
    ("Methylation validation", 8, 12),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def yame_median_mcg_rows(repo_root: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    cache_dir = repo_root / "results" / "figure2" / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    yame = shutil.which("yame")
    rows = []
    for sample in CURRENT_ORDER:
        source, is_fallback = CURRENT_CG_SOURCES[sample]
        suffix = "_legacy" if is_fallback else ""
        cache = cache_dir / f"{sample}_yame_summary{suffix}.tsv"
        if not cache.is_file():
            if source is None or not source.is_file():
                raise FileNotFoundError(f"No YAME summary or CpG source for {sample}")
            if yame is None:
                raise RuntimeError("yame is required to compute median per-pixel mCG")
            temporary = cache.with_suffix(cache.suffix + ".tmp")
            try:
                with temporary.open("w") as handle:
                    subprocess.run([yame, "summary", str(source)], stdout=handle, check=True)
                temporary.replace(cache)
            except BaseException:
                temporary.unlink(missing_ok=True)
                raise
        summary = pd.read_csv(cache, sep="\t")
        required = {"N_query", "Beta"}
        if summary.empty or not required <= set(summary):
            raise ValueError(f"Incomplete YAME summary cache: {cache}")
        beta = pd.to_numeric(summary["Beta"], errors="coerce")
        covered = pd.to_numeric(summary["N_query"], errors="coerce").gt(0)
        value = 100 * float(beta[covered & beta.notna()].median())
        one = metadata.loc[metadata["sample"].eq(sample)].iloc[0]
        rows.append(
            {
                "row_order": one["row_order"],
                "sample": sample,
                "display_name": sample,
                "group": one["group"],
                "species": one["species"],
                "tissue": one["tissue"],
                "unit_type": one["unit_type"],
                "metric": "global_mcg",
                "value": value,
                "unit": "% median pixel beta",
                "source_file": str(source if source is not None else cache),
                "definition": (
                    "legacy_fallback_median_per_pixel_yame_Beta_v1"
                    if is_fallback
                    else "revised_median_per_pixel_yame_Beta_v1"
                ),
            }
        )
    return pd.DataFrame(rows)


def load_current_metrics(repo_root: Path) -> pd.DataFrame:
    path = repo_root / "results" / "figure2" / "figure2_metrics.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame = frame[
        frame["sample"].isin(CURRENT_ORDER)
        & frame["metric"].isin(CURRENT_METRIC_MAP)
        & ~frame["is_fallback"].astype(bool)
    ].copy()
    frame["metric"] = frame["metric"].map(CURRENT_METRIC_MAP)
    frame["display_name"] = frame["sample"]
    frame["group"] = "Current revised workflow"
    frame["unit_type"] = "pixel"
    frame["row_order"] = frame["sample"].map(
        {sample: order for order, sample in enumerate(CURRENT_ORDER, start=1)}
    )
    frame["definition"] = frame["definition_version"]
    frame = frame[
        [
            "row_order",
            "sample",
            "display_name",
            "group",
            "species",
            "tissue",
            "unit_type",
            "metric",
            "value",
            "unit",
            "source_file",
            "definition",
        ]
    ]
    frame = pd.concat([frame, yame_median_mcg_rows(repo_root, frame)], ignore_index=True)
    frame.loc[frame["sample"].eq("SpMETSLE16D"), "display_name"] = "SpMETSLE16D*"
    return frame


def load_external_metrics(repo_root: Path) -> pd.DataFrame:
    path = repo_root / "tools" / "figure2" / "panel_a_external_metrics.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame["metric"] = frame["metric"].replace(EXTERNAL_METRIC_MAP)
    return frame


def format_value(metric: str, value: float) -> str:
    if metric in {
        "current_pass_fraction",
        "spatial_mapping_rate",
        "duplication_rate",
        "conversion_efficiency",
        "global_mcg",
        "consistent_1kb_fraction",
    }:
        return f"{value:.1f}"
    if metric.startswith("correlation_"):
        return f"{value:.3f}"
    if metric == "units_profiled":
        return f"{value / 1e3:.1f}k" if value >= 1e3 else f"{value:.0f}"
    if value >= 1e9:
        return f"{value / 1e9:.2f}B"
    if value >= 1e6:
        return f"{value / 1e6:.1f}M"
    if value >= 1e3:
        return f"{value / 1e3:.0f}k"
    return f"{value:g}"


def bar_fraction(values: np.ndarray, kind: str) -> np.ndarray:
    result = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    if not finite.any():
        return result
    if kind == "percent":
        result[finite] = np.clip(values[finite] / 100.0, 0, 1)
        return result
    if kind == "correlation":
        result[finite] = np.clip(values[finite], 0, 1)
        return result
    maximum = np.nanmax(values)
    if maximum > 0:
        result[finite] = np.log10(values[finite] + 1) / np.log10(maximum + 1)
    else:
        result[finite] = 0
    return result


def context_color(species: str, tissue: str) -> str:
    return CONTEXT_COLORS.get(f"{species} {tissue}", "#777777")


def build_plot(metrics: pd.DataFrame) -> plt.Figure:
    rows = (
        metrics[
            [
                "row_order",
                "sample",
                "display_name",
                "group",
                "species",
                "tissue",
                "unit_type",
            ]
        ]
        .drop_duplicates("sample")
        .sort_values("row_order")
        .reset_index(drop=True)
    )
    metric_names = [item[0] for item in METRICS]
    matrix = (
        metrics.pivot_table(index="sample", columns="metric", values="value", aggfunc="first")
        .reindex(index=rows["sample"], columns=metric_names)
        .to_numpy(dtype=float)
    )
    fractions = np.column_stack(
        [bar_fraction(matrix[:, j], METRICS[j][2]) for j in range(len(METRICS))]
    )

    apply_style()
    fig, ax = plt.subplots(figsize=(17.0, 9.8))
    nrows, ncols = matrix.shape
    label_left = -3.75
    ax.set_xlim(label_left, ncols - 0.45)
    ax.set_ylim(nrows - 0.15, -2.2)
    ax.axis("off")

    ax.text(
        label_left,
        -1.95,
        "A",
        fontsize=12,
        weight="bold",
        ha="left",
        va="center",
    )
    ax.text(
        label_left + 0.35,
        -1.95,
        "Current workflow, SmC-seq QC, and external methylation validation",
        fontsize=10.5,
        weight="bold",
        ha="left",
        va="center",
    )

    for title, start, end in GROUP_HEADERS:
        center = (start + end) / 2
        ax.text(center, -1.48, title, ha="center", va="center", fontsize=7.3, color="#555555")
        ax.plot([start - 0.43, end + 0.43], [-1.28, -1.28], color="#AFAFAF", lw=0.7)
    for j, (_, label, _, _) in enumerate(METRICS):
        ax.text(j, -0.94, label, ha="center", va="center", fontsize=7.0, weight="bold")

    previous_group = None
    for i, row in rows.iterrows():
        if row["group"] != previous_group:
            if i > 0:
                ax.plot(
                    [label_left + 1.75, ncols - 0.5],
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
                fontsize=6.5,
                weight="bold",
                color="#555555",
            )
            previous_group = row["group"]

        ax.add_patch(
            Rectangle(
                (label_left, i - 0.29),
                0.13,
                0.58,
                facecolor=context_color(row["species"], row["tissue"]),
                edgecolor="none",
            )
        )
        ax.text(
            label_left + 0.23,
            i,
            row["display_name"],
            ha="left",
            va="center",
            fontsize=7.1,
            weight="bold" if row["group"] == "Current revised workflow" else "normal",
        )
        ax.text(
            -0.58,
            i,
            row["unit_type"],
            ha="right",
            va="center",
            fontsize=5.8,
            color="#777777",
        )

        for j, (metric, _, _, color) in enumerate(METRICS):
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
            value = matrix[i, j]
            if not np.isfinite(value):
                continue
            ax.add_patch(
                Rectangle(
                    (j - 0.44, i - 0.29),
                    0.88 * fractions[i, j],
                    0.58,
                    facecolor=color,
                    alpha=0.25,
                    edgecolor="none",
                    zorder=-1,
                )
            )
            ax.text(
                j,
                i,
                format_value(metric, value),
                ha="center",
                va="center",
                fontsize=6.6,
                color="#202020",
            )

    ax.text(
        label_left,
        nrows - 0.02,
        "Blank cells are unavailable or not comparable. *16D median mCG uses a legacy .cg fallback.",
        ha="left",
        va="top",
        fontsize=6.1,
        color="#555555",
    )
    ax.text(
        ncols - 0.5,
        nrows - 0.02,
        "Input and CpG-yield cells retain row-specific source units; count bars use log10 scaling within each column.",
        ha="right",
        va="top",
        fontsize=6.1,
        color="#555555",
    )
    return fig


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else repo_root / "results" / "figure2"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    current = load_current_metrics(repo_root)
    external = load_external_metrics(repo_root)
    metrics = pd.concat([current, external], ignore_index=True)
    if metrics.duplicated(["sample", "metric"]).any():
        duplicates = metrics.loc[
            metrics.duplicated(["sample", "metric"], keep=False), ["sample", "metric"]
        ]
        raise ValueError(f"Duplicate sample/metric rows:\n{duplicates}")
    if metrics["value"].isna().any() or (metrics["value"] < 0).any():
        raise ValueError("Plotted metrics must be finite and non-negative")

    data_path = output_dir / "panel_a_current_external_metrics.tsv"
    metrics.sort_values(["row_order", "metric"]).to_csv(data_path, sep="\t", index=False)
    fig = build_plot(metrics)
    stem = output_dir / "panel_a_current_external"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Rows: {metrics['sample'].nunique()} (6 current, 10 published/validation)")
    print(f"Plotted values: {len(metrics)}; missing or incomparable cells: left blank")
    print(f"Wrote: {stem}.pdf/.svg/.png")
    print(f"Data: {data_path}")


if __name__ == "__main__":
    main()
