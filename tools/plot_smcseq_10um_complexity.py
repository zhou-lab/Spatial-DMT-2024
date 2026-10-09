#!/usr/bin/env python3
"""Summarize the SmC-seq 10 um library-complexity occupancy experiment."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DISPLAY = {
    "SmC_E55_10um_rep1": "E5.5 rep. 1",
    "SmC_E55_10um_rep2": "E5.5 rep. 2",
    "SmC_E65_10um": "E6.5",
    "SmC_E75_10um_rep1": "E7.5 rep. 1",
    "SmC_E75_10um_rep2": "E7.5 rep. 2",
}


def component_table(results: Path, summary: pd.DataFrame, kind: str) -> pd.DataFrame:
    rows = []
    for _, sample in summary.iterrows():
        centers = np.fromstring(sample[f"{kind}_component_centers"], sep=",")
        ratios = centers / centers[0]
        assignments = pd.read_csv(
            results / f"{sample['sample']}_library_complexity_assignments.tsv",
            sep="\t",
            usecols=["foreground", f"{kind}_component"],
        )
        foreground = assignments.loc[assignments["foreground"].eq(1)]
        counts = foreground[f"{kind}_component"].value_counts()
        for component, ratio in enumerate(ratios, start=1):
            rows.append(
                {
                    "sample": sample["sample"],
                    "component": component,
                    "ratio": ratio,
                    "pixels": int(counts.get(component, 0)),
                    "fraction": float(counts.get(component, 0) / len(foreground)),
                }
            )
    return pd.DataFrame(rows)


def draw_components(axis, components: pd.DataFrame, order: list[str], title: str) -> None:
    colors = plt.get_cmap("tab10")
    for y, sample in enumerate(order):
        one = components.loc[components["sample"].eq(sample)]
        axis.scatter(
            np.log2(one["ratio"]),
            np.full(len(one), y),
            s=35 + 500 * one["fraction"],
            c=[colors((component - 1) % 10) for component in one["component"]],
            edgecolor="white",
            linewidth=0.8,
            alpha=0.9,
        )
        for _, point in one.loc[one["fraction"].lt(0.01)].iterrows():
            axis.annotate(
                f"{point['ratio']:.0f}x ({point['pixels']} px)",
                (np.log2(point["ratio"]), y),
                xytext=(4, 5),
                textcoords="offset points",
                fontsize=7,
            )
    for ratio, style in [(1, "-"), (2, "--"), (4, "--")]:
        axis.axvline(np.log2(ratio), color="#333333", linestyle=style, linewidth=1)
        axis.text(np.log2(ratio), -0.65, f"{ratio}x", ha="center", va="bottom", fontsize=8)
    axis.set_yticks(range(len(order)), [DISPLAY[sample] for sample in order])
    axis.set_xlim(-0.1, max(5.7, np.log2(components["ratio"].max()) + 0.35))
    ticks = np.arange(0, 6)
    axis.set_xticks(ticks, [f"{2**tick:g}x" for tick in ticks])
    axis.set_xlabel("Component center relative to lowest center")
    axis.set_title(title)
    axis.grid(axis="x", color="#dddddd", linewidth=0.7)
    axis.invert_yaxis()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    summary = pd.read_csv(
        args.results / "library_complexity_occupancy_summary.tsv", sep="\t"
    )
    order = [sample for sample in DISPLAY if sample in set(summary["sample"])]
    summary = summary.set_index("sample").loc[order].reset_index()
    raw = component_table(args.results, summary, "raw")
    normalized = component_table(args.results, summary, "normalized")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    figure, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    draw_components(axes[0, 0], raw, order, "Raw mixture centers")
    draw_components(axes[0, 1], normalized, order, "Row/column-normalized centers")

    y = np.arange(len(summary))
    width = 0.36
    axes[1, 0].barh(
        y - width / 2,
        summary["raw_row_column_adjusted_r2"],
        height=width,
        label="Raw",
        color="#2c7fb8",
    )
    axes[1, 0].barh(
        y + width / 2,
        summary["normalized_row_column_adjusted_r2"],
        height=width,
        label="Normalized",
        color="#41ab5d",
    )
    axes[1, 0].axvline(0, color="#333333", linewidth=0.8)
    axes[1, 0].set_yticks(y, [DISPLAY[sample] for sample in order])
    axes[1, 0].set_xlabel("Adjusted R2 explained by row and column")
    axes[1, 0].set_title("Spatial grid artifact")
    axes[1, 0].legend(frameon=False)
    axes[1, 0].invert_yaxis()

    observed = 100 * summary["median_observed_library_fraction"]
    axes[1, 1].barh(y, observed, color="#756bb1")
    for index, (value, foreground) in enumerate(
        zip(observed, summary["foreground_pixels"])
    ):
        axes[1, 1].text(
            min(value + 1.5, 96),
            index,
            f"{value:.0f}% observed; {foreground:,} px",
            va="center",
            fontsize=8,
        )
    axes[1, 1].set_yticks(y, [DISPLAY[sample] for sample in order])
    axes[1, 1].set_xlim(0, 105)
    axes[1, 1].set_xlabel("Median observed fraction of estimated library (%)")
    axes[1, 1].set_title("Library saturation and analyzed pixels")
    axes[1, 1].invert_yaxis()

    figure.suptitle(
        "Published SmC-seq 10 um pixels: duplicate-extrapolated library complexity",
        fontsize=13,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220, bbox_inches="tight")
    figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
