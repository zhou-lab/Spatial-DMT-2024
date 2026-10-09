#!/usr/bin/env python3
"""Plot validated SmC-seq read counts beside corrected methylation maps.

The source joins are produced by ``compare_smcseq_paper.py``. They use barcode
sequence as the join key, the authors' sample-specific physical A/B positions,
and call-weighted methylation. Only samples whose corrected maps were validated
against the published source data are included.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from plot_style import apply_style


SAMPLES = (
    "E55_10um_rep1",
    "E55_10um_rep2",
    "E75_10um_rep1",
    "E85_20um",
)


def physical_matrix(frame: pd.DataFrame, value: str) -> np.ndarray:
    """Map values to the authors' physical A/B grid and reverse both axes."""
    n_a = int(frame["paper_A"].max())
    n_b = int(frame["paper_B"].max())
    matrix = np.full((n_a, n_b), np.nan, dtype=float)
    matrix[
        frame["paper_A"].to_numpy(dtype=int) - 1,
        frame["paper_B"].to_numpy(dtype=int) - 1,
    ] = pd.to_numeric(frame[value], errors="raise").to_numpy(dtype=float)
    return matrix[::-1, ::-1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--comparison-dir",
        type=Path,
        default=Path("results/figure2/smcseq_paper_comparison"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "results/manuscript_figures/final_plots/"
            "Figure2_smcseq_count_methylation_corrected_v2.png"
        ),
    )
    parser.add_argument("--dpi", type=int, default=360)
    args = parser.parse_args()

    audit = pd.read_csv(args.comparison_dir / "comparison_metrics.tsv", sep="\t")
    validated = audit.loc[
        audit["conclusion"].eq("Strong matched-pixel agreement"), "sample"
    ].tolist()
    if validated != list(SAMPLES):
        raise ValueError(
            "Validated sample set changed; review comparison_metrics.tsv before plotting: "
            f"{validated}"
        )

    frames: list[pd.DataFrame] = []
    count_matrices: list[np.ndarray] = []
    methylation_matrices: list[np.ndarray] = []
    for sample in SAMPLES:
        frame = pd.read_csv(
            args.comparison_dir / f"{sample}_matched.tsv",
            sep="\t",
        )
        frame = frame.loc[frame["retained"].astype(bool)].copy()
        if frame.empty:
            raise ValueError(f"No retained corrected pixels for {sample}")
        frames.append(frame)
        count_matrices.append(np.log10(physical_matrix(frame, "Raw_reads") + 1.0))
        methylation_matrices.append(physical_matrix(frame, "call_weighted_beta"))

    pooled_counts = np.concatenate([m[np.isfinite(m)] for m in count_matrices])
    pooled_methylation = np.concatenate(
        [m[np.isfinite(m)] for m in methylation_matrices]
    )
    count_limits = np.quantile(pooled_counts, [0.02, 0.98])
    methylation_limits = np.quantile(pooled_methylation, [0.02, 0.98])
    count_norm = Normalize(*count_limits, clip=True)
    methylation_norm = Normalize(*methylation_limits, clip=True)

    apply_style()
    fig, axes = plt.subplots(
        len(SAMPLES),
        2,
        figsize=(8.3, 13.3),
        gridspec_kw={
            "left": 0.20,
            "right": 0.96,
            "top": 0.925,
            "bottom": 0.105,
            "hspace": 0.25,
            "wspace": 0.18,
        },
    )
    count_cmap = plt.get_cmap("YlGnBu").copy()
    methylation_cmap = plt.get_cmap("magma").copy()
    count_cmap.set_bad("#E8E8E8")
    methylation_cmap.set_bad("#E8E8E8")

    count_image = methylation_image = None
    for row, (sample, frame, count_matrix, methylation_matrix) in enumerate(
        zip(SAMPLES, frames, count_matrices, methylation_matrices)
    ):
        count_ax, methylation_ax = axes[row]
        count_image = count_ax.imshow(
            np.ma.masked_invalid(count_matrix),
            origin="upper",
            interpolation="nearest",
            cmap=count_cmap,
            norm=count_norm,
            aspect="equal",
        )
        methylation_image = methylation_ax.imshow(
            np.ma.masked_invalid(methylation_matrix),
            origin="upper",
            interpolation="nearest",
            cmap=methylation_cmap,
            norm=methylation_norm,
            aspect="equal",
        )
        count_ax.set_ylabel(
            sample,
            rotation=0,
            ha="right",
            va="center",
            labelpad=13,
            fontsize=8.4,
            fontweight="bold",
        )
        count_ax.text(
            0.0,
            1.025,
            f"n={len(frame):,} validated pixels",
            transform=count_ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=6.3,
            color="#555555",
        )
        for ax in (count_ax, methylation_ax):
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color("#333333")
                spine.set_linewidth(0.55)

    axes[0, 0].set_title(
        "Published raw read count", fontsize=10.3, fontweight="bold", pad=20
    )
    axes[0, 1].set_title(
        "Current CpG methylation", fontsize=10.3, fontweight="bold", pad=20
    )
    fig.suptitle(
        "SmC-seq spatial read depth and methylation after barcode correction",
        x=0.50,
        y=0.972,
        ha="center",
        fontsize=11.7,
        fontweight="bold",
    )

    count_cax = fig.add_axes([0.215, 0.057, 0.31, 0.014])
    methylation_cax = fig.add_axes([0.625, 0.057, 0.31, 0.014])
    count_bar = fig.colorbar(count_image, cax=count_cax, orientation="horizontal")
    methylation_bar = fig.colorbar(
        methylation_image,
        cax=methylation_cax,
        orientation="horizontal",
    )
    count_bar.set_label("log10(raw reads per pixel + 1)", fontsize=7.3)
    methylation_bar.set_label("Call-weighted CpG methylation (β)", fontsize=7.3)
    for bar in (count_bar, methylation_bar):
        bar.ax.tick_params(labelsize=6.3, length=2)

    fig.text(
        0.20,
        0.018,
        (
            "Barcode sequences were joined to the authors’ sample-specific A/B coordinates; "
            "both axes follow the published reversed orientation. Shared scales use pooled "
            "2nd–98th percentiles. Gray denotes pixels absent from the validated join."
        ),
        fontsize=6.2,
        color="#555555",
        wrap=True,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, facecolor="white")
    plt.close(fig)
    print(f"Wrote {args.output}")
    print(
        f"Count range: {count_limits[0]:.4f}–{count_limits[1]:.4f}; "
        f"methylation range: {methylation_limits[0]:.4f}–{methylation_limits[1]:.4f}"
    )


if __name__ == "__main__":
    main()
