#!/usr/bin/env python3
"""Plot Spatial-DMT read counts beside per-pixel CpG methylation."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from plot_style import apply_style
from spatial_utils import spatial_matrix


SAMPLES = (
    "SpMETSLE14DM",
    "SpMETSLE17DM",
    "SpMETSLE21DM",
    "SpMETSLE23DM",
    "SpMETSLB14DM",
)


def methylation_matrix(frame: pd.DataFrame, threshold: int) -> np.ndarray:
    x = pd.to_numeric(frame["x"], errors="raise").astype(int)
    y = pd.to_numeric(frame["y"], errors="raise").astype(int)
    covered = pd.to_numeric(frame["covered_cpgs"], errors="coerce")
    beta = pd.to_numeric(frame["cpg_methylation"], errors="coerce")
    matrix = np.full((int(y.max()), int(x.max())), np.nan, dtype=float)
    keep = covered.ge(threshold) & beta.between(0.0, 1.0)
    matrix[y[keep] - 1, x[keep] - 1] = beta[keep]
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-root",
        type=Path,
        default=Path("/scr1/users/suc3/20260528_runNewSpatialPipe"),
    )
    parser.add_argument(
        "--methylation-dir",
        type=Path,
        default=Path("results/figure2/spatial_methylation_tables"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/manuscript_figures/final_plots/Figure2_spmet_count_methylation_v1.png"),
    )
    parser.add_argument("--min-covered-cpgs", type=int, default=100)
    parser.add_argument("--dpi", type=int, default=360)
    args = parser.parse_args()

    count_matrices: list[np.ndarray] = []
    methylation_matrices: list[np.ndarray] = []
    for sample in SAMPLES:
        count_path = (
            args.run_root
            / "pipeline_output"
            / sample
            / "qc"
            / "table"
            / f"{sample}_spatial_barcode_counts.tsv"
        )
        methylation_path = args.methylation_dir / f"{sample}_spatial_cpg_methylation.tsv"
        count_frame = pd.read_csv(count_path, sep="\t")
        methylation_frame = pd.read_csv(methylation_path, sep="\t")
        count_matrices.append(np.log10(spatial_matrix(count_frame, "dna_reads") + 1.0))
        methylation_matrices.append(
            methylation_matrix(methylation_frame, args.min_covered_cpgs)
        )

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
        figsize=(8.2, 17.0),
        gridspec_kw={"left": 0.18, "right": 0.96, "top": 0.935, "bottom": 0.09,
                     "hspace": 0.18, "wspace": 0.13},
    )
    count_cmap = plt.get_cmap("YlGnBu").copy()
    methylation_cmap = plt.get_cmap("magma").copy()
    methylation_cmap.set_bad("#E6E6E6")

    count_image = methylation_image = None
    for row, (sample, count_matrix, methylation_matrix_) in enumerate(
        zip(SAMPLES, count_matrices, methylation_matrices)
    ):
        count_ax, methylation_ax = axes[row]
        count_image = count_ax.imshow(
            count_matrix,
            origin="upper",
            interpolation="nearest",
            cmap=count_cmap,
            norm=count_norm,
            aspect="equal",
        )
        methylation_image = methylation_ax.imshow(
            methylation_matrix_,
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
        for ax in (count_ax, methylation_ax):
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color("#333333")
                spine.set_linewidth(0.55)

    axes[0, 0].set_title("DNA read count", fontsize=10.5, fontweight="bold", pad=8)
    axes[0, 1].set_title("CpG methylation level", fontsize=10.5, fontweight="bold", pad=8)
    fig.suptitle(
        "Spatial read depth and CpG methylation in Spatial-DMT samples",
        x=0.18,
        y=0.978,
        ha="left",
        fontsize=13,
        fontweight="bold",
    )

    count_cax = fig.add_axes([0.205, 0.045, 0.32, 0.015])
    methylation_cax = fig.add_axes([0.62, 0.045, 0.32, 0.015])
    count_bar = fig.colorbar(count_image, cax=count_cax, orientation="horizontal")
    methylation_bar = fig.colorbar(
        methylation_image, cax=methylation_cax, orientation="horizontal"
    )
    count_bar.set_label("log10(DNA reads per pixel + 1)", fontsize=7.5)
    methylation_bar.set_label("CpG methylation level (β)", fontsize=7.5)
    for bar in (count_bar, methylation_bar):
        bar.ax.tick_params(labelsize=6.5, length=2)

    fig.text(
        0.18,
        0.018,
        (
            "Shared robust scales use pooled 2nd–98th percentiles across the five samples; "
            f"methylation pixels require ≥{args.min_covered_cpgs} covered CpGs."
        ),
        fontsize=6.5,
        color="#555555",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, facecolor="white")
    plt.close(fig)
    print(f"Wrote {args.output}")
    print(
        "Methylation display range: "
        f"{methylation_limits[0]:.4f}–{methylation_limits[1]:.4f}"
    )


if __name__ == "__main__":
    main()
