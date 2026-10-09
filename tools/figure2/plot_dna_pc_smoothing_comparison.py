#!/usr/bin/env python3
"""Compare DNA-PC neighbor smoothing in SLB14DM and E85."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from plot_style import apply_style


ROOT = Path(__file__).resolve().parents[2]


def grid(frame: pd.DataFrame, x: str, y: str, value: str, reverse: bool) -> np.ndarray:
    width, height = int(frame[x].max()), int(frame[y].max())
    matrix = np.full((height, width), np.nan)
    matrix[
        frame[y].to_numpy(dtype=int) - 1,
        frame[x].to_numpy(dtype=int) - 1,
    ] = frame[value].to_numpy(dtype=float)
    return matrix[::-1, ::-1] if reverse else matrix


def scar_metrics(
    frame: pd.DataFrame,
    axis: str,
    index: int,
    other: str,
    sample: str,
    label: str,
) -> dict[str, float | int | str]:
    lookup = frame.set_index([axis, other])
    deltas = []
    for position in sorted(frame[other].unique()):
        keys = [(index + offset, position) for offset in (-1, 0, 1)]
        if not all(key in lookup.index for key in keys):
            continue
        rows = lookup.loc[keys]
        adjacent_count = rows.iloc[[0, 2]].raw_reads.mean()
        adjacent_raw = rows.iloc[[0, 2]].raw_beta.mean()
        adjacent_smooth = rows.iloc[[0, 2]].smoothed_beta.mean()
        deltas.append(
            (
                rows.iloc[1].raw_reads / adjacent_count,
                rows.iloc[1].raw_beta - adjacent_raw,
                rows.iloc[1].smoothed_beta - adjacent_smooth,
            )
        )
    values = np.asarray(deltas)
    rng = np.random.default_rng(20260918 + index)
    bootstrap = np.median(
        values[rng.integers(0, len(values), size=(10_000, len(values)))], axis=1
    )
    count_ci = np.quantile(bootstrap[:, 0], [0.025, 0.975])
    raw_ci = np.quantile(bootstrap[:, 1], [0.025, 0.975])
    smooth_ci = np.quantile(bootstrap[:, 2], [0.025, 0.975])
    return {
        "sample": sample,
        "scar": label,
        "matched_triplets": len(values),
        "median_read_count_ratio": np.median(values[:, 0]),
        "read_count_ratio_ci_low": count_ci[0],
        "read_count_ratio_ci_high": count_ci[1],
        "median_raw_beta_delta": np.median(values[:, 1]),
        "raw_beta_delta_ci_low": raw_ci[0],
        "raw_beta_delta_ci_high": raw_ci[1],
        "median_smoothed_beta_delta": np.median(values[:, 2]),
        "smoothed_beta_delta_ci_low": smooth_ci[0],
        "smoothed_beta_delta_ci_high": smooth_ci[1],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results/figure2/wnn_smoothing_comparison",
    )
    parser.add_argument("--dpi", type=int, default=360)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    slb14 = pd.read_csv(
        args.output_dir / "slb14_dna_knn/pixel_values.tsv", sep="\t"
    )
    slb14["neighbor_count"] = 20

    e85 = pd.read_csv(args.output_dir / "e85_dna_knn/pixel_values.tsv", sep="\t").rename(
        columns={
            "Raw_reads": "raw_reads",
            "call_weighted_beta": "raw_beta",
        }
    )
    e85["neighbor_count"] = 20

    samples = [
        {
            "name": "SpMETSLB14DM",
            "subtitle": "P21 brain, 20 µm · DNA-only 10-PC methylation kNN",
            "frame": slb14,
            "x": "x",
            "y": "y",
            "reverse": False,
        },
        {
            "name": "E85_20um",
            "subtitle": "E8.5, 20 µm · DNA-only 10-PC methylation kNN",
            "frame": e85,
            "x": "paper_B",
            "y": "paper_A",
            "reverse": True,
        },
    ]

    metrics = [
        scar_metrics(slb14, "y", 19, "x", "SpMETSLB14DM", "row 19"),
        scar_metrics(e85, "paper_A", 22, "paper_B", "E85_20um", "row 22"),
        scar_metrics(e85, "paper_B", 14, "paper_A", "E85_20um", "column 14"),
    ]
    pd.DataFrame(metrics).to_csv(
        args.output_dir / "dna_knn_smoothing_scar_metrics.tsv", sep="\t", index=False
    )

    apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(10.7, 7.3), layout="constrained")
    count_cmap = plt.get_cmap("YlGnBu").copy()
    beta_cmap = plt.get_cmap("magma").copy()
    count_cmap.set_bad("#E9E9E9")
    beta_cmap.set_bad("#E9E9E9")
    titles = ["Raw DNA read count", "Raw global CpG β", "Neighbor-smoothed global CpG β"]
    letters = iter("ABCDEF")

    for row, specification in enumerate(samples):
        frame = specification["frame"]
        reverse = specification["reverse"]
        count_matrix = np.log10(
            grid(frame, specification["x"], specification["y"], "raw_reads", reverse) + 1
        )
        raw_matrix = grid(frame, specification["x"], specification["y"], "raw_beta", reverse)
        smooth_matrix = grid(
            frame, specification["x"], specification["y"], "smoothed_beta", reverse
        )
        count_limits = np.nanquantile(count_matrix, [0.02, 0.98])
        beta_values = np.r_[raw_matrix[np.isfinite(raw_matrix)], smooth_matrix[np.isfinite(smooth_matrix)]]
        beta_limits = np.quantile(beta_values, [0.02, 0.98])
        matrices = [count_matrix, raw_matrix, smooth_matrix]
        norms = [Normalize(*count_limits, clip=True), Normalize(*beta_limits, clip=True)]

        for col, (ax, matrix) in enumerate(zip(axes[row], matrices)):
            image = ax.imshow(
                np.ma.masked_invalid(matrix),
                origin="upper",
                interpolation="nearest",
                cmap=count_cmap if col == 0 else beta_cmap,
                norm=norms[0] if col == 0 else norms[1],
                aspect="equal",
            )
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(titles[col], fontsize=9.0, fontweight="bold")
            ax.text(-0.14, 1.03, next(letters), transform=ax.transAxes,
                    fontsize=11, fontweight="bold")
            bar = fig.colorbar(image, ax=ax, orientation="horizontal", fraction=0.05, pad=0.035)
            bar.ax.tick_params(labelsize=6, length=2)
            bar.set_label("log10(reads + 1)" if col == 0 else "methylation β", fontsize=6.8)
        axes[row, 0].set_ylabel(
            f"{specification['name']}\n{specification['subtitle']}",
            rotation=0,
            ha="right",
            va="center",
            labelpad=12,
            fontsize=7.5,
            fontweight="bold",
        )

    fig.suptitle(
        "Global CpG methylation before and after DNA-PC neighbor smoothing",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.012,
        (
            "Scales are robust and shared between raw/smoothed β within each sample. "
            "Both samples use 20-neighbor graphs built only from ten methylation PCs; "
            "E85 uses published raw counts and corrected paper geometry."
        ),
        ha="center",
        fontsize=6.7,
        color="#555555",
        wrap=True,
    )
    output = args.output_dir / "SpMETSLB14DM_E85_neighbor_smoothing.png"
    fig.savefig(output, dpi=args.dpi, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(pd.DataFrame(metrics).to_string(index=False))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
