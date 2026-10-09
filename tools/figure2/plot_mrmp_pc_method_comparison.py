#!/usr/bin/env python3
"""Compare raw beta, prior DNA-PC smoothing, and 5/20-PC MRMP smoothing."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from build_mrmp_smoothing import grid, scar_metric
from plot_style import apply_style


ROOT = Path(__file__).resolve().parents[2]


def load_slb(output_dir: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        output_dir / "SpMETSLB14DM_mrmp_pc_beta_values.tsv",
        sep="\t",
        dtype={"spatial_id": str},
    )
    previous = pd.read_csv(
        ROOT / "results/figure2/wnn_smoothing_comparison/slb14_dna_knn/pixel_values.tsv",
        sep="\t",
        dtype={"spatial_id": str},
    )[["spatial_id", "smoothed_beta"]].rename(columns={"smoothed_beta": "previous_pc_beta"})
    frame["spatial_id"] = frame.spatial_id.str.zfill(4)
    previous["spatial_id"] = previous.spatial_id.str.zfill(4)
    frame = frame.merge(previous, on="spatial_id", validate="one_to_one")
    frame["sample_name"] = "SpMETSLB14DM"
    return frame


def load_e85(output_dir: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        output_dir / "E85_20um_mrmp_pc_beta_values.tsv",
        sep="\t",
        dtype={"spatial_id": str},
    )
    previous = pd.read_csv(
        ROOT / "results/figure2/wnn_smoothing_comparison/e85_dna_knn/pixel_values.tsv",
        sep="\t",
    )[["record", "smoothed_beta"]].rename(columns={"smoothed_beta": "previous_pc_beta"})
    frame = frame.merge(previous, on="record", validate="one_to_one")
    frame["sample_name"] = "E85_20um"
    return frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results/figure2/mrmp_smoothing_comparison",
    )
    parser.add_argument("--dpi", type=int, default=360)
    args = parser.parse_args()

    slb = load_slb(args.output_dir)
    e85 = load_e85(args.output_dir)
    fields = [
        "raw_beta",
        "previous_pc_beta",
        "mrmp_pc5_beta_smoothed",
        "mrmp_pc20_beta_smoothed",
    ]
    titles = [
        "Raw global CpG β",
        "Previous DNA-PC smoothing",
        "MRMP smoothing · 5 PCs",
        "MRMP smoothing · 20 PCs",
    ]
    specifications = [
        {
            "frame": slb,
            "x": "x",
            "y": "y",
            "reverse": False,
            "label": "SpMETSLB14DM\nP21 brain, 20 µm",
        },
        {
            "frame": e85,
            "x": "paper_B",
            "y": "paper_A",
            "reverse": True,
            "label": "E85_20um\nE8.5 embryo, 20 µm",
        },
    ]

    apply_style()
    fig, axes = plt.subplots(2, 4, figsize=(13.2, 7.2), layout="constrained")
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#E9E9E9")
    letters = iter("ABCDEFGH")
    for row, specification in enumerate(specifications):
        frame = specification["frame"]
        matrices = [
            grid(frame, specification["x"], specification["y"], field, specification["reverse"])
            for field in fields
        ]
        pooled = np.concatenate([matrix[np.isfinite(matrix)] for matrix in matrices])
        limits = np.quantile(pooled, [0.02, 0.98])
        norm = Normalize(*limits, clip=True)
        for column, (matrix, title) in enumerate(zip(matrices, titles)):
            ax = axes[row, column]
            image = ax.imshow(
                np.ma.masked_invalid(matrix),
                origin="upper",
                interpolation="nearest",
                cmap=cmap,
                norm=norm,
                aspect="equal",
            )
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(title, fontsize=8.8, fontweight="bold")
            ax.text(-0.14, 1.03, next(letters), transform=ax.transAxes, fontsize=11, fontweight="bold")
            bar = fig.colorbar(image, ax=ax, orientation="horizontal", fraction=0.05, pad=0.035)
            bar.ax.tick_params(labelsize=6, length=2)
            bar.set_label("methylation β", fontsize=6.8)
        axes[row, 0].set_ylabel(
            specification["label"],
            rotation=0,
            ha="right",
            va="center",
            labelpad=12,
            fontsize=7.5,
            fontweight="bold",
        )
    fig.suptitle("Global CpG beta: previous DNA-PC versus MRMP-PC smoothing", fontsize=13, fontweight="bold")
    fig.text(
        0.5,
        -0.012,
        "All panels share one robust beta scale within each sample. MRMP graphs use 20 neighbors and "
        "coverage-qualified beta donors; no scar coordinates enter any smoother.",
        ha="center",
        fontsize=6.7,
        color="#555555",
    )
    output = args.output_dir / "SpMETSLB14DM_E85_raw_previous_mrmp_pc_comparison.png"
    fig.savefig(output, dpi=args.dpi, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    scar_rows = []
    for field in fields:
        scar_rows.append(scar_metric(slb, field, "y", 19, "x", "row 19"))
        scar_rows.append(scar_metric(e85, field, "paper_A", 22, "paper_B", "row 22"))
        scar_rows.append(scar_metric(e85, field, "paper_B", 14, "paper_A", "column 14"))
    pd.DataFrame(scar_rows).to_csv(
        args.output_dir / "raw_previous_mrmp_pc_comparison_scar_metrics.tsv",
        sep="\t",
        index=False,
    )
    print(pd.DataFrame(scar_rows).to_string(index=False))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
