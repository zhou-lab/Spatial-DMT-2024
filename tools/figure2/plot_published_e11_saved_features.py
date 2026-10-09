#!/usr/bin/env python3
"""Plot raw and FigR-smoothed VMRs saved in the published E11 Seurat object."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle

from plot_style import apply_style


ROOT = Path(__file__).resolve().parents[2]


def median_ci(values: np.ndarray, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = np.asarray(values)
    draws = np.median(
        values[rng.integers(0, len(values), size=(10_000, len(values)))], axis=1
    )
    return tuple(np.quantile(draws, [0.025, 0.975]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--audit-dir",
        type=Path,
        default=ROOT / "results/figure2/published_e11_wnn_audit",
    )
    parser.add_argument(
        "--barcode-map",
        type=Path,
        default=ROOT / "barcodes/spatial_barcodes.txt",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "results/figure2/published_e11_wnn_audit/"
        "Published_E11_saved_VMR_smoothing.png",
    )
    parser.add_argument("--scar-x", type=int, default=4)
    args = parser.parse_args()

    values = pd.read_csv(args.audit_dir / "script_feature_values.tsv", sep="\t")
    neighbors = pd.read_csv(args.audit_dir / "weighted_nn_indices.tsv", sep="\t")
    neighbor_columns = [c for c in neighbors if c.startswith("neighbor_")]
    neighbor_indices = neighbors[neighbor_columns].to_numpy(dtype=int) - 1
    barcode = pd.read_csv(
        args.barcode_map,
        sep="\t",
        header=None,
        names=["barcode", "x", "y", "dna_barcode", "legacy_barcode"],
    )
    values["barcode"] = values.cell.str.split("-").str[0].str.split(".").str[0]
    values = values.merge(barcode[["barcode", "x", "y"]], on="barcode", validate="one_to_one")
    raw_columns = [c for c in values if c.endswith("__raw")]
    if not raw_columns:
        raise ValueError("No exported raw VMR columns found")

    apply_style()
    fig, axes = plt.subplots(
        len(raw_columns), 3, figsize=(9.7, 3.35 * len(raw_columns)), layout="constrained"
    )
    axes = np.atleast_2d(axes)
    metrics = []
    letters = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    for row, raw_col in enumerate(raw_columns):
        feature = raw_col.removesuffix("__raw")
        smooth_col = feature + "__smooth"
        raw = np.full((50, 50), np.nan)
        smooth = np.full((50, 50), np.nan)
        for x, y, raw_value, smooth_value in values[
            ["x", "y", raw_col, smooth_col]
        ].itertuples(index=False, name=None):
            iy, ix = int(y) - 1, int(x) - 1
            raw[iy, ix] = raw_value
            smooth[iy, ix] = smooth_value
        change = (smooth - raw) * 100
        pooled = np.r_[raw[np.isfinite(raw)], smooth[np.isfinite(smooth)]]
        low, high = np.quantile(pooled, [0.02, 0.98])
        delta_limit = np.quantile(np.abs(change[np.isfinite(change)]), 0.98)
        panels = [
            (raw, "Raw VMR fraction", "magma", Normalize(low, high), "methylation β"),
            (smooth, "Saved FigR/WNN mean", "magma", Normalize(low, high), "methylation β"),
            (
                change,
                "Smoothed − raw",
                "coolwarm",
                TwoSlopeNorm(vcenter=0, vmin=-delta_limit, vmax=delta_limit),
                "percentage points",
            ),
        ]
        for col, (matrix, title, cmap, norm, cbar_label) in enumerate(panels):
            ax = axes[row, col]
            cmap_obj = plt.get_cmap(cmap).copy()
            cmap_obj.set_bad("#EEEEEE")
            image = ax.imshow(matrix, origin="upper", interpolation="nearest", cmap=cmap_obj, norm=norm)
            ax.add_patch(
                Rectangle(
                    (args.scar_x - 1.5, -0.5), 1, 50, fill=False,
                    edgecolor="#00E5FF", linewidth=1.1
                )
            )
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(title, fontsize=9.5, fontweight="bold")
            ax.text(-0.11, 1.03, next(letters), transform=ax.transAxes,
                    fontsize=11, fontweight="bold")
            bar = fig.colorbar(image, ax=ax, orientation="horizontal", fraction=0.052, pad=0.035)
            bar.set_label(cbar_label, fontsize=6.5)
            bar.ax.tick_params(labelsize=6, length=2)
        chrom, start, end = feature.split(".")
        axes[row, 0].set_ylabel(
            f"{chrom}:{int(start):,}–{int(end):,}", fontsize=8, fontweight="bold"
        )

        triplets = []
        indexed = values.set_index(["x", "y"])
        for y in range(1, 51):
            keys = [(args.scar_x + dx, y) for dx in (-1, 0, 1)]
            if not all(key in indexed.index for key in keys):
                continue
            rr = [indexed.loc[key, raw_col] for key in keys]
            ss = [indexed.loc[key, smooth_col] for key in keys]
            triplets.append((rr[1] - (rr[0] + rr[2]) / 2, ss[1] - (ss[0] + ss[2]) / 2))
        triplets = np.asarray(triplets)
        raw_ci = median_ci(triplets[:, 0], 100 + row)
        smooth_ci = median_ci(triplets[:, 1], 200 + row)
        aligned = values.set_index("cell").loc[neighbors.cell]
        reproduced = aligned[raw_col].to_numpy()[neighbor_indices].mean(axis=1)
        reproduction_error = np.max(
            np.abs(reproduced - aligned[smooth_col].to_numpy())
        )
        metrics.append(
            {
                "feature": feature,
                "scar_triplets": len(triplets),
                "raw_scar_delta": np.median(triplets[:, 0]),
                "raw_ci_low": raw_ci[0],
                "raw_ci_high": raw_ci[1],
                "smooth_scar_delta": np.median(triplets[:, 1]),
                "smooth_ci_low": smooth_ci[0],
                "smooth_ci_high": smooth_ci[1],
                "raw_sd": np.std(values[raw_col]),
                "smooth_sd": np.std(values[smooth_col]),
                "max_abs_neighbor_mean_error": reproduction_error,
            }
        )

    fig.suptitle(
        "Exact VMR smoothing saved in the published E11 object",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.01,
        "Features are the two VMRs named in WNN_DNAmRNA_clustering.R. Cyan outline marks the raw-depth scar at column 4.",
        ha="center",
        fontsize=7,
        color="#555555",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=360, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(metrics).to_csv(
        args.audit_dir / "published_e11_saved_VMR_metrics.tsv", sep="\t", index=False
    )
    print(pd.DataFrame(metrics).to_string(index=False))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
