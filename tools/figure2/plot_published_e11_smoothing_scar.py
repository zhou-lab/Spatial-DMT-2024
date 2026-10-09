#!/usr/bin/env python3
"""Audit the published E11 WNN smoother against the column-4 depth scar.

The WNN neighbor indices come from the GEO-deposited Seurat object.  The
displayed global CpG beta and raw DNA counts come from the current reprocessing
of the same E11#1 library (SpMETSLE14DM).  FigR's smoothScoresNN operation is an
unweighted row mean over the saved neighbor indices; that operation is
implemented directly here for one per-pixel methylation summary.
"""

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


def bootstrap_median_ci(values: np.ndarray, seed: int = 20260918) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    draws = np.median(
        values[rng.integers(0, len(values), size=(10_000, len(values)))], axis=1
    )
    return tuple(np.quantile(draws, [0.025, 0.975]))


def add_scar_outline(ax: plt.Axes, scar_x: int, size: int = 50) -> None:
    ax.add_patch(
        Rectangle(
            (scar_x - 1 - 0.5, -0.5),
            1,
            size,
            fill=False,
            edgecolor="#00E5FF",
            linewidth=1.15,
            zorder=5,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--audit-dir",
        type=Path,
        default=ROOT / "results/figure2/published_e11_wnn_audit",
    )
    parser.add_argument(
        "--counts",
        type=Path,
        default=Path(
            "/scr1/users/suc3/20260528_runNewSpatialPipe/pipeline_output/"
            "SpMETSLE14DM/qc/table/SpMETSLE14DM_spatial_barcode_counts.tsv"
        ),
    )
    parser.add_argument(
        "--methylation",
        type=Path,
        default=ROOT
        / "results/figure2/spatial_methylation_tables/"
        "SpMETSLE14DM_spatial_cpg_methylation.tsv",
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
        "Published_E11_WNN_scar_audit.png",
    )
    parser.add_argument("--scar-x", type=int, default=4)
    args = parser.parse_args()

    nn = pd.read_csv(args.audit_dir / "weighted_nn_indices.tsv", sep="\t")
    neighbor_columns = [c for c in nn if c.startswith("neighbor_")]
    if not neighbor_columns:
        raise ValueError("No neighbor columns in extracted WNN table")

    barcode = pd.read_csv(
        args.barcode_map,
        sep="\t",
        header=None,
        names=["rna_barcode", "x", "y", "dna_barcode", "legacy_barcode"],
        dtype={"rna_barcode": str, "x": int, "y": int},
    )
    barcode_by_cell = barcode.set_index("rna_barcode")

    def normalize_cell(value: str) -> str:
        return str(value).split("-")[0].split(".")[0]

    nn["barcode"] = nn["cell"].map(normalize_cell)
    missing_barcodes = sorted(set(nn["barcode"]) - set(barcode_by_cell.index))
    if missing_barcodes:
        raise ValueError(
            f"{len(missing_barcodes)} deposited cells do not match the project barcode map; "
            f"examples: {missing_barcodes[:5]}"
        )
    nn = nn.join(barcode_by_cell[["x", "y"]], on="barcode")
    metadata = pd.read_csv(args.audit_dir / "cell_metadata.tsv", sep="\t")
    nn = nn.merge(
        metadata[["cell", "SCT.weight", "DNAm.weight"]],
        on="cell",
        how="left",
        validate="one_to_one",
    )

    counts = pd.read_csv(args.counts, sep="\t")[["x", "y", "dna_reads"]]
    methyl = pd.read_csv(args.methylation, sep="\t")[[
        "x", "y", "covered_cpgs", "cpg_methylation"
    ]]
    pixel = counts.merge(methyl, on=["x", "y"], validate="one_to_one")
    pixel_lookup = pixel.set_index(["x", "y"])
    nn = nn.join(pixel_lookup, on=["x", "y"], validate="one_to_one")

    # Seurat stores one-based indices into Neighbor@cell.names.  This is the
    # same indexing used by FigR: rowMeans(mat[, NNmat[x, ]]).
    indices = nn[neighbor_columns].to_numpy(dtype=int) - 1
    if indices.min() < 0 or indices.max() >= len(nn):
        raise ValueError("WNN indices fall outside the extracted cell order")
    beta = nn["cpg_methylation"].to_numpy(dtype=float)
    nn["wnn_smoothed_beta"] = np.mean(beta[indices], axis=1)

    grid_size = 50
    fields = {
        "log10_count": np.full((grid_size, grid_size), np.nan),
        "raw_beta": np.full((grid_size, grid_size), np.nan),
        "smooth_beta": np.full((grid_size, grid_size), np.nan),
        "change": np.full((grid_size, grid_size), np.nan),
    }
    for row in nn.itertuples():
        iy, ix = int(row.y) - 1, int(row.x) - 1
        fields["log10_count"][iy, ix] = np.log10(row.dna_reads + 1)
        fields["raw_beta"][iy, ix] = row.cpg_methylation
        fields["smooth_beta"][iy, ix] = row.wnn_smoothed_beta
        fields["change"][iy, ix] = row.wnn_smoothed_beta - row.cpg_methylation

    triplets = []
    for y in range(1, grid_size + 1):
        rows = nn[(nn.y == y) & nn.x.isin([args.scar_x - 1, args.scar_x, args.scar_x + 1])]
        if len(rows) != 3:
            continue
        rows = rows.set_index("x")
        adjacent_raw = rows.loc[[args.scar_x - 1, args.scar_x + 1], "cpg_methylation"].mean()
        adjacent_smooth = rows.loc[[args.scar_x - 1, args.scar_x + 1], "wnn_smoothed_beta"].mean()
        adjacent_count = rows.loc[[args.scar_x - 1, args.scar_x + 1], "dna_reads"].mean()
        triplets.append(
            {
                "y": y,
                "count_ratio": rows.loc[args.scar_x, "dna_reads"] / adjacent_count,
                "raw_delta": rows.loc[args.scar_x, "cpg_methylation"] - adjacent_raw,
                "smooth_delta": rows.loc[args.scar_x, "wnn_smoothed_beta"] - adjacent_smooth,
                "rna_weight_delta": rows.loc[args.scar_x, "SCT.weight"]
                - rows.loc[[args.scar_x - 1, args.scar_x + 1], "SCT.weight"].mean(),
                "dna_weight_delta": rows.loc[args.scar_x, "DNAm.weight"]
                - rows.loc[[args.scar_x - 1, args.scar_x + 1], "DNAm.weight"].mean(),
            }
        )
    triplets = pd.DataFrame(triplets)

    coords = nn[["x", "y"]].to_numpy(dtype=float)
    source = np.repeat(coords, len(neighbor_columns), axis=0)
    target = coords[indices.ravel()]
    distances = np.sqrt(np.square(source - target).sum(axis=1))
    local_fraction = np.mean(distances <= np.sqrt(2) + 1e-9)
    median_distance = np.median(distances)
    self_fraction = np.mean(distances == 0)

    metrics = {
        "deposited_cells": len(nn),
        "neighbors_per_cell": len(neighbor_columns),
        "matched_scar_triplets": len(triplets),
        "median_count_ratio": triplets.count_ratio.median(),
        "raw_beta_delta": triplets.raw_delta.median(),
        "raw_beta_ci_low": bootstrap_median_ci(triplets.raw_delta)[0],
        "raw_beta_ci_high": bootstrap_median_ci(triplets.raw_delta)[1],
        "smoothed_beta_delta": triplets.smooth_delta.median(),
        "smoothed_beta_ci_low": bootstrap_median_ci(triplets.smooth_delta)[0],
        "smoothed_beta_ci_high": bootstrap_median_ci(triplets.smooth_delta)[1],
        "median_neighbor_grid_distance": median_distance,
        "local_neighbor_fraction": local_fraction,
        "self_neighbor_fraction": self_fraction,
        "scar_rna_weight_delta": triplets.rna_weight_delta.median(),
        "scar_rna_weight_ci_low": bootstrap_median_ci(triplets.rna_weight_delta)[0],
        "scar_rna_weight_ci_high": bootstrap_median_ci(triplets.rna_weight_delta)[1],
        "scar_dna_weight_delta": triplets.dna_weight_delta.median(),
        "scar_dna_weight_ci_low": bootstrap_median_ci(triplets.dna_weight_delta)[0],
        "scar_dna_weight_ci_high": bootstrap_median_ci(triplets.dna_weight_delta)[1],
    }
    pd.DataFrame([metrics]).to_csv(
        args.audit_dir / "published_e11_wnn_scar_metrics.tsv", sep="\t", index=False
    )
    nn[[
        "cell", "barcode", "x", "y", "dna_reads", "covered_cpgs",
        "cpg_methylation", "wnn_smoothed_beta"
    ]].to_csv(args.audit_dir / "published_e11_wnn_pixel_values.tsv", sep="\t", index=False)
    triplets.to_csv(args.audit_dir / "published_e11_wnn_scar_triplets.tsv", sep="\t", index=False)

    apply_style()
    fig = plt.figure(figsize=(12.0, 7.2), layout="constrained")
    gs = fig.add_gridspec(2, 4, height_ratios=[1.05, 0.8])
    axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
    profile_ax = fig.add_subplot(gs[1, :2])
    effect_ax = fig.add_subplot(gs[1, 2])
    distance_ax = fig.add_subplot(gs[1, 3])

    finite_beta = np.r_[
        fields["raw_beta"][np.isfinite(fields["raw_beta"])],
        fields["smooth_beta"][np.isfinite(fields["smooth_beta"])],
    ]
    beta_limits = np.quantile(finite_beta, [0.02, 0.98])
    count_limits = np.nanquantile(fields["log10_count"], [0.02, 0.98])
    change_limit = np.nanquantile(np.abs(fields["change"]), 0.98)
    maps = [
        (fields["log10_count"], "Raw DNA reads", "YlGnBu", Normalize(*count_limits)),
        (fields["raw_beta"], "Global CpG β", "magma", Normalize(*beta_limits)),
        (fields["smooth_beta"], "After WNN mean", "magma", Normalize(*beta_limits)),
        (
            fields["change"] * 100,
            "Smoothing change",
            "coolwarm",
            TwoSlopeNorm(vcenter=0, vmin=-change_limit * 100, vmax=change_limit * 100),
        ),
    ]
    for label, (ax, (matrix, title, cmap, norm)) in zip("ABCD", zip(axes, maps)):
        cmap_obj = plt.get_cmap(cmap).copy()
        cmap_obj.set_bad("#EEEEEE")
        image = ax.imshow(matrix, origin="upper", interpolation="nearest", cmap=cmap_obj, norm=norm)
        add_scar_outline(ax, args.scar_x)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(-0.12, 1.04, label, transform=ax.transAxes, fontsize=12, fontweight="bold")
        cbar = fig.colorbar(image, ax=ax, orientation="horizontal", fraction=0.055, pad=0.04)
        cbar.ax.tick_params(labelsize=6, length=2)
        if label == "A":
            cbar.set_label("log10(reads + 1)", fontsize=7)
        elif label == "D":
            cbar.set_label("percentage points", fontsize=7)
        else:
            cbar.set_label("methylation β", fontsize=7)

    column = nn.groupby("x", as_index=True).agg(
        count=("dna_reads", "median"),
        raw=("cpg_methylation", "median"),
        smooth=("wnn_smoothed_beta", "median"),
    )
    count_scaled = column["count"] / np.median(column["count"])
    profile_ax.plot(column.index, count_scaled, color="#2678A5", lw=1.4, label="DNA reads / median")
    beta_center = column["raw"].median()
    profile_ax.plot(
        column.index,
        (column["raw"] - beta_center) * 100,
        color="#777777",
        lw=1.15,
        label="Global β − median (pp)",
    )
    profile_ax.plot(
        column.index,
        (column["smooth"] - beta_center) * 100,
        color="#B64645",
        lw=1.5,
        label="WNN-smoothed β − median (pp)",
    )
    profile_ax.axvline(args.scar_x, color="#00A8C6", lw=1, ls="--")
    profile_ax.set_xlabel("Array column")
    profile_ax.set_ylabel("Column median (different units)")
    profile_ax.set_title("E  Column profile: depth scar remains", loc="left", fontsize=10, fontweight="bold")
    profile_ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper right")

    effects = [metrics["raw_beta_delta"] * 100, metrics["smoothed_beta_delta"] * 100]
    lows = [metrics["raw_beta_ci_low"] * 100, metrics["smoothed_beta_ci_low"] * 100]
    highs = [metrics["raw_beta_ci_high"] * 100, metrics["smoothed_beta_ci_high"] * 100]
    effect_ax.errorbar(
        effects,
        [0, 1],
        xerr=[np.array(effects) - lows, np.array(highs) - effects],
        fmt="o",
        color="#333333",
        ecolor="#777777",
        capsize=2,
    )
    effect_ax.axvline(0, color="#222222", lw=0.7)
    effect_ax.set_yticks([0, 1], ["Raw β", "WNN mean"])
    effect_ax.invert_yaxis()
    effect_ax.set_xlabel("Scar − adjacent β (pp)\nmedian; bootstrap 95% CI")
    effect_ax.set_title("F  Scar contrast", loc="left", fontsize=10, fontweight="bold")

    bins = np.arange(-0.5, max(20.5, np.ceil(distances.max()) + 0.5), 1)
    distance_ax.hist(distances, bins=bins, color="#6A5ACD", alpha=0.85)
    distance_ax.axvline(np.sqrt(2), color="#333333", lw=0.8, ls="--")
    distance_ax.set_xlim(-0.5, min(35, bins[-1]))
    distance_ax.set_xlabel("WNN edge distance on 50×50 grid")
    distance_ax.set_ylabel("Edges")
    distance_ax.set_title("G  WNN is not spatial", loc="left", fontsize=10, fontweight="bold")
    distance_ax.text(
        0.98,
        0.95,
        f"median = {median_distance:.1f} pixels\nlocal (≤√2) = {local_fraction:.1%}",
        transform=distance_ax.transAxes,
        ha="right",
        va="top",
        fontsize=7,
    )

    fig.suptitle(
        "Published E11 WNN smoothing can soften methylation display, but does not correct the read-depth scar",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.012,
        "Cyan outline: depleted column 4. WNN indices are from GEO GSE270498; counts and global β are from the current reprocessing of the same E11#1 library.",
        ha="center",
        fontsize=7,
        color="#555555",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=360, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(pd.DataFrame([metrics]).to_string(index=False))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
