#!/usr/bin/env python3
"""Plot unsmoothed single-CpG depth and beta for SLB14DM and E85."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from plot_style import apply_style


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/figure2/single_cpg_support"
WNN = ROOT / "results/figure2/wnn_smoothing_comparison"


def make_grid(
    frame: pd.DataFrame, x: str, y: str, value: str, reverse: bool = False
) -> np.ndarray:
    matrix = np.full((int(frame[y].max()), int(frame[x].max())), np.nan)
    matrix[
        frame[y].to_numpy(dtype=int) - 1,
        frame[x].to_numpy(dtype=int) - 1,
    ] = frame[value].to_numpy(dtype=float)
    return matrix[::-1, ::-1] if reverse else matrix


def robust_norm(values: np.ndarray) -> Normalize:
    finite = values[np.isfinite(values)]
    low, high = np.quantile(finite, [0.02, 0.98])
    if high <= low:
        low, high = finite.min(), finite.max()
    if high <= low:
        low, high = low - 0.5, high + 0.5
    return Normalize(low, high, clip=True)


def add_site_values(frame: pd.DataFrame, prefix: str) -> pd.DataFrame:
    frame = frame.copy()
    frame[f"{prefix}_depth"] = frame[f"{prefix}_m"] + frame[f"{prefix}_u"]
    frame[f"{prefix}_beta"] = np.divide(
        frame[f"{prefix}_m"],
        frame[f"{prefix}_depth"],
        out=np.full(len(frame), np.nan),
        where=frame[f"{prefix}_depth"].to_numpy() > 0,
    )
    return frame


def support_metrics(
    sample: str,
    site_class: str,
    coordinate: str,
    ordinal: int,
    depth: np.ndarray,
    beta: np.ndarray,
    neighbors: np.ndarray,
) -> dict[str, object]:
    covered = depth > 0
    covered_depth = depth[covered]
    neighbor_support = covered[neighbors - 1].sum(axis=1)
    return {
        "sample": sample,
        "site_class": site_class,
        "coordinate": coordinate,
        "ordinal0": ordinal,
        "total_pixels": len(depth),
        "covered_pixels": int(covered.sum()),
        "covered_pixel_fraction": float(covered.mean()),
        "pixels_depth_ge_2": int((depth >= 2).sum()),
        "pixels_depth_ge_5": int((depth >= 5).sum()),
        "pixels_depth_ge_10": int((depth >= 10).sum()),
        "median_depth_among_covered": float(np.median(covered_depth)),
        "mean_depth_among_covered": float(np.mean(covered_depth)),
        "aggregate_beta": float(np.nansum(beta * depth) / np.sum(depth)),
        "pixels_with_nonbinary_beta": int(np.sum((beta > 0) & (beta < 1))),
        "median_covered_neighbors_of_20": float(np.median(neighbor_support)),
        "pixels_with_at_least_1_covered_neighbor": int((neighbor_support >= 1).sum()),
        "fraction_with_at_least_1_covered_neighbor": float(
            np.mean(neighbor_support >= 1)
        ),
        "pixels_with_at_least_5_covered_neighbors": int((neighbor_support >= 5).sum()),
        "fraction_with_at_least_5_covered_neighbors": float(
            np.mean(neighbor_support >= 5)
        ),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    coordinates = pd.read_csv(OUT / "selected_site_coordinates.tsv", sep="\t")

    slb_values = pd.read_csv(OUT / "SpMETSLB14DM_pixel_values.tsv", sep="\t")
    slb_pixels = pd.read_csv(
        WNN / "slb14_dna_knn/pixel_values.tsv", sep="\t",
        dtype={"spatial_id": str},
    )
    slb = slb_pixels.merge(slb_values, on="record", validate="one_to_one")
    slb = add_site_values(add_site_values(slb, "top10"), "typical")
    slb_graph = pd.read_csv(
        WNN / "slb14_dna_knn/dna_knn_indices.tsv", sep="\t", dtype={"cell": str}
    )
    slb = slb_graph[["cell"]].merge(
        slb, left_on="cell", right_on="spatial_id", validate="one_to_one"
    )
    slb_neighbors = slb_graph.filter(like="neighbor_").to_numpy(dtype=int)

    e85_values = pd.read_csv(OUT / "E85_20um_pixel_values.tsv", sep="\t")
    e85_pixels = pd.read_csv(WNN / "e85_dna_knn/pixel_values.tsv", sep="\t")
    e85 = e85_pixels.merge(e85_values, on="record", validate="one_to_one")
    e85 = add_site_values(add_site_values(e85, "top10"), "typical")
    e85_graph = pd.read_csv(WNN / "e85_dna_knn/dna_knn_indices.tsv", sep="\t")
    e85 = e85_graph[["cell"]].merge(
        e85, left_on="cell", right_on="Barcodes", validate="one_to_one"
    )
    e85_neighbors = e85_graph.filter(like="neighbor_").to_numpy(dtype=int)

    datasets = [
        {
            "sample": "SpMETSLB14DM",
            "subtitle": "P21 brain, 20 µm",
            "frame": slb,
            "neighbors": slb_neighbors,
            "x": "x",
            "y": "y",
            "reverse": False,
        },
        {
            "sample": "E85_20um",
            "subtitle": "E8.5, 20 µm · corrected geometry",
            "frame": e85,
            "neighbors": e85_neighbors,
            "x": "paper_B",
            "y": "paper_A",
            "reverse": True,
        },
    ]

    metric_rows = []
    for dataset in datasets:
        sample_sites = coordinates.loc[coordinates["sample"] == dataset["sample"]]
        for site_class, prefix in [
            ("top_decile_representative", "top10"),
            ("typical_observed", "typical"),
        ]:
            site = sample_sites.loc[sample_sites.site_class == site_class].iloc[0]
            metric_rows.append(
                support_metrics(
                    dataset["sample"],
                    site_class,
                    site.coordinate,
                    int(site.ordinal0),
                    dataset["frame"][f"{prefix}_depth"].to_numpy(dtype=float),
                    dataset["frame"][f"{prefix}_beta"].to_numpy(dtype=float),
                    dataset["neighbors"],
                )
            )
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "single_cpg_support_metrics.tsv", sep="\t", index=False)

    apply_style()
    fig, axes = plt.subplots(2, 4, figsize=(13.2, 6.8), layout="constrained")
    depth_cmap = plt.get_cmap("viridis").copy()
    beta_cmap = plt.get_cmap("magma").copy()
    depth_cmap.set_bad("#E8E8E8")
    beta_cmap.set_bad("#E8E8E8")
    letters = iter("ABCDEFGH")

    for row, dataset in enumerate(datasets):
        frame = dataset["frame"]
        sample_sites = coordinates.loc[coordinates["sample"] == dataset["sample"]]
        panels = []
        for site_class, prefix, display in [
            (
                "top_decile_representative",
                "top10",
                "Top-decile CpG",
            ),
            ("typical_observed", "typical", "Typical observed CpG"),
        ]:
            site = sample_sites.loc[sample_sites.site_class == site_class].iloc[0]
            depth = make_grid(
                frame, dataset["x"], dataset["y"], f"{prefix}_depth", dataset["reverse"]
            )
            depth[depth == 0] = np.nan
            beta = make_grid(
                frame, dataset["x"], dataset["y"], f"{prefix}_beta", dataset["reverse"]
            )
            support = int(np.isfinite(depth).sum())
            panels.extend(
                [
                    (
                        np.log10(depth),
                        depth_cmap,
                        robust_norm(np.log10(depth)),
                        f"{display}: depth\n{site.coordinate} · {support:,}/{len(frame):,} pixels",
                        "log10(depth)",
                    ),
                    (
                        beta,
                        beta_cmap,
                        Normalize(0, 1),
                        f"{display}: raw β\n{site.coordinate}",
                        "methylation β",
                    ),
                ]
            )
        for col, (matrix, cmap, norm, title, colorbar_label) in enumerate(panels):
            ax = axes[row, col]
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
            ax.set_title(title, fontsize=8.5, fontweight="bold")
            ax.text(
                -0.10, 1.04, next(letters), transform=ax.transAxes,
                fontsize=11, fontweight="bold",
            )
            bar = fig.colorbar(
                image, ax=ax, orientation="horizontal", fraction=0.052, pad=0.035
            )
            bar.ax.tick_params(labelsize=6, length=2)
            bar.set_label(colorbar_label, fontsize=6.8)
        axes[row, 0].set_ylabel(
            f"{dataset['sample']}\n{dataset['subtitle']}",
            rotation=0,
            ha="right",
            va="center",
            labelpad=11,
            fontsize=7.8,
            fontweight="bold",
        )

    fig.suptitle(
        "Single-CpG spatial support before smoothing",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.015,
        (
            "Gray pixels have zero reads at that CpG. Values are direct M/(M+U) estimates; "
            "no spatial or graph smoothing is applied. Typical sites are selected at the median "
            "nonzero pixel breadth within each sample."
        ),
        ha="center",
        fontsize=7,
        color="#555555",
        wrap=True,
    )
    output = OUT / "SpMETSLB14DM_E85_top10_CpG_support.png"
    fig.savefig(output, dpi=360, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(metrics.to_string(index=False))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
