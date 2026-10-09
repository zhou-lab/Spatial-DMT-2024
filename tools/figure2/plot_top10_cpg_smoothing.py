#!/usr/bin/env python3
"""Compare raw and smoothed beta at a representative broadest-decile CpG."""

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


def scar_delta(
    frame: pd.DataFrame, axis: str, index: int, other: str, value: str
) -> tuple[int, float]:
    lookup = frame.set_index([axis, other])
    deltas = []
    for position in sorted(frame[other].unique()):
        keys = [(index + offset, position) for offset in (-1, 0, 1)]
        if not all(key in lookup.index for key in keys):
            continue
        rows = lookup.loc[keys]
        values = rows[value].to_numpy(dtype=float)
        if np.all(np.isfinite(values)):
            deltas.append(values[1] - (values[0] + values[2]) / 2)
    return len(deltas), float(np.median(deltas)) if deltas else np.nan


def prepare_sample(
    sample: str,
    values_path: Path,
    pixels_path: Path,
    graph_path: Path,
    graph_cell: str,
    pixel_cell: str,
) -> tuple[pd.DataFrame, np.ndarray]:
    values = pd.read_csv(values_path, sep="\t")
    pixels = pd.read_csv(pixels_path, sep="\t", dtype={"spatial_id": str})
    pixels = pixels.merge(values, on="record", validate="one_to_one")
    pixels["raw_depth"] = pixels.top10_m + pixels.top10_u
    pixels["raw_beta"] = np.divide(
        pixels.top10_m,
        pixels.raw_depth,
        out=np.full(len(pixels), np.nan),
        where=pixels.raw_depth.to_numpy() > 0,
    )
    graph = pd.read_csv(graph_path, sep="\t", dtype={graph_cell: str})
    pixels = graph[[graph_cell]].merge(
        pixels, left_on=graph_cell, right_on=pixel_cell, validate="one_to_one"
    )
    neighbors = graph.filter(like="neighbor_").to_numpy(dtype=int)
    neighbor_beta = pixels.raw_beta.to_numpy(dtype=float)[neighbors - 1]
    neighbor_support = np.isfinite(neighbor_beta).sum(axis=1)
    pixels["covered_neighbors"] = neighbor_support
    pixels["smoothed_beta"] = np.divide(
        np.nansum(neighbor_beta, axis=1),
        neighbor_support,
        out=np.full(len(pixels), np.nan),
        where=neighbor_support > 0,
    )
    pixels["smoothing_delta"] = pixels.smoothed_beta - pixels.raw_beta
    if len(pixels) != len(graph):
        raise ValueError(f"{sample}: graph/pixel mismatch")
    return pixels, neighbors


def main() -> None:
    coordinates = pd.read_csv(OUT / "selected_site_coordinates.tsv", sep="\t")
    coordinates = coordinates.loc[
        coordinates.site_class == "top_decile_representative"
    ]

    slb, slb_neighbors = prepare_sample(
        "SpMETSLB14DM",
        OUT / "SpMETSLB14DM_pixel_values.tsv",
        WNN / "slb14_dna_knn/pixel_values.tsv",
        WNN / "slb14_dna_knn/dna_knn_indices.tsv",
        "cell",
        "spatial_id",
    )
    e85, e85_neighbors = prepare_sample(
        "E85_20um",
        OUT / "E85_20um_pixel_values.tsv",
        WNN / "e85_dna_knn/pixel_values.tsv",
        WNN / "e85_dna_knn/dna_knn_indices.tsv",
        "cell",
        "Barcodes",
    )

    datasets = [
        {
            "sample": "SpMETSLB14DM",
            "subtitle": "P21 brain · DNA-only methylation kNN",
            "frame": slb,
            "neighbors": slb_neighbors,
            "x": "x",
            "y": "y",
            "reverse": False,
            "scars": [("y", 19, "x", "row 19")],
        },
        {
            "sample": "E85_20um",
            "subtitle": "E8.5 · DNA-only methylation kNN",
            "frame": e85,
            "neighbors": e85_neighbors,
            "x": "paper_B",
            "y": "paper_A",
            "reverse": True,
            "scars": [
                ("paper_A", 22, "paper_B", "row 22"),
                ("paper_B", 14, "paper_A", "column 14"),
            ],
        },
    ]

    metric_rows = []
    scar_rows = []
    for dataset in datasets:
        frame = dataset["frame"]
        raw = frame.raw_beta.to_numpy(dtype=float)
        smooth = frame.smoothed_beta.to_numpy(dtype=float)
        site = coordinates.loc[coordinates["sample"] == dataset["sample"]].iloc[0]
        finite = np.isfinite(raw) & np.isfinite(smooth)
        paired_change = np.abs(smooth[finite] - raw[finite])
        raw_finite = raw[np.isfinite(raw)]
        smooth_finite = smooth[np.isfinite(smooth)]
        metric_rows.append(
            {
                "sample": dataset["sample"],
                "coordinate": site.coordinate,
                "covered_pixels": int(np.isfinite(raw).sum()),
                "total_pixels": len(frame),
                "smoothed_available_pixels": int(np.isfinite(smooth).sum()),
                "raw_and_smoothed_pixels": int(finite.sum()),
                "median_covered_neighbors_when_smoothed": float(
                    np.median(frame.loc[np.isfinite(smooth), "covered_neighbors"])
                ),
                "maximum_covered_neighbors": int(frame.covered_neighbors.max()),
                "pixels_with_at_least_5_covered_neighbors": int(
                    (frame.covered_neighbors >= 5).sum()
                ),
                "median_raw_depth": float(np.median(frame.loc[np.isfinite(raw), "raw_depth"])),
                "raw_smoothed_correlation": (
                    float(np.corrcoef(raw[finite], smooth[finite])[0, 1])
                    if finite.sum() >= 2
                    else np.nan
                ),
                "median_absolute_beta_change": (
                    float(np.median(paired_change)) if len(paired_change) else np.nan
                ),
                "p95_absolute_beta_change": (
                    float(np.quantile(paired_change, 0.95))
                    if len(paired_change)
                    else np.nan
                ),
                "raw_beta_sd": (
                    float(np.std(raw_finite, ddof=1)) if len(raw_finite) >= 2 else np.nan
                ),
                "smoothed_beta_sd": (
                    float(np.std(smooth_finite, ddof=1))
                    if len(smooth_finite) >= 2
                    else np.nan
                ),
                "raw_beta_range": (
                    float(np.ptp(raw_finite)) if len(raw_finite) else np.nan
                ),
                "smoothed_beta_range": (
                    float(np.ptp(smooth_finite)) if len(smooth_finite) else np.nan
                ),
            }
        )
        for axis, index, other, label in dataset["scars"]:
            raw_n, raw_delta = scar_delta(frame, axis, index, other, "raw_beta")
            smooth_n, smooth_delta = scar_delta(
                frame, axis, index, other, "smoothed_beta"
            )
            scar_rows.append(
                {
                    "sample": dataset["sample"],
                    "scar": label,
                    "raw_matched_triplets": raw_n,
                    "raw_beta_delta": raw_delta,
                    "smoothed_matched_triplets": smooth_n,
                    "smoothed_beta_delta": smooth_delta,
                }
            )
    metrics = pd.DataFrame(metric_rows)
    scar_metrics = pd.DataFrame(scar_rows)
    metrics.to_csv(OUT / "top10_cpg_smoothing_metrics.tsv", sep="\t", index=False)
    scar_metrics.to_csv(
        OUT / "top10_cpg_smoothing_scar_metrics.tsv", sep="\t", index=False
    )

    apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(10.6, 7.0), layout="constrained")
    beta_cmap = plt.get_cmap("magma").copy()
    support_cmap = plt.get_cmap("cividis").copy()
    beta_cmap.set_bad("#E8E8E8")
    support_cmap.set_bad("#E8E8E8")
    letters = iter("ABCDEF")

    for row, dataset in enumerate(datasets):
        frame = dataset["frame"]
        site = coordinates.loc[coordinates["sample"] == dataset["sample"]].iloc[0]
        raw = make_grid(
            frame, dataset["x"], dataset["y"], "raw_beta", dataset["reverse"]
        )
        smooth = make_grid(
            frame,
            dataset["x"],
            dataset["y"],
            "smoothed_beta",
            dataset["reverse"],
        )
        support = make_grid(
            frame,
            dataset["x"],
            dataset["y"],
            "covered_neighbors",
            dataset["reverse"],
        )
        beta_values = np.r_[raw[np.isfinite(raw)], smooth[np.isfinite(smooth)]]
        beta_low, beta_high = np.quantile(beta_values, [0.02, 0.98])
        beta_norm = Normalize(beta_low, beta_high, clip=True)
        support_max = max(1.0, float(np.nanmax(support)))
        panels = [
            (raw, beta_cmap, beta_norm, "Raw single-CpG β", "methylation β"),
            (
                smooth,
                beta_cmap,
                beta_norm,
                "20-neighbor mean β",
                "methylation β",
            ),
            (
                support,
                support_cmap,
                Normalize(0, support_max),
                "Covered neighbors",
                "neighbors with a call (of 20)",
            ),
        ]
        for col, (matrix, cmap, norm, title, label) in enumerate(panels):
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
            ax.set_title(title, fontsize=9.4, fontweight="bold")
            ax.text(
                -0.10, 1.03, next(letters), transform=ax.transAxes,
                fontsize=11, fontweight="bold",
            )
            bar = fig.colorbar(
                image, ax=ax, orientation="horizontal", fraction=0.05, pad=0.035
            )
            bar.ax.tick_params(labelsize=6, length=2)
            bar.set_label(label, fontsize=6.8)
        axes[row, 0].set_ylabel(
            f"{dataset['sample']}\n{site.coordinate}\n{dataset['subtitle']}",
            rotation=0,
            ha="right",
            va="center",
            labelpad=11,
            fontsize=7.4,
            fontweight="bold",
        )

    fig.suptitle(
        "Representative broadest-decile CpGs before and after neighbor smoothing",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.012,
        (
            "Smoothing is the unweighted mean of raw single-CpG β over the existing 20-neighbor graph; "
            "the focal pixel is excluded and pixels with zero covered neighbors remain missing. "
            "Raw and smoothed panels share a robust scale within each sample."
        ),
        ha="center",
        fontsize=6.9,
        color="#555555",
        wrap=True,
    )
    output = OUT / "SpMETSLB14DM_E85_top10_CpG_smoothing.png"
    fig.savefig(output, dpi=360, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(metrics.to_string(index=False))
    print(scar_metrics.to_string(index=False))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
