#!/usr/bin/env python3
"""SLB14DM sensitivity analysis with MethScope's mouse-brain lite MRMP."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from build_mrmp_smoothing import (
    grid,
    project_mrmp,
    shared_pattern_knn,
    smooth_mrmp,
    spatial_candidate_pool,
)
from plot_style import apply_style


ROOT = Path(__file__).resolve().parents[2]


def line_delta(frame: pd.DataFrame, field: str) -> tuple[int, float, float, float]:
    lookup = frame.set_index(["y", "x"])
    values = []
    for x in sorted(frame.x.unique()):
        keys = [(y, x) for y in (18, 19, 20)]
        if not all(key in lookup.index for key in keys):
            continue
        observed = lookup.loc[keys, field].to_numpy(dtype=float)
        if np.all(np.isfinite(observed)):
            values.append(observed[1] - observed[[0, 2]].mean())
    values = np.asarray(values)
    rng = np.random.default_rng(20260923)
    bootstrap = np.median(
        values[rng.integers(0, len(values), size=(10_000, len(values)))], axis=1
    )
    interval = np.quantile(bootstrap, [0.025, 0.975])
    return len(values), float(np.median(values)), float(interval[0]), float(interval[1])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results/figure2/mrmp_smoothing_comparison",
    )
    parser.add_argument("--neighbors", type=int, default=20)
    parser.add_argument("--donor-patterns", type=int, default=300)
    parser.add_argument("--min-shared", type=int, default=20)
    parser.add_argument("--scale-floor", type=float, default=0.05)
    parser.add_argument("--dpi", type=int, default=360)
    args = parser.parse_args()

    features = pd.read_csv(
        args.output_dir / "raw_features/SpMETSLB14DM_brain_lite_mrmp_min3.tsv",
        sep="\t",
        dtype={"sample": str},
    )
    columns = [column for column in features if column != "sample"]
    metadata = pd.read_csv(
        ROOT / "results/figure2/wnn_smoothing_comparison/slb14_dna_features/pixel_metadata.tsv",
        sep="\t",
        dtype={"spatial_id": str},
    )
    classification = pd.read_csv(
        args.output_dir / "SpMETSLB14DM_brain_lite_classification.tsv",
        sep="\t",
        dtype={"cell": str},
    )
    frame = (
        metadata.merge(features, left_on="spatial_id", right_on="sample", validate="one_to_one")
        .merge(classification, left_on="spatial_id", right_on="cell", validate="one_to_one")
    )
    values = frame[columns].to_numpy(dtype=np.float32)
    feature_count = np.isfinite(values).sum(axis=1)
    donor_mask = feature_count >= args.donor_patterns
    coordinates = frame[["x", "y"]].to_numpy(dtype=float)
    candidates = spatial_candidate_pool(coordinates, donor_mask, args.neighbors)
    indices, _, _, keep, availability = shared_pattern_knn(
        values,
        donor_mask,
        args.neighbors,
        args.min_shared,
        0.5,
        args.scale_floor,
        256,
        candidates,
    )
    smoothed, support = smooth_mrmp(values, indices, 3)
    raw_pcs, smooth_pcs, explained = project_mrmp(
        values, smoothed, donor_mask, keep, args.scale_floor
    )
    smooth_pcs[np.isfinite(smoothed[:, keep]).sum(axis=1) < args.min_shared] = np.nan
    frame["supported_patterns"] = feature_count
    frame["smoothed_patterns"] = np.isfinite(smoothed).sum(axis=1)
    frame["brain_mrmp_pc1_raw"] = raw_pcs[:, 0]
    frame["brain_mrmp_pc1_smoothed"] = smooth_pcs[:, 0]

    group_rows = []
    groups = {
        "all": np.ones(len(frame), dtype=bool),
        "row_19": frame.y.to_numpy(dtype=int) == 19,
        "adjacent_rows_18_20": frame.y.isin([18, 20]).to_numpy(),
        "other_rows": frame.y.to_numpy(dtype=int) != 19,
    }
    for name, mask in groups.items():
        group_rows.append(
            {
                "group": name,
                "pixels": int(mask.sum()),
                "median_supported_patterns": float(np.median(feature_count[mask])),
                "median_classification_confidence": float(frame.loc[mask, "confidence"].median()),
                "median_classification_certainty": float(frame.loc[mask, "certainty"].median()),
            }
        )
    pd.DataFrame(group_rows).to_csv(
        args.output_dir / "SpMETSLB14DM_brain_mrmp_group_metrics.tsv", sep="\t", index=False
    )
    scar_rows = []
    for field in ["brain_mrmp_pc1_raw", "brain_mrmp_pc1_smoothed"]:
        count, delta, low, high = line_delta(frame, field)
        scar_rows.append(
            {
                "field": field,
                "matched_triplets": count,
                "median_row19_delta": delta,
                "ci_low": low,
                "ci_high": high,
            }
        )
    pd.DataFrame(scar_rows).to_csv(
        args.output_dir / "SpMETSLB14DM_brain_mrmp_scar_metrics.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {
            "component": np.arange(1, len(explained) + 1),
            "explained_variance_fraction": explained,
            "raw_smoothed_correlation": [
                np.corrcoef(raw_pcs[np.isfinite(smooth_pcs[:, i]), i], smooth_pcs[np.isfinite(smooth_pcs[:, i]), i])[0, 1]
                for i in range(len(explained))
            ],
        }
    ).to_csv(
        args.output_dir / "SpMETSLB14DM_brain_mrmp_projection_metrics.tsv",
        sep="\t",
        index=False,
    )
    output_columns = [
        "spatial_id", "x", "y", "covered_cpgs", "raw_beta", "raw_reads",
        "supported_patterns", "smoothed_patterns", "prediction_label", "confidence",
        "certainty", "brain_mrmp_pc1_raw", "brain_mrmp_pc1_smoothed",
    ]
    frame[output_columns].to_csv(
        args.output_dir / "SpMETSLB14DM_brain_mrmp_pixel_values.tsv", sep="\t", index=False
    )
    smooth_frame = pd.DataFrame(smoothed, columns=columns)
    smooth_frame.insert(0, "sample", frame.spatial_id)
    smooth_frame.to_csv(
        args.output_dir / "SpMETSLB14DM_brain_mrmp_smoothed.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
        na_rep="NA",
    )

    apply_style()
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.5), layout="constrained")
    definitions = [
        ("supported_patterns", "Supported classifier patterns", "viridis", Normalize(0, len(columns))),
        ("confidence", "Raw classification confidence", "viridis", Normalize(0, 1)),
        ("brain_mrmp_pc1_raw", "Raw brain-MRMP PC1", "RdBu_r", None),
        ("brain_mrmp_pc1_smoothed", "Smoothed brain-MRMP PC1", "RdBu_r", None),
    ]
    pc_values = np.r_[raw_pcs[:, 0], smooth_pcs[np.isfinite(smooth_pcs[:, 0]), 0]]
    pc_limit = np.quantile(np.abs(pc_values), 0.98)
    for index, (field, title, cmap, norm) in enumerate(definitions):
        if norm is None:
            norm = Normalize(-pc_limit, pc_limit, clip=True)
        matrix = grid(frame, "x", "y", field, False)
        image = axes[index].imshow(
            np.ma.masked_invalid(matrix), origin="upper", interpolation="nearest",
            cmap=cmap, norm=norm, aspect="equal"
        )
        axes[index].set_xticks([])
        axes[index].set_yticks([])
        axes[index].set_title(title, fontsize=8.5, fontweight="bold")
        axes[index].text(-0.14, 1.03, "ABCD"[index], transform=axes[index].transAxes,
                         fontsize=11, fontweight="bold")
        bar = fig.colorbar(image, ax=axes[index], orientation="horizontal", fraction=0.05, pad=0.035)
        bar.ax.tick_params(labelsize=6, length=2)
    fig.suptitle("SLB14DM sensitivity: MethScope mouse-brain lite MRMP", fontsize=12)
    fig.savefig(
        args.output_dir / "SpMETSLB14DM_brain_mrmp_sensitivity.png",
        dpi=args.dpi,
        facecolor="white",
        bbox_inches="tight",
    )
    plt.close(fig)
    print(pd.DataFrame(group_rows).to_string(index=False))
    print(pd.DataFrame(scar_rows).to_string(index=False))


if __name__ == "__main__":
    main()
