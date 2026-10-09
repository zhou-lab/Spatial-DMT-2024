#!/usr/bin/env python3
"""Build a DNA-only methylation kNN graph for SpMETSLB14DM."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prepared-dir",
        type=Path,
        default=ROOT / "results/figure2/wnn_smoothing_comparison/slb14_dna_features",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results/figure2/wnn_smoothing_comparison/slb14_dna_knn",
    )
    parser.add_argument("--neighbors", type=int, default=20)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(
        args.prepared_dir / "pixel_metadata.tsv", sep="\t", dtype={"spatial_id": str}
    )
    pc_frame = pd.read_csv(
        args.prepared_dir / "dna_pc_scores.tsv", sep="\t", dtype={"cell": str}
    )
    if not frame.spatial_id.equals(pc_frame.cell):
        raise ValueError("Pixel metadata and DNA PC scores have different cell ordering")
    pc_columns = [column for column in pc_frame if column.startswith("mPC_")]
    scores = pc_frame[pc_columns].to_numpy(dtype=float)

    search = NearestNeighbors(n_neighbors=args.neighbors + 1, metric="euclidean")
    search.fit(scores)
    all_distances, all_indices = search.kneighbors(scores)
    indices = np.empty((len(frame), args.neighbors), dtype=int)
    distances = np.empty((len(frame), args.neighbors), dtype=float)
    for row in range(len(frame)):
        keep = all_indices[row] != row
        indices[row] = all_indices[row, keep][: args.neighbors]
        distances[row] = all_distances[row, keep][: args.neighbors]

    raw_beta = frame.raw_beta.to_numpy(dtype=float)
    frame["smoothed_beta"] = raw_beta[indices].mean(axis=1)
    graph = pd.DataFrame({"cell": frame.spatial_id})
    distance_frame = pd.DataFrame({"cell": frame.spatial_id})
    for neighbor in range(args.neighbors):
        graph[f"neighbor_{neighbor + 1}"] = indices[:, neighbor] + 1
        distance_frame[f"distance_{neighbor + 1}"] = distances[:, neighbor]
    graph.to_csv(args.output_dir / "dna_knn_indices.tsv", sep="\t", index=False)
    distance_frame.to_csv(
        args.output_dir / "dna_knn_distances.tsv", sep="\t", index=False
    )
    frame.to_csv(args.output_dir / "pixel_values.tsv", sep="\t", index=False)

    coordinates = frame[["x", "y"]].to_numpy(dtype=float)
    source_coordinates = np.repeat(coordinates, args.neighbors, axis=0)
    target_coordinates = coordinates[indices.reshape(-1)]
    grid_distances = np.linalg.norm(source_coordinates - target_coordinates, axis=1)
    scar_mask = frame.y.to_numpy(dtype=int) == 19
    neighbor_is_scar = scar_mask[indices]
    metrics = pd.DataFrame(
        [{
            "pixels": len(frame),
            "dna_pcs": len(pc_columns),
            "neighbors": args.neighbors,
            "median_beta_absolute_change": np.median(
                np.abs(frame.smoothed_beta.to_numpy() - raw_beta)
            ),
            "raw_smoothed_beta_correlation": np.corrcoef(
                raw_beta, frame.smoothed_beta.to_numpy(dtype=float)
            )[0, 1],
            "raw_beta_sd": np.std(raw_beta, ddof=1),
            "smoothed_beta_sd": np.std(frame.smoothed_beta, ddof=1),
            "median_neighbor_grid_distance": np.median(grid_distances),
            "local_neighbor_fraction_le_sqrt2": np.mean(grid_distances <= np.sqrt(2)),
            "self_neighbor_fraction": np.mean(grid_distances == 0),
            "row19_neighbor_fraction_for_row19_sources": np.mean(
                neighbor_is_scar[scar_mask]
            ),
            "row19_neighbor_fraction_for_other_sources": np.mean(
                neighbor_is_scar[~scar_mask]
            ),
            "row19_grid_fraction": np.mean(scar_mask),
            "max_abs_neighbor_mean_error": np.max(
                np.abs(frame.smoothed_beta.to_numpy() - raw_beta[indices].mean(axis=1))
            ),
        }]
    )
    metrics.to_csv(args.output_dir / "dna_knn_metrics.tsv", sep="\t", index=False)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
