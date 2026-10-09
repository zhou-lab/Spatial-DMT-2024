#!/usr/bin/env python3
"""Build a DNA-only methylation kNN graph for E85 and smooth global CpG beta.

E85 is not a genuinely co-profiled RNA+DNA experiment, so WNN is inapplicable.
This reconstruction follows the DNA arm of the Spatial-DMT analysis: select
variable methylation tiles, iteratively impute missing centered values by PCA,
retain 10 PCs, find 20 nearest neighbors, and take an unweighted neighbor mean.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors


MAGIC = 0x314E494247434D53


def read_binned_counts(path: Path) -> tuple[list[int], np.ndarray, np.ndarray, int]:
    with path.open("rb") as handle:
        magic, records, bin_width = np.fromfile(handle, dtype="<u8", count=3)
        if int(magic) != MAGIC:
            raise ValueError("Unexpected binned-CG magic")
        record_ids, methylated, unmethylated = [], [], []
        expected_bins = None
        for _ in range(int(records)):
            record, bins = np.fromfile(handle, dtype="<u8", count=2)
            bins = int(bins)
            if expected_bins is None:
                expected_bins = bins
            if bins != expected_bins:
                raise ValueError(f"Inconsistent bin count: {bins} versus {expected_bins}")
            record_ids.append(int(record))
            methylated.append(np.fromfile(handle, dtype="<u4", count=bins))
            unmethylated.append(np.fromfile(handle, dtype="<u4", count=bins))
        if handle.read(1):
            raise ValueError("Trailing bytes in binned-CG file")
    return record_ids, np.stack(methylated), np.stack(unmethylated), int(bin_width)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matched",
        type=Path,
        default=Path("results/figure2/smcseq_paper_comparison/E85_20um_matched.tsv"),
    )
    parser.add_argument(
        "--binned-cg",
        type=Path,
        default=Path(
            "results/figure2/wnn_smoothing_comparison/E85_20um_cpg_bins_5000.bin"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/figure2/wnn_smoothing_comparison/e85_dna_knn"),
    )
    parser.add_argument("--features", type=int, default=2000)
    parser.add_argument("--pcs", type=int, default=10)
    parser.add_argument("--neighbors", type=int, default=20)
    parser.add_argument("--min-coverage-fraction", type=float, default=0.70)
    parser.add_argument("--max-imputation-iterations", type=int, default=30)
    parser.add_argument("--min-gain", type=float, default=0.01)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(args.matched, sep="\t")
    frame = frame.loc[frame.retained.astype(bool)].copy()
    records, methylated, unmethylated, bin_width = read_binned_counts(args.binned_cg)
    record_index = {record: i for i, record in enumerate(records)}
    order = np.array([record_index[int(record)] for record in frame.record], dtype=int)
    methylated = methylated[order]
    unmethylated = unmethylated[order]
    total = methylated.astype(np.float32) + unmethylated.astype(np.float32)
    beta = np.full(total.shape, np.nan, dtype=np.float32)
    np.divide(methylated, total, out=beta, where=total > 0)

    coverage_fraction = np.mean(total > 0, axis=0)
    variance = np.nanvar(beta, axis=0)
    eligible = np.flatnonzero(
        (coverage_fraction >= args.min_coverage_fraction)
        & np.isfinite(variance)
        & (variance > 0)
    )
    if len(eligible) < args.features:
        raise ValueError(f"Only {len(eligible)} eligible methylation tiles")
    selected = eligible[np.argsort(variance[eligible])[-args.features:]][::-1]
    beta = beta[:, selected]
    feature_mean = np.nanmean(beta, axis=0)
    residual = beta - feature_mean
    missing = np.isnan(residual)
    residual[missing] = 0.0

    convergence = []
    for iteration in range(1, args.max_imputation_iterations + 1):
        previous = residual[missing].copy()
        pca = PCA(n_components=args.pcs, svd_solver="randomized", random_state=20260918)
        scores = pca.fit_transform(residual)
        reconstructed = pca.inverse_transform(scores)
        updated = reconstructed[missing]
        residual[missing] = updated
        mse = float(np.mean(np.square(previous - updated)))
        convergence.append({"iteration": iteration, "missing_value_mse": mse})
        if iteration > 1 and mse / max(x["missing_value_mse"] for x in convergence) < args.min_gain:
            break

    # Refit once after the final missing-value update so graph PCs correspond
    # exactly to the matrix recorded by the convergence loop.
    pca = PCA(n_components=args.pcs, svd_solver="randomized", random_state=20260918)
    scores = pca.fit_transform(residual)
    pc_frame = pd.DataFrame(
        scores, columns=[f"mPC_{component + 1}" for component in range(args.pcs)]
    )
    pc_frame.insert(0, "cell", frame.Barcodes.to_numpy())
    pc_frame.to_csv(args.output_dir / "dna_pc_scores.tsv", sep="\t", index=False)
    search = NearestNeighbors(n_neighbors=args.neighbors + 1, metric="euclidean")
    search.fit(scores)
    distances_all, indices_all = search.kneighbors(scores)
    indices = np.empty((len(frame), args.neighbors), dtype=int)
    distances = np.empty((len(frame), args.neighbors), dtype=float)
    for row in range(len(frame)):
        keep = indices_all[row] != row
        indices[row] = indices_all[row, keep][: args.neighbors]
        distances[row] = distances_all[row, keep][: args.neighbors]
    raw_beta = frame.call_weighted_beta.to_numpy(dtype=float)
    frame["smoothed_beta"] = raw_beta[indices].mean(axis=1)

    # Quantify how strongly the feature-space graph departs from a spatial
    # neighborhood graph.  These diagnostics prevent the smoothed map from
    # being mistaken for a direct correction of spatial scar bands.
    coordinates = frame[["paper_B", "paper_A"]].to_numpy(dtype=float)
    source_coordinates = np.repeat(coordinates, args.neighbors, axis=0)
    target_coordinates = coordinates[indices.reshape(-1)]
    grid_distances = np.linalg.norm(source_coordinates - target_coordinates, axis=1)
    neighbor_mean_error = np.max(
        np.abs(frame.smoothed_beta.to_numpy(dtype=float) - raw_beta[indices].mean(axis=1))
    )

    graph = pd.DataFrame({"cell": frame.Barcodes})
    distance_frame = pd.DataFrame({"cell": frame.Barcodes})
    for j in range(args.neighbors):
        graph[f"neighbor_{j + 1}"] = indices[:, j] + 1
        distance_frame[f"distance_{j + 1}"] = distances[:, j]
    graph.to_csv(args.output_dir / "dna_knn_indices.tsv", sep="\t", index=False)
    distance_frame.to_csv(args.output_dir / "dna_knn_distances.tsv", sep="\t", index=False)
    frame.to_csv(args.output_dir / "pixel_values.tsv", sep="\t", index=False)
    pd.DataFrame(convergence).to_csv(
        args.output_dir / "dna_imputation_convergence.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {
            "feature": [
                f"cpgOrdinal_{i * bin_width + 1}_{(i + 1) * bin_width}" for i in selected
            ],
            "bin_index": selected,
            "coverage_fraction": coverage_fraction[selected],
            "observed_beta_variance": variance[selected],
            "observed_beta_mean": feature_mean,
        }
    ).to_csv(args.output_dir / "dna_tile_features.tsv", sep="\t", index=False)
    metrics = pd.DataFrame(
        [
            {
                "pixels": len(frame),
                "cpg_ordinal_bin_width": bin_width,
                "available_bins": total.shape[1],
                "selected_variable_bins": len(selected),
                "min_feature_coverage_fraction": args.min_coverage_fraction,
                "dna_pcs": args.pcs,
                "neighbors": args.neighbors,
                "imputation_iterations": len(convergence),
                "explained_variance_first_10_pcs": pca.explained_variance_ratio_.sum(),
                "median_beta_absolute_change": np.median(
                    np.abs(frame.smoothed_beta - frame.call_weighted_beta)
                ),
                "raw_smoothed_beta_correlation": np.corrcoef(
                    raw_beta, frame.smoothed_beta.to_numpy(dtype=float)
                )[0, 1],
                "raw_beta_sd": np.std(raw_beta, ddof=1),
                "smoothed_beta_sd": np.std(frame.smoothed_beta, ddof=1),
                "median_neighbor_grid_distance": np.median(grid_distances),
                "local_neighbor_fraction_le_sqrt2": np.mean(grid_distances <= np.sqrt(2)),
                "self_neighbor_fraction": np.mean(grid_distances == 0),
                "max_abs_neighbor_mean_error": neighbor_mean_error,
            }
        ]
    )
    metrics.to_csv(args.output_dir / "dna_knn_metrics.tsv", sep="\t", index=False)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
