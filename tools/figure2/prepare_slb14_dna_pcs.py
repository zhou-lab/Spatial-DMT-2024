#!/usr/bin/env python3
"""Prepare SLB14DM metadata and methylation PCs for DNA-only kNN smoothing."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from build_e85_dna_knn import read_binned_counts


ROOT = Path(__file__).resolve().parents[2]
RUN = Path("/scr1/users/suc3/20260528_runNewSpatialPipe")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binned-cg",
        type=Path,
        default=ROOT / "results/figure2/wnn_smoothing_comparison/SpMETSLB14DM_cpg_bins_5000.bin",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results/figure2/wnn_smoothing_comparison/slb14_dna_features",
    )
    parser.add_argument("--features", type=int, default=2000)
    parser.add_argument("--pcs", type=int, default=10)
    parser.add_argument("--min-coverage-fraction", type=float, default=0.70)
    parser.add_argument("--max-imputation-iterations", type=int, default=30)
    parser.add_argument("--min-gain", type=float, default=0.01)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sample = "SpMETSLB14DM"
    final = RUN / "final_output" / sample
    pipeline = RUN / "pipeline_output" / sample
    index = pd.read_csv(
        final / f"{sample}.cg.idx",
        sep="\t",
        header=None,
        names=["spatial_id", "byte_offset"],
        dtype={"spatial_id": str},
    )
    index["spatial_id"] = index.spatial_id.str.zfill(4)
    index["record"] = np.arange(len(index), dtype=int)
    methylation = pd.read_csv(
        ROOT / "results/figure2/scar_methylation/cache/SpMETSLB14DM_stream.tsv",
        sep="\t",
    )
    index["x"] = index.spatial_id.str[:2].astype(int)
    index["y"] = index.spatial_id.str[2:].astype(int)
    dna = pd.read_csv(
        pipeline / "qc/table/SpMETSLB14DM_spatial_barcode_counts.tsv", sep="\t"
    )
    frame = (
        index.merge(methylation, on="record", validate="one_to_one")
        .merge(dna[["x", "y", "dna_reads"]], on=["x", "y"], validate="one_to_one")
        .sort_values("record")
        .reset_index(drop=True)
    )
    if len(frame) != 2500 or frame.spatial_id.nunique() != len(frame):
        raise ValueError("Expected the complete 50 x 50 SLB14DM DNA grid")

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

    pca = PCA(n_components=args.pcs, svd_solver="randomized", random_state=20260918)
    scores = pca.fit_transform(residual)
    pc_frame = pd.DataFrame(scores, columns=[f"mPC_{i + 1}" for i in range(args.pcs)])
    pc_frame.insert(0, "cell", frame.spatial_id)
    frame.rename(
        columns={"dna_reads": "raw_reads", "call_weighted_beta": "raw_beta"}
    ).to_csv(args.output_dir / "pixel_metadata.tsv", sep="\t", index=False)
    pc_frame.to_csv(args.output_dir / "dna_pc_scores.tsv", sep="\t", index=False)
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
        [{
            "pixels": len(frame),
            "cpg_ordinal_bin_width": bin_width,
            "available_bins": total.shape[1],
            "selected_variable_bins": len(selected),
            "min_feature_coverage_fraction": args.min_coverage_fraction,
            "dna_pcs": args.pcs,
            "imputation_iterations": len(convergence),
            "explained_variance_first_10_dna_pcs": pca.explained_variance_ratio_.sum(),
        }]
    )
    metrics.to_csv(args.output_dir / "dna_preparation_metrics.tsv", sep="\t", index=False)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
