#!/usr/bin/env python3
"""Smooth global CpG beta over a kNN graph built from smoothed MRMP PCs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from sklearn.neighbors import NearestNeighbors

from build_mrmp_smoothing import grid, project_mrmp, scar_metric
from plot_style import apply_style


ROOT = Path(__file__).resolve().parents[2]


def load_sample(
    output_dir: Path, sample: str, feature_file: str
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, list[str]]:
    frame = pd.read_csv(
        output_dir / f"{sample}_pixel_values.tsv", sep="\t", dtype={"spatial_id": str}
    )
    raw = pd.read_csv(
        output_dir / "raw_features" / feature_file,
        sep="\t",
        dtype={"sample": str},
    )
    smoothed = pd.read_csv(
        output_dir / f"{sample}_mrmp_smoothed.tsv.gz",
        sep="\t",
        dtype={"sample": str},
    )
    frame["spatial_id"] = frame.spatial_id.str.zfill(4)
    raw["sample"] = raw["sample"].str.zfill(4)
    smoothed["sample"] = smoothed["sample"].str.zfill(4)
    columns = [column for column in raw if column != "sample"]
    if columns != [column for column in smoothed if column != "sample"]:
        raise ValueError(f"{sample}: raw and smoothed MRMP columns differ")
    if raw["sample"].duplicated().any() or smoothed["sample"].duplicated().any():
        raise ValueError(f"{sample}: duplicate MRMP record names")
    raw = raw.set_index("sample").reindex(frame.spatial_id)
    smoothed = smoothed.set_index("sample").reindex(frame.spatial_id)
    if raw.index.hasnans or smoothed.index.hasnans:
        raise ValueError(f"{sample}: pixel metadata names are absent from an MRMP matrix")
    if raw[columns].isna().all(axis=1).any() or smoothed[columns].isna().all(axis=1).any():
        # Fully missing smoothed records are allowed only when the source was
        # below the graph placement floor; raw records must still exist.
        if raw[columns].isna().all(axis=1).any():
            raise ValueError(f"{sample}: a retained pixel is absent from the raw MRMP export")
    return (
        frame,
        raw[columns].to_numpy(dtype=np.float32),
        smoothed[columns].to_numpy(dtype=np.float32),
        columns,
    )


def beta_knn(
    scores: np.ndarray,
    raw_beta: np.ndarray,
    query_mask: np.ndarray,
    donor_mask: np.ndarray,
    neighbors: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    donor_rows = np.flatnonzero(donor_mask)
    if len(donor_rows) <= neighbors:
        raise ValueError("Not enough eligible MRMP-PC donors")
    search = NearestNeighbors(n_neighbors=neighbors + 1, metric="euclidean")
    search.fit(scores[donor_rows])
    query_rows = np.flatnonzero(query_mask)
    candidate_distances, candidate_positions = search.kneighbors(scores[query_rows])
    candidate_rows = donor_rows[candidate_positions]
    indices = np.full((len(scores), neighbors), -1, dtype=np.int32)
    distances = np.full((len(scores), neighbors), np.nan, dtype=np.float32)
    for local, source in enumerate(query_rows):
        keep = candidate_rows[local] != source
        selected_rows = candidate_rows[local, keep][:neighbors]
        selected_distances = candidate_distances[local, keep][:neighbors]
        if len(selected_rows) < neighbors:
            continue
        indices[source] = selected_rows
        distances[source] = selected_distances
    result = np.full(len(scores), np.nan, dtype=float)
    valid = np.all(indices >= 0, axis=1)
    result[valid] = raw_beta[indices[valid]].mean(axis=1)
    return result, indices, distances


def graph_metrics(
    frame: pd.DataFrame,
    smoothed: np.ndarray,
    indices: np.ndarray,
    distances: np.ndarray,
    method: str,
    x: str,
    y: str,
    line_axis: str,
    line_index: int,
    pcs: int,
) -> dict[str, float | int | str]:
    raw = frame.raw_beta.to_numpy(dtype=float)
    valid = np.isfinite(smoothed)
    edge_valid = indices >= 0
    coordinates = frame[[x, y]].to_numpy(dtype=float)
    source = np.repeat(coordinates, indices.shape[1], axis=0)
    flat_indices = indices.reshape(-1)
    target = np.full_like(source, np.nan)
    target[edge_valid.reshape(-1)] = coordinates[flat_indices[edge_valid.reshape(-1)]]
    grid_distance = np.linalg.norm(source - target, axis=1)
    line = frame[line_axis].to_numpy(dtype=int) == line_index
    source_line = np.repeat(line, indices.shape[1])
    neighbor_line = np.zeros_like(source_line, dtype=bool)
    neighbor_line[edge_valid.reshape(-1)] = line[flat_indices[edge_valid.reshape(-1)]]
    donor_cpgs = np.full(flat_indices.shape, np.nan)
    donor_cpgs[edge_valid.reshape(-1)] = frame.covered_cpgs.to_numpy(dtype=float)[
        flat_indices[edge_valid.reshape(-1)]
    ]
    return {
        "sample": frame.sample_name.iloc[0] if "sample_name" in frame else frame.get("Sample", pd.Series([""])).iloc[0],
        "method": method,
        "pixels": len(frame),
        "pixels_smoothed": int(valid.sum()),
        "neighbors": indices.shape[1],
        "mrmp_pcs": pcs,
        "median_abs_beta_change": float(np.nanmedian(np.abs(smoothed - raw))),
        "raw_smoothed_beta_correlation": float(np.corrcoef(raw[valid], smoothed[valid])[0, 1]),
        "median_pc_distance": float(np.nanmedian(distances)),
        "median_neighbor_grid_distance": float(np.nanmedian(grid_distance)),
        "local_neighbor_fraction_le_sqrt2": float(
            np.mean(grid_distance[np.isfinite(grid_distance)] <= np.sqrt(2))
        ),
        "median_donor_covered_cpgs": float(np.nanmedian(donor_cpgs)),
        "scar_line_neighbor_fraction_for_scar_sources": float(
            np.mean(neighbor_line[source_line & edge_valid.reshape(-1)])
        ),
        "scar_line_neighbor_fraction_for_other_sources": float(
            np.mean(neighbor_line[~source_line & edge_valid.reshape(-1)])
        ),
    }


def save_neighbors(
    frame: pd.DataFrame,
    indices: np.ndarray,
    distances: np.ndarray,
    output: Path,
) -> None:
    result = pd.DataFrame({"sample": frame.spatial_id})
    names = frame.spatial_id.to_numpy()
    for column in range(indices.shape[1]):
        valid = indices[:, column] >= 0
        neighbor_names = np.full(len(frame), "", dtype=object)
        neighbor_names[valid] = names[indices[valid, column]]
        result[f"neighbor_{column + 1}"] = neighbor_names
        result[f"distance_{column + 1}"] = distances[:, column]
    result.to_csv(output, sep="\t", index=False)


def plot(samples: list[dict], output: Path, dpi: int) -> None:
    apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(10.8, 7.4), layout="constrained")
    count_cmap = plt.get_cmap("YlGnBu").copy()
    beta_cmap = plt.get_cmap("magma").copy()
    count_cmap.set_bad("#E9E9E9")
    beta_cmap.set_bad("#E9E9E9")
    titles = ["Raw DNA read count", "Raw global CpG β", "MRMP-PC-neighbor-smoothed CpG β"]
    letters = iter("ABCDEF")
    for row, specification in enumerate(samples):
        frame = specification["frame"]
        count = np.log10(
            grid(frame, specification["x"], specification["y"], "raw_reads", specification["reverse"])
            + 1
        )
        raw = grid(
            frame, specification["x"], specification["y"], "raw_beta", specification["reverse"]
        )
        smoothed = grid(
            frame,
            specification["x"],
            specification["y"],
            "mrmp_pc_beta_smoothed",
            specification["reverse"],
        )
        count_limits = np.nanquantile(count, [0.02, 0.98])
        beta_values = np.r_[raw[np.isfinite(raw)], smoothed[np.isfinite(smoothed)]]
        beta_limits = np.quantile(beta_values, [0.02, 0.98])
        for column, matrix in enumerate([count, raw, smoothed]):
            ax = axes[row, column]
            image = ax.imshow(
                np.ma.masked_invalid(matrix),
                origin="upper",
                interpolation="nearest",
                cmap=count_cmap if column == 0 else beta_cmap,
                norm=Normalize(*(count_limits if column == 0 else beta_limits), clip=True),
                aspect="equal",
            )
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(titles[column], fontsize=9, fontweight="bold")
            ax.text(-0.14, 1.03, next(letters), transform=ax.transAxes, fontsize=11, fontweight="bold")
            bar = fig.colorbar(image, ax=ax, orientation="horizontal", fraction=0.05, pad=0.035)
            bar.ax.tick_params(labelsize=6, length=2)
            bar.set_label("log10(reads + 1)" if column == 0 else "methylation β", fontsize=6.8)
        axes[row, 0].set_ylabel(
            specification["label"], rotation=0, ha="right", va="center", labelpad=12,
            fontsize=7.5, fontweight="bold"
        )
    fig.suptitle("Global CpG beta smoothed over MethScope MRMP-PC neighbors", fontsize=13, fontweight="bold")
    fig.text(
        0.5,
        -0.012,
        "Ten PCs from the locally smoothed 500-pattern MRMP matrix; 20 nearest neighbors; "
        "donors require ≥300 raw supported patterns. No scar coordinates enter the graph.",
        ha="center",
        fontsize=6.7,
        color="#555555",
        wrap=True,
    )
    fig.savefig(output, dpi=dpi, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results/figure2/mrmp_smoothing_comparison",
    )
    parser.add_argument("--neighbors", type=int, default=20)
    parser.add_argument("--pcs", type=int, default=10)
    parser.add_argument("--pc-sensitivity", default="5,10,20")
    parser.add_argument("--donor-patterns", type=int, default=300)
    parser.add_argument("--scale-floor", type=float, default=0.05)
    parser.add_argument("--dpi", type=int, default=360)
    args = parser.parse_args()
    specifications = [
        {
            "sample": "SpMETSLB14DM",
            "feature_file": "SpMETSLB14DM_mrmp_min3.tsv",
            "x": "x",
            "y": "y",
            "line_axis": "y",
            "line_index": 19,
            "reverse": False,
            "label": "SpMETSLB14DM\nP21 brain, 20 µm",
        },
        {
            "sample": "E85_20um",
            "feature_file": "E85_20um_mrmp_all_min3.tsv",
            "x": "paper_B",
            "y": "paper_A",
            "line_axis": "paper_B",
            "line_index": 14,
            "reverse": True,
            "label": "E85_20um\nE8.5 embryo, 20 µm",
        },
    ]
    metric_rows = []
    scar_rows = []
    pc_sensitivity_rows = []
    plot_rows = []
    sensitivity_components = sorted(
        {int(value) for value in args.pc_sensitivity.split(",") if value.strip()}
    )
    if not sensitivity_components or min(sensitivity_components) < 1:
        raise ValueError("--pc-sensitivity must contain positive integers")
    for specification in specifications:
        frame, raw_mrmp, smoothed_mrmp, columns = load_sample(
            args.output_dir, specification["sample"], specification["feature_file"]
        )
        frame["sample_name"] = specification["sample"]
        raw_support = np.isfinite(raw_mrmp).sum(axis=1)
        smoothed_support = np.isfinite(smoothed_mrmp).sum(axis=1)
        donor_mask = raw_support >= args.donor_patterns
        availability = np.isfinite(raw_mrmp[donor_mask]).mean(axis=0)
        keep = availability >= 0.5
        raw_pcs, smoothed_pcs, explained = project_mrmp(
            raw_mrmp, smoothed_mrmp, donor_mask, keep, args.scale_floor,
            components=max(args.pcs, 10, max(sensitivity_components)),
        )
        query_raw = raw_support >= 20
        query_smoothed = smoothed_support >= 20
        methods = [
            ("raw_mrmp_pc_all_donors", raw_pcs[:, : args.pcs], query_raw, query_raw),
            (
                "smoothed_mrmp_pc_all_donors",
                smoothed_pcs[:, : args.pcs],
                query_smoothed,
                query_smoothed,
            ),
            (
                "smoothed_mrmp_pc_qualified_donors",
                smoothed_pcs[:, : args.pcs],
                query_smoothed,
                donor_mask,
            ),
        ]
        for method, scores, query_mask, method_donors in methods:
            beta, indices, distances = beta_knn(
                scores,
                frame.raw_beta.to_numpy(dtype=float),
                query_mask,
                method_donors,
                args.neighbors,
            )
            frame[method] = beta
            metric_rows.append(
                graph_metrics(
                    frame,
                    beta,
                    indices,
                    distances,
                    method,
                    specification["x"],
                    specification["y"],
                    specification["line_axis"],
                    specification["line_index"],
                    args.pcs,
                )
            )
            save_neighbors(
                frame,
                indices,
                distances,
                args.output_dir / f"{specification['sample']}_{method}_neighbors.tsv",
            )
            if method == "smoothed_mrmp_pc_qualified_donors":
                frame["mrmp_pc_beta_smoothed"] = beta
                for component in range(args.pcs):
                    frame[f"smoothed_mrmp_pc{component + 1}"] = smoothed_pcs[:, component]

        for value in [
            "raw_beta",
            "raw_mrmp_pc_all_donors",
            "smoothed_mrmp_pc_all_donors",
            "smoothed_mrmp_pc_qualified_donors",
        ]:
            if specification["sample"] == "SpMETSLB14DM":
                scar_rows.append(scar_metric(frame, value, "y", 19, "x", "row 19"))
            else:
                scar_rows.append(
                    scar_metric(frame, value, "paper_A", 22, "paper_B", "row 22")
                )
                scar_rows.append(
                    scar_metric(frame, value, "paper_B", 14, "paper_A", "column 14")
                )
        for component_count in sensitivity_components:
            beta, _, _ = beta_knn(
                smoothed_pcs[:, :component_count],
                frame.raw_beta.to_numpy(dtype=float),
                query_smoothed,
                donor_mask,
                args.neighbors,
            )
            frame[f"mrmp_pc{component_count}_beta_smoothed"] = beta
            frame["pc_sensitivity_beta"] = beta
            if specification["sample"] == "SpMETSLB14DM":
                row = scar_metric(frame, "pc_sensitivity_beta", "y", 19, "x", "row 19")
                row["method"] = "smoothed_mrmp_pc_qualified_donors"
                row["pcs"] = component_count
                pc_sensitivity_rows.append(row)
            else:
                for axis, index, other, label in [
                    ("paper_A", 22, "paper_B", "row 22"),
                    ("paper_B", 14, "paper_A", "column 14"),
                ]:
                    row = scar_metric(frame, "pc_sensitivity_beta", axis, index, other, label)
                    row["method"] = "smoothed_mrmp_pc_qualified_donors"
                    row["pcs"] = component_count
                    pc_sensitivity_rows.append(row)
            frame.drop(columns="pc_sensitivity_beta", inplace=True)
        output_columns = [
            column
            for column in [
                "spatial_id", "record", "x", "y", "paper_A", "paper_B", "raw_reads",
                "covered_cpgs", "raw_beta", "mrmp_patterns_min3", "mrmp_patterns_smoothed",
                "mrmp_donor_eligible", "raw_mrmp_pc_all_donors",
                "smoothed_mrmp_pc_all_donors", "smoothed_mrmp_pc_qualified_donors",
                "mrmp_pc_beta_smoothed",
            ]
            + [f"mrmp_pc{i}_beta_smoothed" for i in sensitivity_components]
            + [f"smoothed_mrmp_pc{i + 1}" for i in range(args.pcs)]
            if column in frame
        ]
        frame[output_columns].to_csv(
            args.output_dir / f"{specification['sample']}_mrmp_pc_beta_values.tsv",
            sep="\t",
            index=False,
        )
        pd.DataFrame(
            {
                "component": np.arange(1, len(explained) + 1),
                "explained_variance_fraction": explained,
            }
        ).to_csv(
            args.output_dir / f"{specification['sample']}_mrmp_pc_beta_pca_variance.tsv",
            sep="\t",
            index=False,
        )
        plot_rows.append({**specification, "frame": frame})

    pd.DataFrame(metric_rows).to_csv(
        args.output_dir / "mrmp_pc_beta_graph_metrics.tsv", sep="\t", index=False
    )
    pd.DataFrame(scar_rows).to_csv(
        args.output_dir / "mrmp_pc_beta_scar_metrics.tsv", sep="\t", index=False
    )
    pd.DataFrame(pc_sensitivity_rows).to_csv(
        args.output_dir / "mrmp_pc_beta_pcs_sensitivity.tsv", sep="\t", index=False
    )
    plot(
        plot_rows,
        args.output_dir / "SpMETSLB14DM_E85_mrmp_pc_beta_smoothing.png",
        args.dpi,
    )
    print(pd.DataFrame(metric_rows).to_string(index=False))
    print(pd.DataFrame(scar_rows).to_string(index=False))
    print(pd.DataFrame(pc_sensitivity_rows).to_string(index=False))


if __name__ == "__main__":
    main()
