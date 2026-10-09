#!/usr/bin/env python3
"""Build coverage-aware MethScope MRMP neighbor graphs and comparison maps.

The query-to-donor distance is the mean squared standardized beta difference
over MRMP patterns observed in both records.  It therefore does not impute
unobserved patterns or encode the number of missing patterns as a PCA axis.
The production graph further limits donors to records with broad MRMP support.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from plot_style import apply_style


ROOT = Path(__file__).resolve().parents[2]


def load_slb(feature_path: Path) -> tuple[pd.DataFrame, list[str]]:
    features = pd.read_csv(feature_path, sep="\t", dtype={"sample": str})
    columns = [column for column in features if column != "sample"]
    metadata = pd.read_csv(
        ROOT / "results/figure2/wnn_smoothing_comparison/slb14_dna_features/pixel_metadata.tsv",
        sep="\t",
        dtype={"spatial_id": str},
    )
    frame = metadata.merge(
        features, left_on="spatial_id", right_on="sample", validate="one_to_one"
    )
    frame["sample_name"] = "SpMETSLB14DM"
    return frame, columns


def load_e85(feature_path: Path) -> tuple[pd.DataFrame, list[str]]:
    features = pd.read_csv(feature_path, sep="\t", dtype={"sample": str})
    columns = [column for column in features if column != "sample"]
    metadata = pd.read_csv(
        ROOT / "results/figure2/smcseq_paper_comparison/E85_20um_matched.tsv",
        sep="\t",
        dtype={"spatial_id": str},
    )
    metadata["spatial_id"] = metadata.spatial_id.str.zfill(4)
    metadata = metadata.loc[metadata.retained.astype(bool)].copy()
    frame = metadata.merge(
        features, left_on="spatial_id", right_on="sample", validate="one_to_one"
    )
    frame = frame.rename(
        columns={
            "Raw_reads": "raw_reads",
            "call_weighted_beta": "raw_beta",
        }
    )
    frame["sample_name"] = "E85_20um"
    return frame, columns


def shared_pattern_knn(
    values: np.ndarray,
    donor_mask: np.ndarray,
    neighbors: int,
    min_shared: int,
    min_donor_feature_fraction: float,
    scale_floor: float,
    batch_size: int,
    candidate_donor_positions: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Exact kNN under pairwise-complete standardized MRMP distance."""
    finite = np.isfinite(values)
    donor_rows = np.flatnonzero(donor_mask)
    if len(donor_rows) <= neighbors:
        raise ValueError("Not enough coverage-qualified donor records")
    donor_finite = finite[donor_rows]
    feature_fraction = donor_finite.mean(axis=0)
    keep = feature_fraction >= min_donor_feature_fraction
    if keep.sum() < min_shared:
        raise ValueError("Too few patterns meet the donor availability gate")

    selected = values[:, keep].astype(np.float32, copy=True)
    finite = np.isfinite(selected)
    donor_values = selected[donor_rows]
    scale = np.nanstd(donor_values, axis=0, ddof=1).astype(np.float32)
    scale[~np.isfinite(scale)] = scale_floor
    scale = np.maximum(scale, scale_floor)
    selected /= scale
    selected[~finite] = 0.0
    mask = finite.astype(np.float32)

    donor_values = selected[donor_rows]
    donor_mask_values = mask[donor_rows]
    donor_sq = np.square(donor_values)
    n = len(values)
    indices = np.full((n, neighbors), -1, dtype=np.int32)
    distances = np.full((n, neighbors), np.nan, dtype=np.float32)
    shared = np.zeros((n, neighbors), dtype=np.int16)

    donor_position = {int(row): position for position, row in enumerate(donor_rows)}
    for start in range(0, n, batch_size):
        stop = min(start + batch_size, n)
        query = selected[start:stop]
        query_mask = mask[start:stop]
        overlap = query_mask @ donor_mask_values.T
        squared = (
            np.square(query) @ donor_mask_values.T
            + query_mask @ donor_sq.T
            - 2.0 * (query @ donor_values.T)
        )
        np.maximum(squared, 0.0, out=squared)
        with np.errstate(invalid="ignore", divide="ignore"):
            squared /= overlap
        squared[overlap < min_shared] = np.inf
        for local, source in enumerate(range(start, stop)):
            own = donor_position.get(source)
            if own is not None:
                squared[local, own] = np.inf
            candidates = (
                np.arange(len(donor_rows), dtype=int)
                if candidate_donor_positions is None
                else candidate_donor_positions[source]
            )
            candidates = candidates[candidates >= 0]
            candidates = candidates[np.isfinite(squared[local, candidates])]
            if len(candidates) < neighbors:
                continue
            within = np.argpartition(squared[local, candidates], neighbors - 1)[:neighbors]
            chosen = candidates[within]
            chosen = chosen[np.argsort(squared[local, chosen], kind="stable")]
            indices[source] = donor_rows[chosen]
            distances[source] = np.sqrt(squared[local, chosen])
            shared[source] = overlap[local, chosen].astype(np.int16)
    return indices, distances, shared, keep, feature_fraction


def spatial_candidate_pool(
    coordinates: np.ndarray, donor_mask: np.ndarray, pool_size: int
) -> np.ndarray:
    donor_rows = np.flatnonzero(donor_mask)
    size = min(pool_size + 1, len(donor_rows))
    search = NearestNeighbors(n_neighbors=size, metric="euclidean").fit(coordinates[donor_rows])
    candidate_rows = search.kneighbors(coordinates, return_distance=False)
    output = np.full((len(coordinates), min(pool_size, len(donor_rows))), -1, dtype=np.int32)
    donor_position = {int(row): position for position, row in enumerate(donor_rows)}
    for source in range(len(coordinates)):
        candidates = candidate_rows[source]
        own = donor_position.get(source)
        if own is not None:
            candidates = candidates[candidates != own]
        candidates = candidates[: output.shape[1]]
        output[source, : len(candidates)] = candidates
    return output


def smooth(values: np.ndarray, indices: np.ndarray) -> np.ndarray:
    result = np.full(len(values), np.nan, dtype=float)
    valid = np.all(indices >= 0, axis=1)
    result[valid] = values[indices[valid]].mean(axis=1)
    return result


def smooth_mrmp(
    values: np.ndarray, indices: np.ndarray, min_neighbor_support: int
) -> tuple[np.ndarray, np.ndarray]:
    sums = np.zeros(values.shape, dtype=np.float32)
    support = np.zeros(values.shape, dtype=np.int16)
    for column in range(indices.shape[1]):
        valid_source = indices[:, column] >= 0
        donor_values = np.full(values.shape, np.nan, dtype=np.float32)
        donor_values[valid_source] = values[indices[valid_source, column]]
        observed = np.isfinite(donor_values)
        sums += np.where(observed, donor_values, 0.0)
        support += observed
    smoothed = np.full(values.shape, np.nan, dtype=np.float32)
    np.divide(sums, support, out=smoothed, where=support >= min_neighbor_support)
    return smoothed, support


def project_mrmp(
    raw: np.ndarray,
    smoothed: np.ndarray,
    donor_mask: np.ndarray,
    keep: np.ndarray,
    scale_floor: float,
    components: int = 10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    donor = raw[donor_mask][:, keep]
    center = np.nanmean(donor, axis=0)
    scale = np.nanstd(donor, axis=0, ddof=1)
    scale = np.maximum(np.where(np.isfinite(scale), scale, scale_floor), scale_floor)

    def standardize(values: np.ndarray) -> np.ndarray:
        selected = values[:, keep]
        return np.where(np.isfinite(selected), (selected - center) / scale, 0.0)

    raw_standardized = standardize(raw)
    smooth_standardized = standardize(smoothed)
    pca = PCA(n_components=components, svd_solver="randomized", random_state=20260923)
    pca.fit(raw_standardized[donor_mask])
    return (
        pca.transform(raw_standardized),
        pca.transform(smooth_standardized),
        pca.explained_variance_ratio_,
    )


def grid(frame: pd.DataFrame, x: str, y: str, value: str, reverse: bool) -> np.ndarray:
    width, height = int(frame[x].max()), int(frame[y].max())
    matrix = np.full((height, width), np.nan)
    matrix[
        frame[y].to_numpy(dtype=int) - 1,
        frame[x].to_numpy(dtype=int) - 1,
    ] = frame[value].to_numpy(dtype=float)
    return matrix[::-1, ::-1] if reverse else matrix


def graph_metrics(
    frame: pd.DataFrame,
    indices: np.ndarray,
    distances: np.ndarray,
    shared: np.ndarray,
    method: str,
    x: str,
    y: str,
    line_axis: str,
    line_index: int,
) -> dict[str, float | int | str]:
    raw = frame.raw_beta.to_numpy(dtype=float)
    smoothed = frame[method].to_numpy(dtype=float)
    valid = np.isfinite(smoothed)
    edge_valid = indices >= 0
    source_coordinates = frame[[x, y]].to_numpy(dtype=float)
    source = np.repeat(source_coordinates, indices.shape[1], axis=0)
    target_flat = indices.reshape(-1)
    target = np.full_like(source, np.nan)
    target[edge_valid.reshape(-1)] = source_coordinates[target_flat[edge_valid.reshape(-1)]]
    spatial_distance = np.linalg.norm(source - target, axis=1)
    line = frame[line_axis].to_numpy(dtype=int) == line_index
    source_line = np.repeat(line, indices.shape[1])
    neighbor_line = np.zeros_like(source_line, dtype=bool)
    neighbor_line[edge_valid.reshape(-1)] = line[target_flat[edge_valid.reshape(-1)]]
    return {
        "sample": frame.sample_name.iloc[0],
        "method": method,
        "pixels": len(frame),
        "pixels_smoothed": int(valid.sum()),
        "donor_edges": int(edge_valid.sum()),
        "median_abs_beta_change": float(np.nanmedian(np.abs(smoothed - raw))),
        "raw_smoothed_beta_correlation": float(np.corrcoef(raw[valid], smoothed[valid])[0, 1]),
        "median_mrmp_distance": float(np.nanmedian(distances)),
        "median_shared_patterns": float(np.nanmedian(np.where(shared > 0, shared, np.nan))),
        "median_neighbor_grid_distance": float(np.nanmedian(spatial_distance)),
        "local_neighbor_fraction_le_sqrt2": float(np.nanmean(spatial_distance <= np.sqrt(2))),
        "scar_line_neighbor_fraction_for_scar_sources": float(
            np.mean(neighbor_line[source_line & edge_valid.reshape(-1)])
        ),
        "scar_line_neighbor_fraction_for_other_sources": float(
            np.mean(neighbor_line[~source_line & edge_valid.reshape(-1)])
        ),
        "scar_line_grid_fraction": float(line.mean()),
    }


def scar_metric(
    frame: pd.DataFrame,
    value: str,
    axis: str,
    index: int,
    other: str,
    label: str,
) -> dict[str, float | int | str]:
    lookup = frame.set_index([axis, other])
    deltas: list[float] = []
    for position in sorted(frame[other].unique()):
        keys = [(index + offset, position) for offset in (-1, 0, 1)]
        if not all(key in lookup.index for key in keys):
            continue
        rows = lookup.loc[keys]
        observed = rows[value].to_numpy(dtype=float)
        if not np.all(np.isfinite(observed)):
            continue
        deltas.append(float(observed[1] - observed[[0, 2]].mean()))
    values = np.asarray(deltas)
    rng = np.random.default_rng(20260923 + index + len(value))
    bootstrap = np.median(
        values[rng.integers(0, len(values), size=(10_000, len(values)))], axis=1
    )
    interval = np.quantile(bootstrap, [0.025, 0.975])
    return {
        "sample": frame.sample_name.iloc[0],
        "scar": label,
        "method": value,
        "matched_triplets": len(values),
        "median_beta_delta": float(np.median(values)),
        "ci_low": float(interval[0]),
        "ci_high": float(interval[1]),
    }


def save_graph(
    frame: pd.DataFrame,
    indices: np.ndarray,
    distances: np.ndarray,
    shared: np.ndarray,
    path: Path,
) -> None:
    output = pd.DataFrame({"sample": frame.spatial_id})
    for column in range(indices.shape[1]):
        valid = indices[:, column] >= 0
        names = np.full(len(frame), "", dtype=object)
        names[valid] = frame.spatial_id.to_numpy()[indices[valid, column]]
        output[f"neighbor_{column + 1}"] = names
        output[f"distance_{column + 1}"] = distances[:, column]
        output[f"shared_patterns_{column + 1}"] = shared[:, column]
    output.to_csv(path, sep="\t", index=False)


def plot(samples: list[dict], output: Path, dpi: int) -> None:
    apply_style()
    fig, axes = plt.subplots(2, 4, figsize=(13.3, 7.4), layout="constrained")
    count_cmap = plt.get_cmap("YlGnBu").copy()
    support_cmap = plt.get_cmap("viridis").copy()
    pc_cmap = plt.get_cmap("RdBu_r").copy()
    count_cmap.set_bad("#E9E9E9")
    support_cmap.set_bad("#E9E9E9")
    pc_cmap.set_bad("#E9E9E9")
    titles = [
        "Raw DNA read count",
        "Supported MRMP patterns",
        "Raw MRMP PC1",
        "Smoothed MRMP PC1",
    ]
    letters = iter("ABCDEFGH")
    for row, specification in enumerate(samples):
        frame = specification["frame"]
        count = np.log10(
            grid(frame, specification["x"], specification["y"], "raw_reads", specification["reverse"])
            + 1
        )
        supported = grid(
            frame,
            specification["x"],
            specification["y"],
            "mrmp_patterns_min3",
            specification["reverse"],
        )
        raw_pc = grid(
            frame,
            specification["x"],
            specification["y"],
            "mrmp_pc1_raw",
            specification["reverse"],
        )
        smooth_pc = grid(
            frame,
            specification["x"],
            specification["y"],
            "mrmp_pc1_smoothed",
            specification["reverse"],
        )
        count_limits = np.nanquantile(count, [0.02, 0.98])
        pc_values = np.r_[raw_pc[np.isfinite(raw_pc)], smooth_pc[np.isfinite(smooth_pc)]]
        pc_limit = np.quantile(np.abs(pc_values), 0.98)
        matrices = [count, supported, raw_pc, smooth_pc]
        cmaps = [count_cmap, support_cmap, pc_cmap, pc_cmap]
        norms = [
            Normalize(*count_limits, clip=True),
            Normalize(0, 500, clip=True),
            Normalize(-pc_limit, pc_limit, clip=True),
            Normalize(-pc_limit, pc_limit, clip=True),
        ]
        labels = ["log10(reads + 1)", "patterns (of 500)", "PC1 score", "PC1 score"]
        for column, matrix in enumerate(matrices):
            ax = axes[row, column]
            image = ax.imshow(
                np.ma.masked_invalid(matrix),
                origin="upper",
                interpolation="nearest",
                cmap=cmaps[column],
                norm=norms[column],
                aspect="equal",
            )
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(titles[column], fontsize=9, fontweight="bold")
            ax.text(-0.14, 1.03, next(letters), transform=ax.transAxes, fontsize=11, fontweight="bold")
            bar = fig.colorbar(image, ax=ax, orientation="horizontal", fraction=0.05, pad=0.035)
            bar.ax.tick_params(labelsize=6, length=2)
            bar.set_label(labels[column], fontsize=6.8)
        axes[row, 0].set_ylabel(
            specification["label"], rotation=0, ha="right", va="center", labelpad=12,
            fontsize=7.5, fontweight="bold"
        )
    fig.suptitle("Coverage-aware smoothing of MethScope MRMP features", fontsize=13, fontweight="bold")
    fig.text(
        0.5,
        -0.012,
        "500-pattern mm10 whole-genome MRMP; pattern means require ≥3 observed CpGs. "
        "Each pixel averages the same MRMP feature from its 20 nearest coverage-qualified spatial pixels; "
        "no scar coordinate enters the graph.",
        ha="center",
        fontsize=6.7,
        color="#555555",
        wrap=True,
    )
    fig.savefig(output, dpi=dpi, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def plot_pattern_examples(examples: list[dict], output: Path, dpi: int) -> None:
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.8, 7.2), layout="constrained")
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#E9E9E9")
    letters = iter("ABCD")
    for row, specification in enumerate(examples):
        frame = specification["frame"]
        for column, (field, title) in enumerate(
            [("raw_pattern", "Raw MRMP pattern β"), ("smooth_pattern", "Locally smoothed MRMP pattern β")]
        ):
            matrix = grid(
                frame,
                specification["x"],
                specification["y"],
                field,
                specification["reverse"],
            )
            ax = axes[row, column]
            image = ax.imshow(
                np.ma.masked_invalid(matrix),
                origin="upper",
                interpolation="nearest",
                cmap=cmap,
                norm=Normalize(0, 1, clip=True),
                aspect="equal",
            )
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(title, fontsize=9, fontweight="bold")
            ax.text(-0.14, 1.03, next(letters), transform=ax.transAxes, fontsize=11, fontweight="bold")
            bar = fig.colorbar(image, ax=ax, orientation="horizontal", fraction=0.05, pad=0.035)
            bar.ax.tick_params(labelsize=6, length=2)
            bar.set_label("pattern methylation β", fontsize=6.8)
        axes[row, 0].set_ylabel(
            f"{specification['name']}\n{specification['pattern']}",
            rotation=0,
            ha="right",
            va="center",
            labelpad=12,
            fontsize=7.5,
            fontweight="bold",
        )
    fig.suptitle("Representative MethScope MRMP features before and after smoothing", fontsize=12)
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
    parser.add_argument("--donor-patterns", type=int, default=300)
    parser.add_argument("--min-shared", type=int, default=20)
    parser.add_argument("--spatial-pool", type=int, default=100)
    parser.add_argument("--min-neighbor-support", type=int, default=3)
    parser.add_argument("--donor-feature-fraction", type=float, default=0.5)
    parser.add_argument("--scale-floor", type=float, default=0.05)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--dpi", type=int, default=360)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    feature_dir = args.output_dir / "raw_features"

    slb, slb_columns = load_slb(feature_dir / "SpMETSLB14DM_mrmp_min3.tsv")
    e85, e85_columns = load_e85(feature_dir / "E85_20um_mrmp_all_min3.tsv")
    if slb_columns != e85_columns:
        raise ValueError("SLB14DM and E85 were not summarized against the same MRMP")

    graph_metrics_rows: list[dict] = []
    feature_rows: list[dict] = []
    pc_scar_rows: list[dict] = []
    projection_rows: list[dict] = []
    pattern_examples: list[dict] = []
    specifications = [
        (slb, slb_columns, "x", "y", "y", 19, "SpMETSLB14DM"),
        (e85, e85_columns, "paper_B", "paper_A", "paper_B", 14, "E85_20um"),
    ]
    for frame, columns, x, y, line_axis, line_index, name in specifications:
        values = frame[columns].to_numpy(dtype=np.float32)
        feature_count = np.isfinite(values).sum(axis=1)
        frame["mrmp_patterns_min3"] = feature_count
        qualified_donors = feature_count >= args.donor_patterns
        coordinates = frame[[x, y]].to_numpy(dtype=float)
        spatial_pool = spatial_candidate_pool(coordinates, qualified_donors, args.spatial_pool)
        local_pool = spatial_candidate_pool(coordinates, qualified_donors, args.neighbors)
        donor_definitions = [
            ("mrmp_unrestricted", feature_count >= args.min_shared, None),
            ("mrmp_coverage_qualified", qualified_donors, None),
            ("mrmp_spatially_constrained", qualified_donors, spatial_pool),
            ("mrmp_local_spatial", qualified_donors, local_pool),
        ]
        for method, donor_mask, candidates in donor_definitions:
            indices, distances, shared, keep, availability = shared_pattern_knn(
                values,
                donor_mask,
                args.neighbors,
                args.min_shared,
                args.donor_feature_fraction,
                args.scale_floor,
                args.batch_size,
                candidates,
            )
            frame[method] = smooth(frame.raw_beta.to_numpy(dtype=float), indices)
            graph_metrics_rows.append(
                graph_metrics(frame, indices, distances, shared, method, x, y, line_axis, line_index)
            )
            save_graph(frame, indices, distances, shared, args.output_dir / f"{name}_{method}_neighbors.tsv")
            if method == "mrmp_local_spatial":
                frame["mrmp_donor_eligible"] = donor_mask
                smoothed_mrmp, neighbor_support = smooth_mrmp(
                    values, indices, args.min_neighbor_support
                )
                frame["mrmp_patterns_smoothed"] = np.isfinite(smoothed_mrmp).sum(axis=1)
                raw_pcs, smoothed_pcs, explained = project_mrmp(
                    values, smoothed_mrmp, donor_mask, keep, args.scale_floor
                )
                smoothed_pcs[
                    np.isfinite(smoothed_mrmp[:, keep]).sum(axis=1) < args.min_shared
                ] = np.nan
                for component in range(raw_pcs.shape[1]):
                    frame[f"mrmp_pc{component + 1}_raw"] = raw_pcs[:, component]
                    frame[f"mrmp_pc{component + 1}_smoothed"] = smoothed_pcs[:, component]
                    finite = np.isfinite(smoothed_pcs[:, component])
                    projection_rows.append(
                        {
                            "sample": name,
                            "component": component + 1,
                            "explained_variance_fraction": explained[component],
                            "raw_smoothed_correlation_all": np.corrcoef(
                                raw_pcs[finite, component], smoothed_pcs[finite, component]
                            )[0, 1],
                            "raw_smoothed_correlation_donor_eligible": np.corrcoef(
                                raw_pcs[donor_mask, component], smoothed_pcs[donor_mask, component]
                            )[0, 1],
                            "raw_support_spearman": pd.Series(raw_pcs[:, component]).corr(
                                pd.Series(feature_count), method="spearman"
                            ),
                            "smoothed_support_spearman": pd.Series(smoothed_pcs[:, component]).corr(
                                pd.Series(feature_count), method="spearman"
                            ),
                        }
                    )
                smooth_frame = pd.DataFrame(smoothed_mrmp, columns=columns)
                smooth_frame.insert(0, "sample", frame.spatial_id.to_numpy())
                smooth_frame.to_csv(
                    args.output_dir / f"{name}_mrmp_smoothed.tsv.gz",
                    sep="\t",
                    index=False,
                    compression="gzip",
                    na_rep="NA",
                )
                support_frame = pd.DataFrame(neighbor_support, columns=columns)
                support_frame.insert(0, "sample", frame.spatial_id.to_numpy())
                support_frame.to_csv(
                    args.output_dir / f"{name}_mrmp_neighbor_support.tsv.gz",
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
                pd.DataFrame(
                    {
                        "component": np.arange(1, len(explained) + 1),
                        "explained_variance_fraction": explained,
                    }
                ).to_csv(
                    args.output_dir / f"{name}_mrmp_pca_variance.tsv", sep="\t", index=False
                )
                high_availability = availability >= 0.9
                pattern_variance = np.nanvar(values[donor_mask], axis=0)
                candidates_for_plot = np.flatnonzero(keep & high_availability)
                selected_pattern = candidates_for_plot[
                    np.nanargmax(pattern_variance[candidates_for_plot])
                ]
                example_frame = frame.copy()
                example_frame["raw_pattern"] = values[:, selected_pattern]
                example_frame["smooth_pattern"] = smoothed_mrmp[:, selected_pattern]
                pattern_examples.append(
                    {
                        "name": name,
                        "pattern": columns[selected_pattern],
                        "pattern_index": int(selected_pattern),
                        "frame": example_frame,
                        "x": x,
                        "y": y,
                        "reverse": name == "E85_20um",
                        "donor_observed_fraction": float(availability[selected_pattern]),
                        "donor_variance": float(pattern_variance[selected_pattern]),
                    }
                )
                for column, selected, fraction in zip(columns, keep, availability):
                    feature_rows.append(
                        {
                            "sample": name,
                            "pattern": column,
                            "retained_for_distance": bool(selected),
                            "coverage_qualified_donor_fraction_observed": float(fraction),
                        }
                    )

        pixel_columns = [
            column
            for column in [
                "spatial_id", "record", "x", "y", "paper_A", "paper_B", "raw_reads",
                "covered_cpgs", "raw_beta", "mrmp_patterns_min3", "mrmp_donor_eligible",
                "mrmp_patterns_smoothed", "mrmp_unrestricted", "mrmp_coverage_qualified",
                "mrmp_spatially_constrained", "mrmp_local_spatial", "mrmp_pc1_raw",
                "mrmp_pc1_smoothed",
            ]
            if column in frame
        ]
        frame[pixel_columns].to_csv(
            args.output_dir / f"{name}_pixel_values.tsv", sep="\t", index=False
        )

    scar_rows = []
    for value in [
        "raw_beta",
        "mrmp_unrestricted",
        "mrmp_coverage_qualified",
        "mrmp_spatially_constrained",
        "mrmp_local_spatial",
    ]:
        scar_rows.append(scar_metric(slb, value, "y", 19, "x", "row 19"))
        scar_rows.append(scar_metric(e85, value, "paper_A", 22, "paper_B", "row 22"))
        scar_rows.append(scar_metric(e85, value, "paper_B", 14, "paper_A", "column 14"))
    for value in ["mrmp_pc1_raw", "mrmp_pc1_smoothed"]:
        pc_scar_rows.append(scar_metric(slb, value, "y", 19, "x", "row 19"))
        pc_scar_rows.append(scar_metric(e85, value, "paper_A", 22, "paper_B", "row 22"))
        pc_scar_rows.append(scar_metric(e85, value, "paper_B", 14, "paper_A", "column 14"))
    pd.DataFrame(graph_metrics_rows).to_csv(
        args.output_dir / "mrmp_graph_metrics.tsv", sep="\t", index=False
    )
    pd.DataFrame(scar_rows).to_csv(
        args.output_dir / "mrmp_scar_metrics.tsv", sep="\t", index=False
    )
    pd.DataFrame(pc_scar_rows).to_csv(
        args.output_dir / "mrmp_pc1_scar_metrics.tsv", sep="\t", index=False
    )
    pd.DataFrame(projection_rows).to_csv(
        args.output_dir / "mrmp_projection_metrics.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        [
            {key: value for key, value in example.items() if key != "frame"}
            for example in pattern_examples
        ]
    ).to_csv(args.output_dir / "mrmp_pattern_examples.tsv", sep="\t", index=False)
    pd.DataFrame(feature_rows).to_csv(
        args.output_dir / "mrmp_feature_coverage.tsv", sep="\t", index=False
    )

    plot(
        [
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
        ],
        args.output_dir / "SpMETSLB14DM_E85_mrmp_smoothing.png",
        args.dpi,
    )
    plot_pattern_examples(
        pattern_examples,
        args.output_dir / "SpMETSLB14DM_E85_mrmp_pattern_examples.png",
        args.dpi,
    )
    print(pd.DataFrame(graph_metrics_rows).to_string(index=False))
    print(pd.DataFrame(scar_rows).to_string(index=False))
    print(pd.DataFrame(pc_scar_rows).to_string(index=False))


if __name__ == "__main__":
    main()
