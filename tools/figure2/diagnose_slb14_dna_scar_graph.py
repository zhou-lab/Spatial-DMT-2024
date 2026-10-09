#!/usr/bin/env python3
"""Diagnose coverage leakage in the SLB14DM DNA-only smoothing graph."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, pointbiserialr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

from plot_style import apply_style


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results/figure2/wnn_smoothing_comparison"
OUT = ROOT / "results/figure2/slb14_dna_scar_diagnostic"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_e85_dna_knn import read_binned_counts  # noqa: E402


def make_graph(scores: np.ndarray, reference: np.ndarray | None, k: int = 20) -> np.ndarray:
    if reference is None:
        reference = np.ones(len(scores), dtype=bool)
    reference_indices = np.flatnonzero(reference)
    search = NearestNeighbors(n_neighbors=k + 1).fit(scores[reference])
    _, local_indices = search.kneighbors(scores)
    candidates = reference_indices[local_indices]
    indices = np.empty((len(scores), k), dtype=int)
    for query in range(len(scores)):
        indices[query] = candidates[query, candidates[query] != query][:k]
    return indices


def paired_deltas(frame: pd.DataFrame, values: np.ndarray) -> np.ndarray:
    lookup = {
        (int(row.y), int(row.x)): index for index, row in frame.iterrows()
    }
    return np.asarray(
        [
            values[lookup[(19, x)]]
            - (values[lookup[(18, x)]] + values[lookup[(20, x)]]) / 2
            for x in range(1, 51)
        ]
    )


def method_row(
    name: str,
    label: str,
    frame: pd.DataFrame,
    raw_beta: np.ndarray,
    indices: np.ndarray | None,
    scar: np.ndarray,
    rng: np.random.Generator,
    explicit_values: np.ndarray | None = None,
) -> tuple[dict[str, object], np.ndarray]:
    if explicit_values is not None:
        values = explicit_values
        reference_pixels = np.nan
        same_row = np.nan
        donor_cpg = np.nan
        median_grid_distance = np.nan
    elif indices is None:
        values = raw_beta
        reference_pixels = len(frame)
        same_row = np.nan
        donor_cpg = np.nan
        median_grid_distance = np.nan
    else:
        values = raw_beta[indices].mean(axis=1)
        reference_pixels = int(np.unique(indices).size)
        same_row = float((frame.y.to_numpy()[indices[scar]] == 19).mean())
        donor_cpg = float(np.median(frame.covered_cpgs.to_numpy()[indices[scar]]))
        coordinates = frame[["x", "y"]].to_numpy(dtype=float)
        distances = np.linalg.norm(
            np.repeat(coordinates, indices.shape[1], axis=0)
            - coordinates[indices.reshape(-1)],
            axis=1,
        )
        median_grid_distance = float(np.median(distances))
    deltas = paired_deltas(frame, values)
    bootstrap = np.median(
        deltas[rng.integers(0, len(deltas), size=(10_000, len(deltas)))], axis=1
    )
    return (
        {
            "method": name,
            "display_label": label,
            "unique_neighbor_donors_used": reference_pixels,
            "row19_same_row_neighbor_fraction": same_row,
            "median_row19_donor_covered_cpgs": donor_cpg,
            "median_neighbor_grid_distance": median_grid_distance,
            "raw_smoothed_correlation": (
                1.0
                if indices is None and explicit_values is None
                else np.corrcoef(raw_beta, values)[0, 1]
            ),
            "median_absolute_change": (
                0.0
                if indices is None and explicit_values is None
                else np.median(np.abs(raw_beta - values))
            ),
            "row19_median_absolute_change": (
                0.0
                if indices is None and explicit_values is None
                else np.median(np.abs(raw_beta[scar] - values[scar]))
            ),
            "non_scar_median_absolute_change": (
                0.0
                if indices is None and explicit_values is None
                else np.median(np.abs(raw_beta[~scar] - values[~scar]))
            ),
            "row19_beta_contrast": np.median(deltas),
            "row19_beta_contrast_ci_low": np.quantile(bootstrap, 0.025),
            "row19_beta_contrast_ci_high": np.quantile(bootstrap, 0.975),
        },
        values,
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(
        SOURCE / "slb14_dna_knn/pixel_values.tsv", sep="\t", dtype={"spatial_id": str}
    )
    pc_frame = pd.read_csv(
        SOURCE / "slb14_dna_features/dna_pc_scores.tsv", sep="\t", dtype={"cell": str}
    )
    features = pd.read_csv(
        SOURCE / "slb14_dna_features/dna_tile_features.tsv", sep="\t"
    )
    records, methylated, unmethylated, _ = read_binned_counts(
        SOURCE / "SpMETSLB14DM_cpg_bins_5000.bin"
    )
    record_lookup = {record: index for index, record in enumerate(records)}
    order = np.asarray([record_lookup[int(record)] for record in frame.record])
    selected = features.bin_index.to_numpy(dtype=int)
    methylated = methylated[order][:, selected].astype(float)
    unmethylated = unmethylated[order][:, selected].astype(float)
    total = methylated + unmethylated
    tile_beta = np.full(total.shape, np.nan)
    np.divide(methylated, total, out=tile_beta, where=total > 0)
    observed = np.isfinite(tile_beta)
    completeness = observed.mean(axis=1)
    frame["selected_tile_completeness"] = completeness

    scores = pc_frame.filter(like="mPC_").to_numpy(dtype=float)
    raw_beta = frame.raw_beta.to_numpy(dtype=float)
    scar = frame.y.to_numpy(dtype=int) == 19
    adjacent = np.isin(frame.y.to_numpy(dtype=int), [18, 20])
    high_quality = completeness >= 0.90

    group_rows = []
    for group, mask in [
        ("row19", scar),
        ("adjacent_rows18_20", adjacent),
        ("other_rows", ~(scar | adjacent)),
    ]:
        group_rows.append(
            {
                "group": group,
                "pixels": int(mask.sum()),
                "median_raw_reads": np.median(frame.raw_reads[mask]),
                "median_covered_cpgs": np.median(frame.covered_cpgs[mask]),
                "mean_selected_tile_completeness": np.mean(completeness[mask]),
                "median_selected_tile_completeness": np.median(completeness[mask]),
                "median_call_weighted_beta": np.median(raw_beta[mask]),
            }
        )
    pd.DataFrame(group_rows).to_csv(
        OUT / "coverage_group_diagnostics.tsv", sep="\t", index=False
    )

    pc_rows = []
    log_reads = np.log10(frame.raw_reads.to_numpy(dtype=float) + 1)
    log_cpgs = np.log10(frame.covered_cpgs.to_numpy(dtype=float) + 1)
    for number, column in enumerate(pc_frame.filter(like="mPC_"), start=1):
        values = pc_frame[column].to_numpy(dtype=float)
        pooled_sd = np.sqrt((values[scar].var(ddof=1) + values[adjacent].var(ddof=1)) / 2)
        pc_rows.append(
            {
                "pc": number,
                "read_count_correlation": pearsonr(values, log_reads).statistic,
                "covered_cpg_correlation": pearsonr(values, log_cpgs).statistic,
                "tile_completeness_correlation": pearsonr(values, completeness).statistic,
                "row19_indicator_correlation": pointbiserialr(scar, values).statistic,
                "row19_vs_adjacent_standardized_difference": (
                    values[scar].mean() - values[adjacent].mean()
                ) / pooled_sd,
            }
        )
    pd.DataFrame(pc_rows).to_csv(
        OUT / "pc_technical_correlations.tsv", sep="\t", index=False
    )

    current_indices = pd.read_csv(
        SOURCE / "slb14_dna_knn/dna_knn_indices.tsv", sep="\t"
    ).filter(like="neighbor_").to_numpy(dtype=int) - 1
    polynomial = PolynomialFeatures(2, include_bias=True).fit_transform(
        np.c_[log_cpgs, completeness]
    )
    residual_scores = np.empty_like(scores)
    for component in range(scores.shape[1]):
        model = Ridge(alpha=1e-6, fit_intercept=False).fit(
            polynomial, scores[:, component]
        )
        residual_scores[:, component] = scores[:, component] - model.predict(polynomial)
    residual_indices = make_graph(StandardScaler().fit_transform(residual_scores), None)
    drop_pc1_indices = make_graph(scores[:, 1:], None)

    reference_mean = np.nanmean(tile_beta[high_quality], axis=0)
    reference_matrix = np.nan_to_num(tile_beta - reference_mean, nan=0.0)
    reference_pca = PCA(
        n_components=10, svd_solver="randomized", random_state=20260918
    )
    reference_scores = reference_pca.fit(reference_matrix[high_quality]).transform(
        reference_matrix
    )
    reference_indices = make_graph(reference_scores, high_quality)

    coordinate_lookup = {
        (int(row.y), int(row.x)): index for index, row in frame.iterrows()
    }
    local_inpaint = raw_beta.copy()
    for x in range(1, 51):
        target = coordinate_lookup[(19, x)]
        upper = coordinate_lookup[(18, x)]
        lower = coordinate_lookup[(20, x)]
        local_inpaint[target] = (
            frame.methylated_calls.iloc[upper] + frame.methylated_calls.iloc[lower]
        ) / (frame.total_calls.iloc[upper] + frame.total_calls.iloc[lower])

    rng = np.random.default_rng(20260922)
    methods = []
    method_values = {}
    for name, label, indices in [
        ("raw", "Raw β", None),
        ("current_10pc_knn", "Current DNA kNN", current_indices),
        ("drop_pc1", "Drop PC1", drop_pc1_indices),
        ("coverage_residualized", "Residualize coverage", residual_indices),
        ("reference_pca_high_quality_donors", "Reference PCA", reference_indices),
    ]:
        row, values = method_row(name, label, frame, raw_beta, indices, scar, rng)
        methods.append(row)
        method_values[name] = values
    row, values = method_row(
        "local_adjacent_count_inpaint",
        "Local inpaint*",
        frame,
        raw_beta,
        None,
        scar,
        rng,
        explicit_values=local_inpaint,
    )
    methods.append(row)
    method_values["local_adjacent_count_inpaint"] = values
    method_frame = pd.DataFrame(methods)
    method_frame.to_csv(OUT / "method_comparison.tsv", sep="\t", index=False)

    sensitivity = []
    for threshold in [0.70, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99]:
        reference = completeness >= threshold
        mean = np.nanmean(tile_beta[reference], axis=0)
        matrix = np.nan_to_num(tile_beta - mean, nan=0.0)
        pca = PCA(n_components=10, svd_solver="randomized", random_state=20260918)
        threshold_scores = pca.fit(matrix[reference]).transform(matrix)
        indices = make_graph(threshold_scores, reference)
        row, _ = method_row(
            f"reference_pca_{threshold}", f">={threshold:g}",
            frame, raw_beta, indices, scar, rng,
        )
        row["minimum_tile_completeness"] = threshold
        row["eligible_reference_pixels"] = int(reference.sum())
        sensitivity.append(row)
    pd.DataFrame(sensitivity).to_csv(
        OUT / "reference_threshold_sensitivity.tsv", sep="\t", index=False
    )

    donor_is_scar = frame.y.to_numpy()[current_indices[scar]] == 19
    graph_diagnostics = pd.DataFrame(
        [{
            "row19_grid_fraction": scar.mean(),
            "row19_same_row_neighbor_fraction": donor_is_scar.mean(),
            "row19_pixels_with_at_least_one_same_row_neighbor": np.mean(
                donor_is_scar.sum(axis=1) >= 1
            ),
            "median_neighbor_covered_cpgs_for_row19_sources": np.median(
                frame.covered_cpgs.to_numpy()[current_indices[scar]]
            ),
            "median_neighbor_covered_cpgs_for_other_sources": np.median(
                frame.covered_cpgs.to_numpy()[current_indices[~scar]]
            ),
            "median_nonrow19_donor_beta_for_row19_sources": np.median(
                raw_beta[current_indices[scar]][~donor_is_scar]
            ),
        }]
    )
    graph_diagnostics.to_csv(
        OUT / "current_graph_diagnostics.tsv", sep="\t", index=False
    )

    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.3, 7.4), layout="constrained")
    colors = ["#D55E00", "#0072B2", "#A0A0A0"]
    groups = [completeness[scar], completeness[adjacent], completeness[~(scar | adjacent)]]
    box = axes[0, 0].boxplot(groups, patch_artist=True, widths=0.62, showfliers=False)
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    axes[0, 0].set_xticks([1, 2, 3], ["Row 19", "Rows 18/20", "Other"])
    axes[0, 0].set_ylabel("Observed fraction of selected tiles")
    axes[0, 0].set_title("A  Scar pixels lose methylation-feature support", loc="left", fontweight="bold")

    axes[0, 1].scatter(
        completeness[~scar], scores[~scar, 0], s=7, alpha=0.28,
        color="#777777", linewidths=0, label="Other pixels",
    )
    axes[0, 1].scatter(
        completeness[scar], scores[scar, 0], s=16, alpha=0.85,
        color="#D55E00", linewidths=0, label="Row 19",
    )
    axes[0, 1].set_xlabel("Observed fraction of selected tiles")
    axes[0, 1].set_ylabel("DNA PC1 score")
    axes[0, 1].set_title("B  PC1 contains coverage structure", loc="left", fontweight="bold")
    axes[0, 1].legend(frameon=False, fontsize=7)
    axes[0, 1].text(
        0.03, 0.96, f"r = {pc_rows[0]['tile_completeness_correlation']:.2f}",
        transform=axes[0, 1].transAxes, va="top", fontsize=8,
    )

    donor_values = [
        graph_diagnostics.iloc[0].median_neighbor_covered_cpgs_for_row19_sources,
        graph_diagnostics.iloc[0].median_neighbor_covered_cpgs_for_other_sources,
    ]
    axes[1, 0].bar([0, 1], donor_values, color=["#D55E00", "#0072B2"], width=0.62)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xticks([0, 1], ["Row-19\nsources", "Other\nsources"])
    axes[1, 0].set_ylabel("Median covered CpGs among donors")
    axes[1, 0].set_title("C  kNN follows the low-coverage state", loc="left", fontweight="bold")
    axes[1, 0].text(
        0.04, 0.95,
        f"Same-row donors: {100 * donor_is_scar.mean():.1f}%\nGrid baseline: {100 * scar.mean():.1f}%",
        transform=axes[1, 0].transAxes, va="top", fontsize=8,
    )

    positions = np.arange(len(method_frame))
    estimates = 100 * method_frame.row19_beta_contrast.to_numpy()
    lower = estimates - 100 * method_frame.row19_beta_contrast_ci_low.to_numpy()
    upper = 100 * method_frame.row19_beta_contrast_ci_high.to_numpy() - estimates
    axes[1, 1].errorbar(
        positions, estimates, yerr=np.vstack([lower, upper]), fmt="o",
        color="#333333", ecolor="#777777", capsize=3,
    )
    axes[1, 1].axhline(0, color="#999999", linewidth=0.8)
    axes[1, 1].set_xticks(positions, method_frame.display_label, rotation=25, ha="right")
    axes[1, 1].set_ylabel("Row 19 − adjacent rows (percentage points)")
    axes[1, 1].set_title("D  Graph changes the scar contrast", loc="left", fontweight="bold")

    fig.suptitle(
        "Why DNA-only neighbor smoothing does not repair the SLB14DM scar",
        fontsize=13, fontweight="bold",
    )
    output = OUT / "SpMETSLB14DM_DNA_graph_scar_diagnostic.png"
    fig.savefig(output, dpi=360, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    pilot_values = frame[["record", "spatial_id", "x", "y", "raw_reads"]].copy()
    pilot_values["raw_call_weighted_beta"] = raw_beta
    pilot_values["current_dna_knn_beta"] = method_values["current_10pc_knn"]
    pilot_values["reference_pca_beta"] = method_values[
        "reference_pca_high_quality_donors"
    ]
    pilot_values["local_inpaint_beta"] = local_inpaint
    pilot_values.to_csv(OUT / "scar_correction_pilot_pixel_values.tsv", sep="\t", index=False)

    def grid(values: np.ndarray) -> np.ndarray:
        matrix = np.full((50, 50), np.nan)
        matrix[
            frame.y.to_numpy(dtype=int) - 1,
            frame.x.to_numpy(dtype=int) - 1,
        ] = values
        return matrix

    pilot_methods = [
        ("Raw call-weighted β", raw_beta),
        ("Current DNA kNN", method_values["current_10pc_knn"]),
        ("Reference PCA", method_values["reference_pca_high_quality_donors"]),
        ("Local count inpaint*", local_inpaint),
    ]
    all_values = np.concatenate([values for _, values in pilot_methods])
    low, high = np.quantile(all_values[np.isfinite(all_values)], [0.02, 0.98])
    map_fig, map_axes = plt.subplots(1, 4, figsize=(11.3, 3.2), layout="constrained")
    cmap = plt.get_cmap("magma").copy()
    for number, (ax, (title, values)) in enumerate(zip(map_axes, pilot_methods)):
        image = ax.imshow(
            grid(values), origin="upper", interpolation="nearest", cmap=cmap,
            vmin=low, vmax=high, aspect="equal",
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{'ABCD'[number]}  {title}", loc="left", fontsize=9, fontweight="bold")
    bar = map_fig.colorbar(image, ax=map_axes, orientation="horizontal", fraction=0.06, pad=0.04)
    bar.set_label("global CpG methylation β", fontsize=7)
    bar.ax.tick_params(labelsize=6, length=2)
    map_fig.suptitle(
        "SLB14DM scar-correction pilots",
        fontsize=12, fontweight="bold",
    )
    map_fig.text(
        0.5, -0.03,
        "*Local inpainting replaces row 19 with pooled counts from rows 18 and 20; display only, not for inference.",
        ha="center", fontsize=7, color="#555555",
    )
    pilot_output = OUT / "SpMETSLB14DM_scar_correction_pilots.png"
    map_fig.savefig(pilot_output, dpi=360, facecolor="white", bbox_inches="tight")
    plt.close(map_fig)
    print(pd.DataFrame(group_rows).to_string(index=False))
    print(method_frame.to_string(index=False))
    print(f"Wrote {output}")
    print(f"Wrote {pilot_output}")


if __name__ == "__main__":
    main()
