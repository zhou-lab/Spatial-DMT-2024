#!/usr/bin/env python3
"""Probe whether per-pixel covered-CpG counts contain discrete occupancy modes."""

import argparse
import json
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture


def parse_sample(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError("sample must be NAME=TABLE.tsv")
    name, path = value.split("=", 1)
    if not name or not path:
        raise argparse.ArgumentTypeError("sample must be NAME=TABLE.tsv")
    return name, path


def load_grid(name, path, grid_x, grid_y):
    table = pd.read_csv(
        path,
        sep="\t",
        dtype={"spatial_id": str, "Query": str},
    )
    if {"x", "y", "covered_cpgs"}.issubset(table.columns):
        table = table.loc[:, ["x", "y", "covered_cpgs"]].copy()
    elif {"Query", "N_query"}.issubset(table.columns):
        coordinates = table["Query"].str.zfill(4)
        valid = coordinates.str.fullmatch(r"\d{4}")
        table = pd.DataFrame(
            {
                "x": pd.to_numeric(coordinates.loc[valid].str[:2]),
                "y": pd.to_numeric(coordinates.loc[valid].str[2:]),
                "covered_cpgs": pd.to_numeric(table.loc[valid, "N_query"]),
            }
        )
    else:
        raise ValueError(
            f"{path} must contain x/y/covered_cpgs or YAME Query/N_query columns"
        )

    for column in ["x", "y", "covered_cpgs"]:
        table[column] = pd.to_numeric(table[column], errors="raise").astype(int)
    if table.duplicated(["x", "y"]).any():
        raise ValueError(f"{path} contains duplicate spatial coordinates")

    full = pd.MultiIndex.from_product(
        [range(1, grid_x + 1), range(1, grid_y + 1)],
        names=["x", "y"],
    ).to_frame(index=False)
    table = full.merge(table, on=["x", "y"], how="left")
    table["covered_cpgs"] = table["covered_cpgs"].fillna(0).astype(int)
    table.insert(0, "sample", name)
    table["log10_cpgs"] = np.log10(table["covered_cpgs"] + 1)
    return table


def otsu_threshold(values, bins=256):
    positive = np.asarray(values, dtype=float)
    positive = positive[np.isfinite(positive) & (positive > 0)]
    if positive.size < 2 or np.allclose(positive, positive[0]):
        return float(positive[0]) if positive.size else 0.0

    hist, edges = np.histogram(positive, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    probabilities = hist / hist.sum()
    cumulative_weight = np.cumsum(probabilities)
    cumulative_mean = np.cumsum(probabilities * centers)
    total_mean = cumulative_mean[-1]
    denominator = cumulative_weight * (1 - cumulative_weight)
    between = np.zeros_like(denominator)
    valid = denominator > 0
    between[valid] = (
        total_mean * cumulative_weight[valid] - cumulative_mean[valid]
    ) ** 2 / denominator[valid]
    return float(centers[int(np.argmax(between))])


def fit_models(log_counts, max_components, seed):
    values = np.asarray(log_counts, dtype=float).reshape(-1, 1)
    models = {}
    for components in range(1, max_components + 1):
        model = GaussianMixture(
            n_components=components,
            covariance_type="full",
            n_init=20,
            reg_covar=1e-4,
            random_state=seed,
        )
        model.fit(values)
        models[components] = {
            "model": model,
            "bic": float(model.bic(values)),
            "aic": float(model.aic(values)),
        }
    return models


def sorted_component_details(model):
    means = model.means_.ravel()
    variances = model.covariances_.reshape(-1)
    weights = model.weights_.ravel()
    order = np.argsort(means)
    old_to_new = {int(old): new + 1 for new, old in enumerate(order)}
    centers = np.power(10.0, means[order]) - 1
    arithmetic_means = np.power(10.0, means[order]) * np.exp(
        0.5 * variances[order] * math.log(10.0) ** 2
    ) - 1
    ratios = centers / centers[0]
    separations = []
    for left, right in zip(order[:-1], order[1:]):
        pooled_sd = math.sqrt((variances[left] + variances[right]) / 2)
        separations.append(
            float(abs(means[right] - means[left]) / pooled_sd)
            if pooled_sd > 0
            else math.inf
        )
    return {
        "order": order,
        "old_to_new": old_to_new,
        "log10_means": means[order],
        "variances": variances[order],
        "weights": weights[order],
        "centers": centers,
        "arithmetic_means": arithmetic_means,
        "ratios": ratios,
        "separations": separations,
    }


def classify(table, threshold, model, details):
    result = table.copy()
    foreground = result["log10_cpgs"] >= threshold
    result["foreground"] = foreground.astype(int)
    result["coverage_component"] = 0
    result["component_probability"] = np.nan
    if foreground.any():
        values = result.loc[foreground, "log10_cpgs"].to_numpy().reshape(-1, 1)
        raw_labels = model.predict(values)
        probabilities = model.predict_proba(values).max(axis=1)
        labels = np.array([details["old_to_new"][int(label)] for label in raw_labels])
        result.loc[foreground, "coverage_component"] = labels
        result.loc[foreground, "component_probability"] = probabilities
    return result


def grid_effect_r2(assignments):
    foreground = assignments.loc[
        assignments["foreground"] == 1,
        ["x", "y", "log10_cpgs"],
    ].copy()
    response = foreground["log10_cpgs"].to_numpy()
    total_ss = float(np.sum((response - response.mean()) ** 2))
    if total_ss == 0:
        return math.nan, math.nan, math.nan

    row_prediction = foreground.groupby("y")["log10_cpgs"].transform("mean")
    column_prediction = foreground.groupby("x")["log10_cpgs"].transform("mean")
    row_r2 = 1 - float(np.sum((response - row_prediction) ** 2)) / total_ss
    column_r2 = 1 - float(np.sum((response - column_prediction) ** 2)) / total_ss

    row_design = pd.get_dummies(foreground["y"], drop_first=True, dtype=float)
    column_design = pd.get_dummies(foreground["x"], drop_first=True, dtype=float)
    design = np.column_stack(
        [
            np.ones(len(foreground)),
            row_design.to_numpy(),
            column_design.to_numpy(),
        ]
    )
    coefficients, _, _, _ = np.linalg.lstsq(design, response, rcond=None)
    residual_ss = float(np.sum((response - design @ coefficients) ** 2))
    r2 = 1 - residual_ss / total_ss
    predictors = design.shape[1] - 1
    adjusted_r2 = 1 - (1 - r2) * (len(response) - 1) / (
        len(response) - predictors - 1
    )
    return row_r2, column_r2, float(adjusted_r2)


def plot_sample(path, name, assignments, threshold, model, details, bics, grid_x, grid_y):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    foreground = assignments["foreground"].astype(bool)
    x = assignments.loc[foreground, "log10_cpgs"].to_numpy()
    axes[0].hist(x, bins=60, density=True, color="#a7b6c2", edgecolor="white")
    domain = np.linspace(x.min(), x.max(), 500)
    for index, (mean, variance, weight, center) in enumerate(
        zip(
            details["log10_means"],
            details["variances"],
            details["weights"],
            details["centers"],
        ),
        start=1,
    ):
        density = weight * norm.pdf(domain, mean, math.sqrt(variance))
        axes[0].plot(domain, density, linewidth=2, label=f"C{index}: {center:,.0f}")
    axes[0].axvline(threshold, color="#222222", linestyle="--", linewidth=1)
    axes[0].set_xlabel("log10(covered CpGs + 1)")
    axes[0].set_ylabel("Density")
    axes[0].set_title(f"{name}: foreground mixture")
    axes[0].legend(fontsize=8)

    count_matrix = assignments.pivot(index="y", columns="x", values="log10_cpgs")
    count_image = axes[1].imshow(
        count_matrix,
        origin="upper",
        interpolation="nearest",
        cmap="viridis",
        extent=[0.5, grid_x + 0.5, grid_y + 0.5, 0.5],
    )
    fig.colorbar(count_image, ax=axes[1], fraction=0.046, pad=0.04)
    axes[1].set_title("CpG coverage")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")

    component_matrix = assignments.pivot(
        index="y", columns="x", values="coverage_component"
    )
    component_image = axes[2].imshow(
        component_matrix,
        origin="upper",
        interpolation="nearest",
        cmap="tab10",
        vmin=0,
        vmax=max(4, model.n_components),
        extent=[0.5, grid_x + 0.5, grid_y + 0.5, 0.5],
    )
    colorbar = fig.colorbar(component_image, ax=axes[2], fraction=0.046, pad=0.04)
    colorbar.set_label("0 = background; 1..K = coverage component")
    axes[2].set_title("Inferred coverage components")
    axes[2].set_xlabel("X")
    axes[2].set_ylabel("Y")

    bic_text = ", ".join(f"K={key}: {value:,.0f}" for key, value in bics.items())
    fig.suptitle(f"{name} | BIC {bic_text}", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def analyze_sample(name, table, output_dir, max_components, grid_x, grid_y, seed):
    positive = table.loc[table["covered_cpgs"] > 0, "log10_cpgs"].to_numpy()
    threshold = otsu_threshold(positive)
    foreground_values = table.loc[
        table["log10_cpgs"] >= threshold, "log10_cpgs"
    ].to_numpy()
    if foreground_values.size < max(20, max_components * 5):
        raise ValueError(f"{name} has too few foreground pixels for mixture fitting")

    models = fit_models(foreground_values, max_components, seed)
    best_k = min(models, key=lambda key: models[key]["bic"])
    best_model = models[best_k]["model"]
    details = sorted_component_details(best_model)
    assignments = classify(table, threshold, best_model, details)
    assignments.to_csv(
        output_dir / f"{name}_cpg_occupancy_assignments.tsv",
        sep="\t",
        index=False,
    )

    bics = {key: models[key]["bic"] for key in models}
    plot_sample(
        output_dir / f"{name}_cpg_occupancy.png",
        name,
        assignments,
        threshold,
        best_model,
        details,
        bics,
        grid_x,
        grid_y,
    )

    foreground = assignments["foreground"].astype(bool)
    confidence = assignments.loc[foreground, "component_probability"]
    row_r2, column_r2, additive_r2 = grid_effect_r2(assignments)
    expected = np.arange(1, best_k + 1, dtype=float)
    ratio_error = (
        float(np.mean(np.abs(details["ratios"] - expected) / expected))
        if best_k > 1
        else math.nan
    )
    summary = {
        "sample": name,
        "pixels_total": len(assignments),
        "pixels_with_cpgs": int((assignments["covered_cpgs"] > 0).sum()),
        "foreground_pixels": int(foreground.sum()),
        "foreground_threshold_cpgs": float(np.power(10.0, threshold) - 1),
        "median_cpgs_all_positive": float(
            assignments.loc[assignments["covered_cpgs"] > 0, "covered_cpgs"].median()
        ),
        "median_cpgs_foreground": float(
            assignments.loc[foreground, "covered_cpgs"].median()
        ),
        "best_components": best_k,
        "bic_k1": bics[1],
        "bic_improvement_over_k1": float(bics[1] - bics[best_k]),
        "component_centers": ",".join(f"{value:.1f}" for value in details["centers"]),
        "component_ratios": ",".join(f"{value:.3f}" for value in details["ratios"]),
        "component_weights": ",".join(f"{value:.4f}" for value in details["weights"]),
        "minimum_standardized_separation": (
            min(details["separations"]) if details["separations"] else math.nan
        ),
        "integer_ratio_relative_error": ratio_error,
        "median_assignment_probability": float(confidence.median()),
        "fraction_probability_ge_0_8": float((confidence >= 0.8).mean()),
        "row_effect_r2": row_r2,
        "column_effect_r2": column_r2,
        "row_column_additive_adjusted_r2": additive_r2,
    }
    return assignments, summary


def compare_pair(left_name, right_name, assignments, output_dir):
    left = assignments[left_name]
    right = assignments[right_name]
    merged = left.merge(right, on=["x", "y"], suffixes=("_left", "_right"))
    both_foreground = (merged["foreground_left"] == 1) & (
        merged["foreground_right"] == 1
    )
    rho_all = spearmanr(
        merged["log10_cpgs_left"], merged["log10_cpgs_right"]
    ).statistic
    rho_foreground = spearmanr(
        merged.loc[both_foreground, "log10_cpgs_left"],
        merged.loc[both_foreground, "log10_cpgs_right"],
    ).statistic
    ari = adjusted_rand_score(
        merged.loc[both_foreground, "coverage_component_left"],
        merged.loc[both_foreground, "coverage_component_right"],
    )
    binary_agreement = float(
        (merged["foreground_left"] == merged["foreground_right"]).mean()
    )

    fig, ax = plt.subplots(figsize=(6, 5.4))
    ax.hexbin(
        merged["log10_cpgs_left"],
        merged["log10_cpgs_right"],
        gridsize=55,
        bins="log",
        mincnt=1,
        cmap="magma",
    )
    ax.set_xlabel(f"{left_name}: log10(covered CpGs + 1)")
    ax.set_ylabel(f"{right_name}: log10(covered CpGs + 1)")
    ax.set_title(
        f"Matched pixels\nSpearman all={rho_all:.3f}; foreground={rho_foreground:.3f}; "
        f"component ARI={ari:.3f}"
    )
    fig.tight_layout()
    fig.savefig(
        output_dir / f"{left_name}_vs_{right_name}_cpg_occupancy.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)

    return {
        "left_sample": left_name,
        "right_sample": right_name,
        "matched_pixels": len(merged),
        "both_foreground_pixels": int(both_foreground.sum()),
        "spearman_all_pixels": float(rho_all),
        "spearman_foreground": float(rho_foreground),
        "foreground_binary_agreement": binary_agreement,
        "coverage_component_ari": float(ari),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample",
        action="append",
        type=parse_sample,
        required=True,
        metavar="NAME=TABLE.tsv",
    )
    parser.add_argument(
        "--replicate-pair",
        action="append",
        default=[],
        metavar="NAME1,NAME2",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--grid-x", type=int, default=96)
    parser.add_argument("--grid-y", type=int, default=96)
    parser.add_argument("--max-components", type=int, default=4)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()

    if args.max_components < 1:
        parser.error("--max-components must be at least 1")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    assignments = {}
    summaries = []
    inputs = {}
    for name, path in args.sample:
        if name in inputs:
            parser.error(f"duplicate sample name: {name}")
        if not os.path.isfile(path):
            parser.error(f"sample table does not exist: {path}")
        inputs[name] = path
        table = load_grid(name, path, args.grid_x, args.grid_y)
        sample_assignments, summary = analyze_sample(
            name,
            table,
            output_dir,
            args.max_components,
            args.grid_x,
            args.grid_y,
            args.seed,
        )
        assignments[name] = sample_assignments
        summaries.append(summary)

    pd.DataFrame(summaries).to_csv(
        output_dir / "cpg_occupancy_summary.tsv",
        sep="\t",
        index=False,
    )

    comparisons = []
    for value in args.replicate_pair:
        names = value.split(",")
        if len(names) != 2 or any(name not in assignments for name in names):
            parser.error(f"invalid replicate pair: {value}")
        comparisons.append(compare_pair(names[0], names[1], assignments, output_dir))
    if comparisons:
        pd.DataFrame(comparisons).to_csv(
            output_dir / "cpg_occupancy_replicate_comparison.tsv",
            sep="\t",
            index=False,
        )

    metadata = {
        "inputs": inputs,
        "grid": [args.grid_x, args.grid_y],
        "max_components": args.max_components,
        "seed": args.seed,
        "method": (
            "Otsu threshold on positive log10 CpG counts, followed by Gaussian "
            "mixtures selected by BIC. Components are coverage modes, not cell counts."
        ),
    }
    with open(output_dir / "cpg_occupancy_run.json", "w") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")

    print(pd.DataFrame(summaries).to_string(index=False))
    if comparisons:
        print()
        print(pd.DataFrame(comparisons).to_string(index=False))


if __name__ == "__main__":
    main()
