#!/usr/bin/env python3
"""Fit CpG occupancy modes after robust row/column grid normalization."""

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

from cpg_coverage_occupancy import (
    fit_models,
    grid_effect_r2,
    load_grid,
    otsu_threshold,
    parse_sample,
    sorted_component_details,
)


def finite_median(values):
    finite = values[np.isfinite(values)]
    return float(np.median(finite)) if finite.size else 0.0


def median_polish(table, foreground, grid_x, grid_y, iterations, tolerance):
    matrix = (
        table.assign(value=table["log10_cpgs"].where(foreground))
        .pivot(index="y", columns="x", values="value")
        .reindex(index=range(1, grid_y + 1), columns=range(1, grid_x + 1))
        .to_numpy(dtype=float)
    )
    overall = finite_median(matrix.ravel())
    residual = matrix - overall
    row_effects = np.zeros(grid_y)
    column_effects = np.zeros(grid_x)
    converged_at = iterations

    for iteration in range(1, iterations + 1):
        row_updates = np.array([finite_median(row) for row in residual])
        residual -= row_updates[:, np.newaxis]
        row_effects += row_updates

        column_updates = np.array(
            [finite_median(residual[:, index]) for index in range(grid_x)]
        )
        residual -= column_updates[np.newaxis, :]
        column_effects += column_updates

        largest_update = max(
            float(np.max(np.abs(row_updates))),
            float(np.max(np.abs(column_updates))),
        )
        if largest_update < tolerance:
            converged_at = iteration
            break

    normalized = table["log10_cpgs"].to_numpy().copy()
    normalized -= row_effects[table["y"].to_numpy(dtype=int) - 1]
    normalized -= column_effects[table["x"].to_numpy(dtype=int) - 1]
    return normalized, row_effects, column_effects, overall, converged_at


def classify(table, foreground, model, details):
    result = table.copy()
    result["foreground"] = foreground.astype(int)
    result["coverage_component"] = 0
    result["component_probability"] = np.nan
    values = result.loc[foreground, "normalized_log10_cpgs"].to_numpy().reshape(-1, 1)
    raw_labels = model.predict(values)
    labels = np.array([details["old_to_new"][int(label)] for label in raw_labels])
    result.loc[foreground, "coverage_component"] = labels
    result.loc[foreground, "component_probability"] = model.predict_proba(values).max(axis=1)
    return result


def plot_sample(path, name, assignments, model, details, bics, grid_x, grid_y):
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.8))
    foreground = assignments["foreground"].astype(bool)

    raw_matrix = assignments.pivot(index="y", columns="x", values="log10_cpgs")
    normalized_matrix = assignments.pivot(
        index="y", columns="x", values="normalized_log10_cpgs"
    )
    foreground_values = assignments.loc[foreground, "normalized_log10_cpgs"]
    lower, upper = np.quantile(foreground_values, [0.01, 0.99])
    raw_image = axes[0].imshow(
        raw_matrix,
        origin="upper",
        interpolation="nearest",
        cmap="viridis",
        vmin=lower,
        vmax=upper,
        extent=[0.5, grid_x + 0.5, grid_y + 0.5, 0.5],
    )
    fig.colorbar(raw_image, ax=axes[0], fraction=0.046, pad=0.04)
    axes[0].set_title("Raw CpG coverage")

    normalized_image = axes[1].imshow(
        normalized_matrix,
        origin="upper",
        interpolation="nearest",
        cmap="viridis",
        vmin=lower,
        vmax=upper,
        extent=[0.5, grid_x + 0.5, grid_y + 0.5, 0.5],
    )
    fig.colorbar(normalized_image, ax=axes[1], fraction=0.046, pad=0.04)
    axes[1].set_title("Row/column normalized")

    values = foreground_values.to_numpy()
    axes[2].hist(values, bins=60, density=True, color="#a7b6c2", edgecolor="white")
    domain = np.linspace(values.min(), values.max(), 500)
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
        axes[2].plot(domain, density, linewidth=2, label=f"C{index}: {center:,.0f}")
    axes[2].set_xlabel("Normalized log10(covered CpGs + 1)")
    axes[2].set_ylabel("Density")
    axes[2].set_title("Normalized foreground mixture")
    axes[2].legend(frameon=False, fontsize=8)

    component_matrix = assignments.pivot(
        index="y", columns="x", values="coverage_component"
    )
    component_image = axes[3].imshow(
        component_matrix,
        origin="upper",
        interpolation="nearest",
        cmap="tab10",
        vmin=0,
        vmax=max(4, model.n_components),
        extent=[0.5, grid_x + 0.5, grid_y + 0.5, 0.5],
    )
    colorbar = fig.colorbar(component_image, ax=axes[3], fraction=0.046, pad=0.04)
    colorbar.set_label("0 = background; 1..K = normalized component")
    axes[3].set_title("Normalized components")
    for axis in (axes[0], axes[1], axes[3]):
        axis.set_xlabel("X")
        axis.set_ylabel("Y")

    bic_text = ", ".join(f"K={key}: {value:,.0f}" for key, value in bics.items())
    fig.suptitle(f"{name} | median-polish grid normalization | BIC {bic_text}", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def analyze(name, table, args, output_dir):
    threshold = otsu_threshold(
        table.loc[table["covered_cpgs"] > 0, "log10_cpgs"].to_numpy()
    )
    foreground = table["log10_cpgs"] >= threshold
    normalized, row_effects, column_effects, overall, iterations = median_polish(
        table,
        foreground,
        args.grid_x,
        args.grid_y,
        args.polish_iterations,
        args.polish_tolerance,
    )
    table = table.copy()
    table["normalized_log10_cpgs"] = normalized
    table["row_log10_effect"] = row_effects[table["y"].to_numpy(dtype=int) - 1]
    table["column_log10_effect"] = column_effects[
        table["x"].to_numpy(dtype=int) - 1
    ]

    values = table.loc[foreground, "normalized_log10_cpgs"].to_numpy()
    models = fit_models(values, args.max_components, args.seed)
    best_k = min(models, key=lambda key: models[key]["bic"])
    model = models[best_k]["model"]
    details = sorted_component_details(model)
    assignments = classify(table, foreground, model, details)
    assignments.to_csv(
        output_dir / f"{name}_grid_normalized_assignments.tsv", sep="\t", index=False
    )
    bics = {key: value["bic"] for key, value in models.items()}
    plot_sample(
        output_dir / f"{name}_grid_normalized_occupancy.png",
        name,
        assignments,
        model,
        details,
        bics,
        args.grid_x,
        args.grid_y,
    )

    raw_assignments = assignments.copy()
    raw_r2 = grid_effect_r2(raw_assignments)
    normalized_assignments = assignments.copy()
    normalized_assignments["log10_cpgs"] = normalized_assignments[
        "normalized_log10_cpgs"
    ]
    normalized_r2 = grid_effect_r2(normalized_assignments)
    expected = np.arange(1, best_k + 1, dtype=float)
    ratio_error = (
        float(np.mean(np.abs(details["ratios"] - expected) / expected))
        if best_k > 1
        else math.nan
    )
    return {
        "sample": name,
        "foreground_pixels": int(foreground.sum()),
        "foreground_threshold_cpgs": float(np.power(10.0, threshold) - 1),
        "best_components": best_k,
        "bic_improvement_over_k1": float(bics[1] - bics[best_k]),
        "component_centers": ",".join(f"{value:.1f}" for value in details["centers"]),
        "component_ratios": ",".join(f"{value:.3f}" for value in details["ratios"]),
        "component_weights": ",".join(f"{value:.4f}" for value in details["weights"]),
        "integer_ratio_relative_error": ratio_error,
        "raw_row_effect_r2": raw_r2[0],
        "raw_column_effect_r2": raw_r2[1],
        "raw_row_column_adjusted_r2": raw_r2[2],
        "normalized_row_effect_r2": normalized_r2[0],
        "normalized_column_effect_r2": normalized_r2[1],
        "normalized_row_column_adjusted_r2": normalized_r2[2],
        "median_polish_overall_log10": overall,
        "median_polish_iterations": iterations,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", action="append", type=parse_sample, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--grid-x", type=int, default=50)
    parser.add_argument("--grid-y", type=int, default=50)
    parser.add_argument("--max-components", type=int, default=8)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--polish-iterations", type=int, default=50)
    parser.add_argument("--polish-tolerance", type=float, default=1e-6)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    inputs = {}
    for name, path in args.sample:
        if name in inputs:
            parser.error(f"duplicate sample name: {name}")
        inputs[name] = path
        table = load_grid(name, path, args.grid_x, args.grid_y)
        summaries.append(analyze(name, table, args, args.output_dir))
    summary = pd.DataFrame(summaries)
    summary.to_csv(
        args.output_dir / "grid_normalized_occupancy_summary.tsv", sep="\t", index=False
    )
    metadata = {
        "inputs": inputs,
        "grid": [args.grid_x, args.grid_y],
        "max_components": args.max_components,
        "seed": args.seed,
        "normalization": (
            "Tukey-style median polish on original-foreground log10 CpG counts; "
            "the original foreground mask is retained for normalized GMM fitting."
        ),
    }
    (args.output_dir / "grid_normalized_occupancy_run.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
