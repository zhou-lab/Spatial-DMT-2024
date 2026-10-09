#!/usr/bin/env python3
"""Estimate per-pixel library complexity and probe discrete occupancy modes."""

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from scipy.optimize import brentq
from scipy.stats import norm

from cpg_coverage_occupancy import (
    fit_models,
    grid_effect_r2,
    otsu_threshold,
    parse_sample,
    sorted_component_details,
)
from cpg_grid_normalized_occupancy import median_polish


def estimate_library_size(mapped, unique):
    """Maximum-likelihood library size under random sampling with replacement."""
    mapped = int(mapped)
    unique = int(unique)
    if mapped <= 0 or unique <= 0:
        return 0.0
    if unique >= mapped:
        return math.nan

    def residual(library_size):
        return library_size * (-math.expm1(-mapped / library_size)) - unique

    lower = float(unique)
    duplicates = mapped - unique
    upper = max(lower * 2, mapped * mapped / (2.0 * duplicates) * 2)
    while residual(upper) <= 0:
        upper *= 2
    return float(brentq(residual, lower, upper, xtol=1e-6, rtol=1e-12))


def load_complexity(name, path, grid_x, grid_y):
    table = pd.read_csv(path, sep="\t")
    required = {"x", "y", "mapped_reads", "dup_reads"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    table = table.loc[:, ["x", "y", "mapped_reads", "dup_reads"]].copy()
    for column in table.columns:
        table[column] = pd.to_numeric(table[column], errors="raise").astype(int)
    if table.duplicated(["x", "y"]).any():
        raise ValueError(f"{path} contains duplicate spatial coordinates")
    if (table["dup_reads"] > table["mapped_reads"]).any():
        raise ValueError(f"{path} contains duplicate counts above mapped counts")

    full = pd.MultiIndex.from_product(
        [range(1, grid_x + 1), range(1, grid_y + 1)], names=["x", "y"]
    ).to_frame(index=False)
    table = full.merge(table, on=["x", "y"], how="left").fillna(0)
    table[["mapped_reads", "dup_reads"]] = table[
        ["mapped_reads", "dup_reads"]
    ].astype(int)
    table.insert(0, "sample", name)
    table["unique_fragments"] = table["mapped_reads"] - table["dup_reads"]
    table["library_complexity"] = [
        estimate_library_size(mapped, unique)
        for mapped, unique in zip(table["mapped_reads"], table["unique_fragments"])
    ]
    table["estimated_unseen_molecules"] = (
        table["library_complexity"] - table["unique_fragments"]
    )
    table["observed_library_fraction"] = np.divide(
        table["unique_fragments"],
        table["library_complexity"],
        out=np.full(len(table), np.nan),
        where=table["library_complexity"].to_numpy() > 0,
    )
    table["log10_library_complexity"] = np.log10(
        table["library_complexity"] + 1
    )
    return table


def fit_and_classify(values, max_components, seed):
    models = fit_models(values, max_components, seed)
    best_k = min(models, key=lambda key: models[key]["bic"])
    model = models[best_k]["model"]
    details = sorted_component_details(model)
    raw_labels = model.predict(np.asarray(values).reshape(-1, 1))
    labels = np.array(
        [details["old_to_new"][int(label)] for label in raw_labels], dtype=int
    )
    probabilities = model.predict_proba(np.asarray(values).reshape(-1, 1)).max(axis=1)
    return models, best_k, model, details, labels, probabilities


def grid_r2(table, value_column, foreground):
    work = table.loc[:, ["x", "y"]].copy()
    work["log10_cpgs"] = table[value_column]
    work["foreground"] = foreground.astype(int)
    return grid_effect_r2(work)


def draw_histogram(axis, values, model, details, title):
    axis.hist(values, bins=60, density=True, color="#a7b6c2", edgecolor="white")
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
        axis.plot(domain, density, linewidth=2, label=f"C{index}: {center:,.0f}")
    axis.set_xlabel("log10(estimated library molecules + 1)")
    axis.set_ylabel("Density")
    axis.set_title(title)
    axis.legend(frameon=False, fontsize=8)


def draw_grid(axis, table, value_column, title, cmap, vmin=None, vmax=None):
    matrix = table.pivot(index="y", columns="x", values=value_column)
    image = axis.imshow(
        matrix,
        origin="upper",
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extent=[0.5, matrix.shape[1] + 0.5, matrix.shape[0] + 0.5, 0.5],
    )
    axis.set_xlabel("X")
    axis.set_ylabel("Y")
    axis.set_title(title)
    return image


def analyze(name, table, args):
    finite_positive = np.isfinite(table["library_complexity"]) & (
        table["library_complexity"] > 0
    )
    threshold = otsu_threshold(
        table.loc[finite_positive, "log10_library_complexity"].to_numpy()
    )
    foreground = finite_positive & (table["log10_library_complexity"] >= threshold)
    if foreground.sum() < max(20, args.max_components * 5):
        raise ValueError(f"{name} has too few foreground pixels")

    raw_values = table.loc[foreground, "log10_library_complexity"].to_numpy()
    raw = fit_and_classify(raw_values, args.max_components, args.seed)
    raw_models, raw_k, raw_model, raw_details, raw_labels, raw_probabilities = raw

    polish_input = table.copy()
    polish_input["log10_cpgs"] = polish_input["log10_library_complexity"]
    normalized, row_effects, column_effects, overall, iterations = median_polish(
        polish_input,
        foreground,
        args.grid_x,
        args.grid_y,
        args.polish_iterations,
        args.polish_tolerance,
    )
    table = table.copy()
    table["normalized_log10_library_complexity"] = normalized
    table["row_log10_effect"] = row_effects[table["y"].to_numpy() - 1]
    table["column_log10_effect"] = column_effects[table["x"].to_numpy() - 1]
    normalized_values = table.loc[
        foreground, "normalized_log10_library_complexity"
    ].to_numpy()
    normalized_fit = fit_and_classify(
        normalized_values, args.max_components, args.seed
    )
    (
        normalized_models,
        normalized_k,
        normalized_model,
        normalized_details,
        normalized_labels,
        normalized_probabilities,
    ) = normalized_fit

    table["foreground"] = foreground.astype(int)
    table["raw_component"] = 0
    table["raw_component_probability"] = np.nan
    table["normalized_component"] = 0
    table["normalized_component_probability"] = np.nan
    table.loc[foreground, "raw_component"] = raw_labels
    table.loc[foreground, "raw_component_probability"] = raw_probabilities
    table.loc[foreground, "normalized_component"] = normalized_labels
    table.loc[
        foreground, "normalized_component_probability"
    ] = normalized_probabilities

    raw_r2 = grid_r2(table, "log10_library_complexity", foreground)
    normalized_r2 = grid_r2(
        table, "normalized_log10_library_complexity", foreground
    )
    table.to_csv(
        args.output_dir / f"{name}_library_complexity_assignments.tsv",
        sep="\t",
        index=False,
    )

    figure, axes = plt.subplots(2, 3, figsize=(19, 10))
    draw_histogram(
        axes[0, 0], raw_values, raw_model, raw_details, f"Raw mixture (K={raw_k})"
    )
    finite_range = table.loc[foreground, "log10_library_complexity"]
    lower, upper = np.quantile(finite_range, [0.01, 0.99])
    image = draw_grid(
        axes[0, 1],
        table,
        "log10_library_complexity",
        "Raw library complexity",
        "viridis",
        lower,
        upper,
    )
    figure.colorbar(image, ax=axes[0, 1], fraction=0.046, pad=0.04)
    component_cmap = ListedColormap(
        ["#f2f2f2", "#2166ac", "#1b9e77", "#b2182b", "#762a83", "#e08214", "#35978f", "#8073ac", "#4d4d4d"]
    )
    image = draw_grid(
        axes[0, 2], table, "raw_component", "Raw components", component_cmap, 0, 8
    )
    figure.colorbar(image, ax=axes[0, 2], fraction=0.046, pad=0.04)

    draw_histogram(
        axes[1, 0],
        normalized_values,
        normalized_model,
        normalized_details,
        f"Row/column-normalized mixture (K={normalized_k})",
    )
    normalized_range = table.loc[
        foreground, "normalized_log10_library_complexity"
    ]
    lower, upper = np.quantile(normalized_range, [0.01, 0.99])
    image = draw_grid(
        axes[1, 1],
        table,
        "normalized_log10_library_complexity",
        "Normalized library complexity",
        "viridis",
        lower,
        upper,
    )
    figure.colorbar(image, ax=axes[1, 1], fraction=0.046, pad=0.04)
    image = draw_grid(
        axes[1, 2],
        table,
        "normalized_component",
        "Normalized components",
        component_cmap,
        0,
        8,
    )
    figure.colorbar(image, ax=axes[1, 2], fraction=0.046, pad=0.04)
    figure.suptitle(
        f"{name}: per-pixel estimated library complexity | "
        f"grid adjusted R2 {raw_r2[2]:.3f} -> {normalized_r2[2]:.3f}",
        fontsize=12,
    )
    figure.tight_layout()
    figure.savefig(
        args.output_dir / f"{name}_library_complexity_occupancy.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(figure)

    return {
        "sample": name,
        "pixels_total": len(table),
        "finite_complexity_pixels": int(finite_positive.sum()),
        "foreground_pixels": int(foreground.sum()),
        "foreground_threshold_molecules": float(10**threshold - 1),
        "median_mapped_fragments_foreground": float(
            table.loc[foreground, "mapped_reads"].median()
        ),
        "median_unique_fragments_foreground": float(
            table.loc[foreground, "unique_fragments"].median()
        ),
        "median_library_complexity_foreground": float(
            table.loc[foreground, "library_complexity"].median()
        ),
        "median_observed_library_fraction": float(
            table.loc[foreground, "observed_library_fraction"].median()
        ),
        "raw_best_components": raw_k,
        "raw_bic_improvement_over_k1": float(
            raw_models[1]["bic"] - raw_models[raw_k]["bic"]
        ),
        "raw_component_centers": ",".join(
            f"{value:.1f}" for value in raw_details["centers"]
        ),
        "raw_component_ratios": ",".join(
            f"{value:.3f}" for value in raw_details["ratios"]
        ),
        "raw_row_column_adjusted_r2": raw_r2[2],
        "normalized_best_components": normalized_k,
        "normalized_bic_improvement_over_k1": float(
            normalized_models[1]["bic"]
            - normalized_models[normalized_k]["bic"]
        ),
        "normalized_component_centers": ",".join(
            f"{value:.1f}" for value in normalized_details["centers"]
        ),
        "normalized_component_ratios": ",".join(
            f"{value:.3f}" for value in normalized_details["ratios"]
        ),
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
    parser.add_argument(
        "--count-definition",
        default=(
            "Primary nonsupplementary R1 fragments; mapped and duplicate counts "
            "from retained barcode_summary output generated from the merged BAM."
        ),
        help="Provenance definition for the mapped and duplicate input counts.",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    inputs = {}
    for name, path in args.sample:
        if name in inputs:
            parser.error(f"duplicate sample name: {name}")
        inputs[name] = path
        table = load_complexity(name, path, args.grid_x, args.grid_y)
        summaries.append(analyze(name, table, args))

    summary = pd.DataFrame(summaries)
    summary.to_csv(
        args.output_dir / "library_complexity_occupancy_summary.tsv",
        sep="\t",
        index=False,
    )
    metadata = {
        "inputs": inputs,
        "grid": [args.grid_x, args.grid_y],
        "max_components": args.max_components,
        "seed": args.seed,
        "estimator": "unique = L * (1 - exp(-mapped / L))",
        "count_definition": args.count_definition,
        "normalization": "Tukey-style median polish on foreground log10 complexity.",
    }
    (args.output_dir / "library_complexity_occupancy_run.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
