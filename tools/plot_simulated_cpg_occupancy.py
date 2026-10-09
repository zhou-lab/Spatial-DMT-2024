#!/usr/bin/env python3
"""Plot simulated occupancy truth and recovery on a real spatial mask."""

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from scipy.stats import norm
from sklearn.metrics import adjusted_rand_score

from simulate_cpg_occupancy_accessibility import (
    fit_gmms,
    simulate_counts,
    sorted_labels,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--assignments-output", type=Path, required=True)
    parser.add_argument("--baseline-cpgs", type=float, required=True)
    parser.add_argument("--capture-cv", type=float, default=0.1)
    parser.add_argument("--genome-cpgs", type=int, default=21_867_837)
    parser.add_argument("--seed", type=int, default=20260804)
    parser.add_argument("--max-components", type=int, default=6)
    parser.add_argument("--title", default="Occupancy simulation")
    args = parser.parse_args()

    template = pd.read_csv(args.template, sep="\t")
    required = {"x", "y", "foreground"}
    if not required.issubset(template.columns):
        raise ValueError(f"{args.template} must contain {sorted(required)}")
    foreground = template["foreground"].astype(bool).to_numpy()
    occupancies = np.array([1, 2, 4])
    weights = np.array([0.5, 0.3, 0.2])
    rng = np.random.default_rng(args.seed)
    truth = rng.choice(occupancies, foreground.sum(), p=weights)
    counts = simulate_counts(
        rng, args.genome_cpgs, args.baseline_cpgs, truth, args.capture_cv
    )
    log_counts = np.log10(counts + 1)
    models = fit_gmms(log_counts, args.max_components, args.seed)
    values = log_counts.reshape(-1, 1)
    bics = {k: model.bic(values) for k, model in models.items()}
    best_k = min(bics, key=bics.get)
    model = models[best_k]
    predictions = sorted_labels(model, log_counts)

    result = template.loc[:, ["x", "y"]].copy()
    result["simulated_covered_cpgs"] = 0
    result["true_occupancy"] = 0
    result["inferred_component"] = 0
    result.loc[foreground, "simulated_covered_cpgs"] = counts
    result.loc[foreground, "true_occupancy"] = truth
    result.loc[foreground, "inferred_component"] = predictions
    args.assignments_output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.assignments_output, sep="\t", index=False)

    means = model.means_.ravel()
    variances = model.covariances_.reshape(-1)
    component_weights = model.weights_.ravel()
    order = np.argsort(means)
    centers = np.power(10.0, means[order]) - 1
    ari = adjusted_rand_score(truth, predictions)
    accuracy = math.nan
    if best_k == len(occupancies):
        inferred_occupancy = occupancies[predictions - 1]
        accuracy = float(np.mean(inferred_occupancy == truth))

    fig, axes = plt.subplots(1, 4, figsize=(20, 4.8))
    colors = {1: "#2166ac", 2: "#1b9e77", 4: "#b2182b"}
    bins = np.linspace(log_counts.min(), log_counts.max(), 65)
    for occupancy in occupancies:
        axes[0].hist(
            log_counts[truth == occupancy],
            bins=bins,
            density=True,
            histtype="step",
            linewidth=1.8,
            color=colors[occupancy],
            label=f"true {occupancy}x",
        )
    domain = np.linspace(log_counts.min(), log_counts.max(), 500)
    for index in order:
        density = component_weights[index] * norm.pdf(
            domain, means[index], math.sqrt(variances[index])
        )
        axes[0].plot(domain, density, color="#222222", alpha=0.65, linewidth=1)
    axes[0].set_xlabel("log10(covered CpGs + 1)")
    axes[0].set_ylabel("Density")
    axes[0].set_title("Simulated truth and fitted components")
    axes[0].legend(frameon=False, fontsize=8)

    count_matrix = result.pivot(
        index="y", columns="x", values="simulated_covered_cpgs"
    )
    count_cmap = plt.get_cmap("viridis").copy()
    count_cmap.set_bad("#f2f2f2")
    count_image = axes[1].imshow(
        np.log10(count_matrix.where(count_matrix > 0)),
        origin="upper",
        interpolation="nearest",
        cmap=count_cmap,
    )
    fig.colorbar(count_image, ax=axes[1], fraction=0.046, pad=0.04)
    axes[1].set_title("Simulated CpG coverage")

    truth_cmap = ListedColormap(["#2c7fb8", colors[1], colors[2], colors[4]])
    truth_norm = BoundaryNorm([-0.5, 0.5, 1.5, 3, 5], truth_cmap.N)
    truth_matrix = result.pivot(index="y", columns="x", values="true_occupancy")
    truth_image = axes[2].imshow(
        truth_matrix,
        origin="upper",
        interpolation="nearest",
        cmap=truth_cmap,
        norm=truth_norm,
    )
    truth_bar = fig.colorbar(
        truth_image, ax=axes[2], fraction=0.046, pad=0.04, ticks=[0, 1, 2, 4]
    )
    truth_bar.set_label("True occupancy")
    axes[2].set_title("True occupancy")

    inferred_matrix = result.pivot(
        index="y", columns="x", values="inferred_component"
    )
    inferred_image = axes[3].imshow(
        inferred_matrix,
        origin="upper",
        interpolation="nearest",
        cmap="tab10",
        vmin=0,
        vmax=max(4, best_k),
    )
    inferred_bar = fig.colorbar(
        inferred_image, ax=axes[3], fraction=0.046, pad=0.04
    )
    inferred_bar.set_label("0 = background; 1..K = inferred component")
    axes[3].set_title("Inferred components")
    for axis in axes[1:]:
        axis.set_xlabel("X")
        axis.set_ylabel("Y")

    bic_text = ", ".join(f"K={k}: {value:,.0f}" for k, value in bics.items())
    score_text = f"ARI={ari:.3f}"
    if np.isfinite(accuracy):
        score_text += f", accuracy={accuracy:.3f}"
    fig.suptitle(
        f"{args.title} | baseline={args.baseline_cpgs:,.0f}, "
        f"capture CV={args.capture_cv:.2f}, {score_text} | BIC {bic_text}",
        fontsize=10,
    )
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180, bbox_inches="tight")
    plt.close(fig)

    print(
        f"best_k={best_k} centers="
        + ",".join(f"{center:.1f}" for center in centers)
        + f" ari={ari:.4f} accuracy={accuracy:.4f}"
    )


if __name__ == "__main__":
    main()
