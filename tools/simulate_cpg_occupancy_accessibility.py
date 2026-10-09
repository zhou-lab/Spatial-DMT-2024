#!/usr/bin/env python3
"""Simulate when sparse covered-CpG counts identify occupancy states."""

import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture


DEFAULT_DEPTHS = {
    "SpMETSLE14D": (5047.6, 1841),
    "SpMETSLE16D": (5252.8, 2120),
    "SpMETSLE14DM": (217444.2, 1836),
}
DEFAULT_CVS = (0.0, 0.1, 0.2, 0.24, 0.28, 0.35, 0.45, 0.6, 0.8)
DEPTH_SWEEP = (10, 30, 100, 300, 1000, 3000, 5000, 10000, 30000, 100000, 220000)


def parse_scenario(value):
    try:
        name, parameters = value.split("=", 1)
        baseline, pixels = parameters.split(",", 1)
        return name, (float(baseline), int(pixels))
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "scenario must be NAME=BASELINE_CPGS,PIXELS"
        ) from error


def capture_probabilities(genome_cpgs, baseline_cpgs, occupancy, capture_factor):
    draws = baseline_cpgs * occupancy * capture_factor
    return -np.expm1(-draws / genome_cpgs)


def simulate_counts(rng, genome_cpgs, baseline_cpgs, occupancy, capture_cv):
    if capture_cv > 0:
        sigma = math.sqrt(math.log1p(capture_cv**2))
        capture_factor = rng.lognormal(-0.5 * sigma**2, sigma, len(occupancy))
    else:
        capture_factor = np.ones(len(occupancy))
    probabilities = capture_probabilities(
        genome_cpgs, baseline_cpgs, occupancy, capture_factor
    )
    return rng.binomial(genome_cpgs, probabilities)


def fit_gmms(log_counts, max_components, seed):
    values = log_counts.reshape(-1, 1)
    models = {}
    for components in range(1, max_components + 1):
        model = GaussianMixture(
            components,
            covariance_type="full",
            n_init=2,
            reg_covar=1e-4,
            random_state=seed,
        ).fit(values)
        models[components] = model
    return models


def sorted_labels(model, values):
    order = np.argsort(model.means_.ravel())
    mapping = np.empty(len(order), dtype=int)
    mapping[order] = np.arange(1, len(order) + 1)
    return mapping[model.predict(values.reshape(-1, 1))]


def run_replicate(task):
    (
        scenario,
        baseline_cpgs,
        pixels,
        capture_cv,
        replicate,
        genome_cpgs,
        occupancies,
        weights,
        max_components,
        seed,
        fit_mixtures,
    ) = task
    rng = np.random.default_rng(seed + replicate)
    truth = rng.choice(occupancies, pixels, p=weights)
    counts = simulate_counts(
        rng, genome_cpgs, baseline_cpgs, truth, capture_cv
    )
    log_counts = np.log10(counts + 1)
    truth_labels = np.searchsorted(np.asarray(occupancies), truth) + 1
    expected_counts, _ = theoretical_centers(
        genome_cpgs, baseline_cpgs, occupancies
    )
    expected_logs = np.log10(expected_counts + 1)
    oracle_labels = np.argmin(
        np.abs(log_counts[:, np.newaxis] - expected_logs[np.newaxis, :]), axis=1
    ) + 1
    oracle_accuracy = float(np.mean(truth_labels == oracle_labels))

    selected_k = math.nan
    selected_ari = math.nan
    forced_ari = math.nan
    forced_accuracy = math.nan
    bic_margin = math.nan
    if fit_mixtures:
        models = fit_gmms(log_counts, max_components, seed + replicate)
        values = log_counts.reshape(-1, 1)
        bics = {k: model.bic(values) for k, model in models.items()}
        selected_k = min(bics, key=bics.get)
        selected_labels = sorted_labels(models[selected_k], log_counts)
        forced_labels = sorted_labels(models[len(occupancies)], log_counts)
        selected_ari = adjusted_rand_score(truth_labels, selected_labels)
        forced_ari = adjusted_rand_score(truth_labels, forced_labels)
        forced_accuracy = float(np.mean(truth_labels == forced_labels))
        bic_margin = (
            min(bic for k, bic in bics.items() if k != len(occupancies))
            - bics[len(occupancies)]
        )
    return {
        "scenario": scenario,
        "baseline_cpgs": baseline_cpgs,
        "pixels": pixels,
        "capture_cv": capture_cv,
        "replicate": replicate,
        "selected_components": selected_k,
        "selected_ari": selected_ari,
        "forced_3_ari": forced_ari,
        "forced_3_accuracy": forced_accuracy,
        "oracle_accuracy": oracle_accuracy,
        "bic_margin_for_k3": bic_margin,
    }


def summarize(raw):
    return (
        raw.groupby(["scenario", "baseline_cpgs", "pixels", "capture_cv"])
        .agg(
            replicates=("replicate", "size"),
            select_k3_rate=("selected_components", lambda x: np.mean(x == 3)),
            median_selected_k=("selected_components", "median"),
            median_selected_ari=("selected_ari", "median"),
            median_forced_3_ari=("forced_3_ari", "median"),
            median_forced_3_accuracy=("forced_3_accuracy", "median"),
            median_oracle_accuracy=("oracle_accuracy", "median"),
            median_bic_margin_for_k3=("bic_margin_for_k3", "median"),
        )
        .reset_index()
    )


def theoretical_centers(genome_cpgs, baseline_cpgs, occupancies):
    occupancy = np.asarray(occupancies, dtype=float)
    centers = genome_cpgs * (-np.expm1(-baseline_cpgs * occupancy / genome_cpgs))
    return centers, centers / centers[0]


def plot_results(
    path,
    summary,
    genome_cpgs,
    occupancies,
    observed_cv_range,
    scenario_depths,
):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    observed = summary[summary["scenario"].isin(scenario_depths)]
    palette = plt.get_cmap("tab10")
    colors = {
        name: palette(index % 10) for index, name in enumerate(scenario_depths)
    }
    for name, group in observed.groupby("scenario", sort=False):
        group = group.sort_values("capture_cv")
        axes[0].plot(
            group["capture_cv"],
            group["select_k3_rate"],
            marker="o",
            label=name,
            color=colors.get(name),
        )
        axes[1].plot(
            group["capture_cv"],
            group["median_oracle_accuracy"],
            marker="o",
            label=name,
            color=colors.get(name),
        )
    for axis in axes[:2]:
        axis.axvspan(*observed_cv_range, color="#777777", alpha=0.15)
        axis.set_ylim(-0.02, 1.02)
        axis.set_xlabel("Multiplicative capture CV")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Fraction selecting true K=3")
    axes[0].set_title("Unsupervised component recovery")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set_ylabel("Median 1x/2x/4x accuracy")
    axes[1].set_title("Classification with known sparse-count centers")

    depth = summary[summary["scenario"] == "depth_sweep"].copy()
    matrix = depth.pivot(
        index="capture_cv", columns="baseline_cpgs", values="median_oracle_accuracy"
    ).sort_index(ascending=False)
    image = axes[2].imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap="viridis",
        vmin=0,
        vmax=1,
    )
    axes[2].set_xticks(np.arange(len(matrix.columns)))
    axes[2].set_xticklabels(
        [f"{value:g}" for value in matrix.columns], rotation=45, ha="right"
    )
    axes[2].set_yticks(np.arange(len(matrix.index)))
    axes[2].set_yticklabels([f"{value:g}" for value in matrix.index])
    axes[2].set_xlabel("Covered CpGs for a 1x pixel")
    axes[2].set_ylabel("Multiplicative capture CV")
    axes[2].set_title("Depth/CV accessibility boundary")
    fig.colorbar(image, ax=axes[2], label="Median known-center accuracy")

    reference_name = next(iter(scenario_depths))
    centers, ratios = theoretical_centers(
        genome_cpgs, scenario_depths[reference_name][0], occupancies
    )
    fig.suptitle(
        "CpG occupancy accessibility simulation | shallow expected ratios "
        + ", ".join(f"{value:.3f}" for value in ratios),
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--genome-cpgs", type=int, default=21_867_837)
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260804)
    parser.add_argument("--max-components", type=int, default=6)
    parser.add_argument(
        "--observed-cv-range",
        type=float,
        nargs=2,
        metavar=("LOW", "HIGH"),
        default=(0.22, 0.28),
    )
    parser.add_argument(
        "--scenario",
        action="append",
        type=parse_scenario,
        help="Override defaults; repeat NAME=BASELINE_CPGS,PIXELS",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    occupancies = (1, 2, 4)
    weights = (0.5, 0.3, 0.2)
    scenario_depths = dict(args.scenario) if args.scenario else DEFAULT_DEPTHS
    tasks = []
    for name, (baseline, pixels) in scenario_depths.items():
        for cv in DEFAULT_CVS:
            for replicate in range(args.replicates):
                tasks.append(
                    (
                        name,
                        baseline,
                        pixels,
                        cv,
                        replicate,
                        args.genome_cpgs,
                        occupancies,
                        weights,
                        args.max_components,
                        args.seed + len(tasks) * 1000,
                        True,
                    )
                )
    for baseline in DEPTH_SWEEP:
        for cv in DEFAULT_CVS:
            for replicate in range(args.replicates):
                tasks.append(
                    (
                        "depth_sweep",
                        baseline,
                        2000,
                        cv,
                        replicate,
                        args.genome_cpgs,
                        occupancies,
                        weights,
                        args.max_components,
                        args.seed + len(tasks) * 1000,
                        False,
                    )
                )

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        raw = pd.DataFrame(executor.map(run_replicate, tasks, chunksize=10))
    summary = summarize(raw)
    raw.to_csv(args.output_dir / "simulation_replicates.tsv", sep="\t", index=False)
    summary.to_csv(args.output_dir / "simulation_summary.tsv", sep="\t", index=False)

    theory = []
    for name, (baseline, _) in scenario_depths.items():
        centers, ratios = theoretical_centers(
            args.genome_cpgs, baseline, occupancies
        )
        theory.append(
            {
                "scenario": name,
                "baseline_cpgs": baseline,
                "genome_fraction_1x": baseline / args.genome_cpgs,
                "expected_centers": ",".join(f"{value:.1f}" for value in centers),
                "expected_ratios": ",".join(f"{value:.4f}" for value in ratios),
                "four_x_saturation_loss": 1 - centers[-1] / (4 * centers[0]),
            }
        )
    pd.DataFrame(theory).to_csv(
        args.output_dir / "finite_genome_theory.tsv", sep="\t", index=False
    )
    plot_results(
        args.output_dir / "cpg_occupancy_accessibility.png",
        summary,
        args.genome_cpgs,
        occupancies,
        args.observed_cv_range,
        scenario_depths,
    )
    metadata = {
        "genome_cpgs": args.genome_cpgs,
        "occupancies": occupancies,
        "weights": weights,
        "capture_cv_values": DEFAULT_CVS,
        "observed_cv_range": args.observed_cv_range,
        "replicates": args.replicates,
        "seed": args.seed,
        "scenarios": scenario_depths,
        "model": (
            "Finite-genome unique-CpG capture with a mean-one lognormal "
            "pixel-level multiplier; GMMs are fitted to log10(count + 1)."
        ),
    }
    (args.output_dir / "simulation_run.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
