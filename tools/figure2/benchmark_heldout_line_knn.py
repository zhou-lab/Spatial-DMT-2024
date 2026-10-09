#!/usr/bin/env python3
"""Benchmark MRMP- and raw-DNA-PC kNN beta reconstruction on held-out lines.

The target line is selected from read/CpG coverage continuity only.  Its beta
values are then removed from the donor pool, predicted as the unweighted mean
of 20 neighbors, and compared with the untouched call-weighted beta values.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.stats import pearsonr, spearmanr
from sklearn.neighbors import NearestNeighbors

from plot_style import apply_style


ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class SampleSpec:
    sample: str
    label: str
    line_axis: str
    within_axis: str
    known_scar_axis: str
    known_scar_index: int
    reverse_map: bool


SAMPLES = (
    SampleSpec(
        sample="SpMETSLB14DM",
        label="SpMETSLB14DM · P21 brain",
        line_axis="y",
        within_axis="x",
        known_scar_axis="y",
        known_scar_index=19,
        reverse_map=False,
    ),
    SampleSpec(
        sample="E85_20um",
        label="E85_20um · E8.5 embryo",
        line_axis="paper_A",
        within_axis="paper_B",
        known_scar_axis="paper_B",
        known_scar_index=14,
        reverse_map=True,
    ),
)


def load_mrmp_frame(mrmp_dir: Path, spec: SampleSpec) -> pd.DataFrame:
    path = mrmp_dir / f"{spec.sample}_mrmp_pc_beta_values.tsv"
    frame = pd.read_csv(path, sep="\t", dtype={"spatial_id": str})
    frame["spatial_id"] = frame.spatial_id.str.zfill(4)
    required = {
        "spatial_id",
        "raw_beta",
        "raw_reads",
        "covered_cpgs",
        "mrmp_donor_eligible",
        spec.line_axis,
        spec.within_axis,
        *{f"smoothed_mrmp_pc{i}" for i in range(1, 11)},
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if frame.spatial_id.duplicated().any():
        raise ValueError(f"{path}: duplicate spatial IDs")
    return frame


def load_dna_scores(wnn_dir: Path, spec: SampleSpec, order: pd.Series) -> np.ndarray:
    if spec.sample == "SpMETSLB14DM":
        score_path = wnn_dir / "slb14_dna_features/dna_pc_scores.tsv"
        scores = pd.read_csv(score_path, sep="\t", dtype={"cell": str})
        scores["spatial_id"] = scores.cell.str.zfill(4)
    else:
        score_path = wnn_dir / "e85_dna_knn/dna_pc_scores.tsv"
        meta_path = wnn_dir / "e85_dna_knn/pixel_values.tsv"
        scores = pd.read_csv(score_path, sep="\t", dtype={"cell": str})
        meta = pd.read_csv(
            meta_path, sep="\t", dtype={"Barcodes": str, "spatial_id": str}
        )
        meta["spatial_id"] = meta.spatial_id.str.zfill(4)
        scores = scores.merge(
            meta[["Barcodes", "spatial_id"]],
            left_on="cell",
            right_on="Barcodes",
            how="left",
            validate="one_to_one",
        )
    pc_columns = [f"mPC_{i}" for i in range(1, 11)]
    missing = set(pc_columns + ["spatial_id"]).difference(scores.columns)
    if missing:
        raise ValueError(f"{score_path}: missing columns {sorted(missing)}")
    if scores.spatial_id.isna().any() or scores.spatial_id.duplicated().any():
        raise ValueError(f"{score_path}: invalid spatial-ID mapping")
    aligned = scores.set_index("spatial_id").reindex(order)
    if aligned[pc_columns].isna().any().any():
        absent = order[aligned[pc_columns].isna().any(axis=1)].tolist()[:5]
        raise ValueError(f"{score_path}: absent/non-finite scores for {absent}")
    return aligned[pc_columns].to_numpy(dtype=np.float64)


def line_candidates(frame: pd.DataFrame, spec: SampleSpec) -> pd.DataFrame:
    rows: list[dict[str, float | int | bool | str]] = []
    values = sorted(frame[spec.line_axis].dropna().astype(int).unique())
    maximum_pixels = int(frame.groupby(spec.line_axis).size().max())
    minimum_matched = max(30, int(np.ceil(0.60 * maximum_pixels)))
    for line in values:
        focal = frame.loc[frame[spec.line_axis] == line].set_index(spec.within_axis)
        before = frame.loc[frame[spec.line_axis] == line - 1].set_index(spec.within_axis)
        after = frame.loc[frame[spec.line_axis] == line + 1].set_index(spec.within_axis)
        common = focal.index.intersection(before.index).intersection(after.index)
        if not len(common):
            continue
        adjacent = pd.concat([before.loc[common], after.loc[common]])
        read_ratio = float(
            np.nanmedian(focal.loc[common, "raw_reads"])
            / np.nanmedian(adjacent.raw_reads)
        )
        cpg_ratio = float(
            np.nanmedian(focal.loc[common, "covered_cpgs"])
            / np.nanmedian(adjacent.covered_cpgs)
        )
        eligible = (
            len(common) >= minimum_matched
            and 0.8 <= read_ratio <= 1.25
            and 0.8 <= cpg_ratio <= 1.25
            and line != spec.known_scar_index
        )
        rows.append(
            {
                "sample": spec.sample,
                "axis": spec.line_axis,
                "line": int(line),
                "pixels": len(focal),
                "matched_to_both_adjacent_lines": len(common),
                "minimum_matched_required": minimum_matched,
                "median_raw_reads": float(np.nanmedian(focal.raw_reads)),
                "median_covered_cpgs": float(np.nanmedian(focal.covered_cpgs)),
                "read_ratio_to_adjacent": read_ratio,
                "cpg_ratio_to_adjacent": cpg_ratio,
                "coverage_continuity_score": float(
                    abs(np.log(read_ratio)) + abs(np.log(cpg_ratio))
                ),
                "known_scar_line": bool(line == spec.known_scar_index),
                "eligible": eligible,
            }
        )
    result = pd.DataFrame(rows)
    if not result.eligible.any():
        raise ValueError(f"{spec.sample}: no eligible non-scar {spec.line_axis} line")
    selected_index = result.loc[result.eligible, "coverage_continuity_score"].idxmin()
    result["selected"] = False
    result.loc[selected_index, "selected"] = True
    return result


def predict_knn(
    scores: np.ndarray,
    beta: np.ndarray,
    spatial_ids: np.ndarray,
    target_mask: np.ndarray,
    donor_mask: np.ndarray,
    neighbors: int,
    sample: str,
    method: str,
) -> tuple[np.ndarray, pd.DataFrame]:
    if not np.isfinite(scores).all():
        raise ValueError(f"{sample}/{method}: non-finite PC scores")
    donors = np.flatnonzero(donor_mask & ~target_mask & np.isfinite(beta))
    queries = np.flatnonzero(target_mask)
    if len(donors) < neighbors:
        raise ValueError(f"{sample}/{method}: only {len(donors)} donors")
    search = NearestNeighbors(n_neighbors=neighbors, metric="euclidean").fit(scores[donors])
    distances, positions = search.kneighbors(scores[queries])
    selected = donors[positions]
    predictions = beta[selected].mean(axis=1)
    neighbor_rows = []
    for query_position, query_row in enumerate(queries):
        for rank in range(neighbors):
            donor_row = selected[query_position, rank]
            neighbor_rows.append(
                {
                    "sample": sample,
                    "method": method,
                    "query_spatial_id": spatial_ids[query_row],
                    "neighbor_rank": rank + 1,
                    "donor_spatial_id": spatial_ids[donor_row],
                    "pc_distance": distances[query_position, rank],
                    "donor_beta": beta[donor_row],
                }
            )
    return predictions, pd.DataFrame(neighbor_rows)


def metric_values(truth: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    error = prediction - truth
    mse = float(np.mean(np.square(error)))
    variance = float(np.sum(np.square(truth - truth.mean())))
    return {
        "mse": mse,
        "rmse": float(np.sqrt(mse)),
        "mae": float(np.mean(np.abs(error))),
        "mean_error": float(np.mean(error)),
        "pearson_r": float(pearsonr(truth, prediction).statistic),
        "spearman_r": float(spearmanr(truth, prediction).statistic),
        "r_squared": float(1 - np.sum(np.square(error)) / variance) if variance else np.nan,
        "truth_mean": float(np.mean(truth)),
        "prediction_mean": float(np.mean(prediction)),
        "within_1pp_fraction": float(np.mean(np.abs(error) <= 0.01)),
        "within_2pp_fraction": float(np.mean(np.abs(error) <= 0.02)),
    }


def bootstrap_metrics(
    truth: np.ndarray,
    prediction: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> dict[str, tuple[float, float]]:
    draws = rng.integers(0, len(truth), size=(iterations, len(truth)))
    errors = prediction[draws] - truth[draws]
    mse = np.mean(np.square(errors), axis=1)
    distributions = {
        "mse": mse,
        "rmse": np.sqrt(mse),
        "mae": np.mean(np.abs(errors), axis=1),
        "mean_error": np.mean(errors, axis=1),
    }
    return {
        name: tuple(np.quantile(values, [0.025, 0.975]).astype(float))
        for name, values in distributions.items()
    }


def paired_bootstrap_difference(
    truth: np.ndarray,
    mrmp: np.ndarray,
    raw_dna: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> list[dict[str, float | str]]:
    draws = rng.integers(0, len(truth), size=(iterations, len(truth)))
    mrmp_error = mrmp[draws] - truth[draws]
    dna_error = raw_dna[draws] - truth[draws]
    differences = {
        "mse": np.mean(np.square(mrmp_error), axis=1)
        - np.mean(np.square(dna_error), axis=1),
        "mae": np.mean(np.abs(mrmp_error), axis=1)
        - np.mean(np.abs(dna_error), axis=1),
    }
    rows = []
    for metric, values in differences.items():
        low, high = np.quantile(values, [0.025, 0.975])
        rows.append(
            {
                "metric": metric,
                "contrast": "MRMP-PC kNN minus raw-DNA-PC kNN",
                "difference": float(np.mean(values)),
                "ci_low": float(low),
                "ci_high": float(high),
                "probability_mrmp_lower_loss": float(np.mean(values < 0)),
            }
        )
    return rows


def beta_grid(frame: pd.DataFrame, spec: SampleSpec) -> tuple[np.ndarray, int, int]:
    x = frame[spec.within_axis].astype(int)
    y = frame[spec.line_axis].astype(int)
    matrix = np.full((y.max() - y.min() + 1, x.max() - x.min() + 1), np.nan)
    matrix[(y - y.min()).to_numpy(), (x - x.min()).to_numpy()] = frame.raw_beta
    if spec.reverse_map:
        matrix = matrix[::-1, ::-1]
    return matrix, int(x.min()), int(y.min())


def plot_results(
    datasets: list[dict], metrics: pd.DataFrame, output_dir: Path, dpi: int
) -> None:
    apply_style()
    colors = {"MRMP-PC kNN": "#0072B2", "raw-DNA-PC kNN": "#D55E00"}
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.2), layout="constrained")
    letters = iter("ABCDEF")
    for row, data in enumerate(datasets):
        spec: SampleSpec = data["spec"]
        target = data["target"].sort_values(spec.within_axis)
        matrix, x_min, y_min = beta_grid(data["frame"], spec)
        ax = axes[row, 0]
        image = ax.imshow(
            np.ma.masked_invalid(matrix),
            cmap="magma",
            vmin=np.nanquantile(matrix, 0.02),
            vmax=np.nanquantile(matrix, 0.98),
            interpolation="nearest",
            aspect="equal",
        )
        selected_line = data["selected_line"]
        displayed_row = selected_line - y_min
        if spec.reverse_map:
            displayed_row = matrix.shape[0] - 1 - displayed_row
        ax.add_patch(
            Rectangle(
                (-0.5, displayed_row - 0.5),
                matrix.shape[1],
                1,
                fill=False,
                edgecolor="#00E5FF",
                linewidth=1.4,
            )
        )
        ax.set_title("Held-out line selection", fontsize=9, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(image, ax=ax, fraction=0.047, pad=0.03, label="observed β")

        ax = axes[row, 1]
        ax.plot(
            target[spec.within_axis], target.ground_truth_beta,
            color="#222222", linewidth=1.3, marker="o", markersize=2.4,
            label="Ground truth",
        )
        for method in colors:
            field = data["prediction_fields"][method]
            ax.plot(
                target[spec.within_axis], target[field], color=colors[method],
                linewidth=1.0, alpha=0.9, label=method,
            )
        ax.set_xlabel(spec.within_axis)
        ax.set_ylabel("global CpG β")
        ax.set_title(
            f"Held-out {spec.line_axis}={selected_line} profile",
            fontsize=9,
            fontweight="bold",
        )
        ax.legend(fontsize=6.5, frameon=False, ncol=1)

        ax = axes[row, 2]
        eligible = target.evaluation_eligible.to_numpy(dtype=bool)
        truth = target.loc[eligible, "ground_truth_beta"].to_numpy()
        limits = [truth.min() - 0.005, truth.max() + 0.005]
        for method in colors:
            field = data["prediction_fields"][method]
            ax.scatter(
                truth,
                target.loc[eligible, field],
                s=13,
                alpha=0.67,
                color=colors[method],
                edgecolors="none",
                label=method,
            )
        ax.plot(limits, limits, color="#777777", linestyle="--", linewidth=0.8)
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("ground-truth β")
        ax.set_ylabel("predicted β")
        sample_metrics = metrics.loc[metrics["sample"] == spec.sample]
        annotation = "\n".join(
            f"{name.split('-')[0]} RMSE={sample_metrics.loc[sample_metrics.method == name, 'rmse'].iloc[0]:.4f}"
            for name in colors
        )
        ax.text(0.03, 0.97, annotation, transform=ax.transAxes, va="top", fontsize=6.8)
        ax.set_title("Prediction versus truth", fontsize=9, fontweight="bold")

        axes[row, 0].set_ylabel(spec.label, fontsize=8, fontweight="bold")
        for ax in axes[row]:
            ax.text(
                -0.12, 1.04, next(letters), transform=ax.transAxes,
                fontsize=11, fontweight="bold",
            )
    fig.suptitle(
        "Blinded non-scar-line reconstruction with 20-neighbor methylation models",
        fontsize=13,
        fontweight="bold",
    )
    fig.savefig(output_dir / "heldout_line_knn_reconstruction.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(output_dir / "heldout_line_knn_reconstruction.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mrmp-dir",
        type=Path,
        default=ROOT / "results/figure2/mrmp_smoothing_comparison",
    )
    parser.add_argument(
        "--wnn-dir",
        type=Path,
        default=ROOT / "results/figure2/wnn_smoothing_comparison",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results/figure2/heldout_line_reconstruction",
    )
    parser.add_argument("--neighbors", type=int, default=20)
    parser.add_argument("--bootstrap-iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261002)
    parser.add_argument("--dpi", type=int, default=360)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    all_candidates = []
    all_predictions = []
    all_neighbors = []
    all_metrics = []
    all_differences = []
    datasets = []

    for spec in SAMPLES:
        frame = load_mrmp_frame(args.mrmp_dir, spec)
        candidates = line_candidates(frame, spec)
        selected_line = int(candidates.loc[candidates.selected, "line"].iloc[0])
        all_candidates.append(candidates)
        target_mask = frame[spec.line_axis].to_numpy(dtype=int) == selected_line
        spatial_ids = frame.spatial_id.to_numpy()
        beta = frame.raw_beta.to_numpy(dtype=float)

        mrmp_columns = [f"smoothed_mrmp_pc{i}" for i in range(1, 11)]
        mrmp_scores = frame[mrmp_columns].to_numpy(dtype=float)
        dna_scores = load_dna_scores(args.wnn_dir, spec, frame.spatial_id)
        method_inputs = {
            "MRMP-PC kNN": (
                mrmp_scores,
                frame.mrmp_donor_eligible.to_numpy(dtype=bool)
                & np.isfinite(mrmp_scores).all(axis=1),
            ),
            "raw-DNA-PC kNN": (dna_scores, np.isfinite(dna_scores).all(axis=1)),
        }
        target = frame.loc[
            target_mask,
            [
                "spatial_id",
                spec.line_axis,
                spec.within_axis,
                "raw_reads",
                "covered_cpgs",
                "raw_beta",
            ],
        ].copy()
        target.rename(columns={"raw_beta": "ground_truth_beta"}, inplace=True)
        target["sample"] = spec.sample
        target["evaluation_eligible"] = ~(
            (spec.known_scar_axis != spec.line_axis)
            & (target[spec.known_scar_axis] == spec.known_scar_index)
        ) if spec.known_scar_axis in target else True
        prediction_fields = {}
        for method, (scores, donor_mask) in method_inputs.items():
            prediction, neighbor_frame = predict_knn(
                scores,
                beta,
                spatial_ids,
                target_mask,
                donor_mask,
                args.neighbors,
                spec.sample,
                method,
            )
            field = "predicted_beta_" + method.lower().replace("-", "_").replace(" ", "_")
            prediction_fields[method] = field
            target[field] = prediction
            target[f"error_{field.removeprefix('predicted_beta_')}"] = (
                prediction - target.ground_truth_beta.to_numpy()
            )
            all_neighbors.append(neighbor_frame)

        eligible = target.evaluation_eligible.to_numpy(dtype=bool)
        truth = target.loc[eligible, "ground_truth_beta"].to_numpy(dtype=float)
        for method, field in prediction_fields.items():
            prediction = target.loc[eligible, field].to_numpy(dtype=float)
            values = metric_values(truth, prediction)
            intervals = bootstrap_metrics(
                truth, prediction, args.bootstrap_iterations, rng
            )
            row = {
                "sample": spec.sample,
                "heldout_axis": spec.line_axis,
                "heldout_line": selected_line,
                "method": method,
                "neighbors": args.neighbors,
                "pcs": 10,
                "evaluated_pixels": len(truth),
                **values,
            }
            for name, (low, high) in intervals.items():
                row[f"{name}_ci_low"] = low
                row[f"{name}_ci_high"] = high
            all_metrics.append(row)
        difference_rows = paired_bootstrap_difference(
            truth,
            target.loc[eligible, prediction_fields["MRMP-PC kNN"]].to_numpy(),
            target.loc[eligible, prediction_fields["raw-DNA-PC kNN"]].to_numpy(),
            args.bootstrap_iterations,
            rng,
        )
        for row in difference_rows:
            row.update(
                {
                    "sample": spec.sample,
                    "heldout_axis": spec.line_axis,
                    "heldout_line": selected_line,
                    "evaluated_pixels": len(truth),
                }
            )
        all_differences.extend(difference_rows)
        all_predictions.append(target)
        datasets.append(
            {
                "spec": spec,
                "frame": frame,
                "target": target,
                "selected_line": selected_line,
                "prediction_fields": prediction_fields,
            }
        )

    candidate_frame = pd.concat(all_candidates, ignore_index=True)
    prediction_frame = pd.concat(all_predictions, ignore_index=True)
    neighbor_frame = pd.concat(all_neighbors, ignore_index=True)
    metric_frame = pd.DataFrame(all_metrics)
    difference_frame = pd.DataFrame(all_differences)

    neighbor_diagnostics = []
    for data in datasets:
        spec = data["spec"]
        coordinates = data["frame"].set_index("spatial_id")[
            [spec.within_axis, spec.line_axis]
        ]
        sample_neighbors = neighbor_frame.loc[
            neighbor_frame["sample"] == spec.sample
        ].copy()
        sample_neighbors = sample_neighbors.join(
            coordinates.rename(
                columns={spec.within_axis: "query_x", spec.line_axis: "query_y"}
            ),
            on="query_spatial_id",
        ).join(
            coordinates.rename(
                columns={spec.within_axis: "donor_x", spec.line_axis: "donor_y"}
            ),
            on="donor_spatial_id",
        )
        sample_neighbors["grid_distance"] = np.hypot(
            sample_neighbors.query_x - sample_neighbors.donor_x,
            sample_neighbors.query_y - sample_neighbors.donor_y,
        )
        donor_is_known_scar = (
            sample_neighbors.donor_y == spec.known_scar_index
            if spec.known_scar_axis == spec.line_axis
            else sample_neighbors.donor_x == spec.known_scar_index
        )
        sample_neighbors["donor_is_known_scar"] = donor_is_known_scar
        for method, group in sample_neighbors.groupby("method"):
            neighbor_diagnostics.append(
                {
                    "sample": spec.sample,
                    "method": method,
                    "edges": len(group),
                    "median_pc_distance": float(group.pc_distance.median()),
                    "median_grid_distance": float(group.grid_distance.median()),
                    "mean_grid_distance": float(group.grid_distance.mean()),
                    "fraction_grid_distance_le_sqrt2": float(
                        np.mean(group.grid_distance <= np.sqrt(2))
                    ),
                    "known_scar_donor_fraction": float(
                        group.donor_is_known_scar.mean()
                    ),
                    "donor_beta_sd": float(group.donor_beta.std(ddof=1)),
                }
            )
    neighbor_diagnostic_frame = pd.DataFrame(neighbor_diagnostics)

    # Leakage and arithmetic assertions are part of the benchmark contract.
    target_ids = prediction_frame.groupby("sample").spatial_id.agg(set).to_dict()
    leaked = neighbor_frame.apply(
        lambda row: row.donor_spatial_id in target_ids[row["sample"]], axis=1
    )
    if leaked.any():
        raise AssertionError("Held-out-line beta leaked into the kNN donor pool")
    counts = neighbor_frame.groupby(["sample", "method", "query_spatial_id"]).size()
    if not (counts == args.neighbors).all():
        raise AssertionError("A prediction did not use exactly the requested neighbor count")
    recomputed = (
        neighbor_frame.groupby(["sample", "method", "query_spatial_id"]).donor_beta.mean()
    )
    for data in datasets:
        for method, field in data["prediction_fields"].items():
            check = recomputed.loc[data["spec"].sample, method].reindex(
                data["target"].spatial_id
            )
            if not np.allclose(check, data["target"][field], atol=1e-12):
                raise AssertionError(f"{data['spec'].sample}/{method}: neighbor mean mismatch")

    candidate_frame.to_csv(args.output_dir / "line_candidates.tsv", sep="\t", index=False)
    prediction_frame.to_csv(
        args.output_dir / "heldout_line_predictions.tsv", sep="\t", index=False
    )
    neighbor_frame.to_csv(
        args.output_dir / "heldout_line_neighbors.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    metric_frame.to_csv(
        args.output_dir / "heldout_line_reconstruction_metrics.tsv", sep="\t", index=False
    )
    difference_frame.to_csv(
        args.output_dir / "heldout_line_paired_differences.tsv", sep="\t", index=False
    )
    neighbor_diagnostic_frame.to_csv(
        args.output_dir / "heldout_line_neighbor_diagnostics.tsv",
        sep="\t",
        index=False,
    )
    plot_results(datasets, metric_frame, args.output_dir, args.dpi)

    print(candidate_frame.loc[candidate_frame.selected].to_string(index=False))
    print(metric_frame.to_string(index=False))
    print(difference_frame.to_string(index=False))
    print(neighbor_diagnostic_frame.to_string(index=False))


if __name__ == "__main__":
    main()
