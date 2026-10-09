#!/usr/bin/env python3
"""Detect single-line read-count scars and test adjacent methylation similarity."""

from __future__ import annotations

import argparse
import subprocess
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from scipy.stats import ks_2samp, wilcoxon

from plot_style import apply_style
from spatial_utils import detect_paired_scar, detect_scar, score_paired_scars


SPATIAL_SAMPLE = "SpMETSLB14DM"
SMCSEQ_SAMPLE = "E85_20um"


def matrix_from_frame(
    frame: pd.DataFrame,
    row: str,
    column: str,
    value: str,
) -> np.ndarray:
    matrix = np.full(
        (int(frame[row].max()), int(frame[column].max())),
        np.nan,
        dtype=float,
    )
    matrix[
        frame[row].to_numpy(dtype=int) - 1,
        frame[column].to_numpy(dtype=int) - 1,
    ] = pd.to_numeric(frame[value], errors="coerce").to_numpy(dtype=float)
    return matrix


def ensure_spatial_summary(
    sample: str,
    repo_root: Path,
    run_root: Path,
    output_dir: Path,
) -> pd.DataFrame:
    cache_dir = output_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    summary_path = cache_dir / f"{sample}_stream.tsv"
    binary = cache_dir / "summarize_cg_stream"
    source = repo_root / "tools" / "figure2" / "summarize_cg_stream.c"
    cg_path = run_root / "final_output" / sample / f"{sample}.cg"
    if not binary.is_file():
        subprocess.run(
            ["cc", "-O3", str(source), "-lz", "-o", str(binary)],
            check=True,
        )
    if not summary_path.is_file():
        with summary_path.open("w") as handle:
            subprocess.run([str(binary), str(cg_path)], stdout=handle, check=True)
    summary = pd.read_csv(summary_path, sep="\t")
    index = pd.read_csv(f"{cg_path}.idx", sep="\t", header=None, dtype=str)
    if (
        len(summary) != len(index)
        or not index[0].is_unique
        or not index[1].is_unique
        or not pd.to_numeric(index[1]).is_monotonic_increasing
    ):
        raise ValueError("Spatial-DMT CG index does not preserve physical record order")
    summary["spatial_id"] = index[0].to_numpy()
    summary["x"] = summary["spatial_id"].str[:2].astype(int)
    summary["y"] = summary["spatial_id"].str[2:].astype(int)
    return summary


def bootstrap_median_ci(values: np.ndarray, seed: int = 20260918) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(10_000, len(values)), replace=True)
    medians = np.median(draws, axis=1)
    low, high = np.quantile(medians, [0.025, 0.975])
    return float(low), float(high)


def compare_line(
    count_matrix: np.ndarray,
    beta_matrix: np.ndarray,
    coverage_matrix: np.ndarray,
    axis: str,
    index: int,
) -> tuple[dict[str, object], np.ndarray, np.ndarray, np.ndarray]:
    count = count_matrix if axis == "horizontal" else count_matrix.T
    beta = beta_matrix if axis == "horizontal" else beta_matrix.T
    coverage = coverage_matrix if axis == "horizontal" else coverage_matrix.T
    center0 = index - 1

    count_valid = (
        np.isfinite(count[center0])
        & np.isfinite(count[center0 - 1])
        & np.isfinite(count[center0 + 1])
        & (count[center0 - 1] > 0)
        & (count[center0 + 1] > 0)
    )
    count_ratios = (count[center0, count_valid] + 1.0) / (
        np.sqrt(
            count[center0 - 1, count_valid]
            * count[center0 + 1, count_valid]
        )
        + 1.0
    )
    center_count = count[center0, count_valid]
    adjacent_count_pool = np.concatenate(
        [count[center0 - 1, count_valid], count[center0 + 1, count_valid]]
    )

    beta_valid = (
        np.isfinite(beta[center0])
        & np.isfinite(beta[center0 - 1])
        & np.isfinite(beta[center0 + 1])
    )
    center_beta = beta[center0, beta_valid]
    before_beta = beta[center0 - 1, beta_valid]
    after_beta = beta[center0 + 1, beta_valid]
    adjacent_mean = (before_beta + after_beta) / 2.0
    adjacent_pool = np.concatenate([before_beta, after_beta])
    paired_delta = center_beta - adjacent_mean
    ci_low, ci_high = bootstrap_median_ci(paired_delta)
    wilcoxon_result = wilcoxon(paired_delta)
    ks_result = ks_2samp(center_beta, adjacent_pool)
    coverage_ratio = (coverage[center0, beta_valid] + 1.0) / (
        (coverage[center0 - 1, beta_valid] + coverage[center0 + 1, beta_valid])
        / 2.0
        + 1.0
    )
    center_coverage = coverage[center0, beta_valid]
    adjacent_coverage_pool = np.concatenate(
        [coverage[center0 - 1, beta_valid], coverage[center0 + 1, beta_valid]]
    )
    metrics: dict[str, object] = {
        "scar_axis": axis,
        "scar_index": index,
        "count_support_pixels": int(count_valid.sum()),
        "median_count_ratio": float(np.median(count_ratios)),
        "scar_median_read_count": float(np.median(center_count)),
        "adjacent_pooled_median_read_count": float(np.median(adjacent_count_pool)),
        "fraction_count_ratio_below_0.55": float(np.mean(count_ratios < 0.55)),
        "methylation_paired_pixels": len(center_beta),
        "scar_median_beta": float(np.median(center_beta)),
        "adjacent_pooled_median_beta": float(np.median(adjacent_pool)),
        "median_paired_beta_delta": float(np.median(paired_delta)),
        "median_paired_beta_delta_ci_low": ci_low,
        "median_paired_beta_delta_ci_high": ci_high,
        "median_absolute_paired_beta_delta": float(np.median(np.abs(paired_delta))),
        "wilcoxon_p": float(wilcoxon_result.pvalue),
        "ks_statistic": float(ks_result.statistic),
        "ks_p": float(ks_result.pvalue),
        "median_covered_cpg_ratio": float(np.median(coverage_ratio)),
        "scar_median_covered_cpgs": float(np.median(center_coverage)),
        "adjacent_pooled_median_covered_cpgs": float(
            np.median(adjacent_coverage_pool)
        ),
        "practically_equivalent_within_0.02": bool(ci_low >= -0.02 and ci_high <= 0.02),
    }
    return metrics, center_beta, adjacent_pool, paired_delta


def add_scar_box(ax: plt.Axes, matrix: np.ndarray, axis: str, index: int) -> None:
    n_rows, n_cols = matrix.shape
    if axis == "horizontal":
        patch = Rectangle(
            (-0.5, index - 1.5),
            n_cols,
            1.0,
            fill=False,
            edgecolor="#C43C39",
            linewidth=1.3,
        )
    else:
        patch = Rectangle(
            (index - 1.5, -0.5),
            1.0,
            n_rows,
            fill=False,
            edgecolor="#C43C39",
            linewidth=1.3,
        )
    ax.add_patch(patch)


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ordered = np.sort(values)
    return ordered, np.arange(1, len(ordered) + 1) / len(ordered)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument(
        "--run-root",
        type=Path,
        default=Path("/scr1/users/suc3/20260528_runNewSpatialPipe"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/figure2/scar_methylation"),
    )
    parser.add_argument("--dpi", type=int, default=320)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    spatial_counts = pd.read_csv(
        args.run_root
        / "pipeline_output"
        / SPATIAL_SAMPLE
        / "qc"
        / "table"
        / f"{SPATIAL_SAMPLE}_spatial_barcode_counts.tsv",
        sep="\t",
    )
    spatial_beta = ensure_spatial_summary(
        SPATIAL_SAMPLE,
        args.repo_root.resolve(),
        args.run_root.resolve(),
        args.output_dir.resolve(),
    )
    spatial_count_matrix = matrix_from_frame(spatial_counts, "y", "x", "dna_reads")
    spatial_beta_matrix = matrix_from_frame(
        spatial_beta, "y", "x", "call_weighted_beta"
    )
    spatial_coverage_matrix = matrix_from_frame(
        spatial_beta, "y", "x", "covered_cpgs"
    )

    smcseq = pd.read_csv(
        args.repo_root
        / "results"
        / "figure2"
        / "smcseq_paper_comparison"
        / f"{SMCSEQ_SAMPLE}_matched.tsv",
        sep="\t",
    )
    smcseq = smcseq.loc[smcseq["retained"].astype(bool)].copy()
    smcseq_count_matrix = matrix_from_frame(smcseq, "paper_A", "paper_B", "Raw_reads")
    smcseq_beta_matrix = matrix_from_frame(
        smcseq, "paper_A", "paper_B", "call_weighted_beta"
    )
    smcseq_coverage_matrix = matrix_from_frame(
        smcseq, "paper_A", "paper_B", "covered_cpgs"
    )

    datasets = [
        (
            SPATIAL_SAMPLE,
            "Spatial-DMT",
            spatial_count_matrix,
            spatial_beta_matrix,
            spatial_coverage_matrix,
            "current raw DNA reads",
        ),
        (
            SMCSEQ_SAMPLE,
            "SmC-seq (corrected geometry)",
            smcseq_count_matrix,
            smcseq_beta_matrix,
            smcseq_coverage_matrix,
            "published raw reads at validated barcode joins",
        ),
    ]

    results: list[dict[str, object]] = []
    distributions: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    calls = []
    candidate_rows: list[dict[str, object]] = []
    for sample, method, count_matrix, beta_matrix, coverage_matrix, count_source in datasets:
        call = detect_paired_scar(count_matrix)
        # The legacy whole-line profile emits NumPy warnings for the fully
        # missing structural columns in the corrected SmC-seq grid. Retain its
        # original behavior for comparison while keeping the analysis output
        # clean; the paired detector handles those columns explicitly.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="All-NaN slice encountered",
                category=RuntimeWarning,
            )
            old_call = detect_scar(count_matrix)
        metrics, center_beta, adjacent_pool, paired_delta = compare_line(
            count_matrix,
            beta_matrix,
            coverage_matrix,
            call.axis,
            call.index,
        )
        metrics.update(
            {
                "sample": sample,
                "method": method,
                "count_source": count_source,
                "old_detector_axis": old_call.axis,
                "old_detector_start": old_call.start,
                "old_detector_end": old_call.end,
                "old_detector_ratio": old_call.ratio,
                "blank_rows_excluded": int(
                    np.sum(~np.any(np.isfinite(count_matrix) & (count_matrix > 0), axis=1))
                ),
                "blank_columns_excluded": int(
                    np.sum(~np.any(np.isfinite(count_matrix) & (count_matrix > 0), axis=0))
                ),
            }
        )
        results.append(metrics)
        distributions.append((center_beta, adjacent_pool, paired_delta))
        calls.append(call)
        for rank, candidate in enumerate(score_paired_scars(count_matrix), start=1):
            candidate_rows.append(
                {
                    "sample": sample,
                    "rank": rank,
                    **candidate.__dict__,
                }
            )

    metrics_frame = pd.DataFrame(results)
    metrics_frame.to_csv(
        args.output_dir / "scar_artifact_methylation_metrics.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(candidate_rows).to_csv(
        args.output_dir / "scar_candidates.tsv",
        sep="\t",
        index=False,
    )

    apply_style()
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(11.3, 7.8),
        gridspec_kw={
            "left": 0.07,
            "right": 0.98,
            "bottom": 0.105,
            "top": 0.87,
            "wspace": 0.25,
            "hspace": 0.34,
        },
    )
    for row, (dataset, call, result, distribution) in enumerate(
        zip(datasets, calls, results, distributions)
    ):
        sample, method, count_matrix, beta_matrix, _, _ = dataset
        center_beta, adjacent_pool, paired_delta = distribution
        display_count = count_matrix.copy()
        display_beta = beta_matrix.copy()
        display_index = call.index
        if sample == SMCSEQ_SAMPLE:
            display_count = display_count[::-1, ::-1]
            display_beta = display_beta[::-1, ::-1]
            display_index = (
                display_count.shape[0] - call.index + 1
                if call.axis == "horizontal"
                else display_count.shape[1] - call.index + 1
            )

        positive = display_count[np.isfinite(display_count) & (display_count > 0)]
        count_norm = Normalize(*np.quantile(np.log10(positive + 1), [0.02, 0.98]), clip=True)
        beta_values = display_beta[np.isfinite(display_beta)]
        beta_norm = Normalize(*np.quantile(beta_values, [0.02, 0.98]), clip=True)
        count_ax, beta_ax, distribution_ax = axes[row]
        count_ax.imshow(
            np.ma.masked_invalid(np.log10(display_count + 1)),
            cmap="YlGnBu",
            norm=count_norm,
            interpolation="nearest",
            aspect="equal",
        )
        beta_cmap = plt.get_cmap("magma").copy()
        beta_cmap.set_bad("#E8E8E8")
        beta_ax.imshow(
            np.ma.masked_invalid(display_beta),
            cmap=beta_cmap,
            norm=beta_norm,
            interpolation="nearest",
            aspect="equal",
        )
        add_scar_box(count_ax, display_count, call.axis, display_index)
        add_scar_box(beta_ax, display_beta, call.axis, display_index)
        count_ax.set_title(
            f"{sample}: read count\n{call.axis} {call.index}; center/adjacent={call.median_count_ratio:.3f}",
            fontsize=8.2,
        )
        beta_ax.set_title("Call-weighted CpG methylation", fontsize=8.2)
        for ax in (count_ax, beta_ax):
            ax.set_xticks([])
            ax.set_yticks([])
        x_center, y_center = ecdf(center_beta)
        x_adjacent, y_adjacent = ecdf(adjacent_pool)
        distribution_ax.plot(
            x_center,
            y_center,
            color="#C43C39",
            lw=1.8,
            label="scar line",
        )
        distribution_ax.plot(
            x_adjacent,
            y_adjacent,
            color="#246A9A",
            lw=1.8,
            label="adjacent lines",
        )
        distribution_ax.set_xlabel("Per-pixel CpG methylation (β)")
        distribution_ax.set_ylabel("Cumulative fraction")
        distribution_ax.set_title(
            "Methylation distributions\n"
            f"median paired Δ={result['median_paired_beta_delta']:+.3f} "
            f"[{result['median_paired_beta_delta_ci_low']:+.3f}, "
            f"{result['median_paired_beta_delta_ci_high']:+.3f}]",
            fontsize=8.2,
        )
        distribution_ax.legend(frameon=False, fontsize=6.5, loc="lower right")
        distribution_ax.grid(axis="both", color="#E5E5E5", lw=0.45)
        axes[row, 0].text(
            -0.20,
            0.5,
            method,
            transform=axes[row, 0].transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=8.3,
            fontweight="bold",
        )

    fig.suptitle(
        "Automatic paired-line scar detection and local methylation comparison",
        fontsize=13,
        fontweight="bold",
        y=0.965,
    )
    fig.text(
        0.07,
        0.026,
        (
            "Detector: center pixel / geometric mean of its two immediate neighbors; "
            "≥60% must be <0.55 and the median ratio must be ≤0.55.\n"
            "Blank structural lines are excluded. Brackets: 95% bootstrap CI for the "
            "median paired methylation difference."
        ),
        fontsize=6.3,
        color="#555555",
    )
    output = args.output_dir / "Scar_artifact_methylation_comparison.png"
    fig.savefig(output, dpi=args.dpi, facecolor="white")
    plt.close(fig)
    print(f"Wrote {output}")
    print(metrics_frame.to_string(index=False))


if __name__ == "__main__":
    main()
