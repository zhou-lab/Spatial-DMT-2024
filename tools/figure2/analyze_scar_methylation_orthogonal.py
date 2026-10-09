#!/usr/bin/env python3
"""Test column scars in Spatial-DMT and row scars in corrected SmC-seq maps."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from analyze_scar_methylation import (
    add_scar_box,
    compare_line,
    ecdf,
    ensure_spatial_summary,
    matrix_from_frame,
)
from plot_style import apply_style
from spatial_utils import PairedScarCall, score_paired_scars


SPATIAL_SAMPLES = (
    "SpMETSLE14DM",
    "SpMETSLE17DM",
    "SpMETSLE21DM",
)
SMCSEQ_SAMPLES = (
    "E55_10um_rep1",
    "E55_10um_rep2",
    "E85_20um",
)


def requested_call(matrix: np.ndarray, axis: str) -> tuple[PairedScarCall, list[PairedScarCall]]:
    candidates = [call for call in score_paired_scars(matrix) if call.axis == axis]
    if not candidates:
        raise ValueError(f"No supported {axis} scar candidate was detected")
    return candidates[0], candidates


def blank_lines(matrix: np.ndarray) -> tuple[int, int]:
    populated = np.isfinite(matrix) & (matrix > 0)
    return (
        int(np.sum(~np.any(populated, axis=1))),
        int(np.sum(~np.any(populated, axis=0))),
    )


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
        "--comparison-dir",
        type=Path,
        default=Path("results/figure2/smcseq_paper_comparison"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/figure2/scar_methylation_orthogonal"),
    )
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    datasets: list[dict[str, object]] = []
    for sample in SPATIAL_SAMPLES:
        count_frame = pd.read_csv(
            args.run_root
            / "pipeline_output"
            / sample
            / "qc"
            / "table"
            / f"{sample}_spatial_barcode_counts.tsv",
            sep="\t",
        )
        methylation_frame = ensure_spatial_summary(
            sample,
            args.repo_root.resolve(),
            args.run_root.resolve(),
            args.output_dir.resolve(),
        )
        datasets.append(
            {
                "sample": sample,
                "method": "Spatial-DMT",
                "requested_axis": "vertical",
                "count": matrix_from_frame(count_frame, "y", "x", "dna_reads"),
                "beta": matrix_from_frame(
                    methylation_frame, "y", "x", "call_weighted_beta"
                ),
                "coverage": matrix_from_frame(
                    methylation_frame, "y", "x", "covered_cpgs"
                ),
                "count_source": "current raw DNA reads",
                "reverse_display": False,
            }
        )

    for sample in SMCSEQ_SAMPLES:
        frame = pd.read_csv(args.comparison_dir / f"{sample}_matched.tsv", sep="\t")
        frame = frame.loc[frame["retained"].astype(bool)].copy()
        datasets.append(
            {
                "sample": sample,
                "method": "SmC-seq (corrected geometry)",
                "requested_axis": "horizontal",
                "count": matrix_from_frame(frame, "paper_A", "paper_B", "Raw_reads"),
                "beta": matrix_from_frame(
                    frame, "paper_A", "paper_B", "call_weighted_beta"
                ),
                "coverage": matrix_from_frame(
                    frame, "paper_A", "paper_B", "covered_cpgs"
                ),
                "count_source": "published raw reads at validated barcode joins",
                "reverse_display": True,
            }
        )

    results: list[dict[str, object]] = []
    distributions: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    calls: list[PairedScarCall] = []
    candidate_rows: list[dict[str, object]] = []
    for dataset in datasets:
        sample = str(dataset["sample"])
        axis = str(dataset["requested_axis"])
        count = np.asarray(dataset["count"])
        beta = np.asarray(dataset["beta"])
        coverage = np.asarray(dataset["coverage"])
        call, candidates = requested_call(count, axis)
        metrics, center_beta, adjacent_pool, paired_delta = compare_line(
            count, beta, coverage, call.axis, call.index
        )
        blank_rows, blank_columns = blank_lines(count)
        metrics.update(
            {
                "sample": sample,
                "method": dataset["method"],
                "requested_axis": axis,
                "count_source": dataset["count_source"],
                "blank_rows_excluded": blank_rows,
                "blank_columns_excluded": blank_columns,
                "supported_candidates_on_requested_axis": len(candidates),
            }
        )
        results.append(metrics)
        distributions.append((center_beta, adjacent_pool, paired_delta))
        calls.append(call)
        for rank, candidate in enumerate(candidates, start=1):
            candidate_rows.append(
                {"sample": sample, "requested_axis": axis, "rank": rank, **candidate.__dict__}
            )

    metrics_frame = pd.DataFrame(results)
    metrics_frame.to_csv(
        args.output_dir / "orthogonal_scar_methylation_metrics.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(candidate_rows).to_csv(
        args.output_dir / "orthogonal_scar_candidates.tsv",
        sep="\t",
        index=False,
    )

    apply_style()
    n_rows = len(datasets)
    fig, axes = plt.subplots(
        n_rows,
        3,
        figsize=(11.2, 3.25 * n_rows),
        gridspec_kw={
            "left": 0.08,
            "right": 0.98,
            "bottom": 0.07,
            "top": 0.93,
            "wspace": 0.24,
            "hspace": 0.42,
        },
    )
    count_cmap = plt.get_cmap("YlGnBu").copy()
    beta_cmap = plt.get_cmap("magma").copy()
    count_cmap.set_bad("#E8E8E8")
    beta_cmap.set_bad("#E8E8E8")

    for row, (dataset, call, result, distribution) in enumerate(
        zip(datasets, calls, results, distributions)
    ):
        sample = str(dataset["sample"])
        count = np.asarray(dataset["count"]).copy()
        beta = np.asarray(dataset["beta"]).copy()
        display_index = call.index
        if bool(dataset["reverse_display"]):
            count = count[::-1, ::-1]
            beta = beta[::-1, ::-1]
            display_index = (
                count.shape[0] - call.index + 1
                if call.axis == "horizontal"
                else count.shape[1] - call.index + 1
            )

        positive = count[np.isfinite(count) & (count > 0)]
        count_norm = Normalize(
            *np.quantile(np.log10(positive + 1), [0.02, 0.98]), clip=True
        )
        beta_values = beta[np.isfinite(beta)]
        beta_norm = Normalize(*np.quantile(beta_values, [0.02, 0.98]), clip=True)
        count_ax, beta_ax, distribution_ax = axes[row]
        count_ax.imshow(
            np.ma.masked_invalid(np.log10(count + 1)),
            cmap=count_cmap,
            norm=count_norm,
            interpolation="nearest",
            aspect="equal",
        )
        beta_ax.imshow(
            np.ma.masked_invalid(beta),
            cmap=beta_cmap,
            norm=beta_norm,
            interpolation="nearest",
            aspect="equal",
        )
        add_scar_box(count_ax, count, call.axis, display_index)
        add_scar_box(beta_ax, beta, call.axis, display_index)
        line_name = "row" if call.axis == "horizontal" else "column"
        count_ax.set_title(
            f"{sample}: {line_name} {call.index}\n"
            f"center/adjacent={call.median_count_ratio:.3f}; n={call.support_pixels}",
            fontsize=8.0,
        )
        beta_ax.set_title("Call-weighted CpG methylation", fontsize=8.0)
        for ax in (count_ax, beta_ax):
            ax.set_xticks([])
            ax.set_yticks([])

        center_beta, adjacent_pool, _ = distribution
        x_center, y_center = ecdf(center_beta)
        x_adjacent, y_adjacent = ecdf(adjacent_pool)
        distribution_ax.plot(
            x_center, y_center, color="#C43C39", lw=1.7, label="scar line"
        )
        distribution_ax.plot(
            x_adjacent,
            y_adjacent,
            color="#246A9A",
            lw=1.7,
            label="adjacent lines",
        )
        distribution_ax.set_xlabel("Per-pixel CpG methylation (β)")
        distribution_ax.set_ylabel("Cumulative fraction")
        distribution_ax.set_title(
            "Methylation distributions\n"
            f"median paired Δ={result['median_paired_beta_delta']:+.3f} "
            f"[{result['median_paired_beta_delta_ci_low']:+.3f}, "
            f"{result['median_paired_beta_delta_ci_high']:+.3f}]",
            fontsize=8.0,
        )
        distribution_ax.legend(frameon=False, fontsize=6.4, loc="lower right")
        distribution_ax.grid(axis="both", color="#E5E5E5", lw=0.45)
        axes[row, 0].text(
            -0.21,
            0.5,
            str(dataset["method"]),
            transform=axes[row, 0].transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=7.8,
            fontweight="bold",
        )

    fig.suptitle(
        "Direction-swapped scar validation across independent samples",
        fontsize=13,
        fontweight="bold",
        y=0.975,
    )
    fig.text(
        0.08,
        0.025,
        (
            "Spatial-DMT: strongest supported column scar per sample. SmC-seq: strongest "
            "supported row scar per corrected-geometry sample.\n"
            "Brackets show the 95% bootstrap CI for the median matched-pixel methylation difference."
        ),
        fontsize=6.4,
        color="#555555",
    )
    output = args.output_dir / "Orthogonal_scar_methylation_comparison.png"
    fig.savefig(output, dpi=args.dpi, facecolor="white")
    plt.close(fig)
    print(f"Wrote {output}")
    print(metrics_frame.to_string(index=False))


if __name__ == "__main__":
    main()
