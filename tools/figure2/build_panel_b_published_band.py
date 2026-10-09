#!/usr/bin/env python3
"""Build unified-pipeline retention cohorts plus an E6.5 reproduction line.

The plot loads both cohorts from final outputs produced by the shared
Snakefile and writes one canonical PNG for the manuscript figure.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from plot_style import apply_style


LOG = logging.getLogger("panel_b_unified_retention")
RUN_ROOT = Path("/scr1/users/suc3/20260528_runNewSpatialPipe")
SMCSEQ_TRIM_RUN_ROOT = Path("/scr1/users/suc3/20260820_smcseq_retention_speedopt")

METHOD_STYLES = {
    "Spatial-DMT": {"color": "#0072B2", "marker": "o", "linestyle": "-"},
    "SmC-seq": {"color": "#D55E00", "marker": "s", "linestyle": "--"},
    "E6.5 reproduced": {"color": "#666666", "marker": "D", "linestyle": "-."},
}

COHORT_SAMPLES = {
    "Spatial-DMT": (
        "SpMETSLE14DM",
        "SpMETSLE17DM",
        "SpMETSLE21DM",
        "SpMETSLE23DM",
        "SpMETSLB14DM",
    ),
    "SmC-seq": (
        "E55_10um_rep1",
        "E55_10um_rep2",
        "E65_10um_run1",
        "E65_10um_run2",
        "E75_10um_rep1",
        "E75_10um_rep2",
        "E85_20um",
        "E85_placenta_10um",
    ),
}

RETENTION_METRICS = (
    "raw_fraction",
    "trim_demux_fraction",
    "mapped_fraction",
    "qc_retained_mapped_fraction",
)
STAGE_LABELS = (
    "Raw DNA\nreads",
    "Trim passed\nreads",
    "Mapped\nreads",
    "Quality filtered\nreads",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--run-root", type=Path, default=RUN_ROOT)
    parser.add_argument(
        "--smcseq-trim-run-root",
        type=Path,
        default=SMCSEQ_TRIM_RUN_ROOT,
        help="Root containing the optimized SmC-seq trim-only rerun",
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--output-name",
        default="Figure2_experimental_panel_b_spmet_band",
        help="Figure stem and matching metrics/audit prefix",
    )
    parser.add_argument("--dpi", type=int, default=320)
    return parser.parse_args()


def read_curated(path: Path, sample: str) -> dict[str, float]:
    frame = pd.read_csv(path, sep="\t", comment="#")
    row = frame.loc[frame["Sample"].eq(sample)]
    if len(row) != 1:
        raise ValueError(f"Expected one curated row for {sample} in {path}, found {len(row)}")
    result = {}
    for key, value in row.iloc[0].items():
        if key != "Sample":
            result[key] = float(pd.to_numeric(value, errors="coerce"))
    return result


def read_trim_stats(path: Path) -> dict[str, int]:
    frame = pd.read_csv(path, sep="\t", header=None, names=["metric", "value"])
    values = pd.to_numeric(frame["value"], errors="raise")
    return dict(zip(frame["metric"], values.astype(np.int64), strict=True))


def read_raw_read_pairs(path: Path, sample: str) -> int:
    frame = pd.read_csv(path, sep="\t")
    row = frame.loc[frame["Sample"].eq(sample)]
    if len(row) != 1:
        raise ValueError(f"Expected one general-statistics row for {sample} in {path}")
    value = pd.to_numeric(
        row.iloc[0]["custom_content_stats-total_read_pairs"], errors="coerce"
    )
    if not np.isfinite(value):
        raise ValueError(f"Missing raw read-pair count for {sample} in {path}")
    return int(round(float(value)))


def metric_rows(
    method: str,
    sample: str,
    protocol: str,
    values: dict[str, float],
    definitions: dict[str, str],
    sources: dict[str, str],
    source_types: dict[str, str],
) -> list[dict[str, object]]:
    return [
        {
            "method": method,
            "sample": sample,
            "internal_sample": sample,
            "protocol": protocol,
            "metric": metric,
            "value": value,
            "unit": "fraction of raw DNA reads",
            "definition": definitions[metric],
            "source_type": source_types[metric],
            "source_file": sources[metric],
            "local_sources": sources[metric],
            "published_source": "",
        }
        for metric, value in values.items()
    ]


def load_unified_cohorts(
    repo_root: Path,
    run_root: Path = RUN_ROOT,
    smcseq_trim_run_root: Path = SMCSEQ_TRIM_RUN_ROOT,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    snakefile = repo_root / "main" / "Snakefile"

    for method, samples in COHORT_SAMPLES.items():
        protocol = "spatialdmt" if method == "Spatial-DMT" else "smcseq"
        for sample in samples:
            final_curated_source = (
                run_root
                / "final_output"
                / sample
                / "multiqc_report_data"
                / "multiqc_curated_metrics.txt"
            )
            pipeline_curated_source = (
                run_root
                / "pipeline_output"
                / sample
                / "qc"
                / "biscuit_qc_merged"
                / f"{sample}_curated_metrics_mqc.tsv"
            )
            if final_curated_source.is_file():
                curated_source = final_curated_source
                curated_source_type = "unified_pipeline_final"
                curated_source_stage = "final_report"
            elif pipeline_curated_source.is_file():
                curated_source = pipeline_curated_source
                curated_source_type = "unified_pipeline_qc"
                curated_source_stage = "pipeline_qc_before_final_report"
            else:
                raise FileNotFoundError(
                    "Missing unified-pipeline curated metrics; checked "
                    f"{final_curated_source} and {pipeline_curated_source}"
                )
            data = read_curated(curated_source, sample)
            pre_filter_records = data.get("num_reads_merged", np.nan) * 1e6
            mapped_records = data.get("num_mapped_merged", np.nan) * 1e6
            qc_primary_records = data.get("num_reads_qc_bam", np.nan) * 1e6
            if not np.isfinite(
                [pre_filter_records, mapped_records, qc_primary_records]
            ).all():
                raise ValueError(f"Incomplete unified retention metrics for {sample}")

            if protocol == "spatialdmt":
                trim_source = (
                    run_root
                    / "pipeline_output"
                    / sample
                    / "trim"
                    / f"{sample}_trim_stats.txt"
                )
                if not trim_source.is_file():
                    raise FileNotFoundError(f"Missing Spatial-DMT trim statistics: {trim_source}")
                trim_stats = read_trim_stats(trim_source)
                raw_read_pairs = float(trim_stats["Reads_total"])
                raw_reads = 2.0 * raw_read_pairs
                structure_fail_pairs = float(trim_stats["Unmatched_SF"])
                passed_pairs = float(trim_stats["Reads_passed"])
                whitelist_miss_pairs = float(trim_stats.get("Unmatched_WM", 0))
                post_parse_filter_pairs = (
                    raw_read_pairs
                    - passed_pairs
                    - structure_fail_pairs
                    - whitelist_miss_pairs
                )
                if post_parse_filter_pairs < -1:
                    raise ValueError(
                        f"Spatial-DMT trim counters exceed Reads_total for {sample}"
                    )
                if "Reads_accounted" in trim_stats and not np.isclose(
                    float(trim_stats["Reads_accounted"]), raw_read_pairs, rtol=0, atol=0
                ):
                    raise ValueError(
                        f"Spatial-DMT Reads_accounted does not equal Reads_total for {sample}"
                    )
                linker_positive_pairs = raw_read_pairs - structure_fail_pairs
                if linker_positive_pairs < 0:
                    raise ValueError(
                        f"Spatial-DMT linker-positive count is negative for {sample}"
                    )
                if passed_pairs > linker_positive_pairs + 1:
                    raise ValueError(
                        f"Spatial-DMT passed pairs exceed linker-positive pairs for {sample}"
                    )
                trim_demux_reads = 2.0 * linker_positive_pairs
                passed_records = 2.0 * passed_pairs
                if not np.isclose(pre_filter_records, passed_records, rtol=0, atol=1_000):
                    raise ValueError(
                        f"Spatial-DMT trim/BAM read identity failed for {sample}"
                    )
                if not np.isclose(qc_primary_records, pre_filter_records, rtol=0, atol=1):
                    raise ValueError(
                        f"Spatial-DMT QC BAM should equal its pre-filter BAM for {sample}"
                    )
                mapped_reads = mapped_records
                qc_mapped_reads = mapped_reads
                final_definition = (
                    "pre_filter_primary_mapped_reads_divided_by_twice_raw_read_pairs; "
                    "Spatial-DMT QC BAM equals pre-filter BAM because no bsconv rule runs"
                )
                final_source = f"{curated_source};{snakefile}"
                trim_definition = (
                    "twice_(Reads_total_minus_Unmatched_SF)_divided_by_twice_Reads_total; "
                    "linker-positive fraction uses the raw Reads_total denominator"
                )
                trim_source_file = str(trim_source)
                trim_source_type = "unified_pipeline_trim_stats"
            else:
                general_source = (
                    run_root
                    / "final_output"
                    / sample
                    / "multiqc_report_data"
                    / "multiqc_general_stats.txt"
                )
                if not general_source.is_file():
                    raise FileNotFoundError(
                        f"Missing SmC-seq raw read-pair summary: {general_source}"
                    )
                previous_raw_read_pairs = read_raw_read_pairs(general_source, sample)
                trim_source = (
                    smcseq_trim_run_root
                    / "pipeline_output"
                    / sample
                    / "trim"
                    / f"{sample}_trim_stats.txt"
                )
                if not trim_source.is_file():
                    raise FileNotFoundError(
                        f"Missing optimized SmC-seq trim statistics: {trim_source}"
                    )
                trim_stats = read_trim_stats(trim_source)
                raw_read_pairs = float(trim_stats["Reads_total"])
                if raw_read_pairs != previous_raw_read_pairs:
                    raise ValueError(
                        f"SmC-seq raw denominator changed for {sample}: "
                        f"{int(raw_read_pairs)} != {previous_raw_read_pairs}"
                    )
                raw_reads = 2.0 * raw_read_pairs
                passed_pairs = int(trim_stats["Reads_passed"])
                structure_fail_pairs = int(trim_stats.get("Unmatched_SF", 0))
                whitelist_miss_pairs = int(trim_stats.get("Unmatched_WM", 0))
                post_parse_filter_pairs = (
                    raw_read_pairs
                    - passed_pairs
                    - structure_fail_pairs
                    - whitelist_miss_pairs
                )
                if post_parse_filter_pairs < -1:
                    raise ValueError(
                        f"SmC-seq trim counters exceed Reads_total for {sample}"
                    )
                if "Reads_accounted" in trim_stats and not np.isclose(
                    float(trim_stats["Reads_accounted"]), raw_read_pairs, rtol=0, atol=0
                ):
                    raise ValueError(
                        f"SmC-seq Reads_accounted does not equal Reads_total for {sample}"
                    )
                strand_pairs = int(trim_stats["Reads_Watson"]) + int(
                    trim_stats["Reads_Crick"]
                )
                if passed_pairs != strand_pairs:
                    raise ValueError(
                        f"SmC-seq Watson/Crick trim accounting failed for {sample}"
                    )
                # Each passing SmC-seq input pair emits one genomic R1
                # (Watson or Crick). The plotted denominator is 2 x raw
                # read pairs, so this checkpoint carries the single-read
                # penalty just like the mapped and quality-filtered stages.
                trim_demux_reads = float(passed_pairs)
                if not np.isclose(
                    pre_filter_records,
                    passed_pairs,
                    rtol=0,
                    atol=1_000,
                ):
                    raise ValueError(
                        f"SmC-seq trim/pre-filter read identity failed for {sample}"
                    )
                mapped_reads = mapped_records
                if qc_primary_records > mapped_records:
                    raise ValueError(
                        f"SmC-seq post-bsconv primary count exceeds pre-filter mapped count for {sample}"
                    )
                # biscuit bsconv removes unmapped and high-CpH-retention records,
                # so every primary record in this QC BAM is retained mapped input.
                qc_mapped_reads = qc_primary_records
                final_definition = (
                    "post_bsconv_qc_BAM_primary_reads_divided_by_twice_raw_read_pairs; "
                    "SmC-seq aligns genomic R1 only and bsconv removes unmapped and "
                    "high_CpH_retention reads"
                )
                final_source = str(curated_source)
                trim_definition = (
                    "Reads_passed_single_genomic_R1_records_divided_by_twice_Reads_total"
                )
                trim_source_file = str(trim_source)
                trim_source_type = "unified_pipeline_trim_stats"

            values = {
                "raw_fraction": 1.0,
                "trim_demux_fraction": trim_demux_reads / raw_reads,
                "mapped_fraction": mapped_reads / raw_reads,
                "qc_retained_mapped_fraction": qc_mapped_reads / raw_reads,
            }
            finite_values = np.asarray(list(values.values()), dtype=float)
            finite_values = finite_values[np.isfinite(finite_values)]
            if not (
                np.all((0 <= finite_values) & (finite_values <= 1))
                and np.all(np.diff(finite_values) <= 1e-12)
            ):
                raise ValueError(f"Non-monotonic unified retention trajectory for {sample}")

            definitions = {
                "raw_fraction": "twice_raw_DNA_read_pairs_normalized_to_one",
                "trim_demux_fraction": trim_definition,
                "mapped_fraction": (
                    "pre_filter_primary_mapped_reads_divided_by_twice_raw_read_pairs"
                ),
                "qc_retained_mapped_fraction": final_definition,
            }
            raw_source = trim_source_file
            sources = {
                "raw_fraction": raw_source,
                "trim_demux_fraction": trim_source_file,
                "mapped_fraction": str(curated_source),
                "qc_retained_mapped_fraction": final_source,
            }
            source_types = {
                "raw_fraction": "unified_pipeline_trim_stats",
                "trim_demux_fraction": trim_source_type,
                "mapped_fraction": curated_source_type,
                "qc_retained_mapped_fraction": curated_source_type,
            }
            sample_rows = metric_rows(
                method,
                sample,
                protocol,
                values,
                definitions,
                sources,
                source_types,
            )
            rows.extend(sample_rows)
            audit_rows.append(
                {
                    "sample": sample,
                    "method": method,
                    "protocol": protocol,
                    "status": "included",
                    "raw_read_pairs": int(round(raw_read_pairs)),
                    "raw_input_reads": int(round(raw_reads)),
                    "trim_demux_passed_reads": (
                        int(round(trim_demux_reads))
                        if np.isfinite(trim_demux_reads)
                        else np.nan
                    ),
                    # Keep the literal Reads_passed checkpoint beside the
                    # plotted checkpoint. Spatial-DMT's plotted value is
                    # linker-positive by design and can exceed Reads_passed
                    # because later barcode/length exits are retained in the
                    # linker-positive numerator.
                    "reads_passed_pairs": int(round(passed_pairs)),
                    "reads_passed_pair_fraction": passed_pairs / raw_read_pairs,
                    "structure_fail_pairs": int(round(structure_fail_pairs)),
                    "whitelist_miss_pairs": int(round(whitelist_miss_pairs)),
                    "post_parse_filter_pairs": int(round(post_parse_filter_pairs)),
                    "pre_filter_primary_records": int(round(pre_filter_records)),
                    "mapped_reads": int(round(mapped_reads)),
                    "qc_retained_mapped_reads": int(round(qc_mapped_reads)),
                    "trim_checkpoint_status": (
                        "measured" if np.isfinite(trim_demux_reads) else "unavailable"
                    ),
                    "read_unit_check": (
                        "Spatial-DMT: paired BAM primary records = 2 x passed pairs"
                        if protocol == "spatialdmt"
                        else "SmC-seq: one genomic R1 record = one passed pair"
                    ),
                    "checkpoint_unit": (
                        "paired BAM records over 2 x raw read pairs"
                        if protocol == "spatialdmt"
                        else "single genomic R1 records over 2 x raw read pairs"
                    ),
                    "source_file": str(curated_source),
                    "curated_source_stage": curated_source_stage,
                    "trim_source_file": trim_source_file,
                    "raw_validation_source_file": (
                        str(general_source) if protocol == "smcseq" else ""
                    ),
                }
            )

    metrics = pd.DataFrame(rows)
    audit = pd.DataFrame(audit_rows)
    return metrics, audit


def comparison_value(frame: pd.DataFrame, pipeline: str, metric: str) -> float:
    selected = pd.to_numeric(
        frame.loc[(frame["pipeline"] == pipeline) & (frame["metric"] == metric), "value"],
        errors="coerce",
    ).dropna()
    if len(selected) != 1:
        raise ValueError(f"Expected one finite {pipeline} value for {metric!r}")
    return float(selected.iloc[0])


def load_e65_reproduced(repo_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = repo_root / "tools" / "figure2" / "pipeline_comparison.tsv"
    comparison = pd.read_csv(source, sep="\t")
    raw = comparison_value(comparison, "e65_reproduced", "Raw read pairs")
    designed = comparison_value(comparison, "e65_reproduced", "Designed")
    mapped = comparison_value(
        comparison, "e65_reproduced", "Mapped before CpH filter"
    )
    qc_retained = comparison_value(
        comparison, "e65_reproduced", "CpH-filter passed"
    )
    # The reproduced run demultiplexes the spatial barcode from R2 and carries
    # one selected mate forward for the downstream single-end path.  The source
    # table's structure-pass, mapped, and XM-filter-passed counts are therefore
    # single-read records, while Raw read pairs counts both mates.  Normalize
    # this line to 2 x raw pairs just like the live SmC-seq branch.
    raw_reads = 2.0 * raw
    designed_reads = designed
    if not 0 <= qc_retained <= mapped <= designed <= raw:
        raise ValueError("Non-monotonic E6.5 reproduced retention trajectory")

    values = {
        "raw_fraction": 1.0,
        "trim_demux_fraction": designed_reads / raw_reads,
        "mapped_fraction": mapped / raw_reads,
        "qc_retained_mapped_fraction": qc_retained / raw_reads,
    }
    definitions = {
        "raw_fraction": "reproduced_twice_raw_read_pairs_normalized_to_one",
        "trim_demux_fraction": (
            "reproduced_single_mate_structure_pass_records_divided_by_twice_raw_read_pairs"
        ),
        "mapped_fraction": (
            "reproduced_single_mate_mapped_records_divided_by_twice_raw_read_pairs"
        ),
        "qc_retained_mapped_fraction": (
            "reproduced_single_mate_CpH_filter_passed_records_divided_by_twice_raw_read_pairs"
        ),
    }
    sources = {metric: str(source) for metric in values}
    source_types = {
        metric: "reproduced_pipeline_curated_summary" for metric in values
    }
    metrics = pd.DataFrame(
        metric_rows(
            "E6.5 reproduced",
            "E6.5 reproduced",
            "smcseq_reproduced",
            values,
            definitions,
            sources,
            source_types,
        )
    )
    audit = pd.DataFrame(
        [
            {
                "sample": "E6.5 reproduced",
                "method": "E6.5 reproduced",
                "protocol": "smcseq_reproduced",
                "status": "included_single_experiment",
                "raw_read_pairs": int(round(raw)),
                "raw_input_reads": int(round(raw_reads)),
                "trim_demux_passed_reads": int(round(designed_reads)),
                "pre_filter_primary_records": int(round(designed)),
                "mapped_reads": int(round(mapped)),
                "qc_retained_mapped_reads": int(round(qc_retained)),
                "trim_checkpoint_status": "measured",
                "read_unit_check": (
                    "E6.5 reproduced: one selected-mate record per input pair"
                ),
                "checkpoint_unit": "single-mate records over 2 x raw read pairs",
                "source_file": str(source),
            }
        ]
    )
    return metrics, audit


def load_panel_b_metrics(
    repo_root: Path,
    run_root: Path = RUN_ROOT,
    smcseq_trim_run_root: Path = SMCSEQ_TRIM_RUN_ROOT,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cohorts, cohort_audit = load_unified_cohorts(
        repo_root,
        run_root,
        smcseq_trim_run_root,
    )
    reproduced, reproduced_audit = load_e65_reproduced(repo_root)
    metrics = pd.concat([cohorts, reproduced], ignore_index=True)
    audit = pd.concat([cohort_audit, reproduced_audit], ignore_index=True)
    return metrics, audit


def bootstrap_median_ci(
    values: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, np.nan
    draws = rng.choice(values, size=(20_000, len(values)), replace=True)
    medians = np.median(draws, axis=1)
    low, high = np.quantile(medians, [0.025, 0.975])
    return float(np.median(values)), float(low), float(high)


def retention_matrix(metrics: pd.DataFrame, method: str) -> pd.DataFrame:
    subset = metrics.loc[
        (metrics["method"] == method) & metrics["metric"].isin(RETENTION_METRICS),
        ["sample", "metric", "value"],
    ]
    matrix = subset.pivot(index="sample", columns="metric", values="value")
    return matrix.reindex(columns=RETENTION_METRICS)


def plot_retention(ax: plt.Axes, metrics: pd.DataFrame, fontsize: float = 8.0) -> None:
    x = np.arange(len(RETENTION_METRICS), dtype=float)
    rng = np.random.default_rng(20260811)
    legend = []

    for method, style in METHOD_STYLES.items():
        matrix = retention_matrix(metrics, method)
        values = matrix.to_numpy(dtype=float)
        if not len(matrix):
            raise ValueError(f"Missing plotted trajectory for {method}")
        missing = ~np.isfinite(values)
        expected_missing = np.zeros_like(missing, dtype=bool)
        if not np.array_equal(missing, expected_missing):
            raise ValueError(f"Unexpected missing checkpoints for {method}")

        if method == "E6.5 reproduced":
            finite = np.isfinite(values[0])
            ax.plot(
                x[finite],
                values[0, finite],
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=2.2,
                marker=style["marker"],
                markersize=5.8,
                markerfacecolor="white",
                markeredgewidth=1.25,
                zorder=5,
            )
            label = "E6.5 reproduced"
        else:
            for row in values:
                finite = np.isfinite(row)
                ax.plot(
                    x[finite],
                    row[finite],
                    color=style["color"],
                    linestyle=style["linestyle"],
                    linewidth=0.75,
                    marker=style["marker"],
                    markersize=3.0,
                    markerfacecolor="white",
                    markeredgewidth=0.65,
                    alpha=0.18,
                    zorder=2,
                )
            summaries = np.asarray(
                [
                    bootstrap_median_ci(values[:, index], rng)
                    for index in range(values.shape[1])
                ]
            )
            medians, lows, highs = summaries.T
            finite = np.isfinite(medians)
            ax.fill_between(
                x[finite],
                lows[finite],
                highs[finite],
                color=style["color"],
                alpha=0.17,
                linewidth=0,
                zorder=1,
            )
            ax.plot(
                x[finite],
                medians[finite],
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=2.5,
                marker=style["marker"],
                markersize=6.2,
                markerfacecolor="white",
                markeredgewidth=1.4,
                zorder=4,
            )
            label = f"{method}, unified workflow"

        legend.append(
            Line2D(
                [],
                [],
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=2.0,
                marker=style["marker"],
                markerfacecolor="white",
                markeredgewidth=1.1,
                markersize=5.2,
                label=label,
            )
        )

    ax.set_xlim(-0.15, 3.15)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(x, STAGE_LABELS)
    ax.set_ylabel("Fraction of raw DNA reads", fontsize=fontsize)
    ax.tick_params(axis="both", labelsize=fontsize - 0.4)
    ax.grid(axis="y", color="#E4E4E4", linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        handles=legend,
        loc="upper right",
        frameon=False,
        fontsize=fontsize - 0.5,
        handlelength=2.8,
    )


def make_figure(metrics: pd.DataFrame) -> plt.Figure:
    apply_style()
    fig, ax = plt.subplots(figsize=(10.8, 6.4))
    fig.subplots_adjust(left=0.105, right=0.965, bottom=0.18, top=0.84)
    fig.suptitle(
        "Retention from raw reads through alignment and QC checkpoints",
        x=0.105,
        y=0.94,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    plot_retention(ax, metrics)
    return fig


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    repo_root = args.repo_root.resolve()
    run_root = args.run_root.resolve()
    smcseq_trim_run_root = args.smcseq_trim_run_root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else repo_root / "results" / "figure2"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics, audit = load_panel_b_metrics(
        repo_root,
        run_root,
        smcseq_trim_run_root,
    )
    expected_methods = {"Spatial-DMT", "SmC-seq", "E6.5 reproduced"}
    if set(metrics["method"]) != expected_methods:
        raise ValueError(f"Unexpected plotted methods: {sorted(metrics['method'].unique())}")
    if metrics.duplicated(["method", "sample", "metric"]).any():
        raise ValueError("Duplicate method/sample/checkpoint rows")
    missing = metrics["value"].isna()
    if missing.any():
        raise ValueError("Panel B contains missing retention values")

    metrics_path = output_dir / f"{args.output_name}_metrics.tsv"
    audit_path = output_dir / f"{args.output_name}_audit.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False, na_rep="NA")
    audit.to_csv(audit_path, sep="\t", index=False, na_rep="NA")

    fig = make_figure(metrics)
    stem = output_dir / args.output_name
    fig.savefig(
        stem.with_suffix(".png"),
        dpi=args.dpi,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)

    print("Spatial-DMT: unified-workflow cohort loaded")
    print("SmC-seq: unified-workflow cohort loaded")
    print("E6.5 reproduced: one trajectory, no confidence interval")
    print(f"Wrote: {stem}.png")
    print(f"Metrics: {metrics_path}")
    print(f"Audit: {audit_path}")


if __name__ == "__main__":
    main()
