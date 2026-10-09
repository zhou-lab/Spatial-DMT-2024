#!/usr/bin/env python3
"""Build manuscript Figures 2-5 from existing pipeline and audited source data."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter


HERE = Path(__file__).resolve().parent
FIGURE2_TOOLS = HERE.parent / "figure2"
if str(FIGURE2_TOOLS) not in sys.path:
    sys.path.insert(0, str(FIGURE2_TOOLS))

from build_figure2 import (  # noqa: E402
    MATRIX_METRICS,
    _comparison_display,
    _comparison_value,
    _positive_lognorm,
    _spatial_extent,
    compact_tick,
    context_color,
    load_panel_b_cohorts,
    load_pipeline_comparison,
    metric_source,
    metric_value,
    plot_panel_a_multiqc,
    plot_panel_b_cohort_retention,
)
from load_metrics import DEFAULT_SAMPLES, SampleRecord, read_legacy_whitelist  # noqa: E402
from plot_style import (  # noqa: E402
    CONTEXT_COLORS,
    SPATIAL_CMAP,
    WARNING_COLOR,
    apply_style,
    clean_axis,
    panel_title,
)
from spatial_utils import (  # noqa: E402
    detect_scar,
    legacy_order_reconstruction,
    morans_i,
    neighbor_spearman,
    spatial_matrix,
)


LOG = logging.getLogger("manuscript_figures")
RUN_ROOTS = (
    Path("/scr1/users/suc3/20260528_runNewSpatialPipe"),
    Path("/scr1/users/suc3/20260704_runSpMET_latest_mm10"),
    Path("/scr1/users/suc3/20260617_runNewSpatialPipe_hg38"),
    Path("/scr1/users/suc3/20260624_runNewSpatialPipe_hg38"),
    Path("/scr1/users/suc3/20260706_runSpMETSLHSK_more_hg38"),
)
EXCLUDED_RUNS = {
    "SpMETSLE14D_smoke",
    "SpMETSLE16D_ourpipeline_bsconv_y005",
    "SpMETSLHSK11D1",
    "SpMETSLHSK11D2",
}
FAILURE_COLORS = {
    "Structure failure": "#0072B2",
    "Whitelist miss": "#E69F00",
    "Primer": "#CC79A7",
    "Linker 1": "#56B4E9",
    "Linker 2": "#009E73",
    "Tn5 anchor": "#D55E00",
    "Other": "#999999",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path("/home/suc3/Lab/Spatial-DMT-2024"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "/home/suc3/Lab/Spatial-DMT-2024/results/manuscript_figures/data_driven"
        ),
    )
    parser.add_argument("--dpi", type=int, default=320)
    parser.add_argument("--tiff-dpi", type=int, default=600)
    parser.add_argument(
        "--figures",
        nargs="+",
        choices=("2", "3", "4", "5"),
        default=("2", "3", "4", "5"),
    )
    return parser.parse_args()


def save_figure(fig: plt.Figure, stem: Path, dpi: int, tiff_dpi: int) -> None:
    for suffix in ("pdf", "svg"):
        fig.savefig(stem.with_suffix(f".{suffix}"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=tiff_dpi,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def sample_metadata(sample: str) -> tuple[str, str, str]:
    if sample == "E65_20um":
        return "mouse", "embryo", "SmC-seq"
    if "HSK" in sample:
        return "human", "skin", "Spatial-DMT"
    if "SLB" in sample:
        return "mouse", "brain", "Spatial-DMT"
    return "mouse", "embryo", "Spatial-DMT"


def short_sample(sample: str) -> str:
    return sample.replace("SpMETSL", "").replace("SpHSK", "HSK")


def discover_unmatched_diagnostics() -> pd.DataFrame:
    newest: dict[str, Path] = {}
    for root in RUN_ROOTS:
        for path in root.glob("pipeline_output/*/qc/unmatched_diagnostics.tsv"):
            sample = path.parents[1].name
            if sample in EXCLUDED_RUNS or not path.is_file():
                continue
            previous = newest.get(sample)
            if previous is None or path.stat().st_mtime > previous.stat().st_mtime:
                newest[sample] = path

    rows: list[dict] = []
    for sample, path in newest.items():
        source = pd.read_csv(path, sep="\t")
        source["value"] = pd.to_numeric(source["value"], errors="raise")
        values = source.set_index("key")["value"].to_dict()
        required = {
            "structure_fail_n_reads",
            "wl_miss_n_reads",
            "total_n_reads",
            "structure_fail_mapping_rate",
            "wl_miss_mapping_rate",
            "structure_fail_dup_rate",
            "wl_miss_dup_rate",
            "linker1_fail_frac",
            "linker2_fail_frac",
            "tn5_me_fail_frac",
        }
        missing = required - set(values)
        if missing:
            raise ValueError(f"{path} lacks required fields: {sorted(missing)}")
        species, tissue, protocol = sample_metadata(sample)
        common = {
            "sample": sample,
            "species": species,
            "tissue": tissue,
            "protocol": protocol,
            "source_file": str(path),
            "source_type": "revised",
        }
        total = values["total_n_reads"]
        metrics = {
            "unmatched_structure_fraction": values["structure_fail_n_reads"] / total,
            "unmatched_whitelist_fraction": values["wl_miss_n_reads"] / total,
            "structure_mapping_rate": values["structure_fail_mapping_rate"],
            "whitelist_mapping_rate": values["wl_miss_mapping_rate"],
            "structure_duplication_rate": values["structure_fail_dup_rate"],
            "whitelist_duplication_rate": values["wl_miss_dup_rate"],
            "primer_failure_fraction": values.get("primer_fail_frac", np.nan),
            "linker1_failure_fraction": values["linker1_fail_frac"],
            "linker2_failure_fraction": values["linker2_fail_frac"],
            "tn5_failure_fraction": values["tn5_me_fail_frac"],
            "other_structure_failure_fraction": values.get("all_ok_fail_frac", np.nan),
        }
        for metric, value in metrics.items():
            rows.append(
                {
                    **common,
                    "metric": metric,
                    "value": value,
                    "unit": "fraction",
                    "definition_version": "unmatched_diagnostics_v2",
                    "is_fallback": False,
                    "panel": "Figure2",
                }
            )
    result = pd.DataFrame(rows)
    if result["sample"].nunique() < 10:
        raise ValueError("Too few independent unmatched-diagnostic runs were found")
    return result.sort_values(["protocol", "tissue", "sample", "metric"])


def failure_wide(metrics: pd.DataFrame) -> pd.DataFrame:
    metadata = metrics.drop_duplicates("sample").set_index("sample")
    wide = metrics.pivot(index="sample", columns="metric", values="value")
    for column in ("species", "tissue", "protocol", "source_file"):
        wide[column] = metadata[column]
    order = wide.sort_values(["protocol", "tissue", "sample"]).index
    return wide.loc[order]


def plot_failure_pair(
    ax,
    wide: pd.DataFrame,
    first: str,
    second: str,
    first_label: str,
    second_label: str,
    title: str,
    show_labels: bool,
) -> None:
    y = np.arange(len(wide))
    x1 = 100 * wide[first].to_numpy(dtype=float)
    x2 = 100 * wide[second].to_numpy(dtype=float)
    for yi, a, b in zip(y, x1, x2):
        ax.plot([a, b], [yi, yi], color="#BDBDBD", lw=0.7, zorder=1)
    ax.scatter(x1, y, color=FAILURE_COLORS["Structure failure"], marker="o", s=19, zorder=3, label=first_label)
    ax.scatter(x2, y, facecolor="white", edgecolor=FAILURE_COLORS["Whitelist miss"], marker="s", s=20, linewidth=0.9, zorder=3, label=second_label)
    ax.set_yticks(y)
    ax.set_yticklabels([short_sample(sample) for sample in wide.index] if show_labels else [])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Rate (%)")
    ax.legend(loc="lower right", fontsize=6, handletextpad=0.4)
    clean_axis(ax, "x")


def make_figure2(repo_root: Path, output_dir: Path, dpi: int, tiff_dpi: int) -> None:
    metrics = discover_unmatched_diagnostics()
    wide = failure_wide(metrics)
    fig = plt.figure(figsize=(15, 10.5), constrained_layout=False)
    fig.suptitle(
        "Detailed diagnosis of failed spatial-methylation reads",
        fontsize=14,
        weight="bold",
        y=0.985,
    )
    grid = fig.add_gridspec(2, 2, left=0.08, right=0.98, bottom=0.07, top=0.93, hspace=0.28, wspace=0.2)
    axes = [fig.add_subplot(grid[i, j]) for i in range(2) for j in range(2)]

    ax = axes[0]
    panel_title(ax, "A", "Composition of reads excluded from spatial assignment")
    y = np.arange(len(wide))
    structure = 100 * wide["unmatched_structure_fraction"].to_numpy(dtype=float)
    whitelist = 100 * wide["unmatched_whitelist_fraction"].to_numpy(dtype=float)
    ax.barh(y, structure, color=FAILURE_COLORS["Structure failure"], height=0.68, label="Structure failure")
    ax.barh(y, whitelist, left=structure, color=FAILURE_COLORS["Whitelist miss"], height=0.68, label="Whitelist miss")
    ax.set_yticks(y, [short_sample(sample) for sample in wide.index])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Fraction of unmatched reads (%)")
    ax.legend(loc="lower right", fontsize=6)
    clean_axis(ax, "x")

    ax = axes[1]
    panel_title(ax, "B", "First failed structural anchor")
    bottom = np.zeros(len(wide))
    anchor_specs = (
        ("primer_failure_fraction", "Primer"),
        ("linker1_failure_fraction", "Linker 1"),
        ("linker2_failure_fraction", "Linker 2"),
        ("tn5_failure_fraction", "Tn5 anchor"),
        ("other_structure_failure_fraction", "Other"),
    )
    for metric, label in anchor_specs:
        values = 100 * wide[metric].to_numpy(dtype=float)
        finite = np.isfinite(values) & (values > 0)
        ax.barh(y[finite], values[finite], left=bottom[finite], height=0.68, color=FAILURE_COLORS[label], label=label)
        bottom[finite] += values[finite]
    ax.set_yticks(y, [short_sample(sample) for sample in wide.index])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Fraction of structure-fail reads (%)")
    ax.legend(loc="lower right", fontsize=5.8, ncol=2)
    clean_axis(ax, "x")
    ax.text(0.99, -0.11, "Primer stage applies to SmC-seq only", transform=ax.transAxes, ha="right", fontsize=5.8, color="#666666")

    plot_failure_pair(
        axes[2], wide, "structure_mapping_rate", "whitelist_mapping_rate",
        "Structure-fail", "Whitelist-miss", "Mapping of excluded read classes", True,
    )
    panel_title(axes[2], "C", "Mapping of excluded read classes")
    plot_failure_pair(
        axes[3], wide, "structure_duplication_rate", "whitelist_duplication_rate",
        "Structure-fail", "Whitelist-miss", "Duplication among mapped excluded reads", False,
    )
    panel_title(axes[3], "D", "Duplication among mapped excluded reads")

    context_handles = [
        Line2D([], [], color="none", label=f"{wide['protocol'].eq('SmC-seq').sum()} SmC-seq run"),
        Line2D([], [], color="none", label=f"{wide['protocol'].eq('Spatial-DMT').sum()} Spatial-DMT runs"),
    ]
    axes[0].legend(handles=axes[0].get_legend_handles_labels()[0] + context_handles,
                   labels=axes[0].get_legend_handles_labels()[1] + [h.get_label() for h in context_handles],
                   loc="lower right", fontsize=5.8)

    stem = output_dir / "Figure2_failed_read_diagnosis"
    save_figure(fig, stem, dpi, tiff_dpi)
    metrics.to_csv(output_dir / "Figure2_metrics.tsv", sep="\t", index=False, na_rep="NA")
    metrics[["sample", "metric", "source_file", "source_type", "definition_version", "is_fallback"]].to_csv(
        output_dir / "Figure2_provenance.tsv", sep="\t", index=False, na_rep="NA"
    )
    LOG.info("Figure 2: %d independent runs", wide.shape[0])


def representative_records() -> tuple[list[SampleRecord], dict[str, SampleRecord]]:
    records = list(DEFAULT_SAMPLES)
    return records, {record.sample: record for record in records}


def plot_cross_context(container, metrics: pd.DataFrame, records: list[SampleRecord]) -> pd.DataFrame:
    container.axis("off")
    panel_title(container, "C", "QC ranges across tissue contexts and modalities", y=1.08)
    specs = (
        ("pass_fraction", "Preprocessing\npass (%)"),
        ("mapping_rate", "Mapping\nrate (%)"),
        ("duplication_rate", "Duplicate\nrate (%)"),
        ("cpg_per_million", "Unique CpGs /\nM raw pairs"),
    )
    sub = container.get_subplotspec().subgridspec(1, 4, wspace=0.38)
    axes = [container.figure.add_subplot(sub[0, i]) for i in range(4)]
    samples = [record.sample for record in records]
    y = np.arange(len(samples))
    derived = []
    for ax, (metric, title) in zip(axes, specs):
        for yi, record in enumerate(records):
            if metric == "cpg_per_million":
                cpg = metric_value(metrics, record.sample, "unique_cpgs")
                raw = metric_value(metrics, record.sample, "raw_read_pairs")
                value = cpg / (raw / 1e6) if np.isfinite(cpg) and raw > 0 else np.nan
                cpg_source = metric_source(metrics, record.sample, "unique_cpgs")
                raw_source = metric_source(metrics, record.sample, "raw_read_pairs")
                if np.isfinite(value):
                    derived.append(
                        {
                            "sample": record.sample,
                            "species": record.species,
                            "tissue": record.tissue,
                            "modality": record.modality,
                            "pipeline": "revised",
                            "metric": metric,
                            "value": value,
                            "unit": "sites per million raw pairs",
                            "source_file": f"{cpg_source['source_file']}; {raw_source['source_file']}",
                            "source_type": "recomputed_revised",
                            "definition_version": "unique_cpg_sites_per_million_raw_pairs_v1",
                            "is_fallback": False,
                            "panel": "Figure3C",
                        }
                    )
            else:
                value = metric_value(metrics, record.sample, metric)
            if np.isfinite(value):
                marker = "s" if record.modality == "DNA+RNA" else "o"
                ax.scatter(value, yi, s=30, marker=marker, color=context_color(record), edgecolor="white", linewidth=0.45, zorder=3)
        ax.set_title(title, fontsize=7.2, pad=4)
        ax.set_yticks(y)
        ax.set_yticklabels([short_sample(sample) for sample in samples] if ax is axes[0] else [])
        ax.invert_yaxis()
        ax.xaxis.set_major_formatter(FuncFormatter(compact_tick))
        clean_axis(ax, "x")
    axes[0].legend(
        handles=[
            Line2D([], [], marker="o", ls="", color="#555555", label="DNA"),
            Line2D([], [], marker="s", ls="", color="#555555", label="DNA+RNA"),
        ],
        loc="lower right",
        fontsize=5.8,
    )
    return pd.DataFrame(derived)


def plot_rna_reporting(container, metrics: pd.DataFrame, records_by_sample: dict[str, SampleRecord]) -> None:
    container.axis("off")
    panel_title(container, "D", "Coordinated DNA and RNA reporting", y=1.08)
    samples = ["SpMETSLE16D", "SpMETSLB16D"]
    specs = (
        ("rna_unique_mapping_rate", "RNA mapping (%)"),
        ("median_rna_umis", "Median UMIs/pixel"),
        ("median_rna_genes", "Median genes/pixel"),
    )
    sub = container.get_subplotspec().subgridspec(1, 3, wspace=0.35)
    for index, (metric, title) in enumerate(specs):
        ax = container.figure.add_subplot(sub[0, index])
        for yi, sample in enumerate(samples):
            record = records_by_sample[sample]
            value = metric_value(metrics, sample, metric)
            ax.scatter(value, yi, s=38, color=context_color(record), edgecolor="white", linewidth=0.5)
            label = f"{value:.1f}" if metric == "rna_unique_mapping_rate" else f"{value:,.0f}"
            ax.annotate(label, (value, yi), xytext=(5, 0), textcoords="offset points", va="center", fontsize=5.8, color=context_color(record))
        ax.set_yticks(range(len(samples)))
        ax.set_yticklabels([short_sample(s) for s in samples] if index == 0 else [])
        ax.invert_yaxis()
        ax.set_title(title, fontsize=7.2, pad=4)
        ax.margins(x=0.22)
        ax.xaxis.set_major_formatter(FuncFormatter(compact_tick))
        clean_axis(ax, "x")


def make_figure3(repo_root: Path, output_dir: Path, dpi: int, tiff_dpi: int) -> None:
    source = repo_root / "results" / "figure2" / "figure2_metrics.tsv"
    if not source.is_file():
        raise FileNotFoundError(source)
    metrics = pd.read_csv(source, sep="\t")
    records, records_by_sample = representative_records()
    expected = {record.sample for record in records}
    if not expected.issubset(set(metrics["sample"])):
        raise ValueError(f"Representative metrics missing: {sorted(expected - set(metrics['sample']))}")

    cohort, cohort_rows = load_panel_b_cohorts(repo_root, repo_root / "results" / "figure2")
    fig = plt.figure(figsize=(15, 10), constrained_layout=False)
    fig.suptitle(
        "Cohort-level QC across spatial methylation experiments",
        fontsize=14,
        weight="bold",
        y=0.985,
    )
    grid = fig.add_gridspec(2, 2, left=0.04, right=0.985, bottom=0.07, top=0.88, hspace=0.3, wspace=0.16)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])
    plot_panel_a_multiqc(ax_a, metrics, records_by_sample)
    plot_panel_b_cohort_retention(ax_b, cohort)
    derived = plot_cross_context(ax_c, metrics, records)
    plot_rna_reporting(ax_d, metrics, records_by_sample)
    save_figure(fig, output_dir / "Figure3_cohort_qc", dpi, tiff_dpi)

    selected = metrics.loc[metrics["sample"].isin(expected) & metrics["metric"].isin(MATRIX_METRICS)].copy()
    selected["panel"] = "Figure3A/C/D"
    cohort_frame = pd.DataFrame(cohort_rows)
    cohort_frame["panel"] = "Figure3B"
    combined = pd.concat([selected, cohort_frame, derived], ignore_index=True, sort=False)
    combined.to_csv(output_dir / "Figure3_metrics.tsv", sep="\t", index=False, na_rep="NA")
    combined[["sample", "metric", "source_file", "source_type", "definition_version", "is_fallback", "panel"]].to_csv(
        output_dir / "Figure3_provenance.tsv", sep="\t", index=False, na_rep="NA"
    )
    pd.DataFrame([record.__dict__ for record in records]).to_csv(
        output_dir / "Figure3_selected_samples.tsv", sep="\t", index=False
    )
    LOG.info("Figure 3: %d representative samples and two published cohorts", len(records))


def make_figure4(repo_root: Path, output_dir: Path, dpi: int, tiff_dpi: int) -> None:
    comparison, provenance_rows = load_pipeline_comparison(repo_root)
    for row in provenance_rows:
        references = [item.strip() for item in str(row["source_file"]).split(";") if item.strip()]
        existing = [item for item in references if Path(item).is_file()]
        missing = [item for item in references if not Path(item).is_file()]
        row["source_file"] = "; ".join(existing)
        row["historical_upstream_reference"] = "; ".join(missing)
        row["source_status"] = "audited_table_only" if missing else "all_references_found"
    pipelines = ("spmetsle14d", "e65_bsconv", "e65_reproduced")
    labels = {"spmetsle14d": "14D revised\n(context)", "e65_bsconv": "E6.5 ours", "e65_reproduced": "E6.5 reproduced"}
    colors = {"spmetsle14d": "#0072B2", "e65_bsconv": "#E69F00", "e65_reproduced": "#6F6F6F"}
    markers = {"spmetsle14d": "o", "e65_bsconv": "^", "e65_reproduced": "D"}

    fig = plt.figure(figsize=(15, 6.5), constrained_layout=False)
    fig.suptitle(
        "Definition-aware comparison of spatial-methylation processing",
        fontsize=14,
        weight="bold",
        y=0.98,
    )
    grid = fig.add_gridspec(1, 2, width_ratios=(0.9, 1.45), left=0.07, right=0.98, bottom=0.14, top=0.88, wspace=0.22)
    ax = fig.add_subplot(grid[0, 0])
    panel_title(ax, "A", "Stage-wise retention")
    stages = ("Raw", "Designed", "Mapped", "CpH-filter\npassed")
    stage_metrics = ("Raw read pairs", "Designed", "Mapped before CpH filter", "CpH-filter passed")
    for pipeline in pipelines:
        raw = _comparison_value(comparison, pipeline, "Raw read pairs")
        values = np.array([_comparison_value(comparison, pipeline, metric) for metric in stage_metrics])
        ax.plot(range(4), 100 * values / raw, color=colors[pipeline], marker=markers[pipeline], lw=1.5, ms=5, label=labels[pipeline])
        for x, value, count in zip(range(1, 4), 100 * values[1:] / raw, values[1:]):
            if np.isfinite(value):
                ax.annotate(f"{count / 1e6:.1f}M", (x, value), xytext=(3, 5), textcoords="offset points", fontsize=5.5, color=colors[pipeline])
    ax.set_xticks(range(4), stages)
    ax.set_ylim(0, 105)
    ax.set_ylabel("Retained (% of raw read pairs)")
    ax.legend(loc="lower left", fontsize=6.2)
    clean_axis(ax, "y")
    ax.text(0.01, -0.15, "14D is a contextual Spatial-DMT reference; the two E6.5 series use the same raw input.", transform=ax.transAxes, fontsize=5.8, color="#666666")

    ax = fig.add_subplot(grid[0, 1])
    panel_title(ax, "B", "Chemistry and conversion-filter diagnostics")
    diagnostic_metrics = (
        "Mapping rate",
        "CpH-filter passed / mapped",
        "Duplication",
        "CpG methylation",
        "CpH (A/T/G)",
        "Lambda retention",
    )
    ax.set_xlim(-1.9, 2.55)
    ax.set_ylim(len(diagnostic_metrics) - 0.35, -1.2)
    ax.axis("off")
    for j, pipeline in enumerate(pipelines):
        ax.text(j, -0.62, labels[pipeline], ha="center", va="bottom", fontsize=6.6, weight="bold", color=colors[pipeline])
    for i, metric in enumerate(diagnostic_metrics):
        ax.text(-1.82, i, metric.replace(" / ", " /\n") if metric == "CpH-filter passed / mapped" else metric, ha="left", va="center", fontsize=6.5)
        for j, pipeline in enumerate(pipelines):
            comparable = not (
                metric == "Duplication" and pipeline == "e65_reproduced"
            )
            face = "#FFFFFF" if comparable else "#F0F0F0"
            if pipeline == "e65_bsconv" and metric in {"CpH-filter passed / mapped", "CpG methylation", "CpH (A/T/G)"}:
                face = "#FFF2CC"
            ax.add_patch(Rectangle((j - 0.46, i - 0.34), 0.92, 0.68, facecolor=face, edgecolor="#CCCCCC", lw=0.6))
            display = _comparison_display(comparison, pipeline, metric)
            if metric == "CpH (A/T/G)":
                display = display.replace(" / ", " /\n")
            ax.text(j, i, display, ha="center", va="center", fontsize=5.8 if "\n" not in display else 5.2)
    ax.text(
        1, -0.07,
        "E6.5 ours: complete bsconv -y 0.05 result. Reproduced 'designed' denotes structure-pass reads.\n"
        "Blank cells are unavailable; duplication definitions are not compared across incompatible implementations.",
        transform=ax.transAxes, ha="right", va="top", fontsize=5.7, color="#666666",
    )
    save_figure(fig, output_dir / "Figure4_pipeline_comparison", dpi, tiff_dpi)
    provenance = pd.DataFrame(provenance_rows)
    provenance.to_csv(output_dir / "Figure4_metrics.tsv", sep="\t", index=False, na_rep="NA")
    provenance[["sample", "metric", "source_file", "source_type", "definition_version", "is_fallback", "panel", "source_status", "historical_upstream_reference"]].to_csv(
        output_dir / "Figure4_provenance.tsv", sep="\t", index=False, na_rep="NA"
    )
    LOG.info("Figure 4: same-input E6.5 comparison plus contextual 14D reference")


def spatial_sources(repo_root: Path) -> tuple[Path, Path, Path]:
    embryo = RUN_ROOTS[0] / "pipeline_output" / "SpMETSLE14DM" / "qc" / "table" / "SpMETSLE14DM_spatial_barcode_counts.tsv"
    skin = RUN_ROOTS[3] / "pipeline_output" / "SpMETSLHSK09D" / "qc" / "table" / "SpMETSLHSK09D_spatial_barcode_counts.tsv"
    whitelist = repo_root / "barcodes" / "spatial_barcodes.txt"
    for path in (embryo, skin, whitelist):
        if not path.is_file():
            raise FileNotFoundError(path)
    return embryo, skin, whitelist


def plot_spatial_map(ax, matrix: np.ndarray, frame: pd.DataFrame, norm, title: str) -> None:
    ax.imshow(matrix, cmap=SPATIAL_CMAP, norm=norm, interpolation="none", extent=_spatial_extent(frame), aspect="equal", rasterized=True)
    ax.set_title(title, fontsize=8, pad=4)
    ax.set_xlabel("x")
    ax.set_ylabel("y")


def make_figure5(repo_root: Path, output_dir: Path, dpi: int, tiff_dpi: int) -> None:
    embryo_path, skin_path, whitelist_path = spatial_sources(repo_root)
    embryo = pd.read_csv(embryo_path, sep="\t")
    skin = pd.read_csv(skin_path, sep="\t")
    whitelist = read_legacy_whitelist(whitelist_path)
    incorrect = legacy_order_reconstruction(embryo, whitelist)
    good_matrix = spatial_matrix(embryo, "dna_reads")
    bad_matrix = spatial_matrix(incorrect, "dna_reads")
    norm = _positive_lognorm(good_matrix, bad_matrix)
    if not np.array_equal(np.sort(good_matrix.ravel()), np.sort(bad_matrix.ravel())):
        raise ValueError("Controlled coordinate reconstruction did not preserve counts")

    fig = plt.figure(figsize=(15, 9.5), constrained_layout=False)
    fig.suptitle(
        "Spatial diagnostics expose failures invisible to aggregate sequencing QC",
        fontsize=14,
        weight="bold",
        y=0.985,
    )
    grid = fig.add_gridspec(2, 2, left=0.055, right=0.98, bottom=0.07, top=0.93, hspace=0.24, wspace=0.18)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    panel_title(ax_a, "A", "Correct coordinate reconstruction")
    panel_title(ax_b, "B", "Incompatible coordinate reconstruction")
    plot_spatial_map(ax_a, good_matrix, embryo, norm, "SpMETSLE14DM: revised coordinates")
    plot_spatial_map(ax_b, bad_matrix, embryo, norm, "Same counts, legacy barcode order")
    ax_b.text(0.98, 0.02, "Color = raw DNA reads/pixel\n(shared log normalization)", transform=ax_b.transAxes, ha="right", va="bottom", fontsize=5.8, color="#555555")

    ax_c = fig.add_subplot(grid[1, 0])
    panel_title(ax_c, "C", "Aggregate identity versus spatial continuity")
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    ax_c.set_autoscale_on(False)
    ax_c.axis("off")
    total = float(embryo["dna_reads"].sum())
    mapping = 100 * float(embryo["mapped_reads"].sum()) / total
    median = float(embryo.loc[embryo["dna_reads"] > 0, "dna_reads"].median())
    covered = 100 * float((embryo["dna_reads"] > 0).mean())
    aggregate = (
        ("Total DNA fragments", f"{total / 1e9:.3f}B", f"{total / 1e9:.3f}B"),
        ("Mapping rate", f"{mapping:.1f}%", f"{mapping:.1f}%"),
        ("Median reads/covered pixel", f"{median:,.0f}", f"{median:,.0f}"),
        ("Covered pixels", f"{covered:.1f}%", f"{covered:.1f}%"),
    )
    ax_c.text(0.02, 0.88, "Aggregate metric", weight="bold", fontsize=7.2, transform=ax_c.transAxes)
    ax_c.text(0.55, 0.88, "Correct", ha="center", weight="bold", fontsize=7.2, color="#0072B2", transform=ax_c.transAxes)
    ax_c.text(0.79, 0.88, "Incorrect", ha="center", weight="bold", fontsize=7.2, color=WARNING_COLOR, transform=ax_c.transAxes)
    for i, (label, good, bad) in enumerate(aggregate):
        y = 0.78 - i * 0.105
        ax_c.text(0.02, y, label, fontsize=6.6, transform=ax_c.transAxes)
        ax_c.text(0.55, y, good, ha="center", fontsize=6.6, transform=ax_c.transAxes)
        ax_c.text(0.79, y, bad, ha="center", fontsize=6.6, transform=ax_c.transAxes)
        ax_c.plot([0.02, 0.92], [y - 0.045, y - 0.045], color="#E0E0E0", lw=0.5, transform=ax_c.transAxes)
    diagnostics = {
        "Moran's I": (morans_i(good_matrix), morans_i(bad_matrix)),
        "Neighbor Spearman": (neighbor_spearman(good_matrix), neighbor_spearman(bad_matrix)),
    }
    inset = ax_c.inset_axes([0.08, 0.02, 0.82, 0.3])
    y = np.arange(len(diagnostics))
    width = 0.32
    inset.barh(y - width / 2, [v[0] for v in diagnostics.values()], height=width, color="#0072B2", label="Correct")
    inset.barh(y + width / 2, [v[1] for v in diagnostics.values()], height=width, color=WARNING_COLOR, hatch="//", label="Incorrect")
    inset.set_yticks(y, diagnostics.keys())
    inset.set_xlim(0, 1)
    inset.set_xlabel("Spatial continuity coefficient")
    inset.legend(loc="lower right", fontsize=5.8)
    clean_axis(inset, "x")

    container = fig.add_subplot(grid[1, 1])
    container.axis("off")
    panel_title(container, "D", "Localized low-read spatial scar")
    sub = grid[1, 1].subgridspec(1, 2, width_ratios=(1.1, 0.8), wspace=0.24)
    map_ax = fig.add_subplot(sub[0, 0])
    profile_ax = fig.add_subplot(sub[0, 1])
    skin_matrix = spatial_matrix(skin, "dna_reads")
    call = detect_scar(skin_matrix)
    plot_spatial_map(map_ax, skin_matrix, skin, _positive_lognorm(skin_matrix), "SpMETSLHSK09D: raw DNA reads")
    if call.axis == "horizontal":
        map_ax.add_patch(Rectangle((0.5, call.start - 0.5), 50, call.end - call.start + 1, fill=False, edgecolor=WARNING_COLOR, lw=1.6))
        map_ax.annotate("low-read scar", xy=(47, (call.start + call.end) / 2), xytext=(36, max(3, call.start - 8)), color=WARNING_COLOR, fontsize=6.5, arrowprops={"arrowstyle": "->", "color": WARNING_COLOR, "lw": 0.8})
        coords = np.arange(1, len(call.profile) + 1)
        profile_ax.plot(call.profile, coords, color="#333333", lw=1.2)
        profile_ax.axhspan(call.start - 0.5, call.end + 0.5, color=WARNING_COLOR, alpha=0.2)
        profile_ax.invert_yaxis()
        profile_ax.set_xlabel("Median raw reads")
        profile_ax.set_ylabel("Row")
    else:
        map_ax.add_patch(Rectangle((call.start - 0.5, 0.5), call.end - call.start + 1, 50, fill=False, edgecolor=WARNING_COLOR, lw=1.6))
        coords = np.arange(1, len(call.profile) + 1)
        profile_ax.plot(coords, call.profile, color="#333333", lw=1.2)
        profile_ax.axvspan(call.start - 0.5, call.end + 0.5, color=WARNING_COLOR, alpha=0.2)
        profile_ax.set_xlabel("Column")
        profile_ax.set_ylabel("Median raw reads")
    interval = str(call.start) if call.start == call.end else f"{call.start}-{call.end}"
    profile_ax.set_title(f"{call.axis.title()} profile", fontsize=7.2)
    profile_ax.text(0.98, 0.98, f"Highlighted interval: {interval}\ninside / adjacent = {call.inside_median:.0f} / {call.adjacent_median:.0f}\nscar ratio = {call.ratio:.3f}", transform=profile_ax.transAxes, ha="right", va="top", color=WARNING_COLOR, fontsize=6)
    clean_axis(profile_ax, "x" if call.axis == "horizontal" else "y")

    save_figure(fig, output_dir / "Figure5_spatial_failures", dpi, tiff_dpi)
    rows = []
    for metric, value, unit in (
        ("total_dna_fragments", total, "fragments"),
        ("mapping_rate", mapping, "percent"),
        ("median_reads_per_covered_pixel", median, "reads/pixel"),
        ("covered_pixel_fraction", covered, "percent"),
        ("correct_morans_i", diagnostics["Moran's I"][0], "coefficient"),
        ("incorrect_morans_i", diagnostics["Moran's I"][1], "coefficient"),
        ("correct_neighbor_spearman", diagnostics["Neighbor Spearman"][0], "coefficient"),
        ("incorrect_neighbor_spearman", diagnostics["Neighbor Spearman"][1], "coefficient"),
    ):
        rows.append({"sample": "SpMETSLE14DM", "metric": metric, "value": value, "unit": unit, "source_file": f"{embryo_path}; {whitelist_path}", "source_type": "recomputed_revised", "definition_version": "controlled_coordinate_reconstruction_v1", "is_fallback": False, "panel": "Figure5A-C"})
    for metric, value, unit in (
        ("scar_ratio", call.ratio, "ratio"),
        ("scar_inside_median", call.inside_median, "reads/pixel"),
        ("scar_adjacent_median", call.adjacent_median, "reads/pixel"),
        ("scar_start", call.start, "row_or_column"),
        ("scar_end", call.end, "row_or_column"),
    ):
        rows.append({"sample": "SpMETSLHSK09D", "metric": metric, "value": value, "unit": unit, "source_file": str(skin_path), "source_type": "recomputed_revised", "definition_version": "local_median_scar_summary_v1", "is_fallback": False, "panel": "Figure5D"})
    metrics = pd.DataFrame(rows)
    metrics.to_csv(output_dir / "Figure5_metrics.tsv", sep="\t", index=False)
    metrics[["sample", "metric", "source_file", "source_type", "definition_version", "is_fallback", "panel"]].to_csv(
        output_dir / "Figure5_provenance.tsv", sep="\t", index=False
    )
    LOG.info("Figure 5: controlled coordinate reconstruction and %s scar interval %s", call.axis, interval)


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    apply_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    builders = {
        "2": make_figure2,
        "3": make_figure3,
        "4": make_figure4,
        "5": make_figure5,
    }
    for figure in args.figures:
        LOG.info("Building Figure %s", figure)
        builders[figure](args.repo_root.resolve(), args.output_dir.resolve(), args.dpi, args.tiff_dpi)
    LOG.info("Outputs written to %s", args.output_dir.resolve())


if __name__ == "__main__":
    main()
