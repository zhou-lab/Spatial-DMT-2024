#!/usr/bin/env python3
"""Build Figure 2 from revised and explicitly marked legacy QC outputs."""

from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch, Rectangle
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from load_metrics import (  # noqa: E402
    METRIC_META,
    SampleRecord,
    collect_sample_metrics,
    find_spatial_barcode_counts,
    find_trim_stats,
    legacy_pixel_counts,
    load_sample_config,
    read_kv,
    read_legacy_whitelist,
    records_frame,
    run_yame_aggregate,
)
from plot_style import (  # noqa: E402
    CONTEXT_COLORS,
    MATRIX_CMAP,
    MISSING_COLOR,
    PIPELINE_COLORS,
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
from build_panel_b_published_band import (  # noqa: E402
    RETENTION_METRICS,
    load_panel_b_metrics,
    plot_retention,
)


LOG = logging.getLogger("figure2")
TECHNICAL_RUNS = {"SpMETSLE14D_smoke", "SpMETSLE16D_ourpipeline_bsconv_y005"}
MATRIX_METRICS = [
    "raw_read_pairs",
    "pass_fraction",
    "mapping_rate",
    "duplication_rate",
    "conversion_efficiency",
    "unique_cpgs",
    "median_unique_fragments",
    "rna_unique_mapping_rate",
    "median_rna_umis",
    "median_rna_genes",
]

MULTIQC_COLORS = {
    "raw_read_pairs": "#2C7FB8",
    "pass_fraction": "#41B6C4",
    "mapping_rate": "#2CA25F",
    "duplication_rate": "#DE2D26",
    "conversion_efficiency": "#31A354",
    "unique_cpgs": "#756BB1",
    "median_unique_fragments": "#636363",
    "rna_unique_mapping_rate": "#E6550D",
    "median_rna_umis": "#FD8D3C",
    "median_rna_genes": "#FDAE6B",
}


def metric_value(metrics: pd.DataFrame, sample: str, metric: str) -> float:
    selected = metrics.loc[(metrics["sample"] == sample) & (metrics["metric"] == metric)]
    if selected.empty:
        return np.nan
    revised = selected.loc[~selected["is_fallback"].astype(bool)]
    return float((revised if not revised.empty else selected).iloc[0]["value"])


def metric_source(metrics: pd.DataFrame, sample: str, metric: str) -> pd.Series | None:
    selected = metrics.loc[(metrics["sample"] == sample) & (metrics["metric"] == metric)]
    if selected.empty:
        return None
    revised = selected.loc[~selected["is_fallback"].astype(bool)]
    return (revised if not revised.empty else selected).iloc[0]


def apply_revised_metric_snapshot(
    metrics: pd.DataFrame, records: list[SampleRecord], snapshot_path: Path
) -> pd.DataFrame:
    """Restore audited revised cells whose temporary upstream files were cleaned."""
    if not snapshot_path.is_file():
        raise FileNotFoundError(f"Revised metric snapshot not found: {snapshot_path}")
    snapshot = pd.read_csv(snapshot_path, sep="\t")
    required = {
        "sample", "metric", "value", "unit", "definition_version", "original_source"
    }
    missing = required - set(snapshot.columns)
    if missing:
        raise ValueError(f"Revised metric snapshot lacks columns: {sorted(missing)}")
    if snapshot.duplicated(["sample", "metric"]).any():
        raise ValueError("Revised metric snapshot has duplicate sample/metric rows")
    if pd.to_numeric(snapshot["value"], errors="coerce").isna().any():
        raise ValueError("Revised metric snapshot contains a missing or nonnumeric value")

    records_by_sample = {record.sample: record for record in records}
    restored: list[dict] = []
    result = metrics.copy()
    for item in snapshot.itertuples(index=False):
        record = records_by_sample.get(item.sample)
        if record is None:
            continue
        same_cell = result["sample"].eq(item.sample) & result["metric"].eq(item.metric)
        live_revised = same_cell & ~result["is_fallback"].astype(bool)
        if live_revised.any():
            continue
        result = result.loc[~same_cell].copy()
        restored.append(
            {
                "sample": record.sample,
                "species": record.species,
                "tissue": record.tissue,
                "modality": record.modality,
                "pipeline": "revised",
                "metric": item.metric,
                "value": float(item.value),
                "unit": item.unit,
                "source_file": str(snapshot_path),
                "source_type": "revised_snapshot",
                "definition_version": f"{item.definition_version}; upstream={item.original_source}",
                "is_fallback": False,
                "panel": "A",
            }
        )
    if restored:
        LOG.info(
            "Restored %d revised Panel A cells from audited snapshot %s",
            len(restored),
            snapshot_path,
        )
        result = pd.concat([result, pd.DataFrame(restored)], ignore_index=True)
    return result


def pretty_value(metric: str, value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if metric == "raw_read_pairs":
        return f"{value / 1e9:.2f}B" if value >= 1e9 else f"{value / 1e6:.1f}M"
    if metric == "unique_cpgs":
        return f"{value / 1e6:.1f}M"
    if METRIC_META.get(metric, ("", "", ""))[1].startswith("%"):
        return f"{value:.1f}"
    if value >= 1e4:
        return f"{value / 1e3:.0f}k"
    if value >= 100:
        return f"{value:.0f}"
    return f"{value:.1f}"


def compact_tick(value: float, _position: int | None = None) -> str:
    magnitude = abs(value)
    if magnitude >= 1e6:
        return f"{value / 1e6:g}M"
    if magnitude >= 1e3:
        return f"{value / 1e3:g}k"
    return f"{value:g}"


def derived_row(
    record: SampleRecord,
    metric: str,
    value: float,
    unit: str,
    source_file: str | Path,
    source_type: str,
    definition: str,
    panel: str,
    pipeline: str = "revised",
    is_fallback: bool = False,
) -> dict:
    return {
        "sample": record.sample,
        "species": record.species,
        "tissue": record.tissue,
        "modality": record.modality,
        "pipeline": pipeline,
        "metric": metric,
        "value": float(value) if value is not None else np.nan,
        "unit": unit,
        "source_file": str(source_file),
        "source_type": source_type,
        "definition_version": definition,
        "is_fallback": bool(is_fallback),
        "panel": panel,
    }


def context_color(record: SampleRecord) -> str:
    return CONTEXT_COLORS[f"{record.species} {record.tissue}"]


def plot_panel_a_qc_overview(
    ax, metrics: pd.DataFrame, records: list[SampleRecord]
) -> None:
    panel_title(ax, "A", "Representative cohort QC overview")
    samples = [record.sample for record in records]
    raw = np.full((len(samples), len(MATRIX_METRICS)), np.nan)
    fallback = np.zeros_like(raw, dtype=bool)
    for i, sample in enumerate(samples):
        for j, metric in enumerate(MATRIX_METRICS):
            raw[i, j] = metric_value(metrics, sample, metric)
            source = metric_source(metrics, sample, metric)
            fallback[i, j] = bool(source is not None and source["is_fallback"])

    normalized = np.full_like(raw, np.nan)
    for j in range(raw.shape[1]):
        values = raw[:, j]
        finite = np.isfinite(values)
        if not finite.any():
            continue
        lo, hi = np.nanmin(values), np.nanmax(values)
        normalized[finite, j] = 0.5 if hi == lo else (values[finite] - lo) / (hi - lo)
        if MATRIX_METRICS[j] == "duplication_rate":
            normalized[finite, j] = 1 - normalized[finite, j]

    masked = np.ma.masked_invalid(normalized)
    cmap = MATRIX_CMAP.copy()
    cmap.set_bad(MISSING_COLOR)
    ax.imshow(masked, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="none")
    ax.set_xticks(range(len(MATRIX_METRICS)))
    ax.set_xticklabels(
        [METRIC_META[m][0].replace(" ", "\n") for m in MATRIX_METRICS], rotation=42, ha="right"
    )
    labels = [f"{r.sample}  [{r.modality}]" for r in records]
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(labels)
    ax.tick_params(length=0)
    ax.set_xticks(np.arange(-0.5, len(MATRIX_METRICS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(samples), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False)

    for i, record in enumerate(records):
        ax.add_patch(
            Rectangle(
                (-1.04, i - 0.48), 0.14, 0.96, color=context_color(record), clip_on=False, linewidth=0
            )
        )
        for j, metric in enumerate(MATRIX_METRICS):
            if np.isfinite(raw[i, j]):
                label = pretty_value(metric, raw[i, j]) + ("*" if fallback[i, j] else "")
                foreground = "white" if normalized[i, j] > 0.66 else "#202020"
                ax.text(j, i, label, ha="center", va="center", fontsize=5.4, color=foreground)
            else:
                ax.text(j, i, "x", ha="center", va="center", color="#8A8A8A", fontsize=6.5)

    legend = [
        Patch(facecolor=CONTEXT_COLORS["mouse embryo"], label="mouse embryo"),
        Patch(facecolor=CONTEXT_COLORS["mouse brain"], label="mouse brain"),
        Patch(facecolor=CONTEXT_COLORS["human skin"], label="human skin"),
    ]
    ax.legend(handles=legend, ncol=3, loc="lower left", bbox_to_anchor=(0, 1.005), borderaxespad=0)
    ax.text(
        1,
        -0.21,
        "Cell color: within-metric range (duplicates reversed); x, unavailable; *, legacy fallback",
        transform=ax.transAxes,
        ha="right",
        color="#555555",
        fontsize=5.7,
    )
    for spine in ax.spines.values():
        spine.set_visible(False)


def _retention_values(record: SampleRecord, metrics: pd.DataFrame) -> tuple[list[str], np.ndarray, str]:
    trim_path, _, _ = find_trim_stats(record)
    spatial_path, _, _ = find_spatial_barcode_counts(record)
    trim = read_kv(trim_path) if trim_path else {}
    spatial = pd.read_csv(spatial_path, sep="\t") if spatial_path else pd.DataFrame()
    raw = trim.get("Reads_total", trim.get("Reads_detected", np.nan))
    detected = trim.get("Reads_detected", np.nan)
    passed = trim.get("Reads_passed", np.nan)
    assigned = pd.to_numeric(spatial.get("dna_reads"), errors="coerce").sum() if not spatial.empty else np.nan
    mapped = pd.to_numeric(spatial.get("mapped_reads"), errors="coerce").sum() if not spatial.empty else np.nan
    dup = pd.to_numeric(spatial.get("dup_reads"), errors="coerce").sum() if not spatial.empty else np.nan
    labels = ["Input", "Structure\ndetected", "Passed", "Spatially\nassigned", "Mapped", "Unique"]
    values = np.array([raw, detected, passed, assigned, mapped, mapped - dup], dtype=float)
    sources = f"{trim_path}; {spatial_path}"
    return labels, 100 * values / raw, sources


def plot_panel_b_depth_comparison(
    container,
    metrics: pd.DataFrame,
    records_by_sample: dict[str, SampleRecord],
    provenance_rows: list[dict],
) -> None:
    container.axis("off")
    panel_title(container, "B", "Depth-aware QC in a matched shallow and deep sample")
    sub = container.get_subplotspec().subgridspec(1, 2, width_ratios=[1.35, 1], wspace=0.35)
    ax1 = container.figure.add_subplot(sub[0, 0])
    ax2 = container.figure.add_subplot(sub[0, 1])
    colors = {"SpMETSLE14D": "#56B4E9", "SpMETSLE14DM": "#005A8D"}

    for sample, marker in (("SpMETSLE14D", "o"), ("SpMETSLE14DM", "s")):
        record = records_by_sample[sample]
        labels, retention, sources = _retention_values(record, metrics)
        ax1.plot(range(len(labels)), retention, marker=marker, ms=4, lw=1.3, color=colors[sample], label=sample)
        for stage, value in zip(labels, retention):
            provenance_rows.append(
                derived_row(record, f"retention_{stage.replace(chr(10), '_').lower()}", value, "% input", sources, "recomputed_revised", "stage_count/Reads_total_v1", "B")
            )
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, rotation=35, ha="right")
    ax1.set_ylabel("Retained (% of own input)")
    ax1.set_ylim(0, 105)
    ax1.legend(loc="lower left")
    clean_axis(ax1, "y")

    outcome_metrics = ["unique_cpgs", "median_unique_fragments", "covered_pixel_fraction"]
    outcome_labels = ["CpG sites", "Fragments/pixel", "Covered pixels"]
    y = np.arange(len(outcome_metrics))
    for i, metric in enumerate(outcome_metrics):
        shallow = metric_value(metrics, "SpMETSLE14D", metric)
        deep = metric_value(metrics, "SpMETSLE14DM", metric)
        if metric == "unique_cpgs":
            shallow /= 1e6
            deep /= 1e6
        if np.isfinite(shallow) and np.isfinite(deep) and shallow > 0:
            gain = deep / shallow
            ax2.plot([1, gain], [i, i], color="#BDBDBD", lw=1.2, zorder=1)
            ax2.scatter([1, gain], [i, i], c=[colors["SpMETSLE14D"], colors["SpMETSLE14DM"]], s=27, zorder=2)
            ax2.text(gain, i + 0.18, f"{gain:.1f}x", ha="center", va="bottom", fontsize=6.2)
            record = records_by_sample["SpMETSLE14DM"]
            source = metric_source(metrics, record.sample, metric)
            provenance_rows.append(
                derived_row(record, f"deep_to_shallow_{metric}", gain, "fold", source["source_file"] if source is not None else "missing", source["source_type"] if source is not None else "missing", "deep_value/shallow_value_v1", "B", is_fallback=bool(source is not None and source["is_fallback"]))
            )
    ax2.axvline(1, color="#555555", lw=0.7, ls="--")
    ax2.set_xscale("log")
    ax2.set_yticks(y)
    ax2.set_yticklabels(outcome_labels)
    ax2.invert_yaxis()
    ax2.set_xlabel("Deep/shallow gain (fold)")
    ax2.set_title("Outcome gain (deep input: 2.864B pairs)", fontsize=7, pad=3)
    clean_axis(ax2, "x")


def plot_panel_c_cross_context(
    container, metrics: pd.DataFrame, records_by_sample: dict[str, SampleRecord]
) -> None:
    container.axis("off")
    panel_title(container, "C", "QC across tissue contexts and assay modalities", y=1.14)
    samples = ["SpMETSLE16D", "SpMETSLB16D", "SpMETSLHSK03D", "SpMETSLHSK05D"]
    metric_specs = [
        ("mapping_rate", "DNA\nmapped (%)", 1.0),
        ("unique_cpgs", "CpGs per\nM input", None),
        ("median_unique_fragments", "Median unique\nfragments/pixel", 1.0),
        ("covered_pixel_fraction", "Covered\npixels (%)", 1.0),
    ]
    sub = container.get_subplotspec().subgridspec(1, 5, wspace=0.55)
    axes = [container.figure.add_subplot(sub[0, i]) for i in range(5)]
    y = np.arange(len(samples))
    for ax, (metric, title, scale) in zip(axes[:4], metric_specs):
        values = []
        for sample in samples:
            value = metric_value(metrics, sample, metric)
            if metric == "unique_cpgs":
                raw = metric_value(metrics, sample, "raw_read_pairs")
                value = value / (raw / 1e6) if np.isfinite(value) and raw > 0 else np.nan
            values.append(value)
        for yi, sample, value in zip(y, samples, values):
            record = records_by_sample[sample]
            ax.scatter(value, yi, s=27, color=context_color(record), edgecolor="white", linewidth=0.45, zorder=3)
        ax.set_title(title, fontsize=7, pad=3)
        ax.xaxis.set_major_formatter(FuncFormatter(compact_tick))
        ax.set_yticks(y)
        ax.invert_yaxis()
        clean_axis(ax, "x")
    axes[0].set_yticklabels([s.replace("SpMETSL", "") for s in samples])
    for ax in axes[1:4]:
        ax.set_yticklabels([])

    ax = axes[4]
    rna_samples = [s for s in samples if records_by_sample[s].modality == "DNA+RNA"]
    for yi, sample in enumerate(rna_samples):
        record = records_by_sample[sample]
        umis = metric_value(metrics, sample, "median_rna_umis")
        genes = metric_value(metrics, sample, "median_rna_genes")
        ax.scatter(umis, yi, marker="o", s=28, color=context_color(record), label="UMIs" if yi == 0 else None)
        ax.scatter(genes, yi, marker="^", s=30, facecolor="white", edgecolor=context_color(record), linewidth=1, label="Genes" if yi == 0 else None)
    ax.set_title("RNA complexity\nper pixel", fontsize=7, pad=3)
    ax.set_yticks(range(len(rna_samples)))
    ax.set_yticklabels([s.replace("SpMETSL", "") for s in rna_samples])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.legend(loc="lower right", handletextpad=0.3)
    clean_axis(ax, "x")


def _legacy_trim(record: SampleRecord) -> tuple[dict[str, float], Path | None]:
    if record.legacy_root is None:
        return {}, None
    path = record.legacy_root / "trim" / record.sample / f"{record.sample}_stats.txt"
    return (read_kv(path), path) if path.is_file() else ({}, None)


def plot_panel_d_pipeline_comparison(
    container,
    metrics: pd.DataFrame,
    record: SampleRecord,
    revised_spatial: pd.DataFrame,
    legacy_spatial: pd.DataFrame,
    legacy_cache: Path,
    provenance_rows: list[dict],
) -> None:
    container.axis("off")
    panel_title(container, "D", "Diagnostic comparison of revised and legacy processing")
    sub = container.get_subplotspec().subgridspec(1, 2, width_ratios=[1, 1.15], wspace=0.35)
    ax1 = container.figure.add_subplot(sub[0, 0])
    ax2 = container.figure.add_subplot(sub[0, 1])
    trim_path, _, _ = find_trim_stats(record)
    revised_trim = read_kv(trim_path)
    legacy_trim, legacy_trim_path = _legacy_trim(record)
    raw = revised_trim.get("Reads_total", np.nan)
    stages = ["Input", "Structure\ndetected", "Passed", "Mapped", "Unique"]
    revised_mapped = revised_spatial["mapped_reads"].sum()
    revised_unique = (revised_spatial["mapped_reads"] - revised_spatial["dup_reads"]).sum()
    legacy_mapped = legacy_spatial["mapped_reads"].sum()
    legacy_unique = legacy_spatial["unique_fragments"].sum()
    values = {
        "revised": np.array([raw, revised_trim.get("Reads_detected", np.nan), revised_trim.get("Reads_passed", np.nan), revised_mapped, revised_unique]),
        "legacy": np.array([raw, legacy_trim.get("Reads_detected", np.nan), legacy_trim.get("Reads_passed", np.nan), legacy_mapped, legacy_unique]),
    }
    for pipeline, marker, ls in (("revised", "o", "-"), ("legacy", "s", "--")):
        fraction = 100 * values[pipeline] / raw
        ax1.plot(range(len(stages)), fraction, marker=marker, ls=ls, lw=1.3, ms=4, color=PIPELINE_COLORS[pipeline], label=pipeline)
        for stage, value in zip(stages, fraction):
            source = f"{trim_path}; revised spatial table" if pipeline == "revised" else f"{legacy_trim_path}; {legacy_cache}"
            provenance_rows.append(
                derived_row(record, f"{pipeline}_retention_{stage.replace(chr(10), '_').lower()}", value, "% common raw input", source, "recomputed_revised" if pipeline == "revised" else "recomputed_legacy", "primary_mapped_fragment_retention_v1", "D", pipeline=pipeline)
            )
    ax1.set_xticks(range(len(stages)))
    ax1.set_xticklabels(stages, rotation=35, ha="right")
    ax1.set_ylabel("Retained (% common input)")
    ax1.set_ylim(0, 105)
    ax1.legend(loc="lower left")
    clean_axis(ax1, "y")

    revised = revised_spatial.assign(
        revised_unique=revised_spatial["mapped_reads"] - revised_spatial["dup_reads"]
    )[["x", "y", "revised_unique"]]
    legacy = legacy_spatial[["x", "y", "unique_fragments"]].rename(columns={"unique_fragments": "legacy_unique"})
    paired = legacy.merge(revised, on=["x", "y"], how="inner")
    x = paired["legacy_unique"].to_numpy(dtype=float)
    y = paired["revised_unique"].to_numpy(dtype=float)
    rho = float(spearmanr(x, y).statistic)
    ax2.scatter(x + 1, y + 1, s=5, alpha=0.28, color="#0072B2", rasterized=True, linewidth=0)
    limit = max(np.nanmax(x + 1), np.nanmax(y + 1))
    ax2.plot([1, limit], [1, limit], color="#555555", ls=":", lw=0.8)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Legacy unique fragments/pixel + 1")
    ax2.set_ylabel("Revised unique fragments/pixel + 1")
    ax2.text(0.04, 0.96, f"Spearman rho = {rho:.2f}\nn = {len(paired):,} matched pixels", transform=ax2.transAxes, va="top")
    clean_axis(ax2, "both")
    provenance_rows.append(
        derived_row(record, "pixel_unique_fragment_spearman", rho, "rho", f"{legacy_cache}; revised spatial table", "recomputed_legacy", "matched_xy_primary_unique_fragments_v1", "D", pipeline="matched")
    )


def _spatial_extent(frame: pd.DataFrame) -> tuple[float, float, float, float]:
    return (
        float(frame["x"].min()) - 0.5,
        float(frame["x"].max()) + 0.5,
        float(frame["y"].max()) + 0.5,
        float(frame["y"].min()) - 0.5,
    )


def _positive_lognorm(*matrices: np.ndarray) -> LogNorm:
    positive = np.concatenate([m[np.isfinite(m) & (m > 0)] for m in matrices])
    if positive.size == 0:
        raise ValueError("Spatial maps contain no positive values")
    return LogNorm(vmin=max(1, float(np.percentile(positive, 2))), vmax=float(np.percentile(positive, 99.5)))


def plot_panel_e_spatial_failure(
    container,
    record: SampleRecord,
    revised_spatial: pd.DataFrame,
    legacy_whitelist: pd.DataFrame,
    source_path: Path,
    whitelist_path: Path,
    provenance_rows: list[dict],
) -> tuple[LogNorm, LogNorm]:
    container.axis("off")
    panel_title(container, "E", "Conventional QC may pass despite spatial reconstruction failure")
    sub = container.get_subplotspec().subgridspec(2, 3, height_ratios=[0.4, 1], width_ratios=[1, 1, 0.72], hspace=0.13, wspace=0.16)
    header = container.figure.add_subplot(sub[0, :])
    good_ax = container.figure.add_subplot(sub[1, 0])
    bad_ax = container.figure.add_subplot(sub[1, 1])
    diag_ax = container.figure.add_subplot(sub[1, 2])
    header.axis("off")

    bad = legacy_order_reconstruction(revised_spatial, legacy_whitelist)
    good_matrix = spatial_matrix(revised_spatial, "dna_reads")
    bad_matrix = spatial_matrix(bad, "dna_reads")
    norm = _positive_lognorm(good_matrix, bad_matrix)
    extent = _spatial_extent(revised_spatial)
    for ax, matrix, title in (
        (good_ax, good_matrix, "Correct revised coordinates"),
        (bad_ax, bad_matrix, "Controlled legacy-order reconstruction"),
    ):
        image = ax.imshow(matrix, cmap=SPATIAL_CMAP, norm=norm, interpolation="none", extent=extent, aspect="equal")
        ax.set_title(title, fontsize=7.2, pad=3)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.tick_params(labelsize=5.7)
    colorbar = container.figure.colorbar(image, ax=[good_ax, bad_ax], orientation="horizontal", fraction=0.07, pad=0.12, aspect=25)
    colorbar.set_label("DNA fragments/pixel (shared log scale)", fontsize=6.2)
    colorbar.ax.tick_params(labelsize=5.6)

    total = revised_spatial["dna_reads"].sum()
    median = revised_spatial.loc[revised_spatial["dna_reads"] > 0, "dna_reads"].median()
    covered = 100 * (revised_spatial["dna_reads"] > 0).mean()
    mapping = 100 * revised_spatial["mapped_reads"].sum() / revised_spatial["dna_reads"].sum()
    header.text(0, 0.78, "Similar aggregate QC, incorrect spatial reconstruction", color=WARNING_COLOR, weight="bold", fontsize=7.3)
    header.text(0, 0.39, f"Both maps: {total / 1e6:.1f}M fragments  |  {mapping:.1f}% mapped  |  median {median:.0f}/covered pixel  |  {covered:.1f}% pixels covered")
    header.text(0, 0.03, "Controlled failure preserves every per-pixel count but applies incompatible legacy coordinate ordering.", color="#555555", fontsize=5.8)

    diagnostics = [morans_i(good_matrix), morans_i(bad_matrix), neighbor_spearman(good_matrix), neighbor_spearman(bad_matrix)]
    ypos = [3, 2, 1, 0]
    labels = ["Moran's I\ncorrect", "Moran's I\nincorrect", "Neighbor rho\ncorrect", "Neighbor rho\nincorrect"]
    colors = ["#0072B2", WARNING_COLOR, "#0072B2", WARNING_COLOR]
    diag_ax.barh(ypos, diagnostics, color=colors, height=0.62)
    diag_ax.set_yticks(ypos)
    diag_ax.set_yticklabels(labels)
    diag_ax.set_xlabel("Spatial continuity")
    clean_axis(diag_ax, "x")
    for y0, value in zip(ypos, diagnostics):
        diag_ax.text(value, y0, f" {value:.2f}", va="center", fontsize=6)
    for metric, value, definition in (
        ("correct_morans_i", diagnostics[0], "four_neighbor_global_morans_i_v1"),
        ("incorrect_morans_i", diagnostics[1], "four_neighbor_global_morans_i_v1"),
        ("correct_neighbor_spearman", diagnostics[2], "adjacent_pixel_spearman_v1"),
        ("incorrect_neighbor_spearman", diagnostics[3], "adjacent_pixel_spearman_v1"),
    ):
        provenance_rows.append(derived_row(record, metric, value, "coefficient", f"{source_path}; {whitelist_path}", "recomputed_revised", definition, "E"))
    return norm, norm


def plot_panel_f_scar_artifact(
    container,
    record: SampleRecord,
    spatial: pd.DataFrame,
    source_path: Path,
    provenance_rows: list[dict],
) -> object:
    container.axis("off")
    panel_title(container, "F", "Automatic localization of low-read spatial scar artifacts", y=1.08)
    sub = container.get_subplotspec().subgridspec(1, 2, width_ratios=[1.12, 0.82], wspace=0.22)
    map_ax = container.figure.add_subplot(sub[0, 0])
    profile_ax = container.figure.add_subplot(sub[0, 1])
    matrix = spatial_matrix(spatial, "dna_reads")
    call = detect_scar(matrix)
    norm = _positive_lognorm(matrix)
    extent = _spatial_extent(spatial)
    image = map_ax.imshow(matrix, cmap=SPATIAL_CMAP, norm=norm, interpolation="none", extent=extent, aspect="equal")
    if call.axis == "horizontal":
        map_ax.axhspan(call.start - 0.5, call.end + 0.5, facecolor="none", edgecolor=WARNING_COLOR, lw=1.3)
        coordinates = np.arange(1, len(call.profile) + 1)
        profile_ax.plot(call.profile, coordinates, color="#333333", lw=1)
        profile_ax.axhspan(call.start - 0.5, call.end + 0.5, color=WARNING_COLOR, alpha=0.22)
        profile_ax.invert_yaxis()
        profile_ax.set_xlabel("Median DNA fragments")
        profile_ax.set_ylabel("Row")
    else:
        map_ax.axvspan(call.start - 0.5, call.end + 0.5, facecolor="none", edgecolor=WARNING_COLOR, lw=1.3)
        coordinates = np.arange(1, len(call.profile) + 1)
        profile_ax.plot(coordinates, call.profile, color="#333333", lw=1)
        profile_ax.axvspan(call.start - 0.5, call.end + 0.5, color=WARNING_COLOR, alpha=0.22)
        profile_ax.set_ylabel("Median DNA fragments")
        profile_ax.set_xlabel("Column")
    map_ax.set_xlabel("x")
    map_ax.set_ylabel("y")
    map_ax.set_title(f"{record.sample}: detected {call.axis} band", fontsize=7.2, pad=3)
    clean_axis(profile_ax, "x" if call.axis == "horizontal" else "y")
    interval = f"{call.start}" if call.start == call.end else f"{call.start}-{call.end}"
    profile_ax.text(
        0.98,
        0.98,
        f"Detected {call.axis} interval: {interval}\nscar ratio = {call.ratio:.3f}\ninside/adjacent = {call.inside_median:.0f}/{call.adjacent_median:.0f}",
        transform=profile_ax.transAxes,
        ha="right",
        va="top",
        color=WARNING_COLOR,
        fontsize=6.2,
    )
    colorbar = container.figure.colorbar(image, ax=map_ax, orientation="horizontal", fraction=0.06, pad=0.1, aspect=22)
    colorbar.set_label("DNA fragments/pixel (log scale)", fontsize=6.2)
    colorbar.ax.tick_params(labelsize=5.6)
    for metric, value, unit in (
        ("scar_ratio", call.ratio, "ratio"),
        ("scar_inside_median", call.inside_median, "fragments/pixel"),
        ("scar_adjacent_median", call.adjacent_median, "fragments/pixel"),
        ("scar_start", call.start, "row_or_column"),
        ("scar_end", call.end, "row_or_column"),
    ):
        provenance_rows.append(derived_row(record, metric, value, unit, source_path, "recomputed_revised", "local_median_scar_detector_v1", "F"))
    return call


def plot_panel_a_multiqc(
    ax, metrics: pd.DataFrame, records_by_sample: dict[str, SampleRecord]
) -> None:
    """Render a MultiQC general-statistics style table with in-cell bars."""
    panel_title(ax, "A", "Representative cohort QC overview", y=1.08)
    sample_order = [
        "SpMETSLE14D",
        "SpMETSLE16D",
        "SpMETSLE14DM",
        "SpMETSLB16D",
        "SpMETSLHSK03D",
        "SpMETSLHSK05D",
    ]
    records = [records_by_sample[sample] for sample in sample_order]
    headers = [
        "Raw\npairs",
        "Pass\n%",
        "Map\n%",
        "Dup.\n%",
        "Conv.\n%",
        "CpG\nsites",
        "Fragments\n/pixel",
        "RNA map\n%",
        "RNA\nUMIs",
        "RNA\ngenes",
    ]
    nrows, ncols = len(records), len(MATRIX_METRICS)
    values = np.full((nrows, ncols), np.nan)
    fallback = np.zeros((nrows, ncols), dtype=bool)
    for i, record in enumerate(records):
        for j, metric in enumerate(MATRIX_METRICS):
            values[i, j] = metric_value(metrics, record.sample, metric)
            source = metric_source(metrics, record.sample, metric)
            fallback[i, j] = bool(source is not None and source["is_fallback"])

    fractions = np.zeros_like(values)
    for j, metric in enumerate(MATRIX_METRICS):
        finite = np.isfinite(values[:, j])
        if not finite.any():
            continue
        if METRIC_META.get(metric, ("", "", ""))[1].startswith("%"):
            fractions[finite, j] = np.clip(values[finite, j] / 100, 0, 1)
        else:
            maximum = np.nanmax(values[:, j])
            fractions[finite, j] = (
                np.log10(values[finite, j] + 1) / np.log10(maximum + 1)
                if maximum > 0
                else 0
            )

    ax.set_xlim(-2.65, ncols - 0.45)
    ax.set_ylim(nrows - 0.35, -1.25)
    ax.axis("off")
    ax.add_patch(
        Rectangle(
            (-2.53, -0.43),
            ncols + 2.0,
            1.86,
            facecolor="#EFF6FC",
            edgecolor="#6BAED6",
            linewidth=0.9,
            zorder=-3,
        )
    )
    for j, label in enumerate(headers):
        ax.text(j, -0.55, label, ha="center", va="bottom", fontsize=6.2, weight="bold")
    ax.text(3, -1.03, "DNA", ha="center", va="center", fontsize=6.3, color="#555555")
    ax.plot([-0.45, 6.45], [-0.87, -0.87], color="#BDBDBD", lw=0.7)
    ax.text(8, -1.03, "RNA", ha="center", va="center", fontsize=6.3, color="#555555")
    ax.plot([6.55, 9.45], [-0.87, -0.87], color="#BDBDBD", lw=0.7)

    for i, record in enumerate(records):
        ax.add_patch(
            Rectangle(
                (-2.48, i - 0.34),
                0.12,
                0.68,
                color=context_color(record),
                linewidth=0,
            )
        )
        ax.text(
            -2.28,
            i,
            record.sample,
            ha="left",
            va="center",
            fontsize=6.4,
            weight="bold" if i < 2 else "normal",
        )
        for j, metric in enumerate(MATRIX_METRICS):
            ax.add_patch(
                Rectangle(
                    (j - 0.45, i - 0.32),
                    0.9,
                    0.64,
                    facecolor="white",
                    edgecolor="#D9D9D9",
                    linewidth=0.55,
                    zorder=-2,
                )
            )
            value = values[i, j]
            if np.isfinite(value):
                ax.add_patch(
                    Rectangle(
                        (j - 0.45, i - 0.32),
                        0.9 * fractions[i, j],
                        0.64,
                        facecolor=MULTIQC_COLORS[metric],
                        alpha=0.24,
                        linewidth=0,
                        zorder=-1,
                    )
                )
                text_value = pretty_value(metric, value) + (f"*" if fallback[i, j] else "")
            else:
                text_value = "—"
            ax.text(j, i, text_value, ha="center", va="center", fontsize=5.9, color="#202020")

    ax.text(
        1,
        -0.035,
        "Count-bar lengths use log10 scaling; percentage bars use a 0–100% scale; —, unavailable",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.4,
        color="#666666",
    )


def load_pipeline_comparison(repo_root: Path) -> tuple[pd.DataFrame, list[dict]]:
    path = repo_root / "tools" / "figure2" / "pipeline_comparison.tsv"
    frame = pd.read_csv(path, sep="\t")
    run_root = Path("/scr1/users/suc3/20260528_runNewSpatialPipe")
    e65_multiqc = repo_root / "multiqc_smcseq_finished5_data" / "multiqc_curated_metrics.txt"
    reproduced_root = (
        run_root
        / "final_output"
        / "E65_20um"
        / "original_reproduction_20260610"
    )
    rows = []
    for item in frame.to_dict("records"):
        if item["pipeline"] == "spmetsle14d":
            if item["metric"] in {"Raw read pairs", "Designed", "Designed / raw"}:
                upstream = (
                    run_root
                    / "pipeline_output"
                    / "SpMETSLE14D"
                    / "trim"
                    / "SpMETSLE14D_trim_stats.txt"
                )
            elif item["metric"] in {
                "Mapped before CpH filter",
                "Mapping rate",
                "Duplication",
            }:
                upstream = (
                    run_root
                    / "pipeline_output"
                    / "SpMETSLE14D"
                    / "qc"
                    / "table"
                    / "SpMETSLE14D_spatial_barcode_counts.tsv"
                )
            else:
                upstream = (
                    run_root
                    / "pipeline_output"
                    / "SpMETSLE14D"
                    / "qc"
                    / "biscuit_qc_merged"
                    / "SpMETSLE14D_curated_metrics_mqc.tsv"
                )
            sources = f"{path}; {upstream}"
            definition = "spatialdmt_revised_fragment_summary_v1"
        elif item["pipeline"] == "e65_bsconv":
            sources = f"{path}; {e65_multiqc}"
            definition = "smcseq_bsconv_y005_summary_v1"
        else:
            sources = (
                f"{path}; "
                f"{reproduced_root / 'original_vs_current.tsv'}; "
                f"{reproduced_root / 'bismark_xm_filter.tsv'}"
            )
            definition = "reproduced_fragment_unit_summary_v2"
        rows.append(
            {
                "sample": item["pipeline"],
                "species": "mouse",
                "tissue": "embryo",
                "modality": "DNA",
                "pipeline": item["label"],
                "metric": item["metric"],
                "value": item["value"],
                "unit": item["unit"] if pd.notna(item["unit"]) else "",
                "source_file": sources,
                "source_type": "curated_summary",
                "definition_version": definition,
                "is_fallback": False,
                "panel": "B",
            }
        )
    return frame, rows


def _comparison_value(frame: pd.DataFrame, pipeline: str, metric: str) -> float:
    selected = frame.loc[(frame["pipeline"] == pipeline) & (frame["metric"] == metric), "value"]
    return float(selected.iloc[0]) if not selected.empty and pd.notna(selected.iloc[0]) else np.nan


def _comparison_display(frame: pd.DataFrame, pipeline: str, metric: str) -> str:
    selected = frame.loc[(frame["pipeline"] == pipeline) & (frame["metric"] == metric), "display"]
    return str(selected.iloc[0]) if not selected.empty and pd.notna(selected.iloc[0]) else "—"


def load_panel_b_cohorts(
    repo_root: Path, output_dir: Path
) -> tuple[pd.DataFrame, list[dict]]:
    panel_metrics, audit = load_panel_b_metrics(repo_root)
    audit.to_csv(
        output_dir / "panel_b_unified_pipeline_audit.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    plotted = panel_metrics.loc[panel_metrics["metric"].isin(RETENTION_METRICS)]
    expected = {"E6.5 reproduced": 1, "SmC-seq": 8, "Spatial-DMT": 5}
    observed = plotted.groupby("method")["sample"].nunique().to_dict()
    if observed != expected:
        raise ValueError(f"Unexpected Panel B cohort sizes: {observed}")
    missing = plotted["value"].isna()
    if missing.any():
        raise ValueError("Panel B contains unexpected missing retention values")

    rows = []
    for item in plotted.itertuples(index=False):
        sample = item.internal_sample
        tissue = "brain" if str(sample).startswith("SpMETSLB") else "embryo"
        pipeline = (
            "E6.5 reproduced pipeline"
            if item.method == "E6.5 reproduced"
            else f"unified Snakefile; {item.method} branch"
        )
        rows.append(
            {
                "sample": sample,
                "species": "mouse",
                "tissue": tissue,
                "modality": "DNA",
                "pipeline": pipeline,
                "metric": item.metric,
                "value": float(item.value),
                "unit": item.unit,
                "source_file": str(item.source_file).replace(";", "; "),
                "source_type": item.source_type,
                "definition_version": item.definition,
                "is_fallback": False,
                "panel": "B",
            }
        )
    return panel_metrics, rows


def plot_panel_b_cohort_retention(container, metrics: pd.DataFrame) -> None:
    container.axis("off")
    panel_title(container, "B", "Retention from raw reads through alignment and QC", y=1.08)
    ax = container.inset_axes([0.095, 0.17, 0.89, 0.74])
    plot_retention(ax, metrics, fontsize=6.2)
    container.text(
        0.01,
        0.025,
        "Read units: Spatial-DMT aligns both mates; SmC-seq trim/mapped/QC checkpoints use one genomic R1 per raw pair.\n"
        "Spatial-DMT's trim checkpoint is linker-positive; E6.5 reproduced carries one selected mate per pair.\n"
        "Literal Reads_passed and residual drops are in the audit TSV.",
        transform=container.transAxes,
        ha="left",
        va="bottom",
        fontsize=5.15,
        color="#666666",
    )


def plot_panel_b_pipeline_comparison(container, comparison: pd.DataFrame) -> None:
    container.axis("off")
    panel_title(container, "B", "Processing and conversion-filter comparison", y=1.08)
    sub = container.get_subplotspec().subgridspec(
        1, 2, width_ratios=[0.9, 1.45], wspace=0.23
    )
    retention_ax = container.figure.add_subplot(sub[0, 0])
    table_ax = container.figure.add_subplot(sub[0, 1])

    pipelines = ["spmetsle14d", "e65_bsconv", "e65_reproduced"]
    labels = {
        "spmetsle14d": "14D, revised",
        "e65_bsconv": "E6.5, ours",
        "e65_reproduced": "E6.5, reproduced",
    }
    colors = {
        "spmetsle14d": "#0072B2",
        "e65_bsconv": "#E69F00",
        "e65_reproduced": "#7A7A7A",
    }
    markers = {
        "spmetsle14d": "o",
        "e65_bsconv": "^",
        "e65_reproduced": "D",
    }
    stages = ["Raw", "Designed", "Mapped", "CpH-filter\npassed"]
    for pipeline in pipelines:
        raw = _comparison_value(comparison, pipeline, "Raw read pairs")
        values = np.array(
            [
                raw,
                _comparison_value(comparison, pipeline, "Designed"),
                _comparison_value(comparison, pipeline, "Mapped before CpH filter"),
                _comparison_value(comparison, pipeline, "CpH-filter passed"),
            ],
            dtype=float,
        )
        retention_ax.plot(
            range(len(stages)),
            100 * values / raw,
            color=colors[pipeline],
            marker=markers[pipeline],
            ms=4,
            lw=1.25,
            label=labels[pipeline],
        )
    count_annotations = [
        (1, 75.09, "19.92M", "#0072B2", (2, 5)),
        (2, 70.07, "18.59M", "#0072B2", (2, 5)),
        (1, 72.32, "241.02M", "#E69F00", (3, -11)),
        (2, 47.64, "158.76M", "#E69F00", (3, 5)),
        (1, 76.03, "253.38M", "#7A7A7A", (3, 6)),
        (2, 15.62, "52.05M", "#7A7A7A", (3, -11)),
        (3, 31.17, "103.90M", "#E69F00", (3, 5)),
        (3, 14.44, "48.11M", "#7A7A7A", (3, -11)),
    ]
    for x, y, label, color, offset in count_annotations:
        retention_ax.annotate(
            label,
            (x, y),
            xytext=offset,
            textcoords="offset points",
            fontsize=4.8,
            color=color,
        )
    retention_ax.text(
        0.02,
        0.965,
        "Raw: 14D 26.53M; E6.5 333.28M pairs",
        transform=retention_ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.5,
        color="#555555",
    )
    retention_ax.set_xticks(range(len(stages)))
    retention_ax.set_xticklabels(stages, rotation=32, ha="right")
    retention_ax.set_ylabel("Retained (% of raw pairs)")
    retention_ax.set_ylim(0, 105)
    retention_ax.set_title("Stage-wise retention", fontsize=7.2, pad=4)
    retention_ax.legend(
        loc="lower left",
        ncol=1,
        fontsize=5.6,
        handlelength=1.8,
        labelspacing=0.3,
    )
    clean_axis(retention_ax, "y")

    chemistry_metrics = [
        "Mapping rate",
        "CpH-filter passed / mapped",
        "Duplication",
        "CpG methylation",
        "CpH (A/T/G)",
        "Lambda retention",
    ]
    column_labels = [
        "14D\nrevised",
        "E6.5\nours",
        "E6.5\nreproduced",
    ]
    table_ax.set_xlim(-1.75, 2.52)
    table_ax.set_ylim(len(chemistry_metrics) - 0.4, -1.25)
    table_ax.axis("off")
    table_ax.set_title("Chemistry and filtering diagnostics", fontsize=7.2, pad=4)
    for j, label in enumerate(column_labels):
        table_ax.text(j, -0.64, label, ha="center", va="bottom", fontsize=5.8, weight="bold")
    for i, metric in enumerate(chemistry_metrics):
        row_label = (
            "CpH-filter passed\n/ mapped"
            if metric == "CpH-filter passed / mapped"
            else metric
        )
        table_ax.text(-1.68, i, row_label, ha="left", va="center", fontsize=5.9)
        for j, pipeline in enumerate(pipelines):
            face = "#F4F4F4"
            if pipeline == "e65_bsconv" and metric in {
                "CpH-filter passed / mapped",
                "CpG methylation",
                "CpH (A/T/G)",
            }:
                face = "#FFF2CC"
            table_ax.add_patch(
                Rectangle(
                    (j - 0.47, i - 0.34),
                    0.94,
                    0.68,
                    facecolor=face,
                    edgecolor="#D0D0D0",
                    linewidth=0.55,
                )
            )
            display = _comparison_display(comparison, pipeline, metric)
            if metric == "CpH (A/T/G)" and " / " in display:
                display = display.replace(" / ", " /\n")
            table_ax.text(
                j,
                i,
                display,
                ha="center",
                va="center",
                fontsize=5.15 if "\n" in display else 5.55,
            )
    table_ax.text(
        1,
        -0.035,
        "14D: no bsconv filter. E6.5 ours: complete -y 0.05 run; newer -m 2 chemistry incomplete.\n"
        "Reproduced “designed” = structure-pass; duplication definitions differ; λ retention: higher is better.",
        transform=table_ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.25,
        color="#666666",
    )


def plot_panel_c_sample_maps(
    container,
    records: list[SampleRecord],
    spatial_frames: list[pd.DataFrame],
) -> LogNorm:
    container.axis("off")
    panel_title(container, "C", "Spatial raw-read distributions", y=1.08)
    sub = container.get_subplotspec().subgridspec(1, 2, wspace=0.12)
    axes = [container.figure.add_subplot(sub[0, i]) for i in range(2)]
    matrices = [spatial_matrix(frame, "dna_reads") for frame in spatial_frames]
    norm = _positive_lognorm(*matrices)
    for ax, matrix, record, frame in zip(axes, matrices, records, spatial_frames):
        ax.imshow(
            matrix,
            cmap=SPATIAL_CMAP,
            norm=norm,
            interpolation="none",
            extent=_spatial_extent(frame),
            aspect="equal",
            rasterized=True,
        )
        ax.add_patch(
            Rectangle(
                (0.31, 1.025),
                0.025,
                0.048,
                transform=ax.transAxes,
                facecolor=CONTEXT_COLORS["mouse embryo"],
                edgecolor="none",
                clip_on=False,
            )
        )
        ax.text(
            0.35,
            1.049,
            record.sample,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=7.4,
            weight="bold",
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.tick_params(labelsize=5.7)
    container.text(
        0.99,
        1.01,
        "Heatmap color = raw DNA read count per pixel (shared scale)",
        transform=container.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.8,
        color="#555555",
    )
    return norm


def load_panel_c_spatial_frame(path: Path, whitelist_path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    revised_columns = {"x", "y", "dna_reads"}
    if revised_columns.issubset(frame.columns):
        return frame
    if not {"barcode", "count"}.issubset(frame.columns):
        raise ValueError(f"Unsupported spatial count schema in {path}: {frame.columns.tolist()}")

    whitelist = read_legacy_whitelist(whitelist_path)
    if whitelist["barcode"].duplicated().any():
        raise ValueError(f"Duplicate barcodes in spatial whitelist: {whitelist_path}")
    merged = frame.merge(whitelist, on="barcode", how="left", validate="one_to_one")
    if merged[["x", "y"]].isna().any().any():
        missing = merged.loc[merged["x"].isna() | merged["y"].isna(), "barcode"].tolist()
        raise ValueError(f"Legacy spatial barcodes missing coordinates: {missing[:5]}")
    merged["dna_reads"] = pd.to_numeric(merged["count"], errors="raise")
    return merged


def plot_panel_d_scar_map(container, record: SampleRecord, spatial: pd.DataFrame) -> LogNorm:
    container.axis("off")
    panel_title(container, "D", "Localized low-read spatial scar artifact", y=1.08)
    map_ax = container.inset_axes([0.13, 0.02, 0.74, 0.9])
    matrix = spatial_matrix(spatial, "dna_reads")
    norm = _positive_lognorm(matrix)
    map_ax.imshow(
        matrix,
        cmap=SPATIAL_CMAP,
        norm=norm,
        interpolation="none",
        extent=_spatial_extent(spatial),
        aspect="equal",
        rasterized=True,
    )
    map_ax.set_title(f"{record.sample}: total DNA fragments", fontsize=7.4, pad=4)
    map_ax.add_patch(
        Rectangle(
            (0.5, 18.5),
            50.0,
            1.0,
            fill=False,
            edgecolor=WARNING_COLOR,
            linewidth=1.5,
        )
    )
    map_ax.annotate(
        "low-read scar",
        xy=(47, 19),
        xytext=(39, 12),
        color=WARNING_COLOR,
        fontsize=6.3,
        weight="bold",
        arrowprops={"arrowstyle": "->", "color": WARNING_COLOR, "lw": 0.9},
    )
    map_ax.set_xlabel("x")
    map_ax.set_ylabel("y")
    map_ax.tick_params(labelsize=5.8)
    return norm


def validate_inputs(
    records: list[SampleRecord],
    metrics: pd.DataFrame,
    generated_panels: set[str],
    sample_map_norm: LogNorm,
) -> None:
    if generated_panels != set("ABCD"):
        raise AssertionError(f"Required panels missing: {sorted(set('ABCD') - generated_panels)}")
    samples = {record.sample for record in records}
    forbidden = samples & TECHNICAL_RUNS
    if forbidden:
        raise AssertionError(f"Technical/control samples included: {sorted(forbidden)}")
    if metrics["sample"].isna().any() or metrics["metric"].isna().any():
        raise AssertionError("Metric identifiers contain missing labels")
    if np.isinf(pd.to_numeric(metrics["value"], errors="coerce")).any():
        raise AssertionError("Infinite plotted metric detected")
    duplicate_units = metrics.groupby("metric")["unit"].nunique()
    if (duplicate_units > 1).any():
        bad = duplicate_units[duplicate_units > 1].index.tolist()
        raise AssertionError(f"Inconsistent units for metrics: {bad}")
    if not isinstance(sample_map_norm, LogNorm):
        raise AssertionError("Panel C sample maps do not share logarithmic normalization")
    missing_samples = samples - set(metrics["sample"])
    if missing_samples:
        LOG.warning("Selected samples explicitly missing all metrics: %s", sorted(missing_samples))
    missing_sources = []
    for source in metrics["source_file"].dropna().astype(str).unique():
        for part in source.split("; "):
            if part in {"missing", "revised spatial table"} or "*.bam" in part:
                continue
            if not Path(part).exists():
                missing_sources.append(part)
    if missing_sources:
        LOG.warning("Derived provenance includes non-path labels or unavailable cache inputs: %s", missing_sources[:5])


def make_figure2(
    repo_root: Path,
    output_dir: Path,
    records: list[SampleRecord],
    force: bool,
    skip_expensive: bool,
) -> tuple[plt.Figure, pd.DataFrame]:
    apply_style()
    cache_dir = output_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    for record in records:
        rows, _ = collect_sample_metrics(record, cache_dir, force, skip_expensive, include_yame=False)
        if not rows:
            LOG.warning("No revised or legacy metrics found for %s", record.sample)
        all_rows.extend(rows)
    deep_record = next(record for record in records if record.sample == "SpMETSLE14DM")
    if not any(row["sample"] == deep_record.sample and row["metric"] == "unique_cpgs" for row in all_rows):
        aggregate, cg_path, cg_type, cg_fallback = run_yame_aggregate(
            deep_record, cache_dir, force, skip_expensive
        )
        if aggregate is not None and not aggregate.empty:
            all_rows.append(
                derived_row(
                    deep_record,
                    "unique_cpgs",
                    float(aggregate.iloc[0]["N_query"]),
                    "sites",
                    cg_path,
                    cg_type,
                    "yame_musum_union_N_query_v1",
                    "A",
                    is_fallback=cg_fallback,
                )
            )
        else:
            LOG.warning("Deep-sample aggregate unique CpGs are unavailable")
    metrics = pd.DataFrame(all_rows)
    metrics = apply_revised_metric_snapshot(
        metrics, records, repo_root / "tools" / "figure2" / "revised_shallow_metrics_snapshot.tsv"
    )
    metrics["panel"] = "A"
    records_by_sample = {record.sample: record for record in records}
    panel_b_metrics, panel_b_rows = load_panel_b_cohorts(repo_root, output_dir)

    fig = plt.figure(figsize=(15, 10), constrained_layout=False)
    outer = fig.add_gridspec(2, 2, hspace=0.23, wspace=0.16)
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.055, top=0.925)
    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    ax_c = fig.add_subplot(outer[1, 0])
    ax_d = fig.add_subplot(outer[1, 1])
    generated: set[str] = set()
    derived: list[dict] = []

    plot_panel_a_multiqc(ax_a, metrics, records_by_sample)
    generated.add("A")
    plot_panel_b_cohort_retention(ax_b, panel_b_metrics)
    derived.extend(panel_b_rows)
    generated.add("B")

    c_records = [
        records_by_sample["SpMETSLE14D"],
        records_by_sample["SpMETSLE16D"],
    ]
    c_sources = [find_spatial_barcode_counts(record) for record in c_records]
    c_paths = [source[0] for source in c_sources]
    if any(path is None for path in c_paths):
        raise FileNotFoundError("Panel C requires revised 14D and 16D spatial counts")
    whitelist_path = repo_root / "barcodes" / "spatial_barcodes.txt"
    c_frames = [load_panel_c_spatial_frame(path, whitelist_path) for path in c_paths]
    sample_map_norm = plot_panel_c_sample_maps(ax_c, c_records, c_frames)
    for record, source, frame in zip(c_records, c_sources, c_frames):
        path, source_type, is_fallback = source
        source_file = f"{path}; {whitelist_path}" if is_fallback else path
        definition = (
            "legacy_barcode_count_join_spatial_whitelist_v1"
            if is_fallback
            else "revised_spatial_raw_dna_reads_v1"
        )
        derived.append(
            derived_row(
                record,
                "raw_read_spatial_map_total",
                float(frame["dna_reads"].sum()),
                "fragments",
                source_file,
                source_type,
                definition,
                "C",
                is_fallback=is_fallback,
            )
        )
    generated.add("C")

    d_record = SampleRecord(
        "SpMETSLHSK09D", "human", "skin", "DNA",
        Path("/scr1/users/suc3/20260624_runNewSpatialPipe_hg38"), Path("/scr1/users/suc3/tmp")
    )
    d_path, _, _ = find_spatial_barcode_counts(d_record)
    if d_path is None:
        raise FileNotFoundError("Panel D revised SpMETSLHSK09D spatial counts not found")
    d_spatial = pd.read_csv(d_path, sep="\t")
    plot_panel_d_scar_map(ax_d, d_record, d_spatial)
    derived.append(
        derived_row(
            d_record,
            "scar_map_total_fragments",
            float(d_spatial["dna_reads"].sum()),
            "fragments",
            d_path,
            "recomputed_revised",
            "revised_spatial_dna_fragments_v1",
            "D",
        )
    )
    generated.add("D")

    metrics = pd.concat([metrics, pd.DataFrame(derived)], ignore_index=True)
    validate_inputs(records, metrics, generated, sample_map_norm)
    fig.suptitle(
        "Spatially aware quality-control outputs across samples, cohorts, modalities, and failure modes",
        x=0.5,
        y=0.985,
        fontsize=11,
        weight="bold",
    )
    return fig, metrics


def save_outputs(
    fig: plt.Figure,
    metrics: pd.DataFrame,
    records: list[SampleRecord],
    output_dir: Path,
    formats: list[str],
    dpi: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics = metrics.sort_values(["panel", "sample", "metric", "pipeline"], kind="stable")
    metrics.to_csv(output_dir / "figure2_metrics.tsv", sep="\t", index=False, na_rep="NA")
    provenance_columns = [
        "panel", "sample", "pipeline", "metric", "unit", "source_file", "source_type",
        "definition_version", "is_fallback",
    ]
    metrics[provenance_columns].drop_duplicates().to_csv(
        output_dir / "figure2_provenance.tsv", sep="\t", index=False
    )
    selected = records_frame(records)
    selected["figure_role"] = "representative cohort"
    scar_record = {
        "sample": "SpMETSLHSK09D",
        "species": "human",
        "tissue": "skin",
        "modality": "DNA",
        "revised_root": "/scr1/users/suc3/20260624_runNewSpatialPipe_hg38",
        "legacy_root": "/scr1/users/suc3/tmp",
        "figure_role": "scar artifact",
    }
    panel_b_selected = metrics.loc[
        metrics["panel"].eq("B"),
        ["sample", "species", "tissue", "modality", "pipeline"],
    ].drop_duplicates()
    panel_b_selected["revised_root"] = np.where(
        panel_b_selected["pipeline"].str.startswith("revised Spatial-DMT"),
        "/scr1/users/suc3/20260528_runNewSpatialPipe",
        "",
    )
    panel_b_selected["legacy_root"] = ""
    panel_b_selected["figure_role"] = "Panel B retention cohort: " + panel_b_selected.pop(
        "pipeline"
    )
    selected = pd.concat(
        [selected, pd.DataFrame([scar_record]), panel_b_selected], ignore_index=True
    )
    selected.to_csv(output_dir / "figure2_selected_samples.tsv", sep="\t", index=False)

    for fmt in formats:
        fmt = fmt.lower().lstrip(".")
        path = output_dir / f"Figure2.{fmt}"
        kwargs = {"bbox_inches": "tight"}
        if fmt == "png":
            kwargs["dpi"] = max(300, dpi)
        elif fmt in {"tif", "tiff"}:
            kwargs["dpi"] = max(600, dpi)
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        else:
            kwargs["dpi"] = dpi
        fig.savefig(path, **kwargs)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--sample-config", type=Path)
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--skip-expensive", action="store_true")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--format", default="pdf,svg,png,tiff", help="Comma-separated output formats")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = (args.output_dir or repo_root / "results" / "figure2").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    records = load_sample_config(args.sample_config)
    fig, metrics = make_figure2(
        repo_root, output_dir, records, args.force_recompute, args.skip_expensive
    )
    save_outputs(fig, metrics, records, output_dir, args.format.split(","), args.dpi)
    plt.close(fig)

    fallbacks = metrics.loc[metrics["is_fallback"].astype(bool), ["sample", "metric"]].drop_duplicates()
    print(f"Panel A: {len(records)} samples, {len(MATRIX_METRICS)} metrics")
    print(
        "Panel B: unified Spatial-DMT and SmC-seq cohorts at shared alignment/QC "
        "checkpoints + E6.5 reproduced"
    )
    print("Panel C: 14D and 16D raw-read maps with shared hidden normalization")
    print("Panel D: SpMETSLHSK09D low-read scar map with manual annotation")
    print(f"Legacy fallbacks recorded: {len(fallbacks)}")


if __name__ == "__main__":
    main()
