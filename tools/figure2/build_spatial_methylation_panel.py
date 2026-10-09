#!/usr/bin/env python3
"""Build a shared-scale 2D per-pixel CpG-methylation panel.

The retained QC reports contain aggregate methylation metrics but no spatial
methylation heatmap. This builder reads the indexed per-pixel YAME ``.cg``
files, masks pixels with too little CpG coverage, and writes one comparable
panel for the current Spatial-DMT and SmC-seq cohorts.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

TOOLS_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_ROOT))
from plot_spatial_methylation import summarize_cells, summarize_cells_bulk  # noqa: E402
from figure2.plot_style import apply_style  # noqa: E402


RUN_ROOT = Path("/scr1/users/suc3/20260528_runNewSpatialPipe")
SAMPLES = (
    ("Spatial-DMT", "SpMETSLE14DM"),
    ("Spatial-DMT", "SpMETSLE17DM"),
    ("Spatial-DMT", "SpMETSLE21DM"),
    ("Spatial-DMT", "SpMETSLE23DM"),
    ("Spatial-DMT", "SpMETSLB14DM"),
    ("SmC-seq", "E55_10um_rep1"),
    ("SmC-seq", "E55_10um_rep2"),
    ("SmC-seq", "E65_10um_run1"),
    ("SmC-seq", "E65_10um_run2"),
    ("SmC-seq", "E75_10um_rep1"),
    ("SmC-seq", "E75_10um_rep2"),
    ("SmC-seq", "E85_20um"),
    ("SmC-seq", "E85_placenta_10um"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--run-root", type=Path, default=RUN_ROOT)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--output-name", default="Figure2_spatial_methylation_v1")
    parser.add_argument("--min-covered-cpgs", type=int, default=100)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--yame", default="yame")
    parser.add_argument("--reuse-tsv", action="store_true")
    parser.add_argument(
        "--per-cell-subset",
        action="store_true",
        help="Use the legacy one-subset-process-per-pixel extraction",
    )
    parser.add_argument("--dpi", type=int, default=320)
    return parser.parse_args()


def write_rows(path: Path, sample: str, rows: list[dict], threshold: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "spatial_id",
        "x",
        "y",
        "covered_cpgs",
        "cpg_methylation",
        "depth",
        "passes_threshold",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "sample": sample,
                    "spatial_id": row["spatial_id"],
                    "x": row["x"],
                    "y": row["y"],
                    "covered_cpgs": row["covered_cpgs"],
                    "cpg_methylation": (
                        f'{row["beta"]:.6f}' if np.isfinite(row["beta"]) else "NA"
                    ),
                    "depth": f'{row["depth"]:.6f}' if np.isfinite(row["depth"]) else "NA",
                    "passes_threshold": int(row["covered_cpgs"] >= threshold),
                }
            )


def load_or_extract(
    sample: str,
    cg_path: Path,
    table_path: Path,
    yame: str,
    threads: int,
    threshold: int,
    reuse: bool,
    per_cell_subset: bool,
) -> pd.DataFrame:
    if reuse and table_path.is_file():
        frame = pd.read_csv(table_path, sep="\t")
        logging.info("Reusing %s", table_path)
        return frame
    if not cg_path.is_file() or not Path(f"{cg_path}.idx").is_file():
        raise FileNotFoundError(f"Missing indexed .cg source for {sample}: {cg_path}")
    logging.info("Summarizing per-pixel methylation: %s", sample)
    rows = (
        summarize_cells(str(cg_path), yame, threads)
        if per_cell_subset
        else summarize_cells_bulk(str(cg_path), yame)
    )
    write_rows(table_path, sample, rows, threshold)
    return pd.DataFrame(
        [
            {
                "sample": sample,
                "spatial_id": row["spatial_id"],
                "x": row["x"],
                "y": row["y"],
                "covered_cpgs": row["covered_cpgs"],
                "cpg_methylation": row["beta"],
                "depth": row["depth"],
                "passes_threshold": int(row["covered_cpgs"] >= threshold),
            }
            for row in rows
        ]
    )


def matrix_for(frame: pd.DataFrame, threshold: int) -> tuple[np.ndarray, int, int, int, float]:
    x = pd.to_numeric(frame["x"], errors="raise").astype(int)
    y = pd.to_numeric(frame["y"], errors="raise").astype(int)
    covered = pd.to_numeric(frame["covered_cpgs"], errors="coerce")
    beta = pd.to_numeric(frame["cpg_methylation"], errors="coerce")
    grid_x = int(x.max())
    grid_y = int(y.max())
    matrix = np.full((grid_y, grid_x), np.nan, dtype=float)
    keep = covered.ge(threshold) & beta.notna() & beta.between(0, 1)
    matrix[y[keep] - 1, x[keep] - 1] = beta[keep]
    values = beta[keep].to_numpy(dtype=float)
    return matrix, grid_x, grid_y, int(keep.sum()), float(np.median(values)) if values.size else np.nan


def build_figure(
    frames: list[tuple[str, str, pd.DataFrame]],
    threshold: int,
    output_stem: Path,
    dpi: int,
) -> pd.DataFrame:
    apply_style()
    ncols = 4
    nrows = int(np.ceil(len(frames) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(12.0, 3.25 * nrows),
        squeeze=False,
        constrained_layout=False,
    )
    fig.subplots_adjust(left=0.045, right=0.90, bottom=0.055, top=0.925, wspace=0.18, hspace=0.42)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#EEEEEE")
    norm = Normalize(vmin=0.0, vmax=1.0)
    provenance: list[dict[str, object]] = []
    for index, (method, sample, frame) in enumerate(frames):
        ax = axes.flat[index]
        matrix, grid_x, grid_y, n_kept, median = matrix_for(frame, threshold)
        # Use one vector quadrilateral per pixel so the SVG/PDF remain
        # editable; imshow would embed each heatmap as a raster image.
        x_edges = np.arange(0.5, grid_x + 1.5)
        y_edges = np.arange(0.5, grid_y + 1.5)
        image = ax.pcolormesh(
            x_edges,
            y_edges,
            np.ma.masked_invalid(matrix),
            cmap=cmap,
            norm=norm,
            shading="flat",
            rasterized=False,
        )
        ax.set_title(
            f"{sample}\n{method}  |  median β={median:.3f}  |  n={n_kept:,}",
            loc="left",
            fontsize=7.2,
            fontweight="bold",
            linespacing=1.15,
            pad=3,
        )
        ax.set_xlim(0.5, grid_x + 0.5)
        ax.set_ylim(grid_y + 0.5, 0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.55)
            spine.set_color("#333333")
        provenance.append(
            {
                "sample": sample,
                "method": method,
                "grid_x": grid_x,
                "grid_y": grid_y,
                "min_covered_cpgs": threshold,
                "pixels_indexed": len(frame),
                "pixels_plotted": n_kept,
                "median_cpg_methylation": median,
            }
        )
    for ax in axes.flat[len(frames) :]:
        ax.axis("off")
    fig.suptitle(
        "Spatial CpG methylation level across current Spatial-DMT and SmC-seq samples",
        x=0.045,
        y=0.968,
        ha="left",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.045,
        0.018,
        f"Per-pixel YAME β; pixels with < {threshold} covered CpGs are gray. Shared color scale: 0–1.",
        fontsize=6.6,
        color="#555555",
    )
    # Matplotlib's default colorbar is an embedded gradient image in SVG.
    # Draw the gradient as vector quads so the complete deliverable stays
    # editable, including its legend.
    cbar_ax = fig.add_axes([0.925, 0.16, 0.022, 0.64])
    gradient = np.linspace(0.0, 1.0, 256)[:, None]
    cbar_ax.pcolormesh(
        [0.0, 1.0],
        np.linspace(0.0, 1.0, 257),
        gradient,
        cmap=cmap,
        norm=norm,
        shading="flat",
        rasterized=False,
    )
    cbar_ax.set_ylim(0.0, 1.0)
    cbar_ax.set_xlim(0.0, 1.0)
    cbar_ax.set_xticks([])
    cbar_ax.set_yticks(np.linspace(0.0, 1.0, 6))
    cbar_ax.yaxis.tick_right()
    cbar_ax.tick_params(labelsize=6.5)
    cbar_ax.set_ylabel("CpG methylation level (β)", fontsize=8, labelpad=8)
    for fmt in ("png", "svg", "pdf"):
        kwargs = {"bbox_inches": "tight"}
        if fmt == "png":
            kwargs["dpi"] = dpi
        fig.savefig(output_stem.with_suffix(f".{fmt}"), **kwargs)
    plt.close(fig)
    return pd.DataFrame(provenance)


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    if args.min_covered_cpgs < 0 or args.threads < 1:
        raise ValueError("--min-covered-cpgs must be nonnegative and --threads >= 1")
    repo_root = args.repo_root.resolve()
    run_root = args.run_root.resolve()
    output_dir = (args.output_dir or repo_root / "results" / "figure2").resolve()
    table_dir = output_dir / "spatial_methylation_tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    frames: list[tuple[str, str, pd.DataFrame]] = []
    for method, sample in SAMPLES:
        cg_path = run_root / "final_output" / sample / f"{sample}.cg"
        table_path = table_dir / f"{sample}_spatial_cpg_methylation.tsv"
        frame = load_or_extract(
            sample,
            cg_path,
            table_path,
            args.yame,
            args.threads,
            args.min_covered_cpgs,
            args.reuse_tsv,
            args.per_cell_subset,
        )
        frames.append((method, sample, frame))
    output_stem = output_dir / args.output_name
    provenance = build_figure(frames, args.min_covered_cpgs, output_stem, args.dpi)
    provenance.insert(0, "source_type", "indexed YAME per-pixel .cg summary")
    provenance.insert(1, "qc_report_spatial_methylation_map", False)
    provenance.insert(2, "source_root", str(run_root))
    provenance.to_csv(f"{output_stem}_provenance.tsv", sep="\t", index=False, na_rep="NA")
    print(f"Wrote {output_stem}.png/.svg/.pdf")
    print(f"Wrote {output_stem}_provenance.tsv")


if __name__ == "__main__":
    main()
