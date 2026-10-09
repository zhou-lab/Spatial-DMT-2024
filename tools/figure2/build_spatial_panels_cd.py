#!/usr/bin/env python3
"""Build standalone Figure 2 spatial panels C and D."""

from __future__ import annotations

import argparse
import base64
import io
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from PIL import Image

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from build_figure2 import load_panel_c_spatial_frame  # noqa: E402
from load_metrics import SampleRecord, find_spatial_barcode_counts  # noqa: E402
from plot_style import (  # noqa: E402
    CONTEXT_COLORS,
    SPATIAL_CMAP,
    WARNING_COLOR,
    apply_style,
)
from spatial_utils import spatial_matrix  # noqa: E402


SPATIALDMT_SCAR_REPORT = Path(
    "/scr1/users/suc3/20260624_runNewSpatialPipe_hg38/"
    "final_output/SpMETSLHSK09D/qc_report.html"
)
SMCSEQ_SCAR_REPORT = Path(
    "/scr1/users/suc3/20260528_runNewSpatialPipe/"
    "final_output/E65_10um_run1/qc_report.html"
)
EMBEDDED_PLOT_CAPTION = "DNA reads per barcode (log10)"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--dpi", type=int, default=320)
    parser.add_argument("--tiff-dpi", type=int, default=600)
    return parser.parse_args()


def extract_embedded_plot(
    report: Path,
    caption: str,
    cache_path: Path,
) -> np.ndarray:
    if not report.is_file():
        raise FileNotFoundError(report)
    html = report.read_text()
    matches = re.findall(
        r'<img src="data:image/png;base64,([^"]+)"[^>]*>'
        r"<figcaption>([^<]+)</figcaption>",
        html,
        flags=re.DOTALL,
    )
    encoded = [payload for payload, label in matches if label == caption]
    if len(encoded) != 1:
        raise ValueError(
            f"Expected one embedded plot {caption!r} in {report}, found {len(encoded)}"
        )
    payload = base64.b64decode(encoded[0])
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_bytes(payload)
    return np.asarray(Image.open(io.BytesIO(payload)).convert("RGB"))


def crop_pipeline_map(image: np.ndarray) -> np.ndarray:
    """Crop the square map area from barcode_summary.py's saved 7-inch plot."""
    height, width = image.shape[:2]
    top = int(round(0.048 * height))
    left = int(round(0.081 * width))
    side = int(round(0.872 * height))
    if top + side > height or left + side > width:
        raise ValueError(f"Unexpected embedded plot dimensions: {image.shape}")
    return image[top : top + side, left : left + side]


def add_interval_boxes(
    ax: plt.Axes,
    grid_x: int,
    grid_y: int,
    intervals: tuple[tuple[str, int, int], ...],
) -> None:
    for axis, start, end in intervals:
        if axis == "horizontal":
            patch = Rectangle(
                (0.5, start - 0.5),
                grid_x,
                end - start + 1,
                fill=False,
                edgecolor=WARNING_COLOR,
                linewidth=1.35,
            )
        elif axis == "vertical":
            patch = Rectangle(
                (start - 0.5, 0.5),
                end - start + 1,
                grid_y,
                fill=False,
                edgecolor=WARNING_COLOR,
                linewidth=1.35,
            )
        else:
            raise ValueError(f"Unsupported interval axis: {axis}")
        ax.add_patch(patch)


def plot_embedded_scar_map(
    ax: plt.Axes,
    image: np.ndarray,
    grid_x: int,
    grid_y: int,
    sample: str,
    method_color: str,
    intervals: tuple[tuple[str, int, int], ...],
) -> None:
    ax.imshow(
        crop_pipeline_map(image),
        extent=(0.5, grid_x + 0.5, grid_y + 0.5, 0.5),
        interpolation="none",
        aspect="equal",
        rasterized=True,
    )
    add_interval_boxes(ax, grid_x, grid_y, intervals)
    ax.add_patch(
        Rectangle(
            (0.0, 1.012),
            0.026,
            0.038,
            transform=ax.transAxes,
            facecolor=method_color,
            edgecolor="none",
            clip_on=False,
        )
    )
    ax.text(
        0.04,
        1.031,
        sample,
        transform=ax.transAxes,
        va="center",
        fontsize=8.0,
        weight="bold",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.tick_params(labelsize=6.3)


def positive_shared_norm(*matrices: np.ndarray) -> LogNorm:
    positive = np.concatenate(
        [matrix[np.isfinite(matrix) & (matrix > 0)] for matrix in matrices]
    )
    if not positive.size:
        raise ValueError("Spatial maps contain no positive values")
    return LogNorm(
        vmin=max(1.0, float(np.percentile(positive, 2))),
        vmax=float(np.percentile(positive, 99.5)),
    )


def load_reconstruction_maps(
    repo_root: Path,
) -> tuple[list[SampleRecord], list[pd.DataFrame], list[dict[str, object]]]:
    records = [
        SampleRecord(
            "SpMETSLE14D",
            "mouse",
            "embryo",
            "DNA",
            Path("/scr1/users/suc3/20260528_runNewSpatialPipe"),
            Path("/scr1/users/suc3/tmp"),
        ),
        SampleRecord(
            "SpMETSLE16D",
            "mouse",
            "embryo",
            "DNA+RNA",
            Path("/scr1/users/suc3/20260528_runNewSpatialPipe"),
            Path("/scr1/users/suc3/tmp"),
        ),
    ]
    whitelist = repo_root / "barcodes" / "spatial_barcodes.txt"
    frames: list[pd.DataFrame] = []
    provenance: list[dict[str, object]] = []
    statuses = {
        "SpMETSLE14D": "spatially sound",
        "SpMETSLE16D": "incorrect spatial reconstruction",
    }
    for record in records:
        path, source_type, is_fallback = find_spatial_barcode_counts(record)
        if path is None:
            raise FileNotFoundError(f"Spatial counts not found for {record.sample}")
        frame = load_panel_c_spatial_frame(path, whitelist)
        frames.append(frame)
        provenance.append(
            {
                "panel": "D",
                "sample": record.sample,
                "method": "Spatial-DMT",
                "map_role": statuses[record.sample],
                "metric": "raw DNA reads per spatial barcode",
                "source_file": str(path),
                "source_type": source_type,
                "is_fallback": is_fallback,
                "embedded_plot": "",
                "highlighted_intervals": "",
                "coordinate_transform": (
                    f"barcode joined to {whitelist}"
                    if is_fallback
                    else "pipeline x/y coordinates"
                )
                + (
                    "; display rotated 90 degrees counterclockwise: "
                    "display_x=source_y, display_y=51-source_x"
                    if record.sample == "SpMETSLE16D"
                    else "; no display rotation"
                ),
            }
        )
    return records, frames, provenance


def plot_reconstruction_map(
    ax: plt.Axes,
    record: SampleRecord,
    frame: pd.DataFrame,
    norm: LogNorm,
    warning: bool,
    rotation_k: int = 0,
) -> None:
    matrix = spatial_matrix(frame, "dna_reads")
    if rotation_k:
        matrix = np.rot90(matrix, k=rotation_k)
    extent = (
        float(frame["x"].min()) - 0.5,
        float(frame["x"].max()) + 0.5,
        float(frame["y"].max()) + 0.5,
        float(frame["y"].min()) - 0.5,
    )
    ax.imshow(
        matrix,
        cmap=SPATIAL_CMAP,
        norm=norm,
        interpolation="none",
        extent=extent,
        aspect="equal",
        rasterized=True,
    )
    ax.add_patch(
        Rectangle(
            (0.0, 1.012),
            0.026,
            0.038,
            transform=ax.transAxes,
            facecolor=CONTEXT_COLORS["mouse embryo"],
            edgecolor="none",
            clip_on=False,
        )
    )
    ax.text(
        0.04,
        1.031,
        record.sample,
        transform=ax.transAxes,
        va="center",
        fontsize=8.0,
        weight="bold",
    )
    if warning:
        for spine in ax.spines.values():
            spine.set_color(WARNING_COLOR)
            spine.set_linewidth(1.0)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.tick_params(labelsize=6.3)


def make_figure(
    repo_root: Path,
    output_dir: Path,
) -> tuple[plt.Figure, pd.DataFrame]:
    apply_style()
    cache_dir = output_dir / "cache" / "spatial_panels_cd"
    spatialdmt_image = extract_embedded_plot(
        SPATIALDMT_SCAR_REPORT,
        EMBEDDED_PLOT_CAPTION,
        cache_dir / "SpMETSLHSK09D_spatial_log_dna_reads.png",
    )
    smcseq_image = extract_embedded_plot(
        SMCSEQ_SCAR_REPORT,
        EMBEDDED_PLOT_CAPTION,
        cache_dir / "E65_10um_run1_spatial_log_dna_reads.png",
    )

    scar_specs = (
        {
            "sample": "SpMETSLHSK09D",
            "method": "Spatial-DMT",
            "color": "#0072B2",
            "image": spatialdmt_image,
            "grid": (50, 50),
            "intervals": (
                ("horizontal", 19, 19),
                ("horizontal", 39, 39),
                ("vertical", 17, 17),
            ),
            "report": SPATIALDMT_SCAR_REPORT,
        },
        {
            "sample": "E65_10um_run1",
            "method": "SmC-seq",
            "color": "#D55E00",
            "image": smcseq_image,
            "grid": (96, 96),
            "intervals": (
                ("horizontal", 24, 24),
                ("horizontal", 58, 60),
                ("horizontal", 69, 69),
            ),
            "report": SMCSEQ_SCAR_REPORT,
        },
    )

    records, frames, provenance = load_reconstruction_maps(repo_root)
    matrices = [spatial_matrix(frame, "dna_reads") for frame in frames]
    shared_norm = positive_shared_norm(*matrices)

    fig = plt.figure(figsize=(11.5, 8.2), constrained_layout=False)
    fig.subplots_adjust(left=0.055, right=0.98, bottom=0.06, top=0.915)
    outer = fig.add_gridspec(2, 1, hspace=0.13)

    c_grid = outer[0].subgridspec(
        2, 2, height_ratios=(0.10, 0.90), hspace=0.06, wspace=0.11
    )
    c_header = fig.add_subplot(c_grid[0, :])
    c_header.axis("off")
    c_header.text(
        0.0,
        0.52,
        r"$\bf{C}$  Low-read scar artifacts in both workflow branches",
        transform=c_header.transAxes,
        ha="left",
        va="center",
        fontsize=8.5,
    )
    c_axes = [fig.add_subplot(c_grid[1, index]) for index in range(2)]
    for ax, spec in zip(c_axes, scar_specs):
        grid_x, grid_y = spec["grid"]
        plot_embedded_scar_map(
            ax,
            spec["image"],
            grid_x,
            grid_y,
            spec["sample"],
            spec["color"],
            spec["intervals"],
        )
        provenance.append(
            {
                "panel": "C",
                "sample": spec["sample"],
                "method": spec["method"],
                "map_role": "observed low-read scar example",
                "metric": "log10(DNA reads per spatial barcode + 1)",
                "source_file": str(spec["report"]),
                "source_type": "revised_embedded_qc_plot",
                "is_fallback": False,
                "embedded_plot": EMBEDDED_PLOT_CAPTION,
                "highlighted_intervals": "; ".join(
                    f"{axis}:{start}-{end}"
                    for axis, start, end in spec["intervals"]
                ),
                "coordinate_transform": "none; pipeline x/y orientation retained",
            }
        )
    d_grid = outer[1].subgridspec(
        2, 2, height_ratios=(0.10, 0.90), hspace=0.06, wspace=0.11
    )
    d_header = fig.add_subplot(d_grid[0, :])
    d_header.axis("off")
    d_header.text(
        0.0,
        0.52,
        r"$\bf{D}$  Spatial reconstruction outcome",
        transform=d_header.transAxes,
        ha="left",
        va="center",
        fontsize=8.5,
    )
    d_axes = [fig.add_subplot(d_grid[1, index]) for index in range(2)]
    for ax, record, frame, warning, rotation_k in zip(
        d_axes, records, frames, (False, True), (0, 1)
    ):
        plot_reconstruction_map(
            ax, record, frame, shared_norm, warning, rotation_k=rotation_k
        )
    fig.suptitle(
        "Spatial QC examples across scar artifacts and reconstruction outcomes",
        fontsize=13,
        weight="bold",
        y=0.985,
    )
    return fig, pd.DataFrame(provenance)


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else repo_root / "results" / "figure2"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, provenance = make_figure(repo_root, output_dir)
    stem = output_dir / "Figure2_spatial_panels_CD"
    fig.savefig(stem.with_suffix(".png"), dpi=args.dpi, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=args.tiff_dpi,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    provenance_path = output_dir / "figure2_spatial_panels_CD_provenance.tsv"
    provenance.to_csv(provenance_path, sep="\t", index=False, na_rep="NA")
    print("Panel C: Spatial-DMT and SmC-seq scar maps with selected bands highlighted")
    print("Panel D: SpMETSLE14D sound and SpMETSLE16D incorrect spatial reconstruction")
    print(f"Wrote: {stem}.png/.pdf/.svg/.tiff")
    print(f"Provenance: {provenance_path}")


if __name__ == "__main__":
    main()
