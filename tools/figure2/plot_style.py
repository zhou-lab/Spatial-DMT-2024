"""Shared plotting style for Figure 2."""

from __future__ import annotations

import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap


CONTEXT_COLORS = {
    "mouse embryo": "#0072B2",
    "mouse brain": "#009E73",
    "human skin": "#D55E00",
}
PIPELINE_COLORS = {"revised": "#222222", "legacy": "#8A8A8A"}
WARNING_COLOR = "#B2182B"
MISSING_COLOR = "#E6E6E6"
SPATIAL_CMAP = mpl.colormaps["viridis"].copy()
SPATIAL_CMAP.set_bad("#F1F1F1")
MATRIX_CMAP = LinearSegmentedColormap.from_list(
    "qc_matrix", ["#F4F4F4", "#9ECAE1", "#2166AC"]
)


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.2,
            "axes.titlesize": 8.5,
            "axes.labelsize": 7.2,
            "axes.linewidth": 0.65,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "legend.fontsize": 6.3,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def clean_axis(ax, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color="#D9D9D9", linewidth=0.45, alpha=0.7)
        ax.set_axisbelow(True)


def panel_title(ax, letter: str, title: str, y: float = 1.04) -> None:
    ax.set_title(rf"$\bf{{{letter}}}$  {title}", loc="left", y=y, pad=2)
