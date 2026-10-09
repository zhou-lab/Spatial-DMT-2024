"""Spatial-table transformations and diagnostics used by Figure 2."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


@dataclass(frozen=True)
class ScarCall:
    axis: str
    start: int
    end: int
    ratio: float
    inside_median: float
    adjacent_median: float
    profile: np.ndarray


@dataclass(frozen=True)
class PairedScarCall:
    """A single depleted line scored against its two immediate neighbors."""

    axis: str
    index: int
    median_count_ratio: float
    depleted_fraction: float
    support_pixels: int
    center_occupancy: float


def score_paired_scars(
    matrix: np.ndarray,
    depletion_threshold: float = 0.55,
    max_median_ratio: float = 0.55,
    min_depleted_fraction: float = 0.60,
    min_support_fraction: float = 0.50,
    min_center_occupancy: float = 0.50,
    pseudocount: float = 1.0,
) -> list[PairedScarCall]:
    """Rank one-row/column scars using pixel-matched immediate neighbors.

    Each center pixel is divided by the geometric mean of the pixels directly
    across the candidate line on its two sides. A candidate must have adequate
    populated flanking support and a non-trivial center occupancy. The latter
    rule excludes completely blank structural rows/columns rather than calling
    them sequencing scars.
    """
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("Scar detection requires a two-dimensional matrix")
    calls: list[PairedScarCall] = []
    for axis, oriented in (("horizontal", matrix), ("vertical", matrix.T)):
        n_lines, line_length = oriented.shape
        required = max(5, int(np.ceil(line_length * min_support_fraction)))
        for center0 in range(1, n_lines - 1):
            before = oriented[center0 - 1]
            center = oriented[center0]
            after = oriented[center0 + 1]
            flank_support = (
                np.isfinite(before)
                & np.isfinite(after)
                & (before > 0)
                & (after > 0)
            )
            if not flank_support.any():
                continue
            center_present = np.isfinite(center) & (center > 0)
            center_occupancy = float(center_present[flank_support].mean())
            if center_occupancy < min_center_occupancy:
                continue
            valid = flank_support & np.isfinite(center)
            support_pixels = int(valid.sum())
            if support_pixels < required:
                continue
            flank_reference = np.sqrt(before[valid] * after[valid])
            ratios = (center[valid] + pseudocount) / (
                flank_reference + pseudocount
            )
            median_ratio = float(np.median(ratios))
            depleted_fraction = float(np.mean(ratios < depletion_threshold))
            if (
                median_ratio <= max_median_ratio
                and depleted_fraction >= min_depleted_fraction
            ):
                calls.append(
                    PairedScarCall(
                        axis=axis,
                        index=center0 + 1,
                        median_count_ratio=median_ratio,
                        depleted_fraction=depleted_fraction,
                        support_pixels=support_pixels,
                        center_occupancy=center_occupancy,
                    )
                )
    return sorted(
        calls,
        key=lambda call: (
            call.median_count_ratio,
            -call.depleted_fraction,
            -call.support_pixels,
        ),
    )


def detect_paired_scar(matrix: np.ndarray, **kwargs) -> PairedScarCall:
    """Return the strongest supported single-line paired scar call."""
    calls = score_paired_scars(matrix, **kwargs)
    if not calls:
        raise ValueError("No supported single-line scar candidate was detected")
    return calls[0]


def spatial_matrix(df: pd.DataFrame, value: str = "dna_reads") -> np.ndarray:
    """Return a Y-by-X matrix in increasing pipeline coordinate order."""
    return (
        df.pivot(index="y", columns="x", values=value)
        .sort_index()
        .sort_index(axis=1)
        .to_numpy(dtype=float)
    )


def morans_i(matrix: np.ndarray) -> float:
    """Four-neighbor global Moran's I, ignoring missing pixels."""
    z = np.asarray(matrix, dtype=float)
    valid = np.isfinite(z)
    if valid.sum() < 3:
        return float("nan")
    mean = np.nanmean(z)
    centered = z - mean
    numerator = 0.0
    weight = 0
    for a, b, va, vb in (
        (centered[:, :-1], centered[:, 1:], valid[:, :-1], valid[:, 1:]),
        (centered[:-1, :], centered[1:, :], valid[:-1, :], valid[1:, :]),
    ):
        keep = va & vb
        numerator += float(np.sum(a[keep] * b[keep]) * 2.0)
        weight += int(keep.sum() * 2)
    denominator = float(np.nansum(centered[valid] ** 2))
    if not weight or not denominator:
        return float("nan")
    return float(valid.sum() / weight * numerator / denominator)


def legacy_order_reconstruction(
    revised: pd.DataFrame, legacy_whitelist: pd.DataFrame
) -> pd.DataFrame:
    """Controlled failure: assign encoded-CB counts using legacy barcode order.

    This preserves the observed values exactly while emulating application of
    the old whitelist ordering to the revised 8-nt coordinate codec.
    """
    values = revised.sort_values("cb")[
        ["dna_reads", "mapped_reads", "dup_reads"]
    ].reset_index(drop=True)
    coords = legacy_whitelist.sort_values("barcode")[["x", "y"]].reset_index(drop=True)
    if len(values) != len(coords):
        raise ValueError(
            f"Coordinate reconstruction requires equal grids: {len(values)} != {len(coords)}"
        )
    return pd.concat([coords, values], axis=1)


def neighbor_spearman(matrix: np.ndarray) -> float:
    """Spearman correlation across all horizontal and vertical neighbors."""
    a = np.r_[matrix[:, :-1].ravel(), matrix[:-1, :].ravel()]
    b = np.r_[matrix[:, 1:].ravel(), matrix[1:, :].ravel()]
    keep = np.isfinite(a) & np.isfinite(b)
    return float(spearmanr(a[keep], b[keep]).statistic) if keep.sum() > 2 else np.nan


def detect_scar(
    matrix: np.ndarray,
    max_width: int = 3,
    edge_margin: int = 3,
    support_percentile: float = 35.0,
) -> ScarCall:
    """Detect the strongest locally depleted row/column band.

    Candidate medians are compared with equally local flanks. Candidates whose
    flanks are below the profile's 35th percentile are excluded to avoid
    calling empty tissue boundaries as scars.
    """
    best: ScarCall | None = None
    for axis, name in ((0, "vertical"), (1, "horizontal")):
        profile = np.nanmedian(matrix, axis=axis)
        support = np.nanpercentile(profile, support_percentile)
        n = len(profile)
        for width in range(1, min(max_width, n // 5) + 1):
            flank = max(edge_margin, width)
            for start0 in range(flank, n - width - flank + 1):
                inside = float(np.nanmedian(profile[start0 : start0 + width]))
                adjacent = float(
                    np.nanmedian(
                        np.r_[
                            profile[start0 - flank : start0],
                            profile[start0 + width : start0 + width + flank],
                        ]
                    )
                )
                if adjacent <= 0 or adjacent < support:
                    continue
                ratio = inside / adjacent
                call = ScarCall(
                    name,
                    start0 + 1,
                    start0 + width,
                    ratio,
                    inside,
                    adjacent,
                    profile,
                )
                if best is None or call.ratio < best.ratio:
                    best = call
    if best is None:
        raise ValueError("No supported scar candidate was detected")
    return best
