#!/usr/bin/env python3
"""Plot per-cell CpG methylation from an indexed YAME .cg file."""

import argparse
import csv
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COORD_RE = re.compile(r"^(\d{2})(\d{2})$")


def parse_float(value):
    return float(value) if value not in {"NA", "NaN", "nan", ""} else np.nan


def indexed_cells(cg_path):
    cells = []
    with open(f"{cg_path}.idx") as handle:
        for line in handle:
            cell = line.split("\t", 1)[0]
            if COORD_RE.fullmatch(cell):
                cells.append(cell)
    if not cells:
        raise RuntimeError(f"No spatial coordinates found in {cg_path}.idx")
    return cells


def summarize_cell(cg_path, yame, cell):
    subset = subprocess.Popen(
        [yame, "subset", cg_path, cell],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert subset.stdout is not None
    summary = subprocess.Popen(
        [yame, "summary", "-H", "-q", cell, "-"],
        stdin=subset.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    subset.stdout.close()
    stdout, summary_stderr = summary.communicate()
    subset_stderr = subset.stderr.read().decode() if subset.stderr else ""
    subset_code = subset.wait()
    if subset_code or summary.returncode:
        raise RuntimeError(
            f"YAME failed for cell {cell}: {subset_stderr}{summary_stderr}"
        )

    fields = stdout.rstrip("\n").split("\t")
    if len(fields) != 11:
        raise RuntimeError(f"Unexpected YAME summary for cell {cell}: {stdout!r}")
    match = COORD_RE.fullmatch(cell)
    assert match is not None
    return {
        "spatial_id": cell,
        "x": int(match.group(1)),
        "y": int(match.group(2)),
        "covered_cpgs": int(fields[5]),
        "beta": parse_float(fields[9]),
        "depth": parse_float(fields[10]),
    }


def summarize_cells(cg_path, yame, threads):
    cells = indexed_cells(cg_path)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        rows = list(
            executor.map(
                lambda cell: summarize_cell(cg_path, yame, cell),
                cells,
            )
        )
    return rows


def summarize_cells_bulk(cg_path, yame):
    """Summarize every indexed record with one sequential YAME pass.

    This is substantially cheaper than launching one ``subset`` process per
    pixel for large spatial files.  YAME emits one summary row per record;
    format-3 ``N_query`` is the number of covered CpGs and ``Beta`` is the
    methylated fraction used by the per-pixel map.
    """
    completed = subprocess.run(
        [yame, "summary", "-H", cg_path],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    rows = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 11:
            raise RuntimeError(f"Unexpected YAME summary row: {line!r}")
        cell = fields[1]
        match = COORD_RE.fullmatch(cell)
        if match is None:
            continue
        rows.append(
            {
                "spatial_id": cell,
                "x": int(match.group(1)),
                "y": int(match.group(2)),
                "covered_cpgs": int(fields[5]),
                "beta": parse_float(fields[9]),
                "depth": parse_float(fields[10]),
            }
        )
    expected = len(indexed_cells(cg_path))
    if len(rows) != expected:
        raise RuntimeError(
            f"YAME summary returned {len(rows)} spatial records; index has {expected}"
        )
    return rows


def write_table(path, sample, rows, min_covered_cpgs):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "spatial_id",
                "x",
                "y",
                "covered_cpgs",
                "cpg_methylation",
                "depth",
                "passes_threshold",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    sample,
                    row["spatial_id"],
                    row["x"],
                    row["y"],
                    row["covered_cpgs"],
                    f'{row["beta"]:.6f}' if np.isfinite(row["beta"]) else "NA",
                    f'{row["depth"]:.6f}' if np.isfinite(row["depth"]) else "NA",
                    int(row["covered_cpgs"] >= min_covered_cpgs),
                ]
            )


def plot_map(
    path,
    sample,
    rows,
    grid_x,
    grid_y,
    min_covered_cpgs,
    vmin,
    vmax,
):
    matrix = np.full((grid_y, grid_x), np.nan)
    retained = []
    for row in rows:
        x = row["x"]
        y = row["y"]
        if not (1 <= x <= grid_x and 1 <= y <= grid_y):
            continue
        if row["covered_cpgs"] < min_covered_cpgs:
            continue
        matrix[y - 1, x - 1] = row["beta"]
        retained.append(row["beta"])

    if not retained:
        raise RuntimeError(
            f"No cells have at least {min_covered_cpgs} covered CpGs"
        )

    cmap = plt.get_cmap("RdYlBu_r").copy()
    cmap.set_bad("#eeeeee")
    fig, ax = plt.subplots(figsize=(7.4, 7.1))
    image = ax.imshow(
        matrix,
        origin="upper",
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extent=[0.5, grid_x + 0.5, grid_y + 0.5, 0.5],
        aspect="equal",
    )
    colorbar = fig.colorbar(image, ax=ax, fraction=0.042, pad=0.03)
    colorbar.set_label("CpG methylation level")
    ax.set_title(f"{sample}\nSpatial CpG methylation", fontsize=13)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_xticks(np.arange(10, grid_x + 1, 10))
    ax.set_yticks(np.arange(10, grid_y + 1, 10))
    ax.text(
        0.0,
        -0.10,
        (
            f"{len(retained):,} cells with >= {min_covered_cpgs:,} covered CpGs; "
            f"median beta = {np.median(retained):.3f}"
        ),
        transform=ax.transAxes,
        fontsize=9,
        color="#444444",
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cg", required=True, help="Indexed multi-cell YAME .cg file")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True, help="Output PNG path")
    parser.add_argument(
        "--table",
        default="",
        help="Per-cell TSV path (default: output path with .tsv suffix)",
    )
    parser.add_argument("--grid-x", type=int, default=96)
    parser.add_argument("--grid-y", type=int, default=96)
    parser.add_argument("--min-covered-cpgs", type=int, default=100)
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=1.0)
    parser.add_argument("--yame", default="yame")
    parser.add_argument("--threads", type=int, default=24)
    args = parser.parse_args()

    if not os.path.isfile(args.cg):
        parser.error(f"CG file does not exist: {args.cg}")
    if not os.path.isfile(f"{args.cg}.idx"):
        parser.error(f"CG index does not exist: {args.cg}.idx")
    if args.vmin >= args.vmax:
        parser.error("--vmin must be smaller than --vmax")
    if args.threads < 1:
        parser.error("--threads must be at least 1")

    output = os.path.abspath(args.output)
    table = os.path.abspath(
        args.table or f"{os.path.splitext(output)[0]}.tsv"
    )
    os.makedirs(os.path.dirname(output), exist_ok=True)
    os.makedirs(os.path.dirname(table), exist_ok=True)

    rows = summarize_cells(args.cg, args.yame, args.threads)
    write_table(table, args.sample, rows, args.min_covered_cpgs)
    plot_map(
        output,
        args.sample,
        rows,
        args.grid_x,
        args.grid_y,
        args.min_covered_cpgs,
        args.vmin,
        args.vmax,
    )
    passing = sum(
        row["covered_cpgs"] >= args.min_covered_cpgs for row in rows
    )
    print(f"Wrote {output}")
    print(f"Wrote {table}")
    print(f"{passing:,}/{len(rows):,} indexed cells passed the coverage threshold")


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)
