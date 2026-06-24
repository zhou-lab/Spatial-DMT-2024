#!/usr/bin/env python3
"""Summarize a STARsolo raw matrix by spatial barcode and render RNA maps."""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


BLUE_YELLOW = LinearSegmentedColormap.from_list(
    "blue_yellow", ["#000033", "#2255aa", "#faf9cf"]
)
GREEN_YELLOW = LinearSegmentedColormap.from_list(
    "green_yellow", ["#002b20", "#16856b", "#fff2a8"]
)


def read_whitelist(path):
    coords = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                coords[fields[0]] = (int(fields[1]), int(fields[2]))
    return coords


def read_barcodes(path):
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def aggregate_matrix(path, n_barcodes):
    umis = np.zeros(n_barcodes, dtype=np.int64)
    genes = np.zeros(n_barcodes, dtype=np.int64)
    dimensions_seen = False
    with open(path) as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if len(fields) != 3:
                continue
            if not dimensions_seen:
                dimensions_seen = True
                continue
            _, barcode_index, value = fields
            index = int(barcode_index) - 1
            umis[index] += int(value)
            genes[index] += 1
    return umis, genes


def plot_metric(df, col, title, cmap, output_path, grid_x, grid_y):
    fig, ax = plt.subplots(figsize=(7, 7))
    vals = df[col].values
    scatter = ax.scatter(
        df["x"],
        df["y"],
        c=vals,
        cmap=cmap,
        s=max(10, 4000 / max(grid_x, grid_y)),
        alpha=0.7,
        vmin=np.nanpercentile(vals, 2),
        vmax=np.nanpercentile(vals, 98),
    )
    plt.colorbar(scatter, ax=ax, fraction=0.03, pad=0.02)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_xlim(0.5, grid_x + 0.5)
    ax.set_ylim(grid_y + 0.5, 0.5)
    ax.set_aspect("equal")
    ax.grid(False)
    fig.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_barcode_rank(df, sample, output_path):
    counts = df["rna_umis"].sort_values(ascending=False).values
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(np.arange(1, len(counts) + 1), counts, s=5, alpha=0.6, c="#16856b")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Barcode rank (log)")
    ax.set_ylabel("RNA UMIs (log)")
    ax.set_title(f"{sample} - Barcode rank vs RNA UMI count")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--whitelist", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--plots-dir", required=True)
    args = parser.parse_args()

    coords = read_whitelist(args.whitelist)
    barcodes = read_barcodes(args.barcodes)
    umis, genes = aggregate_matrix(args.matrix, len(barcodes))

    rows = []
    missing = []
    for barcode, n_umis, n_genes in zip(barcodes, umis, genes):
        if barcode not in coords:
            missing.append(barcode)
            continue
        x, y = coords[barcode]
        rows.append((barcode, x, y, int(n_umis), int(n_genes)))
    if missing:
        raise ValueError(f"{len(missing)} STARsolo barcodes are absent from the spatial whitelist")

    os.makedirs(os.path.dirname(os.path.abspath(args.table)), exist_ok=True)
    os.makedirs(args.plots_dir, exist_ok=True)
    df = pd.DataFrame(rows, columns=["barcode", "x", "y", "rna_umis", "rna_genes"])
    df.to_csv(args.table, sep="\t", index=False)

    df["log_rna_umis"] = np.log10(df["rna_umis"] + 1)
    df["log_rna_genes"] = np.log10(df["rna_genes"] + 1)
    grid_x = int(df["x"].max())
    grid_y = int(df["y"].max())
    for col, title, cmap in [
        ("rna_umis", "RNA UMIs per barcode", GREEN_YELLOW),
        ("log_rna_umis", "RNA UMIs per barcode (log10)", GREEN_YELLOW),
        ("rna_genes", "Detected RNA genes per barcode", BLUE_YELLOW),
        ("log_rna_genes", "Detected RNA genes per barcode (log10)", BLUE_YELLOW),
    ]:
        plot_metric(
            df,
            col,
            title,
            cmap,
            os.path.join(args.plots_dir, f"spatial_{col}.png"),
            grid_x,
            grid_y,
        )
    plot_barcode_rank(df, args.sample, os.path.join(args.plots_dir, "rna_barcode_rank.png"))
    print(f"Wrote {args.table}")


if __name__ == "__main__":
    main()
