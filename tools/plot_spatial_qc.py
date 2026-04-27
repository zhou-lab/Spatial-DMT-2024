#!/usr/bin/env python3
"""
Generate spatial heatmap PNGs from per-barcode QC metrics.

Input: {sample}_spatial_barcode_counts.tsv  (barcode, x, y, dna_reads, mapped_reads, dup_reads)
Output: {output_dir}/spatial_*.png
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


BLUE_YELLOW = LinearSegmentedColormap.from_list('blue_yellow', ['#000033', '#2255aa', '#faf9cf'])
FIRE        = LinearSegmentedColormap.from_list('fire',        ['#220022', '#dd4400', '#ffff00'])
RED_GRAY    = LinearSegmentedColormap.from_list('red_gray',    ['#e0e0e0', '#ee3322', '#880000'])


def load_data(spatial_counts_path):
    df = pd.read_csv(spatial_counts_path, sep='\t')
    df['dna_reads']    = pd.to_numeric(df['dna_reads'],    errors='coerce').fillna(0)
    df['mapped_reads'] = pd.to_numeric(df.get('mapped_reads', 0), errors='coerce').fillna(0)
    df['dup_reads']    = pd.to_numeric(df.get('dup_reads', 0),    errors='coerce').fillna(0)

    df['log_dna_reads']  = np.log10(df['dna_reads'] + 1)
    df['mapping_rate']   = np.where(df['dna_reads'] > 0, df['mapped_reads'] / df['dna_reads'], np.nan)
    df['dup_rate']       = np.where(df['mapped_reads'] > 0, df['dup_reads'] / df['mapped_reads'], np.nan)
    return df


def plot_metric(df, col, title, cmap, output_path, img_path=None, vmin=None, vmax=None):
    fig, ax = plt.subplots(figsize=(7, 7))

    if img_path and os.path.exists(img_path):
        img = mpimg.imread(img_path)
        ax.imshow(img, extent=[0.5, 50.5, 0.5, 50.5], aspect='auto')

    vals = df[col].values
    sc = ax.scatter(
        df['x'], df['y'],
        c=vals, cmap=cmap, s=80, alpha=0.7,
        vmin=vmin if vmin is not None else np.nanpercentile(vals, 2),
        vmax=vmax if vmax is not None else np.nanpercentile(vals, 98),
    )
    plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel('Lane (x)')
    ax.set_ylabel('Position (y)')
    ax.set_xlim(0.5, 50.5)
    ax.set_ylim(50.5, 0.5)
    ax.set_aspect('equal')
    ax.grid(False)
    fig.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_barcode_rank(df, sample, output_path):
    counts = df['dna_reads'].sort_values(ascending=False).values
    in_spatial = np.ones(len(counts), dtype=bool)  # all rows are designed spatial barcodes
    ranks = np.arange(1, len(counts) + 1)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(ranks, counts, s=5, alpha=0.6, c='steelblue')
    ax.axvline(x=len(counts), color='k', linestyle='--', alpha=0.5,
               label=f'{len(counts)} designed barcodes')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Barcode rank (log)')
    ax.set_ylabel('DNA read count (log)')
    ax.set_title(f'{sample} — Barcode rank vs read count')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Generate spatial QC heatmaps')
    parser.add_argument('--spatial-counts', required=True,
                        help='{sample}_spatial_barcode_counts.tsv')
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--img', default='', help='Optional tissue image path')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = load_data(args.spatial_counts)
    img = args.img if args.img else None

    metrics = [
        ('dna_reads',     'DNA reads per barcode',          BLUE_YELLOW, None, None),
        ('log_dna_reads', 'DNA reads per barcode (log10)',  BLUE_YELLOW, None, None),
        ('mapping_rate',  'Mapping rate',                   FIRE,        0.0,  1.0),
        ('dup_rate',      'Duplication rate',               RED_GRAY,    0.0,  1.0),
    ]

    for col, title, cmap, vmin, vmax in metrics:
        if col not in df.columns:
            continue
        out = os.path.join(args.output_dir, f'spatial_{col}.png')
        plot_metric(df, col, title, cmap, out, img_path=img, vmin=vmin, vmax=vmax)
        print(f'  {out}')

    rank_out = os.path.join(args.output_dir, 'barcode_rank.png')
    plot_barcode_rank(df, args.sample, rank_out)
    print(f'  {rank_out}')


if __name__ == '__main__':
    main()
