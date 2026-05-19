#!/usr/bin/env python3
"""
Count reads per barcode from demuxed FASTQs and dupsifter stats, map to
spatial (x, y) positions using the spatial whitelist, and generate QC plots.

Whitelist format (tab-separated, no header):
  col0: barcode_A  col1: x  col2: y  col3: barcode (demux sequence)

Outputs:
  {table_dir}/{sample}_spatial_barcode_counts.tsv
  {plots_dir}/spatial_{metric}.png  (dna_reads, log_dna_reads, mapping_rate, dup_rate)
  {plots_dir}/barcode_rank.png
"""
import argparse
import gzip
import os
import re

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


def count_fastq_reads(fq_path):
    opener = gzip.open if fq_path.endswith('.gz') else open
    n = 0
    with opener(fq_path, 'rt') as fh:
        for _ in fh:
            n += 1
    return n // 4


def parse_dupsifter(stat_path):
    """Return (both_mapped_pairs, dup_pairs) from a dupsifter .stat file."""
    both_mapped = dup = 0
    with open(stat_path) as fh:
        for line in fh:
            m = re.search(r'number of reads with both reads mapped:\s+(\d+)', line)
            if m:
                both_mapped = int(m.group(1))
            m = re.search(r'number of reads with both reads marked as duplicates:\s+(\d+)', line)
            if m:
                dup = int(m.group(1))
    return both_mapped, dup


def load_spatial_whitelist(path, barcode_col=3):
    """Return list of (x, y, barcode); x=col1, y=col2 in whitelist."""
    entries = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) > barcode_col:
                entries.append((int(parts[1]), int(parts[2]), parts[barcode_col]))
    return entries


def plot_metric(df, col, title, cmap, output_path, img_path=None, vmin=None, vmax=None,
                spot_size=40, spot_alpha=0.50):
    fig, ax = plt.subplots(figsize=(7, 7))
    if img_path and os.path.exists(img_path):
        img = mpimg.imread(img_path)
        ax.imshow(img, extent=[0.5, 50.5, 0.5, 50.5], aspect='auto', zorder=0)
    vals = df[col].values
    sc = ax.scatter(
        df['x'], df['y'], c=vals, cmap=cmap, s=spot_size, alpha=spot_alpha,
        edgecolors='none', linewidths=0, zorder=2,
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
    print(f'  {output_path}')


def plot_barcode_rank(df, sample, output_path):
    counts = df['dna_reads'].sort_values(ascending=False).values
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(np.arange(1, len(counts) + 1), counts, s=5, alpha=0.6, c='steelblue')
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
    print(f'  {output_path}')


def main():
    parser = argparse.ArgumentParser(
        description='Per-barcode read counts, spatial mapping, and QC plots'
    )
    parser.add_argument('--dmux-dir',    required=True, help='Dir containing {barcode}_R1.fq.gz')
    parser.add_argument('--bam-dir',     required=True, help='Dir containing {barcode}.dupsifter.stat')
    parser.add_argument('--whitelist',   required=True, help='Spatial barcode whitelist (tab-separated)')
    parser.add_argument('--barcode-col', type=int, default=3,
                        help='0-indexed barcode column in whitelist (default: 3)')
    parser.add_argument('--sample',      required=True, help='Sample name')
    parser.add_argument('--table-dir',   required=True, help='Output directory for TSV')
    parser.add_argument('--plots-dir',   required=True, help='Output directory for PNGs')
    parser.add_argument('--img',         default='',    help='Optional tissue image path')
    parser.add_argument('--spot-size',   type=float, default=40,
                        help='Scatter marker area in points^2 for spatial overlays (default: 40)')
    parser.add_argument('--spot-alpha',  type=float, default=0.50,
                        help='Scatter marker opacity for spatial overlays (default: 0.50)')
    args = parser.parse_args()

    os.makedirs(args.table_dir, exist_ok=True)
    os.makedirs(args.plots_dir, exist_ok=True)

    # Count reads per barcode from dmux FASTQs and dupsifter stats
    barcodes = sorted(
        f.replace('_R1.fq.gz', '')
        for f in os.listdir(args.dmux_dir)
        if f.endswith('_R1.fq.gz')
    )
    counts = {}
    for bc in barcodes:
        fq   = os.path.join(args.dmux_dir, f'{bc}_R1.fq.gz')
        stat = os.path.join(args.bam_dir,  f'{bc}.dupsifter.stat')
        dna_reads   = count_fastq_reads(fq)  if os.path.exists(fq)   else 0
        mapped, dup = parse_dupsifter(stat)  if os.path.exists(stat) else (0, 0)
        counts[bc] = (dna_reads, mapped, dup)

    # Write spatial counts TSV
    spatial  = load_spatial_whitelist(args.whitelist, args.barcode_col)
    tsv_path = os.path.join(args.table_dir, f'{args.sample}_spatial_barcode_counts.tsv')
    with open(tsv_path, 'w') as out:
        out.write('barcode\tx\ty\tdna_reads\tmapped_reads\tdup_reads\n')
        for x, y, bc in spatial:
            dna_reads, mapped, dup = counts.get(bc, (0, 0, 0))
            out.write(f'{bc}\t{x}\t{y}\t{dna_reads}\t{mapped}\t{dup}\n')
    n_matched = sum(1 for _, _, bc in spatial if bc in counts and counts[bc][0] > 0)
    print(f'Matched {n_matched}/{len(spatial)} spatial barcodes with reads')
    print(f'Wrote {tsv_path}')

    # Generate QC plots
    df = pd.read_csv(tsv_path, sep='\t')
    df['log_dna_reads'] = np.log10(df['dna_reads'] + 1)
    df['mapping_rate']  = np.where(df['dna_reads'] > 0, df['mapped_reads'] / df['dna_reads'], np.nan)
    df['dup_rate']      = np.where(df['mapped_reads'] > 0, df['dup_reads'] / df['mapped_reads'], np.nan)
    img = args.img if args.img else None

    for col, title, cmap, vmin, vmax in [
        ('dna_reads',     'DNA reads per barcode',         BLUE_YELLOW, None, None),
        ('log_dna_reads', 'DNA reads per barcode (log10)', BLUE_YELLOW, None, None),
        ('mapping_rate',  'Mapping rate',                  FIRE,        0.0,  1.0),
        ('dup_rate',      'Duplication rate',              RED_GRAY,    0.0,  1.0),
    ]:
        plot_metric(df, col, title, cmap,
                    os.path.join(args.plots_dir, f'spatial_{col}.png'),
                    img_path=img, vmin=vmin, vmax=vmax,
                    spot_size=args.spot_size, spot_alpha=args.spot_alpha)

    plot_barcode_rank(df, args.sample, os.path.join(args.plots_dir, 'barcode_rank.png'))


if __name__ == '__main__':
    main()
