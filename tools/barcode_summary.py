#!/usr/bin/env python3
"""Per-cell read counts (total / mapped / dup) from the merged deduped BAM,
mapped onto the (X, Y) spatial grid and rendered as heatmaps + a barcode-rank
plot.

Inputs: one merged coord-sorted BAM produced by `dna_biscuit_align`. Every
alignment carries CB:Z:<8-char ACGT> set by biscuit -9 from the
<origid>_<CB>_<UMI> read-name convention. We iterate the merged BAM once,
aggregate per-CB counts, and decode CB back to (X, Y) via cb_codec.

Outputs:
  {table_dir}/{sample}_spatial_barcode_counts.tsv     per-cell counts
  {plots_dir}/spatial_{metric}.png                    heatmaps over the grid
  {plots_dir}/barcode_rank.png                        rank vs read count

Replaces the earlier per-cell-FASTQ + per-cell-dupsifter-stat path which is
no longer produced after the demux-after-align rewrite.
"""
import argparse, os, sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import pysam

from cb_codec import coord_to_acgt8, UNMATCHED_CB


BLUE_YELLOW = LinearSegmentedColormap.from_list('blue_yellow', ['#000033', '#2255aa', '#faf9cf'])
FIRE        = LinearSegmentedColormap.from_list('fire',        ['#220022', '#dd4400', '#ffff00'])
RED_GRAY    = LinearSegmentedColormap.from_list('red_gray',    ['#e0e0e0', '#ee3322', '#880000'])


def aggregate_merged_bam(bam_path):
    """One pass over the merged BAM; returns dict[cb_str] -> [total, mapped, dup].
    Counts R1 only (one bucket per pair) so totals match per-pair semantics.
    Skips secondary/supplementary alignments."""
    counts = defaultdict(lambda: [0, 0, 0])
    n_missing_cb = 0
    bf = pysam.AlignmentFile(bam_path, "rb")
    for read in bf:
        if read.is_secondary or read.is_supplementary:
            continue
        if not read.is_read1:
            continue
        try:
            cb = read.get_tag("CB")
        except KeyError:
            n_missing_cb += 1
            continue
        bucket = counts[cb]
        bucket[0] += 1
        if not read.is_unmapped:
            bucket[1] += 1
            if read.is_duplicate:
                bucket[2] += 1
    bf.close()
    if n_missing_cb:
        sys.stderr.write(f"WARN: {n_missing_cb} primary R1 reads had no CB tag\n")
    return counts


def plot_metric(df, col, title, cmap, output_path, grid_x, grid_y,
                img_path=None, vmin=None, vmax=None):
    fig, ax = plt.subplots(figsize=(7, 7))
    if img_path and os.path.exists(img_path):
        img = mpimg.imread(img_path)
        ax.imshow(img, extent=[0.5, grid_x + 0.5, 0.5, grid_y + 0.5], aspect='auto')
    vals = df[col].values
    sc = ax.scatter(
        df['x'], df['y'], c=vals, cmap=cmap,
        s=max(10, 4000 / max(grid_x, grid_y)),
        alpha=0.7,
        vmin=vmin if vmin is not None else np.nanpercentile(vals, 2),
        vmax=vmax if vmax is not None else np.nanpercentile(vals, 98),
    )
    plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_xlim(0.5, grid_x + 0.5)
    ax.set_ylim(grid_y + 0.5, 0.5)
    ax.set_aspect('equal')
    ax.grid(False)
    fig.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  {output_path}')


def plot_barcode_rank(df, sample, output_path, grid_x, grid_y):
    counts = df['dna_reads'].sort_values(ascending=False).values
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(np.arange(1, len(counts) + 1), counts, s=5, alpha=0.6, c='steelblue')
    designed = grid_x * grid_y
    ax.axvline(x=designed, color='k', linestyle='--', alpha=0.5,
               label=f'{designed} designed barcodes ({grid_x}x{grid_y})')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Barcode rank (log)')
    ax.set_ylabel('DNA reads (log)')
    ax.set_title(f'{sample} — Barcode rank vs read count')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  {output_path}')


def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--merged-bam', required=True,
                   help='coord-sorted deduped BAM with CB:Z tags (from biscuit -9)')
    p.add_argument('--sample',     required=True)
    p.add_argument('--table-dir',  required=True)
    p.add_argument('--plots-dir',  required=True)
    p.add_argument('--grid-x',     type=int, default=96, help='grid X dim (default 96 = smcseq)')
    p.add_argument('--grid-y',     type=int, default=96, help='grid Y dim (default 96 = smcseq)')
    p.add_argument('--img',        default='', help='optional tissue image overlay path')
    args = p.parse_args()

    os.makedirs(args.table_dir, exist_ok=True)
    os.makedirs(args.plots_dir, exist_ok=True)

    counts = aggregate_merged_bam(args.merged_bam)

    tsv_path = os.path.join(args.table_dir, f'{args.sample}_spatial_barcode_counts.tsv')
    rows = []
    n_cells_with_reads = 0
    for x in range(1, args.grid_x + 1):
        for y in range(1, args.grid_y + 1):
            cb = coord_to_acgt8(x, y)
            total, mapped, dup = counts.get(cb, (0, 0, 0))
            if total > 0:
                n_cells_with_reads += 1
            rows.append((cb, x, y, total, mapped, dup))
    unm_total, unm_mapped, unm_dup = counts.get(UNMATCHED_CB, (0, 0, 0))

    with open(tsv_path, 'w') as out:
        out.write('cb\tx\ty\tdna_reads\tmapped_reads\tdup_reads\n')
        for cb, x, y, t, m, d in rows:
            out.write(f'{cb}\t{x}\t{y}\t{t}\t{m}\t{d}\n')
        out.write(f'{UNMATCHED_CB}\t0\t0\t{unm_total}\t{unm_mapped}\t{unm_dup}\n')
    print(f'Wrote {tsv_path}')
    print(f'  {n_cells_with_reads}/{args.grid_x * args.grid_y} cells with reads')
    print(f'  UNMATCHED: {unm_total} R1 reads ({unm_mapped} mapped, {unm_dup} dup)')

    df = pd.read_csv(tsv_path, sep='\t')
    df_grid = df[df['cb'] != UNMATCHED_CB].copy()
    df_grid['log_dna_reads'] = np.log10(df_grid['dna_reads'] + 1)
    df_grid['mapping_rate']  = np.where(df_grid['dna_reads'] > 0,
                                        df_grid['mapped_reads'] / df_grid['dna_reads'], np.nan)
    df_grid['dup_rate']      = np.where(df_grid['mapped_reads'] > 0,
                                        df_grid['dup_reads'] / df_grid['mapped_reads'], np.nan)

    img = args.img if args.img else None
    for col, title, cmap, vmin, vmax in [
        ('dna_reads',     'DNA reads per cell',          BLUE_YELLOW, None, None),
        ('log_dna_reads', 'DNA reads per cell (log10)',  BLUE_YELLOW, None, None),
        ('mapping_rate',  'Mapping rate',                FIRE,        0.0,  1.0),
        ('dup_rate',      'Duplication rate',            RED_GRAY,    0.0,  1.0),
    ]:
        plot_metric(df_grid, col, title, cmap,
                    os.path.join(args.plots_dir, f'spatial_{col}.png'),
                    args.grid_x, args.grid_y, img_path=img, vmin=vmin, vmax=vmax)

    plot_barcode_rank(df_grid, args.sample,
                      os.path.join(args.plots_dir, 'barcode_rank.png'),
                      args.grid_x, args.grid_y)


if __name__ == '__main__':
    main()
