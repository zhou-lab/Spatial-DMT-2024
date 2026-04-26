#!/usr/bin/env python3
"""
Map per-barcode counts to spatial (x, y) positions using the spatial whitelist.

Whitelist format (tab-separated, no header):
  col0: barcode_A  col1: x  col2: y  col3: barcode (used for demux)  col4: barcode_combined

Outputs:
  {output_dir}/{sample}_spatial_barcode_counts.tsv   x, y, barcode, dna_reads, mapped_reads, dup_reads
  {output_dir}/{sample}_barcodes_bool_index.tsv       barcode, count, in_spatial
"""
import argparse
import os


def load_spatial_barcodes(whitelist_path):
    """Return list of (x, y, barcode) from col1, col2, col3 (0-indexed)."""
    entries = []
    with open(whitelist_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4:
                entries.append((int(parts[1]), int(parts[2]), parts[3]))
    return entries


def load_barcode_counts(counts_path):
    """Return dict: barcode -> dict of metric columns."""
    counts = {}
    with open(counts_path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if parts:
                bc = parts[0]
                counts[bc] = {header[i]: parts[i] for i in range(1, len(header))}
    return counts, header[1:]


def main():
    parser = argparse.ArgumentParser(description='Match barcode counts to spatial positions')
    parser.add_argument('--barcode-counts', required=True,
                        help='TSV from count_barcodes.py (barcode + metric columns)')
    parser.add_argument('--whitelist', required=True,
                        help='barcodes/spatial_barcodes.txt')
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    spatial = load_spatial_barcodes(args.whitelist)
    counts, metric_cols = load_barcode_counts(args.barcode_counts)
    spatial_set = {bc for _, _, bc in spatial}

    # Spatial counts: one row per designed pixel
    spatial_out = os.path.join(args.output_dir, f'{args.sample}_spatial_barcode_counts.tsv')
    with open(spatial_out, 'w') as out:
        out.write('barcode\tx\ty\t' + '\t'.join(metric_cols) + '\n')
        for x, y, bc in spatial:
            row = counts.get(bc, {col: '0' for col in metric_cols})
            vals = '\t'.join(row.get(col, '0') for col in metric_cols)
            out.write(f'{bc}\t{x}\t{y}\t{vals}\n')

    # Bool index: all observed barcodes flagged True/False
    bool_out = os.path.join(args.output_dir, f'{args.sample}_barcodes_bool_index.tsv')
    all_barcodes = sorted(counts.items(), key=lambda kv: -int(kv[1].get('dna_reads', 0)))
    with open(bool_out, 'w') as out:
        out.write('barcode\tdna_reads\tin_spatial\n')
        for bc, row in all_barcodes:
            flag = bc in spatial_set
            out.write(f'{bc}\t{row.get("dna_reads", 0)}\t{flag}\n')

    n_matched = sum(1 for _, _, bc in spatial if bc in counts and int(counts[bc].get('dna_reads', 0)) > 0)
    print(f'Matched {n_matched}/{len(spatial)} spatial barcodes with reads')
    print(f'Wrote {spatial_out}')
    print(f'Wrote {bool_out}')


if __name__ == '__main__':
    main()
