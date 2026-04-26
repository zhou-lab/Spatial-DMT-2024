#!/usr/bin/env python3
"""
Count reads per barcode from demuxed FASTQ files and per-barcode dupsifter stats.

Outputs:
  {output_dir}/{sample}_barcode_counts.tsv   barcode <TAB> dna_reads <TAB> mapped <TAB> dup_rate
"""
import argparse
import gzip
import os
import re


def count_fastq_reads(fq_path):
    opener = gzip.open if fq_path.endswith('.gz') else open
    n = 0
    with opener(fq_path, 'rt') as fh:
        for _ in fh:
            n += 1
    return n // 4


def parse_dupsifter(stat_path):
    """Return (both_mapped_pairs, dup_pairs) from a dupsifter .stat file.

    'both reads mapped'            = pairs where R1+R2 both aligned
    'both reads marked as dup'     = pairs where both reads are duplicates
    """
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


def main():
    parser = argparse.ArgumentParser(description='Count reads per barcode from demux FASTQs')
    parser.add_argument('--dmux-dir', required=True, help='Dir containing {barcode}_R1.fq.gz')
    parser.add_argument('--bam-dir', required=True, help='Dir containing {barcode}.dupsifter.stat')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, f'{args.sample}_barcode_counts.tsv')

    barcodes = sorted(
        os.path.basename(f).replace('_R1.fq.gz', '')
        for f in os.listdir(args.dmux_dir)
        if f.endswith('_R1.fq.gz')
    )

    with open(out_path, 'w') as out:
        out.write('barcode\tdna_reads\tmapped_reads\tdup_reads\n')
        for bc in barcodes:
            fq = os.path.join(args.dmux_dir, f'{bc}_R1.fq.gz')
            stat = os.path.join(args.bam_dir, f'{bc}.dupsifter.stat')

            dna_reads = count_fastq_reads(fq) if os.path.exists(fq) else 0
            mapped, dup = parse_dupsifter(stat) if os.path.exists(stat) else (0, 0)
            out.write(f'{bc}\t{dna_reads}\t{mapped}\t{dup}\n')

    print(f'Wrote {len(barcodes)} barcodes to {out_path}')


if __name__ == '__main__':
    main()
