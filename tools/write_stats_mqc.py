#!/usr/bin/env python3
"""
Generate a MultiQC generalstats TSV from trim stats, spatial barcode counts,
and bulk ChromHMM methylation levels (TssA, Tx) derived from a yame summary file.

Reads:
  - {sample}_trim_stats.txt      (key-value: Reads_total, Reads_passed, WL_unmatched, ...)
  - {sample}_spatial_barcode_counts.tsv (legacy fallback for the designed count)
  - {sample}_chromhmm_summary.tsv (yame summary output)
  - optional lambda BAM flagstat; if supplied, genomic read pairs are
    Reads_total minus lambda pairs with both mates mapped.

Writes:
  - {output_dir}/{sample}_stats_mqc.tsv  (MultiQC generalstats)
"""
import argparse
import os
import re


def lambda_both_mapped_pairs(flagstat_path):
    """Return lambda read pairs where both mates mapped from samtools flagstat."""
    if not flagstat_path:
        return None
    try:
        with open(flagstat_path) as handle:
            for line in handle:
                m = re.match(r'^(\d+)\s+\+\s+\d+\s+with itself and mate mapped$', line.rstrip())
                if m:
                    return int(m.group(1)) // 2
    except FileNotFoundError:
        return None
    return None


def genomic_total(stats, lambda_flagstat=None):
    """DNA read pairs excluding the barcode85/lambda bucket."""
    n_total = stats.get('Reads_total')
    lambda_pairs = lambda_both_mapped_pairs(lambda_flagstat)
    if n_total is not None and lambda_pairs is not None:
        return max(0, n_total - lambda_pairs)
    if stats.get('Reads_detected') is not None:
        return stats['Reads_detected']
    if n_total is None:
        return None
    for key in ('Reads_barcode85', 'Reads_written_barcode85'):
        if stats.get(key) is not None:
            return n_total - stats[key]
    return n_total


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--trim-stats',       required=True)
    p.add_argument('--spatial-counts',   required=True)
    p.add_argument('--chromhmm-summary', required=True)
    p.add_argument('--lambda-flagstat',
                   help='samtools flagstat for lambda BAM; uses pairs with both mates mapped as temporary control subtraction')
    p.add_argument('--sample',           required=True)
    p.add_argument('--output-dir',       required=True)
    args = p.parse_args()

    stats = {}
    with open(args.trim_stats) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                try:
                    stats[parts[0]] = int(parts[1])
                except ValueError:
                    pass

    n_designed_from_spatial = 0
    with open(args.spatial_counts) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                try:
                    n_designed_from_spatial += int(parts[3])  # dna_reads column
                except ValueError:
                    pass

    n_total = genomic_total(stats, args.lambda_flagstat)
    n_assigned = stats.get('Reads_passed')
    # Reads_passed is the pre-alignment whitelist-assigned count. Spatial
    # counts may be post-bsconv and therefore must not define barcode yield.
    n_designed = n_assigned if n_assigned is not None else n_designed_from_spatial
    n_wl_miss = stats.get('WL_unmatched', stats.get('Unmatched_WM', 0))
    n_linkers = None if n_assigned is None else n_assigned + n_wl_miss

    # Read bulk ChromHMM methylation levels from summary file
    chromhmm = {}
    if os.path.exists(args.chromhmm_summary):
        with open(args.chromhmm_summary) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 10 and parts[0] != 'QFile':
                    state = parts[3]   # Mask column (state name)
                    beta  = parts[9]   # Beta column (methylation level)
                    if state in ('TssA', 'Tss', 'Tx'):
                        try:
                            chromhmm[state] = float(beta)
                        except ValueError:
                            pass

    def pct(num, denom):
        if num is None or not denom:
            return 'NA'
        return f"{100.0 * num / denom:.2f}"

    def fmt_beta(val):
        if val is None:
            return 'NA'
        return f"{val:.4f}"

    # hg38 names the active TSS state TssA; the mm10 ChromHMM mask names it Tss.
    tssa = chromhmm.get('TssA', chromhmm.get('Tss'))
    tx   = chromhmm.get('Tx')

    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, f"{args.sample}_stats_mqc.tsv")
    with open(out_path, 'w') as f:
        f.write('# id: "stats"\n')
        f.write('# section_name: "Spatial Methyl QC Stats"\n')
        f.write('# plot_type: "generalstats"\n')
        f.write('# pconfig:\n')
        f.write('#    - total read pairs:       {title: "Genomic Read Pairs",     format: "{:,.0f}"}\n')
        f.write('#    - both linkers:           {title: "Both Linkers",           format: "{:,.0f}"}\n')
        f.write('#    - with designed barcodes: {title: "Designed Barcodes",      format: "{:,.0f}"}\n')
        f.write('#    - both linkers / total:   {title: "Linkers / Total (%)",    min: 0, max: 100, suffix: "%"}\n')
        f.write('#    - designed / linkers:     {title: "Designed / Linkers (%)", min: 0, max: 100, suffix: "%"}\n')
        f.write('#    - designed / total:       {title: "Designed / Total (%)",   min: 0, max: 100, suffix: "%"}\n')
        f.write('#    - TssA meth:              {title: "TssA Meth",              min: 0, max: 1,   format: "{:.4f}"}\n')
        f.write('#    - Tx meth:                {title: "Tx Meth",                min: 0, max: 1,   format: "{:.4f}"}\n')
        f.write('Sample\ttotal read pairs\tboth linkers\twith designed barcodes\t'
                'both linkers / total\tdesigned / linkers\tdesigned / total\t'
                'TssA meth\tTx meth\n')
        f.write(f"{args.sample}\t"
                f"{n_total if n_total is not None else 'NA'}\t"
                f"{n_linkers if n_linkers is not None else 'NA'}\t"
                f"{n_designed}\t"
                f"{pct(n_linkers, n_total)}\t"
                f"{pct(n_designed, n_linkers)}\t"
                f"{pct(n_designed, n_total)}\t"
                f"{fmt_beta(tssa)}\t"
                f"{fmt_beta(tx)}\n")
    print(f"Written: {out_path}")


if __name__ == '__main__':
    main()
