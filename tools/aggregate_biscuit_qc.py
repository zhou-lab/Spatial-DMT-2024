#!/usr/bin/env python3
"""
Merge per-barcode BISCUIT QC tables into bulk tables for MultiQC.

Reads per-barcode files from --input-dir (named {sample}_{barcode}_*.txt),
merges them across all barcodes in --align-txt, and writes bulk merged
files to --output-dir for consumption by MultiQC's BISCUIT module.

Per-spatial-barcode files (from --spatial-counts) are used for median
metric calculations only.
"""
import argparse
import os
import statistics
from collections import defaultdict


# ─────────────────────────────────────────────────────────────────────────────
# Parsers
# ─────────────────────────────────────────────────────────────────────────────

def parse_covdist_table(filepath):
    counts = defaultdict(int)
    with open(filepath) as f:
        for line in list(f)[2:]:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                counts[int(parts[0])] = int(parts[1])
    return counts

def parse_mapq_table(filepath):
    counts = defaultdict(int)
    with open(filepath) as f:
        for line in list(f)[2:]:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                counts[parts[0]] = int(parts[1])
    return counts

def parse_retention_table(filepath):
    counts = defaultdict(int)
    with open(filepath) as f:
        for line in list(f)[2:]:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                counts[(parts[0], parts[1], parts[2])] = int(parts[3])
    return counts

def parse_dup_report(filepath):
    counts = {}
    with open(filepath) as f:
        for line in list(f)[1:]:
            if ':' in line:
                parts = line.strip().split(':\t')
                if len(parts) >= 2:
                    counts[parts[0]] = int(parts[1])
    return counts

def parse_strand_table(filepath):
    values = {k: 0 for k in ('R1_f_BSW','R1_f_BSC','R1_r_BSW','R1_r_BSC',
                              'R2_f_BSW','R2_f_BSC','R2_r_BSW','R2_r_BSC')}
    with open(filepath) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        for tag, key_w, key_c in [
            ('R1 (f):', 'R1_f_BSW', 'R1_f_BSC'),
            ('R1 (r):', 'R1_r_BSW', 'R1_r_BSC'),
            ('R2 (f):', 'R2_f_BSW', 'R2_f_BSC'),
            ('R2 (r):', 'R2_r_BSW', 'R2_r_BSC'),
        ]:
            if tag in line:
                parts = line.replace(tag, '').strip().split()
                if parts:
                    values[key_w] = int(parts[0])
                if i + 1 < len(lines):
                    nxt = lines[i + 1].strip().split()
                    if nxt and nxt[0].isdigit():
                        values[key_c] = int(nxt[0])
    return values

def parse_total_base_conversion(filepath):
    if not os.path.exists(filepath):
        return {}
    with open(filepath) as f:
        lines = f.readlines()
    if len(lines) < 3:
        return {}
    headers = lines[1].strip().split('\t')
    values  = lines[2].strip().split('\t')
    result = {}
    for h, v in zip(headers, values):
        try:
            result[h.strip()] = float(v.strip())
        except ValueError:
            pass
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Mergers
# ─────────────────────────────────────────────────────────────────────────────

def _merge_sum(dicts):
    merged = defaultdict(int)
    for d in dicts:
        for k, v in d.items():
            merged[k] += v
    return merged


# ─────────────────────────────────────────────────────────────────────────────
# Writers
# ─────────────────────────────────────────────────────────────────────────────

def write_covdist_table(filepath, header, counts):
    with open(filepath, 'w') as f:
        f.write(f"{header}\n")
        f.write("depth\tcount\n")
        for depth in sorted(counts):
            f.write(f"{depth}\t{counts[depth]}\n")

def write_mapq_table(filepath, counts):
    with open(filepath, 'w') as f:
        f.write("BISCUITqc Mapping Quality Table\n")
        f.write("MapQ\tCount\n")
        if 'unmapped' in counts:
            f.write(f"unmapped\t{counts['unmapped']}\n")
        for k in sorted((k for k in counts if k != 'unmapped'), key=int):
            f.write(f"{k}\t{counts[k]}\n")

def write_retention_table(filepath, header, counts):
    with open(filepath, 'w') as f:
        f.write(f"{header}\n")
        f.write("ReadInPair\tPosition\tConversion/Retention\tCount\n")
        for key in sorted(counts, key=lambda x: (int(x[0]), int(x[1]), x[2])):
            f.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{counts[key]}\n")

def write_dup_report(filepath, counts):
    fields = ["Number of duplicate reads", "Number of reads",
              "Number of duplicate q40-reads", "Number of q40-reads"]
    with open(filepath, 'w') as f:
        f.write("BISCUITqc Read Duplication Table\n")
        for field in fields:
            if field in counts:
                f.write(f"{field}:\t{counts[field]}\n")

def write_strand_table(filepath, values):
    with open(filepath, 'w') as f:
        f.write("BISCUITqc Strand Table\n")
        f.write("Strand Distribution:\n")
        f.write("strand\\BS      BSW (f)      BSC (r)\n")
        for tag, key_w, key_c in [
            ('R1 (f)', 'R1_f_BSW', 'R1_f_BSC'),
            ('R1 (r)', 'R1_r_BSW', 'R1_r_BSC'),
            ('R2 (f)', 'R2_f_BSW', 'R2_f_BSC'),
            ('R2 (r)', 'R2_r_BSW', 'R2_r_BSC'),
        ]:
            f.write(f"     {tag}:   {values[key_w]}\n")
            f.write(f"{values[key_c]}\n")


# ─────────────────────────────────────────────────────────────────────────────
# Metric helpers
# ─────────────────────────────────────────────────────────────────────────────

def _conversion_rate(retention_dict):
    rates = [retention_dict[ctx] for ctx in ('CA', 'CC', 'CT') if ctx in retention_dict]
    return (1.0 - sum(rates) / len(rates)) if rates else None

def _mappability_rate(mapq):
    total = sum(mapq.values())
    return (total - mapq.get('unmapped', 0)) / total if total else None

def _cvg_rate(covdist):
    total = sum(covdist.values())
    return (total - covdist.get(0, 0)) / total if total else None

def _retention_rate(ret):
    r = sum(v for (_, _, cr), v in ret.items() if cr == 'R')
    c = sum(v for (_, _, cr), v in ret.items() if cr == 'C')
    return r / (r + c) if (r + c) else None

def _dup_rate(dup):
    total = dup.get('Number of reads', 0)
    return dup.get('Number of duplicate reads', 0) / total if total else None

def _median(values):
    valid = [v for v in values if v is not None]
    return statistics.median(valid) if valid else None


# ─────────────────────────────────────────────────────────────────────────────
# Main merge
# ─────────────────────────────────────────────────────────────────────────────

def merge_qc_files(input_dir, output_dir, sample, barcodes, spatial_barcodes):
    os.makedirs(output_dir, exist_ok=True)

    covdist_types = [
        ('covdist_all_cpg_table.txt',         'BISCUITqc Depth Distribution - All CpGs'),
        ('covdist_all_cpg_topgc_table.txt',   'BISCUITqc Depth Distribution - All Top GC CpGs'),
        ('covdist_all_cpg_botgc_table.txt',   'BISCUITqc Depth Distribution - All Bot GC CpGs'),
    ]

    def _load(bcs, suffix, parser):
        result = []
        for bc in bcs:
            fp = os.path.join(input_dir, f"{sample}_{bc}_{suffix}")
            if os.path.exists(fp):
                result.append(parser(fp))
        return result

    # Covdist tables
    for suffix, header in covdist_types:
        counts_list = _load(barcodes, suffix, parse_covdist_table)
        if counts_list:
            merged = _merge_sum(counts_list)
            write_covdist_table(os.path.join(output_dir, f"{sample}_{suffix}"), header, merged)
            print(f"  covdist {suffix}: merged {len(counts_list)} barcodes")

    # MapQ
    mapq_list = _load(barcodes, 'mapq_table.txt', parse_mapq_table)
    if mapq_list:
        write_mapq_table(os.path.join(output_dir, f"{sample}_mapq_table.txt"), _merge_sum(mapq_list))
        print(f"  mapq: merged {len(mapq_list)} barcodes")

    # CpG / CpH retention
    for tag, suffix, header in [
        ('CpG', 'CpGRetentionByReadPos.txt', 'BISCUITqc CpG Retention by Read Position Table'),
        ('CpH', 'CpHRetentionByReadPos.txt', 'BISCUITqc CpH Retention by Read Position Table'),
    ]:
        lst = _load(barcodes, suffix, parse_retention_table)
        if lst:
            write_retention_table(os.path.join(output_dir, f"{sample}_{suffix}"),
                                  header, _merge_sum(lst))
            print(f"  {tag} retention: merged {len(lst)} barcodes")

    # Dup report
    dup_list = _load(barcodes, 'dup_report.txt', parse_dup_report)
    if dup_list:
        write_dup_report(os.path.join(output_dir, f"{sample}_dup_report.txt"), _merge_sum(dup_list))
        print(f"  dup report: merged {len(dup_list)} barcodes")

    # Strand table
    strand_list = _load(barcodes, 'strand_table.txt', parse_strand_table)
    if strand_list:
        write_strand_table(os.path.join(output_dir, f"{sample}_strand_table.txt"),
                           _merge_sum(strand_list))
        print(f"  strand table: merged {len(strand_list)} barcodes")

    # Median metrics across spatial barcodes
    sp_mapq  = _load(spatial_barcodes, 'mapq_table.txt',          parse_mapq_table)
    sp_covd  = _load(spatial_barcodes, 'covdist_all_cpg_table.txt', parse_covdist_table)
    sp_cpg   = _load(spatial_barcodes, 'CpGRetentionByReadPos.txt',  parse_retention_table)
    sp_cph   = _load(spatial_barcodes, 'CpHRetentionByReadPos.txt',  parse_retention_table)
    sp_dup   = _load(spatial_barcodes, 'dup_report.txt',             parse_dup_report)
    sp_conv  = _load(spatial_barcodes, 'totalBaseConversionRate.txt', parse_total_base_conversion)
    all_conv = _load(barcodes,         'totalBaseConversionRate.txt', parse_total_base_conversion)

    median_mapq = _median([_mappability_rate(d) for d in sp_mapq])
    median_cvg  = _median([_cvg_rate(d)         for d in sp_covd])
    median_cpg  = _median([_retention_rate(d)   for d in sp_cpg])
    median_cph  = _median([_retention_rate(d)   for d in sp_cph])
    median_dup  = _median([_dup_rate(d)          for d in sp_dup])
    median_conv = _median([_conversion_rate(d)   for d in sp_conv])
    all_conv_rates = [_conversion_rate(d) for d in all_conv if _conversion_rate(d) is not None]
    mean_conv   = statistics.mean(all_conv_rates) if all_conv_rates else None

    def pct(v): return f"{v*100:.2f}%" if v is not None else "N/A"
    print(f"\nMedian metrics (spatial barcodes, n={len(spatial_barcodes)}):")
    print(f"  Mapping rate    : {pct(median_mapq)}")
    print(f"  Coverage rate   : {pct(median_cvg)}")
    print(f"  CpG retention   : {pct(median_cpg)}")
    print(f"  CpH retention   : {pct(median_cph)}")
    print(f"  Duplication rate: {pct(median_dup)}")
    print(f"  Conversion rate : {pct(median_conv)}  (mean all: {pct(mean_conv)})")


def main():
    parser = argparse.ArgumentParser(
        description='Merge per-barcode BISCUIT QC tables into bulk files for MultiQC'
    )
    parser.add_argument('--input-dir',      required=True, help='Dir with per-barcode BISCUITqc files')
    parser.add_argument('--output-dir',     required=True, help='Dir to write merged files')
    parser.add_argument('--sample',         required=True, help='Sample ID')
    parser.add_argument('--align-txt',      required=True, help='align_list.txt (whitelisted barcodes)')
    parser.add_argument('--spatial-counts', required=True,
                        help='spatial_barcode_counts.tsv (spatial barcodes for median metrics)')
    args = parser.parse_args()

    barcodes = [l.strip() for l in open(args.align_txt) if l.strip()]
    with open(args.spatial_counts) as f:
        f.readline()  # skip header
        spatial_barcodes = [l.split('\t')[0] for l in f if l.strip()]

    print(f"Merging {len(barcodes)} barcodes ({len(spatial_barcodes)} spatial)...")
    merge_qc_files(args.input_dir, args.output_dir, args.sample, barcodes, spatial_barcodes)
    print("Done.")


if __name__ == '__main__':
    main()
