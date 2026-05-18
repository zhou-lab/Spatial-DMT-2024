#!/usr/bin/env python3
"""Emit a sample-level curated-metrics MultiQC custom-content _mqc.tsv.

Mirrors the Ultima WGBS pattern (labpipelines/main/20240320_UltimaWGBS_singleEnd.smk
final_report rule) -- one row per sample with rich tooltip-style headers.

For Spatial-5mC (smcseq + spatialdmt) the row is the sample-aggregate after
cell-aware dedup. Per-cell metrics live in the spatial-map qc_report.html,
not here.

Columns (suffix encodes source so the tooltip is self-documenting):
  num_reads_merged              merged-BAM flagstat "in total"
  num_mapped_merged             merged-BAM flagstat "mapped" (primary)
  mapping_rate_merged           derived
  dup_rate_dupsifter            dupsifter -B stats (cell-aware dup signature)
  num_mapped_lambda             lambda flagstat mapped
  lambda_retention_yame         yame summary lambda .cg -- mean beta over
                                lambda CpGs (smcseq: full-fastq aligned to
                                lambda, so this is near 0 if no spike-in;
                                spatialdmt: from dedicated lambda-barcode
                                reads, near 0 = good bisulfite conversion)
  cpg_coverage_yame             covered CpGs (M+U > 0) on {sample}.cg
  cpg_depth_yame                mean M+U over covered CpGs
  CA_retention_biscuitqc        mean beta in CA context -- residual non-CpG;
                                near 0 = complete bisulfite conversion
  CC_retention_biscuitqc        mean beta in CC context -- residual non-CpG
  CG_retention_biscuitqc        mean beta in CG context -- the genuine CpG
                                methylation signal (not residual!)
  CT_retention_biscuitqc        mean beta in CT context -- residual non-CpG

The CA/CC/CG/CT columns all come from the same per-cell VCF aggregation done
in dna_biscuit_pileup (vcf2bed -e -t c streamed to awk; no second pileup).

Usage:
    python write_curated_metrics_mqc.py \\
        --sample E65_20um \\
        --merged-flagstat ... \\
        --merged-dupstat ... \\
        --biscuit-conv ... \\
        --cg ... \\
        --lambda-flagstat ... \\
        --lambda-cg ... \\
        --output {sample}_curated_metrics_mqc.tsv
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

HEADER = """\
# id: curated_metrics
# section_name: 'Curated metrics'
# description: 'Per-sample summary built directly from QC source files. Column-name suffixes name the source (_merged = merged BAM flagstat, _dupsifter = dupsifter.stat, _yame = yame summary on .cg, _biscuitqc = vcf2bed -e -t c aggregated over per-cell pileup VCFs, _lambda = lambda alignment). Hover any header for the exact source and semantics. CA/CC/CG/CT are all derived from the SAME per-cell pileup VCFs (no second sample-level pileup) -- CG is the genuine CpG methylation signal, CA/CC/CT are residual non-CpG retention (~0 = complete bisulfite conversion).'
# plot_type: 'table'
# format: 'tsv'
# headers:
#     num_reads_merged:
#         title: 'Reads (merged)'
#         description: 'Total reads in the merged dedup BAM (samtools flagstat "in total"). Includes both R1 and R2 of paired reads.'
#         format: '{:,.1f}'
#         suffix: ' M'
#     num_mapped_merged:
#         title: 'Mapped'
#         description: 'Primary mapped reads in the merged dedup BAM (samtools flagstat "mapped"). Excludes secondary/supplementary.'
#         format: '{:,.1f}'
#         suffix: ' M'
#     mapping_rate_merged:
#         title: 'Mapping rate'
#         description: 'num_mapped_merged / num_reads_merged.'
#         format: '{:.1%}'
#     dup_rate_dupsifter:
#         title: 'Duplication rate'
#         description: 'Cell-aware (dupsifter -B) duplication rate from {sample}_merged_dedup.dupsifter.stat: (fwd_dup + rev_dup) / total_reads_processed. Reflects PCR duplicates collapsed per (CB, position) signature.'
#         format: '{:.1%}'
#     num_mapped_lambda:
#         title: 'Mapped (λ)'
#         description: 'Primary mapped reads in the lambda BAM (samtools flagstat). For smcseq the full tagged fastq is aligned to lambda; for spatialdmt only barcode-85 reads are.'
#         format: '{:,.1f}'
#         suffix: ' M'
#     lambda_retention_yame:
#         title: 'λ retention'
#         description: 'Mean beta over lambda CpGs (yame summary on {sample}_Lambda.cg). Near 0 = good bisulfite conversion (unmethylated lambda fully converted); high values suggest under-conversion or contamination.'
#         format: '{:.3f}'
#     cpg_coverage_yame:
#         title: 'CpG coverage'
#         description: 'Number of reference CpGs with M+U > 0 from yame summary on {sample}.cg (sample-aggregate across all cells in the merged dedup BAM). Universe = the CpG set the .cg was packed against.'
#         format: '{:,.1f}'
#         suffix: ' M'
#     cpg_depth_yame:
#         title: 'CpG depth'
#         description: 'Mean M+U over covered CpGs from yame summary on {sample}.cg. Sample-aggregate.'
#         format: '{:.2f}'
#         suffix: 'x'
#     CA_retention_biscuitqc:
#         title: 'CA retention'
#         description: 'Mean beta in CA context, record-pooled across per-cell pileup VCFs (vcf2bed -e -t c -> awk in dna_biscuit_pileup). Residual non-CpG methylation; near 0 = complete bisulfite conversion.'
#         format: '{:.2e}'
#     CC_retention_biscuitqc:
#         title: 'CC retention'
#         description: 'Mean beta in CC context, record-pooled across per-cell pileup VCFs. Residual non-CpG methylation; near 0 = complete conversion.'
#         format: '{:.2e}'
#     CG_retention_biscuitqc:
#         title: 'CG methylation'
#         description: 'Mean beta in CG context, record-pooled across per-cell pileup VCFs. This is the genuine CpG methylation signal -- NOT a residual. Typical mammal values 0.6-0.8 genome-wide.'
#         format: '{:.3f}'
#     CT_retention_biscuitqc:
#         title: 'CT retention'
#         description: 'Mean beta in CT context, record-pooled across per-cell pileup VCFs. Residual non-CpG methylation; near 0 = complete conversion.'
#         format: '{:.2e}'
"""


def parse_flagstat(path):
    """Return (total_reads, mapped_primary) from samtools flagstat output."""
    total = mapped = 0
    with open(path) as f:
        for line in f:
            m = re.match(r'^(\d+)\s+\+\s+\d+\s+(.*)', line)
            if not m:
                continue
            n, label = int(m.group(1)), m.group(2)
            if label.startswith('in total'):
                total = n
            elif label.startswith('mapped (') or label == 'mapped':
                ## first "mapped" line (primary); flagstat lists primary mapped first
                if mapped == 0:
                    mapped = n
    return total, mapped


def parse_dupsifter(path):
    """Return dup_rate from dupsifter -B stats: (fwd_dup + rev_dup) / total."""
    total = dupf = dupr = 0
    with open(path) as f:
        for line in f:
            if ':' not in line:
                continue
            key, _, val = line.partition(':')
            try:
                v = int(val.strip())
            except ValueError:
                continue
            k = key.strip().lower()
            if 'individual reads processed' in k:
                total = v
            elif 'forward strand marked as duplicates' in k:
                dupf = v
            elif 'reverse strand marked as duplicates' in k:
                dupr = v
    return (dupf + dupr) / total if total else 0.0


def parse_biscuit_conv(path):
    """Return (CA, CC, CG, CT) from BISCUITqc totalBaseConversionRate table.

    File layout (from dna_biscuit_pileup):
        BISCUITqc Conversion Rate by Base Average Table
        CA<TAB>CC<TAB>CG<TAB>CT
        <ca><TAB><cc><TAB><cg><TAB><ct>      (-1 if <20 records for context)
    """
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f if l.strip()]
    cols = lines[-2].split('\t')
    vals = lines[-1].split('\t')
    d = dict(zip(cols, vals))
    def f(x):
        try: return float(x)
        except (TypeError, ValueError): return float('nan')
    return f(d.get('CA', 'nan')), f(d.get('CC', 'nan')), \
           f(d.get('CG', 'nan')), f(d.get('CT', 'nan'))


def yame_summary(cg_path):
    """Return (n_query, depth, beta) from `yame summary {cg}` first data row.

    Header columns (yame summary, 2026): QFile, Query, MFile, Mask, N_univ,
    N_query, N_mask, N_overlap, Log2OddsRatio, Beta, Depth. Parse by name to
    survive future yame column re-ordering.
    """
    out = subprocess.check_output(['yame', 'summary', cg_path], text=True)
    lines = [l for l in out.splitlines() if l.strip()]
    if len(lines) < 2:
        return float('nan'), float('nan'), float('nan')
    hdr = lines[0].split('\t')
    row = lines[1].split('\t')
    d = dict(zip(hdr, row))
    def f(k):
        try: return float(d.get(k, 'nan'))
        except (TypeError, ValueError): return float('nan')
    return f('N_query'), f('Depth'), f('Beta')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sample',          required=True)
    p.add_argument('--merged-flagstat', required=True)
    p.add_argument('--merged-dupstat',  required=True)
    p.add_argument('--biscuit-conv',    required=True)
    p.add_argument('--cg',              required=True)
    p.add_argument('--lambda-flagstat', required=True)
    p.add_argument('--lambda-cg',       required=True)
    p.add_argument('--output',          required=True)
    args = p.parse_args()

    total, mapped       = parse_flagstat(args.merged_flagstat)
    dup_rate            = parse_dupsifter(args.merged_dupstat)
    ca, cc, cg, ct      = parse_biscuit_conv(args.biscuit_conv)
    lam_total, lam_map  = parse_flagstat(args.lambda_flagstat)
    nq,  depth, _       = yame_summary(args.cg)
    _,   _,     lam_beta = yame_summary(args.lambda_cg)

    map_rate = (mapped / total) if total else 0.0

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        f.write(HEADER)
        f.write('Sample\tnum_reads_merged\tnum_mapped_merged\tmapping_rate_merged\t'
                'dup_rate_dupsifter\tnum_mapped_lambda\tlambda_retention_yame\t'
                'cpg_coverage_yame\tcpg_depth_yame\t'
                'CA_retention_biscuitqc\tCC_retention_biscuitqc\t'
                'CG_retention_biscuitqc\tCT_retention_biscuitqc\n')
        f.write(f'{args.sample}\t{total/1e6:.3f}\t{mapped/1e6:.3f}\t{map_rate:.6f}\t'
                f'{dup_rate:.6f}\t{lam_map/1e6:.3f}\t{lam_beta:.6f}\t'
                f'{nq/1e6:.3f}\t{depth:.3f}\t'
                f'{ca:.6e}\t{cc:.6e}\t{cg:.6f}\t{ct:.6e}\n')


if __name__ == '__main__':
    main()
