#!/usr/bin/env python3
"""
Generate a self-contained HTML QC report embedding all spatial heatmap PNGs as base64.

Inputs:
  --spatial-counts   TSV with per-barcode spatial metrics
  --plots-dir        Directory of spatial_*.png and barcode_rank.png
  --rna-counts       Optional TSV with per-barcode RNA UMI/gene counts
  --rna-summary      Optional STARsolo GeneFull/Summary.csv
  --multiqc-html     Path to multiqc_report.html (linked, not embedded)
  --curated-metrics  Optional curated metrics _mqc.tsv used for pre-filter mapping
  --sample           Sample name
  --output           Output HTML path
"""
import argparse
import base64
import csv
import os

import pandas as pd

from generate_dna_r2_diagnostic_html import (
    DIAGNOSTIC_CSS,
    build_sections as build_diagnostic_sections,
    parse_barcode_stats,
    parse_bbduk,
    parse_kv,
    parse_unmatched,
)


PLOT_ROWS = [
    (
        ('barcode_rank.png', 'DNA barcode rank vs read count'),
        ('rna_barcode_rank.png', 'RNA barcode rank vs UMI count'),
    ),
    (
        ('spatial_dna_reads.png', 'DNA reads per barcode'),
        ('spatial_rna_umis.png', 'RNA UMIs per barcode'),
    ),
    (
        ('spatial_log_dna_reads.png', 'DNA reads per barcode (log10)'),
        ('spatial_log_rna_umis.png', 'RNA UMIs per barcode (log10)'),
    ),
    (
        ('spatial_mapping_rate.png', 'DNA mapping rate'),
        ('spatial_rna_genes.png', 'Detected RNA genes per barcode'),
    ),
    (
        ('spatial_dup_rate.png', 'DNA duplication rate'),
        ('spatial_log_rna_genes.png', 'Detected RNA genes per barcode (log10)'),
    ),
]


def img_to_b64(path):
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode()


def read_curated_metrics(path):
    header = None
    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            if header is None:
                header = line.split('\t')
                continue
            return dict(zip(header, line.split('\t')))
    return {}


def metric_float(metrics, key):
    try:
        return float(metrics.get(key, ''))
    except (TypeError, ValueError):
        return None


def compute_dna_summary(spatial_counts_path, curated_metrics_path=None):
    df = pd.read_csv(spatial_counts_path, sep='\t')
    df['dna_reads']    = pd.to_numeric(df['dna_reads'],    errors='coerce').fillna(0)
    df['mapped_reads'] = pd.to_numeric(df.get('mapped_reads', 0), errors='coerce').fillna(0)
    df['dup_reads']    = pd.to_numeric(df.get('dup_reads', 0),    errors='coerce').fillna(0)

    total_reads   = int(df['dna_reads'].sum())
    n_covered     = int((df['dna_reads'] > 0).sum())
    n_designed    = len(df)
    total_mapped  = int(df['mapped_reads'].sum())
    total_dup     = int(df['dup_reads'].sum())
    med_reads     = float(df.loc[df['dna_reads'] > 0, 'dna_reads'].median()) if n_covered else 0.0
    map_rate      = total_mapped / total_reads if total_reads else 0.0
    dup_rate      = total_dup / total_mapped if total_mapped else 0.0

    summary = {
        'DNA reads in QC BAM': f'{total_reads:,}',
        'Designed spatial barcodes': f'{n_designed:,}',
        'Barcodes with ≥1 read': f'{n_covered:,} ({n_covered/n_designed*100:.1f}%)',
        'Median reads (covered barcodes)': f'{med_reads:.0f}',
        'Overall mapping rate': f'{map_rate*100:.1f}%',
        'Overall duplication rate': f'{dup_rate*100:.1f}%',
    }
    if curated_metrics_path and os.path.exists(curated_metrics_path):
        curated = read_curated_metrics(curated_metrics_path)
        pre_reads = metric_float(curated, 'num_reads_merged')
        pre_mapped = metric_float(curated, 'num_mapped_merged')
        pre_map_rate = metric_float(curated, 'mapping_rate_merged')
        qc_reads = metric_float(curated, 'num_reads_qc_bam')
        bsconv_retained = metric_float(curated, 'bsconv_retained_fraction')
        if pre_reads is not None:
            summary['Pre-filter DNA reads'] = f'{pre_reads:.3f} M'
        if pre_mapped is not None:
            summary['Pre-filter mapped reads'] = f'{pre_mapped:.3f} M'
        if pre_map_rate is not None:
            summary['Overall mapping rate'] = f'{pre_map_rate*100:.1f}%'
        if qc_reads is not None:
            summary['DNA reads in QC BAM'] = f'{qc_reads:.3f} M'
        if bsconv_retained is not None:
            summary['BSCONV retained'] = f'{bsconv_retained*100:.1f}%'
    return summary


def read_star_summary(path):
    metrics = {}
    with open(path, newline='') as handle:
        for row in csv.reader(handle):
            if len(row) >= 2:
                metrics[row[0]] = row[1]
    return metrics


def compute_rna_summary(rna_counts_path, star_summary_path):
    df = pd.read_csv(rna_counts_path, sep='\t')
    df['rna_umis'] = pd.to_numeric(df['rna_umis'], errors='coerce').fillna(0)
    df['rna_genes'] = pd.to_numeric(df['rna_genes'], errors='coerce').fillna(0)
    covered = df[df['rna_umis'] > 0]
    star = read_star_summary(star_summary_path)

    def integer_metric(name):
        try:
            return f"{int(float(star[name])):,}"
        except (KeyError, TypeError, ValueError):
            return 'NA'

    def percent_metric(name):
        try:
            return f"{100 * float(star[name]):.1f}%"
        except (KeyError, TypeError, ValueError):
            return 'NA'

    n_covered = len(covered)
    n_designed = len(df)
    return {
        'Total RNA reads': integer_metric('Number of Reads'),
        'Reads with valid barcodes': percent_metric('Reads With Valid Barcodes'),
        'Unique genome mapping rate': percent_metric('Reads Mapped to Genome: Unique'),
        'Estimated cells': integer_metric('Estimated Number of Cells'),
        'Barcodes with >=1 raw UMI': f'{n_covered:,} ({n_covered/n_designed*100:.1f}%)',
        'Total raw RNA UMIs': f"{int(df['rna_umis'].sum()):,}",
        'Median raw UMIs (covered barcodes)': f"{covered['rna_umis'].median():.0f}" if n_covered else '0',
        'Median raw genes (covered barcodes)': f"{covered['rna_genes'].median():.0f}" if n_covered else '0',
    }


def summary_table(title, summary):
    rows = ''.join(
        f'<tr><td>{key}</td><td>{value}</td></tr>' for key, value in summary.items()
    )
    return (
        f'<section class="summary-table"><h3>{title}</h3><table>'
        f'<tr><th>Metric</th><th>Value</th></tr>{rows}</table></section>'
    )


def plot_figure(plots_dir, fname, caption):
    if not fname:
        return ''
    p = os.path.join(plots_dir, fname)
    if not os.path.exists(p):
        return ''
    b64 = img_to_b64(p)
    return (
        f'<figure><img src="data:image/png;base64,{b64}" alt="{caption}">'
        f'<figcaption>{caption}</figcaption></figure>'
    )


def build_html(sample, dna_summary, plots_dir, multiqc_html, diagnostic_section, rna_summary=None):
    layout_class = '' if rna_summary else ' single'
    summary_html = summary_table('DNA', dna_summary)
    if rna_summary:
        summary_html += summary_table('RNA', rna_summary)

    plot_rows = []
    for dna_plot, rna_plot in PLOT_ROWS:
        figures = [plot_figure(plots_dir, *dna_plot)]
        if rna_summary:
            figures.append(plot_figure(plots_dir, *rna_plot))
        plot_rows.append(f'<div class="plot-row{layout_class}">{"".join(figures)}</div>')

    extra_plots = []
    known_plots = {plot[0] for row in PLOT_ROWS for plot in row}
    for fname in sorted(os.listdir(plots_dir)):
        if fname.endswith('.png') and fname not in known_plots:
            extra_plots.append(plot_figure(plots_dir, fname, os.path.splitext(fname)[0]))

    multiqc_link = ''
    if multiqc_html and os.path.exists(multiqc_html):
        multiqc_link = f'<p><a href="{os.path.basename(multiqc_html)}" target="_blank">Open MultiQC report</a></p>'

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Spatial-DMT QC: {sample}</title>
<style>
  body {{ font-family: sans-serif; max-width: 1200px; margin: auto; padding: 1em; }}
  h1   {{ border-bottom: 2px solid #333; }}
  h2   {{ border-bottom: 1px solid #aaa; }}
  h3   {{ margin: 0 0 0.6em; }}
  table {{ border-collapse: collapse; width: 100%; margin-bottom: 1em; }}
  td, th {{ border: 1px solid #ccc; padding: 6px 12px; }}
  th {{ background: #f0f0f0; }}
  .summary-grid {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 1.5em; }}
  .plot-row {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 1em; margin-bottom: 1em; }}
  .summary-grid.single, .plot-row.single {{ grid-template-columns: 1fr; }}
  .extra-gallery {{ display: flex; flex-wrap: wrap; gap: 1em; }}
  figure {{ margin: 0; text-align: center; min-width: 0; }}
  figure img {{ width: 100%; max-width: 560px; border: 1px solid #ddd; }}
  figcaption {{ font-size: 0.85em; color: #555; }}
{DIAGNOSTIC_CSS}
  @media (max-width: 800px) {{
    .summary-grid, .plot-row {{ grid-template-columns: 1fr; }}
  }}
</style>
</head>
<body>
<h1>Spatial-DMT QC Report: {sample}</h1>

<h2>Summary</h2>
<div class="summary-grid{layout_class}">
{summary_html}
</div>
{multiqc_link}

{diagnostic_section}

<h2>Spatial Maps</h2>
{''.join(plot_rows)}
	<div class="extra-gallery">
	{''.join(extra_plots)}
	</div>
	</body>
	</html>
	"""


def main():
    parser = argparse.ArgumentParser(description='Generate self-contained HTML QC report')
    parser.add_argument('--spatial-counts', required=True)
    parser.add_argument('--plots-dir', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--multiqc-html', default='')
    parser.add_argument('--curated-metrics', default='')
    parser.add_argument('--protocol', choices=['spatialdmt', 'smcseq'], default='spatialdmt')
    parser.add_argument('--trim-stats', required=True)
    parser.add_argument('--unmatched-diagnostics', required=True)
    parser.add_argument('--rna-primer-stats')
    parser.add_argument('--rna-linker1-stats')
    parser.add_argument('--rna-linker2-stats')
    parser.add_argument('--rna-barcode-stats')
    parser.add_argument('--rna-counts')
    parser.add_argument('--rna-summary')
    args = parser.parse_args()

    if bool(args.rna_counts) != bool(args.rna_summary):
        parser.error('--rna-counts and --rna-summary must be supplied together')
    rna_diagnostic_paths = [
        args.rna_primer_stats,
        args.rna_linker1_stats,
        args.rna_linker2_stats,
        args.rna_barcode_stats,
    ]
    if any(rna_diagnostic_paths) and not all(rna_diagnostic_paths):
        parser.error('All RNA diagnostic input arguments must be supplied together')
    dna_summary = compute_dna_summary(args.spatial_counts, args.curated_metrics)
    rna_summary = (
        compute_rna_summary(args.rna_counts, args.rna_summary)
        if args.rna_counts
        else None
    )
    diagnostic_rna = None
    if all(rna_diagnostic_paths):
        diagnostic_rna = (
            parse_bbduk(args.rna_primer_stats),
            parse_bbduk(args.rna_linker1_stats),
            parse_bbduk(args.rna_linker2_stats),
            parse_barcode_stats(args.rna_barcode_stats),
        )
    diagnostic_section = build_diagnostic_sections(
        args.sample,
        parse_kv(args.trim_stats),
        parse_unmatched(args.unmatched_diagnostics),
        diagnostic_rna,
        args.protocol,
    )
    html_doc = build_html(
        args.sample,
        dna_summary,
        args.plots_dir,
        args.multiqc_html,
        diagnostic_section,
        rna_summary,
    )

    with open(args.output, 'w') as f:
        f.write(html_doc)
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
