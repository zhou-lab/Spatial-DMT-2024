#!/usr/bin/env python3
"""
Generate a self-contained HTML QC report embedding all spatial heatmap PNGs as base64.

Inputs:
  --spatial-counts   TSV with per-barcode spatial metrics
  --plots-dir        Directory of spatial_*.png and barcode_rank.png
  --multiqc-html     Path to multiqc_report.html (linked, not embedded)
  --sample           Sample name
  --output           Output HTML path
"""
import argparse
import base64
import os

import pandas as pd

from generate_dna_r2_diagnostic_html import (
    DIAGNOSTIC_CSS,
    build_sections as build_diagnostic_sections,
    parse_kv,
    parse_unmatched,
)


PLOT_ORDER = [
    ('barcode_rank.png', 'DNA barcode rank vs read count'),
    ('spatial_dna_reads.png', 'DNA reads per barcode'),
    ('spatial_log_dna_reads.png', 'DNA reads per barcode (log10)'),
    ('spatial_mapping_rate.png', 'DNA mapping rate'),
    ('spatial_dup_rate.png', 'DNA duplication rate'),
]


def img_to_b64(path):
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode()


def compute_dna_summary(spatial_counts_path):
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

    return {
        'Total DNA reads (demux)': f'{total_reads:,}',
        'Designed spatial barcodes': f'{n_designed:,}',
        'Barcodes with ≥1 read': f'{n_covered:,} ({n_covered/n_designed*100:.1f}%)',
        'Median reads (covered barcodes)': f'{med_reads:.0f}',
        'Overall mapping rate': f'{map_rate*100:.1f}%',
        'Overall duplication rate': f'{dup_rate*100:.1f}%',
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


def build_html(sample, dna_summary, plots_dir, multiqc_html, diagnostic_section):
    layout_class = ' single'
    summary_html = summary_table('DNA', dna_summary)

    plot_rows = []
    for plot in PLOT_ORDER:
        plot_rows.append(f'<div class="plot-row{layout_class}">{plot_figure(plots_dir, *plot)}</div>')

    extra_plots = []
    known_plots = {plot[0] for plot in PLOT_ORDER}
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
    parser.add_argument('--trim-stats', required=True)
    parser.add_argument('--unmatched-diagnostics', required=True)
    args = parser.parse_args()

    dna_summary = compute_dna_summary(args.spatial_counts)
    diagnostic_section = build_diagnostic_sections(
        args.sample,
        parse_kv(args.trim_stats),
        parse_unmatched(args.unmatched_diagnostics),
    )
    html_doc = build_html(
        args.sample,
        dna_summary,
        args.plots_dir,
        args.multiqc_html,
        diagnostic_section,
    )

    with open(args.output, 'w') as f:
        f.write(html_doc)
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
