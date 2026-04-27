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


PLOT_ORDER = [
    ('barcode_rank.png',      'Barcode rank vs DNA read count'),
    ('spatial_dna_reads.png', 'DNA reads per barcode'),
    ('spatial_log_dna_reads.png', 'DNA reads per barcode (log10)'),
    ('spatial_mapping_rate.png',  'Mapping rate'),
    ('spatial_dup_rate.png',      'Duplication rate'),
]


def img_to_b64(path):
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode()


def compute_summary(spatial_counts_path):
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


def build_html(sample, summary, plots_dir, multiqc_html):
    rows = ''.join(
        f'<tr><td>{k}</td><td>{v}</td></tr>' for k, v in summary.items()
    )

    imgs_html = []
    for fname, caption in PLOT_ORDER:
        p = os.path.join(plots_dir, fname)
        if os.path.exists(p):
            b64 = img_to_b64(p)
            imgs_html.append(
                f'<figure><img src="data:image/png;base64,{b64}" alt="{caption}">'
                f'<figcaption>{caption}</figcaption></figure>'
            )

    multiqc_link = ''
    if multiqc_html and os.path.exists(multiqc_html):
        rel = os.path.relpath(multiqc_html, os.path.dirname(multiqc_html))
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
  table {{ border-collapse: collapse; margin-bottom: 1em; }}
  td, th {{ border: 1px solid #ccc; padding: 6px 12px; }}
  th {{ background: #f0f0f0; }}
  .gallery {{ display: flex; flex-wrap: wrap; gap: 1em; }}
  figure {{ margin: 0; text-align: center; }}
  figure img {{ max-width: 500px; border: 1px solid #ddd; }}
  figcaption {{ font-size: 0.85em; color: #555; }}
</style>
</head>
<body>
<h1>Spatial-DMT QC Report: {sample}</h1>

<h2>DNA Summary</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
{rows}
</table>
{multiqc_link}

<h2>Spatial Maps</h2>
<div class="gallery">
{''.join(imgs_html)}
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
    args = parser.parse_args()

    summary = compute_summary(args.spatial_counts)
    html = build_html(args.sample, summary, args.plots_dir, args.multiqc_html)

    with open(args.output, 'w') as f:
        f.write(html)
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
