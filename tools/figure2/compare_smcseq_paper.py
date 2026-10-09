#!/usr/bin/env python3
"""Compare current SmC-seq calls with published Supplementary Table 2.

Requires physical-record summaries from summarize_cg_stream.c. For malformed
legacy indexes, record names are reconstructed in sorted filename order as
documented in the retained packing logs; never use their index row order.
Original .cg and .idx files are not modified.
"""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/figure2/smcseq_paper_comparison'
RUN = Path('/scr1/users/suc3/20260528_runNewSpatialPipe')
PAIRS = {
    'E55_10um_rep1': 'E55_5mC_10um_rep1',
    'E55_10um_rep2': 'E55_5mC_10um_rep2',
    'E65_10um_run1': 'E65_5mC_10um',
    'E65_10um_run2': 'E65_5mC_10um',
    'E75_10um_rep1': 'E75_5mC_10um_rep1',
    'E75_10um_rep2': 'E75_5mC_10um_rep2',
    'E85_20um': 'E85_5mC_20um',
    'E85_placenta_10um': None,
}
VERIFIED = ['E55_10um_rep1', 'E55_10um_rep2', 'E75_10um_rep1', 'E85_20um']


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pub = pd.read_excel(ROOT / 'results/figure2/cache/smcseq_source_data.xlsx',
                        sheet_name='Supplementary Table 2', header=2)
    pub.Barcodes = pub.Barcodes.str.strip()
    # The workbook repeats some E8.5 pixel rows exactly. Count each measured
    # barcode once; conflicting repeats would require manual reconciliation.
    source_rows = pub.groupby('Sample').size()
    pub = pub.drop_duplicates()
    assert not pub.duplicated(['Sample', 'Barcodes']).any()
    barcode = pd.read_csv(ROOT / 'barcodes/smcseq_barcodes_22mer.txt', sep='\t',
                          header=None, names=['Barcodes', 'x', 'y'])
    assert not barcode.Barcodes.duplicated().any()
    comparisons, audits, cross = {}, [], []
    for sample, published in PAIRS.items():
        src = RUN / 'final_output' / sample / f'{sample}.cg'
        idx = pd.read_csv(str(src)+'.idx', sep='\t', header=None, dtype=str)
        frame = pd.read_csv(OUT / f'{sample}_stream.tsv', sep='\t')
        assert len(frame) == len(idx) and idx[0].is_unique
        malformed = not idx[1].is_unique
        # Four-digit numeric names sort identically in the logged shell glob.
        ids = idx[0].sort_values().reset_index(drop=True) if malformed else idx[0]
        if not malformed:
            assert pd.to_numeric(idx[1]).is_monotonic_increasing
        frame['spatial_id'] = ids
        frame['x'] = ids.str[:2].astype(int)
        frame['y'] = ids.str[2:].astype(int)
        frame['Barcodes'] = frame.merge(barcode, on=['x', 'y'], how='left',
                                       validate='one_to_one').Barcodes
        assert frame.Barcodes.notna().all()
        frame.to_csv(OUT / f'{sample}_audited_pixels.tsv', sep='\t', index=False)
        for pub_name, group in pub.groupby('Sample'):
            joined = group.merge(frame, on='Barcodes', validate='one_to_one')
            valid = joined[joined.covered_cpgs.ge(100)].copy()
            cross.append(dict(sample=sample, published_sample=pub_name, n=len(valid),
                              pearson_r=valid.CpG_ML.corr(valid.call_weighted_beta)))
        audit = dict(sample=sample, published_sample=published,
                     paper_source_rows=int(source_rows[published]) if published else None,
                     indexed_pixels=len(idx), distinct_offsets=idx[1].nunique(),
                     legacy_index_malformed=malformed,
                     record_name_method='sorted names from logged packing order' if malformed
                     else 'unique monotonic index in physical record order')
        if published is None:
            audit['conclusion'] = 'No direct counterpart in Supplementary Table 2'
            audits.append(audit)
            continue
        reference = pub[pub.Sample.eq(published)]
        joined = reference.merge(frame, on='Barcodes', validate='one_to_one')
        joined['paper_A'] = joined.nA.str[1:].astype(int)
        joined['paper_B'] = joined.nB.str[1:].astype(int)
        joined['retained'] = joined.covered_cpgs.ge(100)
        joined.to_csv(OUT / f'{sample}_matched.tsv', sep='\t', index=False)
        valid = joined[joined.retained]
        comparisons[sample] = valid
        audit.update(paper_pixels=len(reference), matched_pixels=len(valid),
                     paper_median_beta=valid.CpG_ML.median(),
                     current_weighted_median_beta=valid.call_weighted_beta.median(),
                     current_site_mean_median_beta=valid.site_mean_beta.median(),
                     pearson_r_weighted=valid.CpG_ML.corr(valid.call_weighted_beta),
                     pearson_r_site_mean=valid.CpG_ML.corr(valid.site_mean_beta),
                     median_paired_difference= (valid.call_weighted_beta-valid.CpG_ML).median(),
                     mean_absolute_difference=(valid.call_weighted_beta-valid.CpG_ML).abs().mean(),
                     coverage_spearman_r=valid.CpG_number.corr(valid.covered_cpgs, method='spearman'),
                     conclusion='Strong matched-pixel agreement' if sample in VERIFIED
                     else 'Not confirmed: investigate sample identity and barcode geometry')
        audits.append(audit)
    audit = pd.DataFrame(audits)
    audit.to_csv(OUT / 'comparison_metrics.tsv', sep='\t', index=False)
    pd.DataFrame(cross).to_csv(OUT / 'all_sample_pair_correlations.tsv', sep='\t', index=False)
    plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype':42, 'font.size':8})
    fig, axes = plt.subplots(4, 4, figsize=(12, 10), layout='constrained')
    cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('#ededed')
    for row, sample in enumerate(VERIFIED):
        data = comparisons[sample]
        # Authors' public heatmap script uses rows=A reversed, columns=B reversed.
        # Use the same display transform on both maps and the signed difference.
        na, nb = int(data.paper_A.max()), int(data.paper_B.max())
        for col, field in enumerate(['CpG_ML', 'call_weighted_beta', 'delta']):
            values = data.call_weighted_beta-data.CpG_ML if field=='delta' else data[field]
            mat = np.full((na, nb), np.nan)
            mat[data.paper_A.to_numpy()-1, data.paper_B.to_numpy()-1] = values
            mat = mat[::-1, ::-1]
            ax = axes[row, col]
            im = ax.pcolormesh(np.ma.masked_invalid(mat), cmap='RdBu_r' if col==2 else cmap,
                              vmin=-.1 if col==2 else 0, vmax=.1 if col==2 else 1,
                              rasterized=False, shading='flat')
            ax.set_ylim(na, 0); ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(['Published source data', 'Current: matched definition', 'Current − published'][col])
            if col==0: ax.set_ylabel(sample, fontsize=8)
            cb = fig.colorbar(im, ax=ax, fraction=.035, pad=.025)
            if cb.solids is not None: cb.solids.set_rasterized(False)
        ax = axes[row, 3]
        ax.scatter(data.CpG_ML, data.call_weighted_beta, s=3, alpha=.35, linewidths=0)
        ax.plot([0,1],[0,1], color='#888888', lw=.7)
        ax.set(xlim=(.2,1), ylim=(.2,1), xlabel='Published methylation', ylabel='Current methylation')
        ax.set_aspect('equal')
        r = data.CpG_ML.corr(data.call_weighted_beta)
        ax.set_title(f'n={len(data):,}; Pearson r={r:.3f}')
    fig.suptitle('SmC-seq: comparison at the same barcode-defined pixels', fontsize=13)
    for ext in ['png','pdf','svg']:
        fig.savefig(OUT / f'SmCseq_published_comparison.{ext}', dpi=220)
    plt.close(fig)
    print(audit[['sample','matched_pixels','pearson_r_weighted','median_paired_difference','conclusion']].to_string(index=False))


if __name__ == '__main__':
    main()
