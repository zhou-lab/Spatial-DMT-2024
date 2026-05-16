#!/usr/bin/env python3
"""Diagnose the Unmatched bucket: alignment stats by category + SF-anchor-fail
breakdown.

spatialmeth_trimtag.py routes everything not heading to a position cell (or to
the spatialdmt lambda bucket) into the Unmatched fastq pair, suffixing each
qname with:
  _SF — structure-fail: R2 grammar couldn't be parsed
  _WM — whitelist-miss: grammar parsed cleanly but BC1+BC2 not in whitelist
        even within Hamming-1.

This script reads the post-align Unmatched BAM (for per-category alignment
stats) AND the original Unmatched R2 fastq (for the structural-anchor probe on
SF reads), and emits a single long-format TSV. The R2 probe re-runs each
protocol anchor independently and reports the FIRST stage that fails for each
SF read -- tells you whether the unmatched bucket is dominated by primer
failures (wrong chemistry / contamination), linker drift, or short-insert
Tn5-readthrough.

Usage:
  unmatched_diagnostics.py --protocol smcseq \\
      --bam      bam_unmatched/{sample}_unmatched.bam \\
      --r2-fastq trim/{sample}_Unmatched_R2.fq.gz \\
      --out      qc/unmatched_diagnostics.tsv
"""
import argparse, gzip, sys
from collections import defaultdict
import pysam
import regex as _regex


## --- Anchor patterns (must match spatialmeth_trimtag.py) -----------------
SMCSEQ = {
    "primer":  _regex.compile(r'(GGTGTAGTGGGTTTGGAGG){s<=3}'),
    "linker1_W": _regex.compile(r'(..CC.C...C.....C.C.C..C...C...){s<=4}'),
    "linker1_C": _regex.compile(r'(..TT.T...T.....T.T.T..T...T...){s<=4}'),
    "linker2_W": _regex.compile(r'(CCC.....C..CC....C...C...CC...){s<=4}'),
    "linker2_C": _regex.compile(r'(TTT.....T..TT....T...T...TT...){s<=4}'),
    "tn5_me":  _regex.compile(r'(AGATGTGT){s<=2}'),
}
SMCSEQ_BC_LEN = 11   # bp per side
SMCSEQ_PRIMER_LEN = 19
SMCSEQ_LINKER_LEN = 30

SPATIALDMT = {
    "linker1": _regex.compile(r'(GTGGTTGATGTTTTGTATTGGTGTATGATT){e<=2}'),
    "linker2": _regex.compile(r'(ATTTATGTGTTTGAGAGGTTAGAGTATTTG){e<=2}'),
    "tn5_me":  _regex.compile(r'(AGATGTGTATAAGAGATAG){e<=1}'),
}


def categorize_sf_smcseq(s):
    """Return the FIRST anchor that fails: 'primer' | 'linker1' | 'linker2' |
    'tn5_me'.  Mirrors parse_r2_smcseq() in spatialmeth_trimtag.py."""
    m = SMCSEQ["primer"].match(s)
    if not m: return "primer"
    pos = m.end() + SMCSEQ_BC_LEN  # past primer + BC1
    m = SMCSEQ["linker1_W"].match(s, pos=pos) or SMCSEQ["linker1_C"].match(s, pos=pos)
    if not m: return "linker1"
    pos = m.end() + SMCSEQ_BC_LEN  # past linker1 + BC2
    m = SMCSEQ["linker2_W"].match(s, pos=pos) or SMCSEQ["linker2_C"].match(s, pos=pos)
    if not m: return "linker2"
    m_te = SMCSEQ["tn5_me"].match(s, pos=m.end())
    if not m_te: return "tn5_me"
    return "all_ok"  # shouldn't happen if read was tagged SF


def categorize_sf_spatialdmt(s):
    """For spatialdmt the trimtag check requires len(lk1)==1, len(lk2)==1,
    len(te2)==1 ALL true. Report which one(s) failed -- bitmask-style."""
    lk1_ok = len(SPATIALDMT["linker1"].findall(s)) == 1
    lk2_ok = len(SPATIALDMT["linker2"].findall(s)) == 1
    te2_ok = len(SPATIALDMT["tn5_me"].findall(s)) == 1
    if not lk1_ok and not lk2_ok and not te2_ok: return "all_three"
    if not lk1_ok: return "linker1"
    if not lk2_ok: return "linker2"
    if not te2_ok: return "tn5_me"
    return "all_ok"  # shouldn't happen


def bam_align_stats(bam_path):
    """Per-suffix alignment stats from the BAM. Counts R1-primary only."""
    counts = defaultdict(lambda: {"n_reads": 0, "n_mapped": 0, "n_dup": 0})
    n_other = 0
    bf = pysam.AlignmentFile(bam_path, "rb")
    for read in bf:
        if read.is_secondary or read.is_supplementary: continue
        if not read.is_read1: continue
        suffix = read.query_name.rsplit("_", 1)[-1]
        if suffix not in ("SF", "WM"):
            n_other += 1
            continue
        b = counts[suffix]
        b["n_reads"] += 1
        if not read.is_unmapped:
            b["n_mapped"] += 1
            if read.is_duplicate:
                b["n_dup"] += 1
    bf.close()
    if n_other:
        sys.stderr.write(f"WARN: {n_other} R1-primary reads had unexpected "
                         f"qname suffix; excluded from align stats\n")
    return counts


def sf_anchor_breakdown(r2_fastq, protocol):
    """For each SF read in the unmatched R2 fastq, classify which structural
    anchor was the FIRST to fail. Returns dict[stage] -> count."""
    categorize = categorize_sf_smcseq if protocol == "smcseq" else categorize_sf_spatialdmt
    counts = defaultdict(int)
    opener = gzip.open if r2_fastq.endswith(".gz") else open
    with opener(r2_fastq, "rt") as f:
        while True:
            h = f.readline()
            if not h: break
            seq = f.readline().rstrip("\n")
            f.readline(); f.readline()
            ## suffix is at the end of qname (before any whitespace)
            qname = h[1:].split(None, 1)[0]
            if qname.rsplit("_", 1)[-1] != "SF":
                continue
            counts[categorize(seq)] += 1
    return counts


def write_mqc_summary(path, sample, align, sf_breakdown, protocol):
    """MultiQC custom-content: single table combining per-category reads +
    mapping_rate AND the SF anchor-fail breakdown (counts) for the sample.
    One row per sample so the report aggregates cleanly across runs."""
    sf = align.get("SF", {"n_reads": 0, "n_mapped": 0})
    wm = align.get("WM", {"n_reads": 0, "n_mapped": 0})
    sf_rate = (sf["n_mapped"] / sf["n_reads"]) if sf["n_reads"] else 0.0
    wm_rate = (wm["n_mapped"] / wm["n_reads"]) if wm["n_reads"] else 0.0
    if protocol == "smcseq":
        anchor_stages = ("primer", "linker1", "linker2", "tn5_me")
    else:
        anchor_stages = ("linker1", "linker2", "tn5_me", "all_three")
    cols = ["sf_n_reads", "sf_mapping_rate", "wm_n_reads", "wm_mapping_rate"] \
         + [f"sf_fail_{s}" for s in anchor_stages]
    vals = [sf["n_reads"], f"{sf_rate:.4f}", wm["n_reads"], f"{wm_rate:.4f}"] \
         + [sf_breakdown.get(s, 0) for s in anchor_stages]
    with open(path, "w") as f:
        f.write("# id: 'unmatched_summary'\n")
        f.write("# section_name: 'Unmatched bucket'\n")
        f.write("# description: 'Reads not routed to a position cell (or to lambda). Columns: SF=structure-fail (R2 grammar parse failed), WM=whitelist-miss (parsed but BC1+BC2 not in whitelist within Hamming-1). sf_fail_* columns count SF reads by the first R2 anchor stage that failed (sequential parse for smcseq, parallel for spatialdmt).'\n")
        f.write("# plot_type: 'table'\n")
        f.write("# pconfig:\n")
        f.write("#     id: 'unmatched_summary_table'\n")
        f.write("#     title: 'Unmatched bucket'\n")
        f.write("Sample\t" + "\t".join(cols) + "\n")
        f.write(f"{sample}\t" + "\t".join(str(v) for v in vals) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam",      required=True)
    ap.add_argument("--r2-fastq", required=True)
    ap.add_argument("--protocol", required=True, choices=["smcseq", "spatialdmt"])
    ap.add_argument("--sample",   required=True, help="sample name for MultiQC")
    ap.add_argument("--out",      required=True, help="long-format raw TSV")
    ap.add_argument("--mqc-summary", required=True,
                    help="MultiQC custom-content TSV: combined unmatched summary")
    args = ap.parse_args()

    align = bam_align_stats(args.bam)
    sf_breakdown = sf_anchor_breakdown(args.r2_fastq, args.protocol)

    with open(args.out, "w") as f:
        f.write("section\tkey\tvalue\n")
        ## align stats per category (long-format raw -- keeps mapped + dup for
        ## the curious; the MultiQC view only surfaces reads + mapping_rate).
        for cat_name, sx in [("structure_fail", "SF"), ("wl_miss", "WM")]:
            c = align.get(sx, {"n_reads": 0, "n_mapped": 0, "n_dup": 0})
            n, m, d = c["n_reads"], c["n_mapped"], c["n_dup"]
            f.write(f"align_stats\t{cat_name}_n_reads\t{n}\n")
            f.write(f"align_stats\t{cat_name}_n_mapped\t{m}\n")
            f.write(f"align_stats\t{cat_name}_mapping_rate\t{m/n if n else 0:.4f}\n")
            f.write(f"align_stats\t{cat_name}_n_dup\t{d}\n")
            f.write(f"align_stats\t{cat_name}_dup_rate\t{d/m if m else 0:.4f}\n")
        total_reads = sum(align[s]["n_reads"] for s in align)
        total_mapped = sum(align[s]["n_mapped"] for s in align)
        f.write(f"align_stats\ttotal_n_reads\t{total_reads}\n")
        f.write(f"align_stats\ttotal_mapping_rate\t{total_mapped/total_reads if total_reads else 0:.4f}\n")
        ## SF anchor-fail breakdown
        total_sf = sum(sf_breakdown.values())
        for stage in ("primer", "linker1", "linker2", "tn5_me", "all_three", "all_ok"):
            n = sf_breakdown.get(stage, 0)
            if stage == "primer" and args.protocol == "spatialdmt":
                continue
            if stage == "all_three" and args.protocol == "smcseq":
                continue
            f.write(f"sf_anchor_fail\t{stage}_fail_n\t{n}\n")
            f.write(f"sf_anchor_fail\t{stage}_fail_frac\t{n/total_sf if total_sf else 0:.4f}\n")
        f.write(f"sf_anchor_fail\ttotal_sf\t{total_sf}\n")

    write_mqc_summary(args.mqc_summary, args.sample, align, sf_breakdown, args.protocol)

    sys.stderr.write(
        f"wrote {args.out}, {args.mqc_summary}\n"
        f"  SF={align['SF']['n_reads']} reads, WM={align['WM']['n_reads']} reads\n"
        f"  SF anchor-fail: {dict(sf_breakdown)}\n")


if __name__ == "__main__":
    main()
