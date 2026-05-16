#!/usr/bin/env python3
"""Bucket the unmatched BAM by qname suffix to report per-category alignment
stats. spatialmeth_trimtag.py appends `_SF` (structure-fail) or `_WM`
(whitelist-miss) to every read going into the Unmatched fastq pair; this
script reads the post-align BAM, parses the suffix back out, and emits a
single TSV summarizing per-category counts.

Categories:
  SF — structure-fail: R2 primer/linker/Tn5-ME couldn't be located. R2 still
       contains the structural prefix, so most of these are expected to be
       unmapped.
  WM — whitelist-miss: structure parsed cleanly but the 16/22-mer barcode
       isn't in the whitelist even within Hamming-1. R2 was trimmed to the
       genomic start, so mapping rate should be comparable to matched reads.

Each row counts R1-primary reads (one row per pair). Secondary/supplementary
alignments are skipped.
"""
import argparse, sys
from collections import defaultdict
import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("in_bam")
    ap.add_argument("out_tsv")
    args = ap.parse_args()

    counts = defaultdict(lambda: {"n_reads": 0, "n_mapped": 0, "n_dup": 0})
    n_other_suffix = 0

    bf = pysam.AlignmentFile(args.in_bam, "rb")
    for read in bf:
        if read.is_secondary or read.is_supplementary:
            continue
        if not read.is_read1:
            continue
        suffix = read.query_name.rsplit("_", 1)[-1]
        if suffix not in ("SF", "WM"):
            n_other_suffix += 1
            continue
        bucket = counts[suffix]
        bucket["n_reads"] += 1
        if not read.is_unmapped:
            bucket["n_mapped"] += 1
            if read.is_duplicate:
                bucket["n_dup"] += 1
    bf.close()

    if n_other_suffix:
        sys.stderr.write(f"WARN: {n_other_suffix} reads had unexpected qname "
                         f"suffix (not _SF or _WM); excluded from breakdown\n")

    rows = []
    for cat_name, cat in [("structure_fail", "SF"), ("wl_miss", "WM")]:
        c = counts.get(cat, {"n_reads": 0, "n_mapped": 0, "n_dup": 0})
        rows.append((cat_name, c["n_reads"], c["n_mapped"], c["n_dup"]))
    total_reads = sum(r[1] for r in rows)
    total_mapped = sum(r[2] for r in rows)
    total_dup = sum(r[3] for r in rows)
    rows.append(("total", total_reads, total_mapped, total_dup))

    with open(args.out_tsv, "w") as f:
        f.write("category\tn_reads\tn_mapped\tmapping_rate\tn_dup\tdup_rate\n")
        for name, n, m, d in rows:
            map_rate = m / n if n else 0.0
            dup_rate = d / m if m else 0.0
            f.write(f"{name}\t{n}\t{m}\t{map_rate:.4f}\t{d}\t{dup_rate:.4f}\n")
    sys.stderr.write(f"wrote {args.out_tsv} ({total_reads} R1-primary reads)\n")


if __name__ == "__main__":
    main()
