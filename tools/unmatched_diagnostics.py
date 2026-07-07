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
SF reads), and emits a single long-format TSV. The R2 probe follows the
protocol grammar in order and reports the FIRST stage that fails for each SF
read -- tells you whether the unmatched bucket is dominated by primer
failures (wrong chemistry / contamination), linker drift, or short-insert
Tn5-readthrough.

Usage:
  unmatched_diagnostics.py --protocol smcseq \\
      --bam      bam_unmatched/{sample}_unmatched.bam \\
      --r2-fastq trim/{sample}_Unmatched_R2.fq.gz \\
      --out      qc/unmatched_diagnostics.tsv
"""
import argparse, gzip, multiprocessing as mp, sys
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
    """Return the first failed Spatial-DMT grammar stage.

    R2 = [8 BC1][30 linker1][8 BC2][30 linker2][19 Tn5 ME][genomic].
    A downstream stage is evaluated only after every preceding stage passed.
    Match anchors using the same unique-anywhere criterion as trimtag.
    """
    if len(SPATIALDMT["linker1"].findall(s)) != 1:
        return "linker1"
    if len(SPATIALDMT["linker2"].findall(s)) != 1:
        return "linker2"
    if len(SPATIALDMT["tn5_me"].findall(s)) != 1:
        return "tn5_me"
    return "all_ok"  # SF caused by a non-grammar constraint or legacy output


def bam_align_stats(bam_path, threads=1):
    """Per-suffix alignment stats from the BAM. Counts R1-primary only."""
    counts = defaultdict(lambda: {"n_reads": 0, "n_mapped": 0, "n_dup": 0})
    n_other = 0
    bf = pysam.AlignmentFile(bam_path, "rb", threads=max(1, threads))
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


def _merge_counts(dst, src):
    for k, v in src.items():
        dst[k] += v


def _classify_sf_batch(args):
    protocol, seqs = args
    counts = defaultdict(int)
    overlaps = defaultdict(int)
    if protocol == "smcseq":
        for seq in seqs:
            counts[categorize_sf_smcseq(seq)] += 1
        return dict(counts), dict(overlaps)

    for seq in seqs:
        match_counts = {
            name: len(SPATIALDMT[name].findall(seq))
            for name in ("linker1", "linker2", "tn5_me")
        }
        if all(n == 0 for n in match_counts.values()):
            overlaps["all_three_absent"] += 1
        for name in ("linker1", "linker2", "tn5_me"):
            if match_counts[name] != 1:
                counts[name] += 1
                break
        else:
            counts["all_ok"] += 1
    return dict(counts), dict(overlaps)


def _sf_seq_batches(r2_fastq, batch_size):
    opener = gzip.open if r2_fastq.endswith(".gz") else open
    batch = []
    with opener(r2_fastq, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().rstrip("\n")
            f.readline(); f.readline()
            ## suffix is at the end of qname (before any whitespace)
            qname = h[1:].split(None, 1)[0]
            if qname.rsplit("_", 1)[-1] != "SF":
                continue
            batch.append(seq)
            if len(batch) >= batch_size:
                yield batch
                batch = []
    if batch:
        yield batch


def sf_anchor_breakdown(r2_fastq, protocol, threads=1, batch_size=100000):
    """For each SF read in the unmatched R2 fastq, classify which structural
    anchor was the FIRST to fail. Also count Spatial-DMT reads where all three
    anchors have zero matches. Returns (first_failure_counts, overlap_counts)."""
    counts = defaultdict(int)
    overlaps = defaultdict(int)

    threads = max(1, threads)
    batches = ((protocol, batch) for batch in _sf_seq_batches(r2_fastq, batch_size))
    pool = None
    if threads == 1:
        results = map(_classify_sf_batch, batches)
    else:
        pool = mp.Pool(processes=threads)
        results = pool.imap_unordered(_classify_sf_batch, batches, chunksize=1)

    try:
        for batch_counts, batch_overlaps in results:
            _merge_counts(counts, batch_counts)
            _merge_counts(overlaps, batch_overlaps)
    except Exception:
        if pool is not None:
            pool.terminate()
            pool.join()
            pool = None
        raise
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    return counts, overlaps


## Per-column MultiQC headers. `title` = short column label shown in the table
## header; `description` = the longer tooltip shown on hover. Keep these
## paragraph-length -- this is the only in-report documentation a user gets
## when they wonder what a column means.
_ANCHOR_HOVER = {
    "primer": (
        "smcseq R2 starts with a 19-bp GGTGTAGTGGGTTTGGAGG primer. If the "
        "primer can't be located at R2 position 0 within 3 substitutions, "
        "parsing aborts here. High counts here usually mean the wrong "
        "chemistry was sequenced, the library is contaminated, or the cycle-1 "
        "base call is broken."),
    "linker1": (
        "After the primer + 11-bp BC1, the parser expects a 30-bp linker1 "
        "(Watson C-preserved or Crick T-preserved variant) within 4 "
        "substitutions. Failure here usually means BC1 has insertions/deletions "
        "that shifted the linker out of position, or the linker chemistry is "
        "drifting."),
    "linker2": (
        "After linker1 + 11-bp BC2, expects 30-bp linker2 (same C/T duality). "
        "Failure here usually means BC2 has indels."),
    "tn5_me": (
        "Final anchor: 8-bp AGATGTGT (Tn5 ME) after linker2. Failure here "
        "usually means the insert is shorter than ~110 bp, so R2 reads past "
        "the genomic end into adapter -- the Tn5 ME isn't in R2 anymore. "
        "Common for size-selected libraries; not a quality problem."),
}


def _yq(s):
    """YAML single-quoted scalar with apostrophe-escaping (doubled '')."""
    return "'" + s.replace("'", "''") + "'"


def _format_headers(anchor_stages):
    """Build the `headers:` YAML block for MultiQC custom-content table."""
    lines = ["# headers:"]
    def col(name, title, desc, fmt=None, min_=None, max_=None, scale=None, suffix=None):
        lines.append(f"#     {name}:")
        lines.append(f"#         title: {_yq(title)}")
        lines.append(f"#         description: {_yq(desc)}")
        if fmt    is not None: lines.append(f"#         format: {_yq(fmt)}")
        if min_   is not None: lines.append(f"#         min: {min_}")
        if max_   is not None: lines.append(f"#         max: {max_}")
        if scale  is not None: lines.append(f"#         scale: {_yq(scale)}")
        if suffix is not None: lines.append(f"#         suffix: {_yq(suffix)}")
    col("sf_n_reads", "SF reads",
        "Number of R1-primary reads in the structure-fail bucket. These are "
        "reads whose R2 sequence didn't pass the protocol's structural-anchor "
        "parse (see sf_fail_* columns for which anchor failed). Their R2 still "
        "carries the structural prefix, so mapping rate against the genome is "
        "expected to be low.",
        fmt="{:,.0f}", min_=0, scale="Blues")
    col("sf_mapping_rate", "SF map rate",
        "Fraction of SF reads that the aligner placed against the main genome "
        "(R1-primary, any quality). Typically near zero because the unmapped "
        "structural prefix dominates R2. Elevated rate here would suggest the "
        "structural parse is too strict and is dropping real library reads.",
        fmt="{:.2%}", min_=0, max_=1, scale="Reds")
    col("wm_n_reads", "WM reads",
        "Number of R1-primary reads in the whitelist-miss bucket. Structure "
        "parsed cleanly, but the 16/22-mer BC1+BC2 isn't in the whitelist "
        "even within Hamming-1. R2 was genomic-trimmed before alignment, so "
        "mapping rate should be comparable to matched cells.",
        fmt="{:,.0f}", min_=0, scale="Blues")
    col("wm_mapping_rate", "WM map rate",
        "Fraction of WM reads that mapped. If close to the matched-cell rate, "
        "these are real library reads with sequencing errors >1 bp in the "
        "barcode -- could be salvageable with a more permissive whitelist "
        "(Hamming-2). If much lower, they're more likely off-target or "
        "contamination.",
        fmt="{:.2%}", min_=0, max_=1, scale="Reds")
    anchor_descriptions = {s: _ANCHOR_HOVER[s] for s in anchor_stages}
    for s in anchor_stages:
        col(f"sf_fail_{s}",
            {"primer":"SF: primer", "linker1":"SF: linker1",
             "linker2":"SF: linker2", "tn5_me":"SF: Tn5 ME"}[s],
            f"Number of SF reads whose first failing anchor was {s}. "
            f"{anchor_descriptions[s]}",
            fmt="{:,.0f}", min_=0, scale="Oranges")
    return "\n".join(lines)


def write_mqc_summary(path, sample, align, sf_breakdown, protocol):
    """MultiQC custom-content: single table combining per-category reads +
    mapping_rate AND the SF anchor-fail breakdown (counts) for the sample.
    One row per sample so the report aggregates cleanly across runs.
    Per-column hover descriptions explain each metric in-report."""
    sf = align.get("SF", {"n_reads": 0, "n_mapped": 0})
    wm = align.get("WM", {"n_reads": 0, "n_mapped": 0})
    sf_rate = (sf["n_mapped"] / sf["n_reads"]) if sf["n_reads"] else 0.0
    wm_rate = (wm["n_mapped"] / wm["n_reads"]) if wm["n_reads"] else 0.0
    if protocol == "smcseq":
        anchor_stages = ("primer", "linker1", "linker2", "tn5_me")
    else:
        anchor_stages = ("linker1", "linker2", "tn5_me")
    cols = ["sf_n_reads", "sf_mapping_rate", "wm_n_reads", "wm_mapping_rate"] \
         + [f"sf_fail_{s}" for s in anchor_stages]
    vals = [sf["n_reads"], f"{sf_rate:.4f}", wm["n_reads"], f"{wm_rate:.4f}"] \
         + [sf_breakdown.get(s, 0) for s in anchor_stages]
    section_description = (
        "Per-sample summary of the Unmatched track -- the bucket of read pairs "
        "that didn't make it onto a spatial-position cell or into the lambda "
        "spike-in. Two source paths: <strong>SF</strong> (structure-fail, "
        "the R2 protocol grammar -- primer + BC + linker + ME -- couldn't be "
        "parsed) and <strong>WM</strong> (whitelist-miss, the grammar parsed "
        "cleanly but BC1+BC2 was not in the spatial whitelist even within "
        "Hamming-1). The sf_fail_* columns count SF reads by the FIRST anchor "
        "stage that failed during the sequential parse, telling you whether "
        "the unmatched bucket is primer/linker noise or just short-insert "
        "fragments where R2 reads past the genomic end.")
    with open(path, "w") as f:
        f.write(f"# id: {_yq('unmatched_summary')}\n")
        f.write(f"# section_name: {_yq('Unmatched bucket')}\n")
        f.write(f"# description: {_yq(section_description)}\n")
        f.write(f"# plot_type: {_yq('table')}\n")
        f.write("# pconfig:\n")
        f.write(f"#     id: {_yq('unmatched_summary_table')}\n")
        f.write(f"#     title: {_yq('Unmatched bucket: read disposition and SF anchor-fail breakdown')}\n")
        f.write(_format_headers(anchor_stages) + "\n")
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
    ap.add_argument("--threads", type=int, default=1,
                    help="worker processes for SF anchor classification")
    ap.add_argument("--batch-size", type=int, default=100000,
                    help="SF reads per worker batch")
    args = ap.parse_args()

    align = bam_align_stats(args.bam, threads=args.threads)
    sf_breakdown, sf_overlaps = sf_anchor_breakdown(
        args.r2_fastq, args.protocol, threads=args.threads,
        batch_size=args.batch_size)

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
        for stage in ("primer", "linker1", "linker2", "tn5_me", "all_ok"):
            n = sf_breakdown.get(stage, 0)
            if stage == "primer" and args.protocol == "spatialdmt":
                continue
            f.write(f"sf_anchor_fail\t{stage}_fail_n\t{n}\n")
            f.write(f"sf_anchor_fail\t{stage}_fail_frac\t{n/total_sf if total_sf else 0:.4f}\n")
        f.write(f"sf_anchor_fail\ttotal_sf\t{total_sf}\n")
        if args.protocol == "spatialdmt":
            all_three_absent = sf_overlaps.get("all_three_absent", 0)
            f.write(f"sf_anchor_overlap\tall_three_absent_n\t{all_three_absent}\n")
            f.write(
                "sf_anchor_overlap\tall_three_absent_frac\t"
                f"{all_three_absent/total_sf if total_sf else 0:.4f}\n"
            )

    write_mqc_summary(args.mqc_summary, args.sample, align, sf_breakdown, args.protocol)

    sys.stderr.write(
        f"wrote {args.out}, {args.mqc_summary}\n"
        f"  SF={align['SF']['n_reads']} reads, WM={align['WM']['n_reads']} reads\n"
        f"  SF anchor-fail: {dict(sf_breakdown)}\n")


if __name__ == "__main__":
    main()
