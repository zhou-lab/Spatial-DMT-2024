#!/usr/bin/env python3
"""Trim adapters/linkers, assign spatial barcode, route reads by whitelist hit.

For each structure-passing read, look up its 22-mer (smcseq) or 16-mer
(spatialdmt) barcode in the whitelist (exact, or Hamming-1 with --bc-mismatch).

Whitelist HIT  -> <prefix>_R{1,2}.fq.gz with read name rewritten as
                  `<origid>_<CB>_<UMI>` where:
                    CB  = 8-char ACGT encoding of (X, Y) (see tools/cb_codec.py;
                          ACGT-only is required by dupsifter -B).
                    UMI = 1-bp R2 UMI for smcseq, "N" placeholder for spatialdmt.
                  biscuit align -9 parses these names as second-to-last=BC,
                  last=UMI, and emits CB:Z + RX:Z SAM tags. dupsifter -B
                  dedups per cell on the CB tag. split_bam_by_cb.py decodes
                  the CB back to "XXYY" for the per-cell BAM name.

Whitelist MISS -> <prefix>_Unmatched_R{1,2}.fq.gz with the ORIGINAL read name
                  (no CB tag). These reads went through structural parse (so
                  primer/linker/Tn5 ME all matched) but the BC1+BC2 sequence
                  isn't in the whitelist even within Hamming-1. They're aligned
                  by dna_biscuit_align_unmatched and QC'd separately so the
                  per-cell stats stay clean.

Structure-FAIL reads (primer / linker / Tn5 ME not found) are dropped.

Two protocols (--protocol):
  spatialdmt — 50x50 grid, 8+8 bp barcodes (default)
  smcseq     — 96x96 grid, 11+11 bp + GGTGTAGTGGGTTTGGAGG primer prefix

Outputs (one set per chunk; concatenate across chunks in the calling rule):
  <prefix>_R1.fq.gz / <prefix>_R2.fq.gz             matched, tagged fastqs
  <prefix>_Unmatched_R1.fq.gz / _R2.fq.gz           structure-pass + WL-miss
  <prefix>_Lambda_R1.fq.gz / _R2.fq.gz              spatialdmt only (smcseq skips)
  <prefix>_stats.txt                                read/whitelist counters
  <prefix>_barcodes_all.txt                         observed-barcode frequency
"""
import argparse, gzip, os, subprocess, sys
import regex as _regex
from fuzzysearch import find_near_matches
from cb_codec import coord_to_acgt8


def open_fastq_file(fp):
    return gzip.open(fp, "rt") if fp.endswith(".gz") else open(fp, "r")


def iter_fastq(fh):
    """Yield (rid, seq_str, qual_ascii_str) from a text-mode fastq fh.
    rid is the first whitespace-separated token of the header (matching
    Biopython's SeqRecord.id semantics)."""
    while True:
        h = fh.readline()
        if not h: return
        rid = h[1:].split(None, 1)[0]
        seq = fh.readline().rstrip("\n")
        fh.readline()  # '+'
        qual = fh.readline().rstrip("\n")
        yield rid, seq, qual


class PigzWriter:
    """Write text records to a pigz subprocess (parallel BGZ-compatible gzip)."""
    def __init__(self, out_path, threads=2):
        self._fh = open(out_path, "wb")
        self._proc = subprocess.Popen(
            ["pigz", "-c", "-p", str(threads)],
            stdin=subprocess.PIPE, stdout=self._fh,
        )
    def write(self, s):
        self._proc.stdin.write(s.encode())
    def close(self):
        self._proc.stdin.close()
        self._proc.wait()
        self._fh.close()


def write_fastq_record(h, rid, seq, qual_str):
    h.write("@{}\n{}\n+\n{}\n".format(rid, seq, qual_str))


def load_whitelist_with_coords(path, col):
    """barcode (col `col`) -> (X, Y) from cols 2, 3 (1-indexed)."""
    if not path:
        return None
    wl_map = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= max(col, 3):
                wl_map[parts[col - 1]] = (int(parts[1]), int(parts[2]))
    return wl_map


def assign_barcode_hamming1(bc, wl_map, ambig_policy="drop"):
    """Return (matched_bc, (X, Y)) or (None, None).
    Unambiguous single hits are returned regardless of ambig_policy; only
    multi-hit (genuinely ambiguous) cases are subject to drop-vs-first."""
    if bc in wl_map:
        return bc, wl_map[bc]
    hits = []
    for i in range(len(bc)):
        orig = bc[i]
        for b in "ACGT":
            if b == orig:
                continue
            cand = bc[:i] + b + bc[i+1:]
            if cand in wl_map:
                hits.append(cand)
    if len(hits) == 0:
        return None, None
    if len(hits) == 1:
        return hits[0], wl_map[hits[0]]
    if ambig_policy == "first":
        chosen = sorted(hits)[0]
        return chosen, wl_map[chosen]
    return None, None


parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("read1_fastq")
parser.add_argument("read2_fastq")
parser.add_argument("-o", "--output_prefix", default="out")
parser.add_argument("-a", "--adapters", nargs='+',
    default=["CTATCTCTTATA", "AGATGCGAGAAGCCAACGCTTG"])
parser.add_argument("--protocol", choices=["spatialdmt", "smcseq"], default="spatialdmt",
    help="spatialdmt = 50x50 grid, 8+8 bp barcodes (default); "
         "smcseq = 96x96 grid, 11+11 bp with GGTGTAGTGGGTTTGGAGG primer")
parser.add_argument("-l1", "--linker1", default=None)
parser.add_argument("-l2", "--linker2", default=None)
parser.add_argument("-t",  "--trim_end2", default=None)
parser.add_argument("--min_len", type=int, default=10, help="min length after trimming")
parser.add_argument("--whitelist", required=True, help="whitelist file (cols: bc, X, Y, ...)")
parser.add_argument("--whitelist_col", type=int, default=None,
    help="column with the matchable barcode (auto: 4 for spatialdmt, 1 for smcseq)")
parser.add_argument("--bc-side-len", type=int, default=None,
    help="bp per side; auto: 8 for spatialdmt, 11 for smcseq")
parser.add_argument("--barcode_len", type=int, default=None,
    help="total barcode length; auto: 2 * bc-side-len")
parser.add_argument("--bc-mismatch", type=int, default=1,
    help="allowed barcode mismatches (Hamming); 0 or 1 recommended")
parser.add_argument("--bc-ambig", choices=["drop", "first"], default="drop")

args = parser.parse_args()

## Protocol-specific defaults
if args.protocol == "smcseq":
    ## R2 = [19 primer][11 BC1][30 linker1][11 BC2][30 linker2][8 AGATGTGT][1 UMI][genomic]
    if args.bc_side_len is None: args.bc_side_len = 11
    if args.trim_end2 is None: args.trim_end2 = "AGATGTGT"
    if args.whitelist_col is None: args.whitelist_col = 1
    SMCSEQ_PRIMER  = _regex.compile(r'(GGTGTAGTGGGTTTGGAGG){s<=3}')
    SMCSEQ_LK1_W   = _regex.compile(r'(..CC.C...C.....C.C.C..C...C...){s<=4}')
    SMCSEQ_LK1_C   = _regex.compile(r'(..TT.T...T.....T.T.T..T...T...){s<=4}')
    SMCSEQ_LK2_W   = _regex.compile(r'(CCC.....C..CC....C...C...CC...){s<=4}')
    SMCSEQ_LK2_C   = _regex.compile(r'(TTT.....T..TT....T...T...TT...){s<=4}')
    SMCSEQ_TRIM    = _regex.compile(r'(AGATGTGT){s<=2}')
    SMCSEQ_UMI_LEN = 1
else:
    ## R2 = [8 BC1][30 linker1][8 BC2][30 linker2][19 AGATGTGTATAAGAGATAG][genomic]
    if args.bc_side_len is None: args.bc_side_len = 8
    if args.linker1 is None: args.linker1 = "GTGGTTGATGTTTTGTATTGGTGTATGATT"
    if args.linker2 is None: args.linker2 = "ATTTATGTGTTTGAGAGGTTAGAGTATTTG"
    if args.trim_end2 is None: args.trim_end2 = "AGATGTGTATAAGAGATAG"
    if args.whitelist_col is None: args.whitelist_col = 4

if args.barcode_len is None:
    args.barcode_len = 2 * args.bc_side_len


def parse_r2_smcseq(s):
    """Return (bc1_s, bc1_e, bc2_s, bc2_e, umi_start, genomic_start) or None."""
    m = SMCSEQ_PRIMER.match(s)
    if not m: return None
    bc1_s = m.end()
    bc1_e = bc1_s + args.bc_side_len
    m = SMCSEQ_LK1_W.match(s, pos=bc1_e) or SMCSEQ_LK1_C.match(s, pos=bc1_e)
    if not m: return None
    bc2_s = m.end()
    bc2_e = bc2_s + args.bc_side_len
    m = SMCSEQ_LK2_W.match(s, pos=bc2_e) or SMCSEQ_LK2_C.match(s, pos=bc2_e)
    if not m: return None
    m_te = SMCSEQ_TRIM.match(s, pos=m.end())
    if not m_te: return None
    return (bc1_s, bc1_e, bc2_s, bc2_e, m_te.end(), m_te.end() + SMCSEQ_UMI_LEN)


wl_map = load_whitelist_with_coords(args.whitelist, args.whitelist_col)
if not wl_map:
    sys.exit(f"ERROR: empty whitelist from {args.whitelist} col {args.whitelist_col}")

## MATCHED reads (whitelist hit, exact or Hamming-1) get a real CB tag and go
## to <prefix>_R{1,2}.fq.gz for the per-cell pipeline. Structure-fail (_SF)
## and whitelist-miss (_WM) reads go to <prefix>_Unmatched_R{1,2}.fq.gz with
## a suffix tag in the read name -- aligned separately to bam_unmatched/ and
## summarised by mqc_unmatched_diagnostics. Keeps the merged BAM clean
## (every CB is a real coord; no sentinel cell, no dupsifter over-dedup risk).
r1_out = PigzWriter(args.output_prefix + "_R1.fq.gz")
r2_out = PigzWriter(args.output_prefix + "_R2.fq.gz")
unmatched_r1_out = PigzWriter(args.output_prefix + "_Unmatched_R1.fq.gz")
unmatched_r2_out = PigzWriter(args.output_prefix + "_Unmatched_R2.fq.gz")
## Lambda spike-in detection (barcode-85 motif) only fires for spatialdmt.
## smcseq has no dedicated lambda barcode -- skip creating the empty files.
WRITE_LAMBDA = (args.protocol == "spatialdmt")
if WRITE_LAMBDA:
    lambda_r1_out = PigzWriter(args.output_prefix + "_Lambda_R1.fq.gz")
    lambda_r2_out = PigzWriter(args.output_prefix + "_Lambda_R2.fq.gz")
else:
    lambda_r1_out = lambda_r2_out = None

BARCODE85_MOTIF = "TTATTTTT"

n_reads = n_pass = n_trim = n_trim_bases = 0
n_reads_barcode85 = n_written_barcode85 = 0
n_wl_exact = n_wl_hamming1 = n_wl_unmatched = 0
n_unmatched_sf = n_unmatched_wm = 0
barcode_cnt = {}


def write_unmatched(rid_suffix, rid1, seq1, qual1, rid2, seq2, qual2):
    """Send a read pair to the unmatched fastqs with a category suffix on the
    read name. rid_suffix = 'SF' (structure-fail) or 'WM' (whitelist-miss).
    The mqc_unmatched_diagnostics rule parses the suffix back out."""
    write_fastq_record(unmatched_r1_out, f"{rid1}_{rid_suffix}", seq1, qual1)
    write_fastq_record(unmatched_r2_out, f"{rid2}_{rid_suffix}", seq2, qual2)

with open_fastq_file(args.read1_fastq) as r1h, open_fastq_file(args.read2_fastq) as r2h:
    for (rid1_raw, r1seq, r1qual), (rid2_raw, r2seq, r2qual) in zip(iter_fastq(r1h), iter_fastq(r2h)):
        n_reads += 1
        is_barcode85 = BARCODE85_MOTIF in r1seq or BARCODE85_MOTIF in r2seq
        if is_barcode85:
            n_reads_barcode85 += 1

        ## R1 adapter trim (fuzzysearch unchanged: byte-identical to v0)
        lambda_seq1 = r1seq
        lambda_qual1 = r1qual
        hits = []
        for a in args.adapters:
            hits.extend(find_near_matches(a, r1seq, max_l_dist=1))
        if hits:
            trim_i = min(m.start for m in hits)
            if trim_i < len(lambda_seq1):
                n_trim += 1
                n_trim_bases += (len(lambda_seq1) - trim_i)
                lambda_seq1 = lambda_seq1[:trim_i]
                lambda_qual1 = lambda_qual1[:trim_i]

        ## R2 demux: protocol-specific. Structure-fail reads go to the Unmatched
        ## fastq pair (with _SF suffix on read name) instead of being dropped --
        ## the unmatched alignment + diagnostics rule then bucket on the suffix.
        if args.protocol == "smcseq":
            parsed = parse_r2_smcseq(r2seq)
            if parsed is None:
                write_unmatched("SF", rid1_raw, lambda_seq1, lambda_qual1,
                                      rid2_raw, r2seq, r2qual)
                n_unmatched_sf += 1
                continue
            b1s, l1s, b2s, l2s, umi_start, genomic_start = parsed
            umi_seq = r2seq[umi_start:genomic_start]
            new_seq2 = r2seq[genomic_start:]
            new_qual2 = r2qual[genomic_start:]
        else:
            te2 = find_near_matches(args.trim_end2, r2seq, max_l_dist=1)
            if te2:
                trim2_i = te2[0].end
                lambda_seq2 = r2seq[trim2_i:]
                lambda_qual2 = r2qual[trim2_i:]
            else:
                lambda_seq2 = r2seq
                lambda_qual2 = r2qual
            if is_barcode85 and len(lambda_seq1) >= args.min_len and len(lambda_seq2) >= args.min_len:
                write_fastq_record(lambda_r1_out, rid1_raw, lambda_seq1, lambda_qual1)
                write_fastq_record(lambda_r2_out, rid2_raw, lambda_seq2, lambda_qual2)
                n_written_barcode85 += 1
            lk1 = find_near_matches(args.linker1, r2seq, max_l_dist=2)
            lk2 = find_near_matches(args.linker2, r2seq, max_l_dist=2)
            if not (len(lk1)==1 and len(lk2)==1 and len(te2)==1):
                write_unmatched("SF", rid1_raw, lambda_seq1, lambda_qual1,
                                      rid2_raw, r2seq, r2qual)
                n_unmatched_sf += 1
                continue
            l1s = lk1[0].start; l2s = lk2[0].start
            b1s = max(0, l1s - args.bc_side_len)
            b2s = max(0, l2s - args.bc_side_len)
            umi_seq = "N"
            new_seq2 = lambda_seq2
            new_qual2 = lambda_qual2

        barcode = r2seq[b1s:l1s] + r2seq[b2s:l2s]
        if len(barcode) != args.barcode_len:
            continue
        if len(lambda_seq1) < args.min_len or len(new_seq2) < args.min_len:
            continue

        barcode_cnt[barcode] = barcode_cnt.get(barcode, 0) + 1

        ## Whitelist lookup -> coord on hit. Unmatched reads route to a separate
        ## fastq pair (no CB tag in name) for the dna_biscuit_align_unmatched
        ## rule -- avoids putting fake-cell reads in the main per-cell pipeline.
        coord = None
        if barcode in wl_map:
            coord = wl_map[barcode]
            n_wl_exact += 1
        elif args.bc_mismatch >= 1:
            m, c = assign_barcode_hamming1(barcode, wl_map, args.bc_ambig)
            if m is not None:
                coord = c
                n_wl_hamming1 += 1
            else:
                n_wl_unmatched += 1
        else:
            n_wl_unmatched += 1

        if coord is not None:
            cb_tag = coord_to_acgt8(*coord)
            rid1 = f"{rid1_raw}_{cb_tag}_{umi_seq}"
            rid2 = f"{rid2_raw}_{cb_tag}_{umi_seq}"
            write_fastq_record(r1_out, rid1, lambda_seq1, lambda_qual1)
            write_fastq_record(r2_out, rid2, new_seq2, new_qual2)
            n_pass += 1
        else:
            ## structure-pass + whitelist-miss -> unmatched track with _WM suffix.
            ## Keep the genomic-trim of R1 and R2; no CB tag.
            write_unmatched("WM", rid1_raw, lambda_seq1, lambda_qual1,
                                  rid2_raw, new_seq2, new_qual2)
            n_unmatched_wm += 1


r1_out.close()
r2_out.close()
unmatched_r1_out.close()
unmatched_r2_out.close()
if WRITE_LAMBDA:
    lambda_r1_out.close()
    lambda_r2_out.close()

with open(args.output_prefix + "_stats.txt", "w") as fs:
    fs.write(f"Reads_total\t{n_reads}\n")
    fs.write(f"Reads_detected\t{n_reads - n_reads_barcode85}\n")
    fs.write(f"Reads_passed\t{n_pass}\n")
    fs.write(f"Reads_trimmed\t{n_trim}\n")
    fs.write(f"Bases_trimmed\t{n_trim_bases}\n")
    fs.write(f"Reads_written_barcode85\t{n_written_barcode85}\n")
    fs.write(f"WL_exact\t{n_wl_exact}\n")
    fs.write(f"WL_hamming1\t{n_wl_hamming1}\n")
    fs.write(f"WL_unmatched\t{n_wl_unmatched}\n")
    fs.write(f"Unmatched_SF\t{n_unmatched_sf}\n")
    fs.write(f"Unmatched_WM\t{n_unmatched_wm}\n")
    fs.write(f"Unmatched_written\t{n_unmatched_sf + n_unmatched_wm}\n")

sorted_barcodes = sorted(barcode_cnt.items(), key=lambda x: x[1], reverse=True)
with open(args.output_prefix + "_barcodes_all.txt", "w") as fbc:
    if sorted_barcodes:
        cnt_max = sorted_barcodes[0][1]
        for bc, cnt in sorted_barcodes:
            fbc.write("{}\t{}\t{}\n".format(bc, cnt, "T" if cnt > cnt_max/100 else "F"))
