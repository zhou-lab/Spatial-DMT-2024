#!/usr/bin/env python3
import argparse, gzip, os, sys
import regex as _regex
from fuzzysearch import find_near_matches
from Bio import SeqIO
from Bio.Seq import Seq

def open_fastq_file(fp):
    return gzip.open(fp, "rt") if fp.endswith(".gz") else open(fp, "r")

def write_fastq_record(h, rid, seq, qual_ints):
    h.write("@{}\n".format(rid))
    h.write(str(seq) + "\n")
    h.write("+\n")
    h.write("".join(chr(q + 33) for q in qual_ints) + "\n")

def load_whitelist(path, col=4):
    if not path: return None
    wl = set()
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= col:
                wl.add(parts[col-1])
    return wl

parser = argparse.ArgumentParser(
    description="Trim adapters/linkers and (optionally) demultiplex in one pass."
)
parser.add_argument("read1_fastq")
parser.add_argument("read2_fastq")
parser.add_argument("-o","--output_prefix", default="out")
parser.add_argument("-a","--adapters", nargs='+',
    default=["CTATCTCTTATA","AGATGCGAGAAGCCAACGCTTG"])
parser.add_argument("--protocol", choices=["spatialdmt","smcseq"], default="spatialdmt",
    help="spatialdmt = Fan-lab 50x50 8+8bp barcodes (default); "
         "smcseq = Liu-lab 96x96 11+11bp barcodes with GGTGTAGTGGGTTTGGAGG primer")
parser.add_argument("-l1","--linker1", default=None)
parser.add_argument("-l2","--linker2", default=None)
parser.add_argument("-t","--trim_end2", default=None)
parser.add_argument("--min_len", type=int, default=10, help="min length after trimming")
parser.add_argument("--whitelist", help="whitelist file (barcode in column --whitelist_col)")
parser.add_argument("--whitelist_col", type=int, default=None)
parser.add_argument("--dmux_dir", help="dir to write per-barcode FASTQs for this chunk")
parser.add_argument("--chunk_id", default="", help="chunk ID suffix")
parser.add_argument("--bc-side-len", type=int, default=None,
    help="bp per side barcode; auto: 8 for spatialdmt, 11 for smcseq")
parser.add_argument("--barcode_len", type=int, default=None,
    help="total barcode length; auto: 2 * bc-side-len")
parser.add_argument("--bc-mismatch", type=int, default=0,
                    help="allowed barcode mismatches (Hamming); 0 or 1 recommended")
parser.add_argument("--bc-ambig", choices=["drop","first"], default="drop",
                    help="if multiple 1-mismatch matches: drop or pick first")

args = parser.parse_args()

## Protocol-specific defaults
if args.protocol == "smcseq":
    ## Liu-lab spatial 5mC: [19bp primer][11bp BC1][30bp linker1][11bp BC2][30bp linker2][8bp AGATGTGT][1bp UMI][genomic]
    ## linker1/linker2 each have two variants (Watson C-preserved, Crick C->T) and we try both per read
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
    ## spatialdmt (Fan-lab DBiT chemistry):
    ##   R2 = [8bp BC1][30bp linker1][8bp BC2][30bp linker2][19bp AGATGTGTATAAGAGATAG][genomic]
    if args.bc_side_len is None: args.bc_side_len = 8
    if args.linker1 is None: args.linker1 = "GTGGTTGATGTTTTGTATTGGTGTATGATT"
    if args.linker2 is None: args.linker2 = "ATTTATGTGTTTGAGAGGTTAGAGTATTTG"
    if args.trim_end2 is None: args.trim_end2 = "AGATGTGTATAAGAGATAG"
    if args.whitelist_col is None: args.whitelist_col = 4

if args.barcode_len is None:
    args.barcode_len = 2 * args.bc_side_len


def parse_r2_smcseq(s):
    """Return (bc1_s, bc1_e, bc2_s, bc2_e, genomic_start) or None.
    Anchors at R2 start; tries both Watson and Crick linker variants per pair."""
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
    return (bc1_s, bc1_e, bc2_s, bc2_e, m_te.end() + SMCSEQ_UMI_LEN)

dmux_handles_R1 = {}
dmux_handles_R2 = {}

OPEN_HANDLE_CAP = 7000

def get_dmux_handles(barcode):
    if barcode not in dmux_handles_R1:
        # If you want LRU safety, call _maybe_evict_lru() here
        r1p = os.path.join(args.dmux_dir, f"{barcode}_R1.{args.chunk_id}.fq")
        r2p = os.path.join(args.dmux_dir, f"{barcode}_R2.{args.chunk_id}.fq")
        # Use large buffers to cut syscalls (1 MiB)
        dmux_handles_R1[barcode] = open(r1p, "a", buffering=1024*1024)
        dmux_handles_R2[barcode] = open(r2p, "a", buffering=1024*1024)
    return dmux_handles_R1[barcode], dmux_handles_R2[barcode]

wl = load_whitelist(args.whitelist, args.whitelist_col) if args.whitelist else None
do_dmux = wl is not None and args.dmux_dir is not None
if do_dmux:
    os.makedirs(args.dmux_dir, exist_ok=True)

r1_out = gzip.open(args.output_prefix + "_R1.fq.gz", "wt")
r2_out = gzip.open(args.output_prefix + "_R2.fq.gz", "wt")
lambda_r1_out = gzip.open(args.output_prefix + "_Lambda_R1.fq.gz", "wt")
lambda_r2_out = gzip.open(args.output_prefix + "_Lambda_R2.fq.gz", "wt")

BARCODE85_MOTIF = "TTATTTTT"

n_reads = n_pass = n_trim = n_trim_bases = n_reads_barcode85 = n_written_barcode85 = 0
dmux_counts_R1 = {}
dmux_counts_R2 = {}
barcode_cnt = {}

def assign_barcode_hamming1(bc, wl, ambig_policy="drop"):
    if bc in wl:
        return bc
    hits = []
    bases = ("A","C","G","T")
    bc_list = list(bc)
    for i, orig in enumerate(bc_list):
        for b in bases:
            if b == orig: continue
            cand = bc[:i] + b + bc[i+1:]
            if cand in wl:
                hits.append(cand)
                if ambig_policy == "first":
                    return cand
    if not hits:
        return None
    if ambig_policy == "first":
        return sorted(hits)[0]
    return None

with open_fastq_file(args.read1_fastq) as r1h, open_fastq_file(args.read2_fastq) as r2h:
    for rec1, rec2 in zip(SeqIO.parse(r1h, "fastq"), SeqIO.parse(r2h, "fastq")):
        n_reads += 1
        r1seq = rec1.seq
        r2seq = rec2.seq
        is_barcode85 = BARCODE85_MOTIF in str(r1seq) or BARCODE85_MOTIF in str(r2seq)
        if is_barcode85:
            n_reads_barcode85 += 1

        lambda_seq1 = r1seq
        lambda_qual1 = rec1.letter_annotations["phred_quality"]
        hits = []
        for a in args.adapters:
            hits.extend(find_near_matches(a, str(r1seq), max_l_dist=1))
        if hits:
            trim_i = min(m.start for m in hits)
            if trim_i < len(lambda_seq1):
                n_trim += 1
                n_trim_bases += (len(lambda_seq1) - trim_i)
                lambda_seq1 = lambda_seq1[:trim_i]
                lambda_qual1 = lambda_qual1[:trim_i]

        if args.protocol == "smcseq":
            parsed = parse_r2_smcseq(str(r2seq))
            if parsed is None:
                continue
            b1s, l1s, b2s, l2s, genomic_start = parsed
            lambda_seq2 = r2seq[genomic_start:]
            lambda_qual2 = rec2.letter_annotations["phred_quality"][genomic_start:]
        else:
            te2 = find_near_matches(args.trim_end2, str(r2seq), max_l_dist=1)
            if te2:
                trim2_i = te2[0].end
                lambda_seq2 = r2seq[trim2_i:]
                lambda_qual2 = rec2.letter_annotations["phred_quality"][trim2_i:]
            else:
                lambda_seq2 = r2seq
                lambda_qual2 = rec2.letter_annotations["phred_quality"]
            ## lambda write may precede demux fail (kept for spatialdmt parity)
            if is_barcode85 and len(lambda_seq1) >= args.min_len and len(lambda_seq2) >= args.min_len:
                write_fastq_record(lambda_r1_out, rec1.id, lambda_seq1, lambda_qual1)
                write_fastq_record(lambda_r2_out, rec2.id, lambda_seq2, lambda_qual2)
                n_written_barcode85 += 1
            lk1 = find_near_matches(args.linker1, str(r2seq), max_l_dist=2)
            lk2 = find_near_matches(args.linker2, str(r2seq), max_l_dist=2)
            if not (len(lk1)==1 and len(lk2)==1 and len(te2)==1):
                continue
            l1s = lk1[0].start; l2s = lk2[0].start
            b1s = max(0, l1s - args.bc_side_len)
            b2s = max(0, l2s - args.bc_side_len)

        barcode = str(r2seq[b1s:l1s] + r2seq[b2s:l2s])
        
        if barcode in barcode_cnt:
            barcode_cnt[barcode] += 1
        else:
            barcode_cnt[barcode] = 1
                
        if len(barcode) != args.barcode_len:
            continue

        new_seq1 = lambda_seq1
        new_qual1 = lambda_qual1
        new_seq2 = lambda_seq2
        new_qual2 = lambda_qual2

        if len(new_seq1) < args.min_len or len(new_seq2) < args.min_len:
            continue

        # Write trimmed chunk-level outputs (keep sample-wide merge compatible)
        rid1 = f"{barcode}_{rec1.id}"
        rid2 = f"{barcode}_{rec2.id}"
        write_fastq_record(r1_out, rid1, new_seq1, new_qual1)
        write_fastq_record(r2_out, rid2, new_seq2, new_qual2)
        n_pass += 1
        
        assigned = None
        if do_dmux and wl is not None:
            if args.bc_mismatch <= 0:
                if barcode in wl:
                    assigned = barcode
            else:
            # allow 1 mismatch (you can allow >1 but not recommended)
                if barcode in wl:
                    assigned = barcode
                elif args.bc_mismatch >= 1:
                    assigned = assign_barcode_hamming1(barcode, wl, args.bc_ambig)

        if assigned is not None:
            h1, h2 = get_dmux_handles(assigned)
            h1.write(f"@{rid1}\n{str(new_seq1)}\n+\n{''.join(chr(q+33) for q in new_qual1)}\n")
            h2.write(f"@{rid2}\n{str(new_seq2)}\n+\n{''.join(chr(q+33) for q in new_qual2)}\n")
            dmux_counts_R1[assigned] = dmux_counts_R1.get(assigned, 0) + 1
            dmux_counts_R2[assigned] = dmux_counts_R2.get(assigned, 0) + 1
        elif do_dmux:
            h1, h2 = get_dmux_handles("UNMATCHED")
            h1.write(f"@{rid1}\n{str(new_seq1)}\n+\n{''.join(chr(q+33) for q in new_qual1)}\n")
            h2.write(f"@{rid2}\n{str(new_seq2)}\n+\n{''.join(chr(q+33) for q in new_qual2)}\n")
            dmux_counts_R1["UNMATCHED"] = dmux_counts_R1.get("UNMATCHED", 0) + 1
            dmux_counts_R2["UNMATCHED"] = dmux_counts_R2.get("UNMATCHED", 0) + 1



r1_out.close()
r2_out.close()
lambda_r1_out.close()
lambda_r2_out.close()

for h in dmux_handles_R1.values():
    try: h.close()
    except: pass
for h in dmux_handles_R2.values():
    try: h.close()
    except: pass

# Write per-chunk barcode list/counts (for later merging)
# barcodes file (sorted desc)
if dmux_counts_R1 or dmux_counts_R2:
    _outdir = os.path.dirname(args.output_prefix)
    if _outdir:
        os.makedirs(_outdir, exist_ok=True)

# Barcodes summary across all passed reads (not only whitelisted)
# Recompute from dmux_counts_R1 if present, else empty
if dmux_counts_R1:
    sorted_bc = sorted(dmux_counts_R1.items(), key=lambda x: x[1], reverse=True)
    with open(args.output_prefix + "_barcodes.txt","w") as fbc:
        mx = sorted_bc[0][1]
        for bc, cnt in sorted_bc:
            fbc.write(f"{bc}\t{cnt}\t{'T' if cnt > mx/100.0 else 'F'}\n")
else:
    # still write an empty file to avoid missing targets
    open(args.output_prefix + "_barcodes.txt","w").close()

with open(args.output_prefix + "_stats.txt","w") as fs:
    fs.write(f"Reads_total\t{n_reads}\n")
    fs.write(f"Reads_written_barcode85\t{n_written_barcode85}\n")
    fs.write(f"Reads_detected\t{n_reads - n_reads_barcode85}\n")
    fs.write(f"Reads_passed\t{n_pass}\n")
    fs.write(f"Reads_trimmed\t{n_trim}\n")
    fs.write(f"Bases_trimmed\t{n_trim_bases}\n")

# Per-chunk dmux counts (mates separate) to be merged later
if do_dmux:
    with open(args.output_prefix + f"_R1_dmux.{args.chunk_id}.txt","w") as f1:
        for bc,c in sorted(dmux_counts_R1.items()):
            f1.write(f"{bc}\t{c}\n")
    with open(args.output_prefix + f"_R2_dmux.{args.chunk_id}.txt","w") as f2:
        for bc,c in sorted(dmux_counts_R2.items()):
            f2.write(f"{bc}\t{c}\n")


sorted_barcodes = sorted(barcode_cnt.items(), key=lambda x: x[1], reverse=True)
output_bc_handle = open(args.output_prefix + "_barcodes_all.txt", "w")
cnt_max = sorted_barcodes[0][1]
for barcode, cnt in sorted_barcodes:
    output_bc_handle.write("{}\t{}\t{}\n".format(barcode, cnt, "T" if cnt > cnt_max/100 else "F"))
output_bc_handle.close()

