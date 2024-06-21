#!/usr/bin/env python3

import argparse
import gzip
from fuzzysearch import find_near_matches
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Create the argument parser
parser = argparse.ArgumentParser(description="Process paired-end FASTQ files with adapter filtering. Example usage: python spatialmeth_trimadapters.py reads1.fastq.gz reads2.fastq.gz -o output")

# Add positional arguments
parser.add_argument("read1_fastq", help="Path to the read 1 FASTQ file.")
parser.add_argument("read2_fastq", help="Path to the read 2 FASTQ file.")

# Add optional flag argument
parser.add_argument("-o", "--output_prefix", default="out", help="Output file prefix.")
parser.add_argument("-a", "--adapters", nargs='+', default=["CTATCTCTTATA", "AGATGCGAGAAGCCAACGCTTG"], help="Read 1 adapter sequence to trim.")
parser.add_argument("-l1", "--linker1", default="GTGGTTGATGTTTTGTATTGGTGTATGATT", help="First linker sequence.")
parser.add_argument("-l2", "--linker2", default="ATTTATGTGTTTGAGAGGTTAGAGTATTTG", help="Second linker sequence.")
parser.add_argument("-t", "--trim_end2", default="AGATGTGTATAAGAGATAG", help="Trim read 2 until this subsequence.")

# Parse the command-line arguments
args = parser.parse_args()

# Access the arguments
read1_fastq_path = args.read1_fastq
read2_fastq_path = args.read2_fastq

# Function to open the FASTQ file, handling gzip compression if necessary
def open_fastq_file(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")
    
def write_seqrecord_to_fastq(id, seq, qual, f):
    f.write("@{}\n".format(id))
    f.write(str(seq) + "\n")
    f.write("+\n")
    f.write("".join(chr(q + 33) for q in qual) + "\n")

# Create output file handles
output_read1_handle = gzip.open(args.output_prefix + "_R1.fq.gz", "wt")
output_read2_handle = gzip.open(args.output_prefix + "_R2.fq.gz", "wt")

# Read the paired-end FASTQ files simultaneously
n_reads = 0
n_reads_passed = 0
n_reads_trimmed = 0
n_bases_trimmed = 0

barcode_cnt = {}
with open_fastq_file(args.read1_fastq) as read1_handle, open_fastq_file(args.read2_fastq) as read2_handle:
    for read1_record, read2_record in zip(SeqIO.parse(read1_handle, "fastq"), SeqIO.parse(read2_handle, "fastq")):
        n_reads += 1
        
        # Access the sequences of read 1 and read 2
        read1_seq = read1_record.seq
        read2_seq = read2_record.seq

        matches_linker1 = find_near_matches(args.linker1, str(read2_seq), max_l_dist=2)
        matches_linker2 = find_near_matches(args.linker2, str(read2_seq), max_l_dist=2)
        matches_trim2 = find_near_matches(args.trim_end2, str(read2_seq), max_l_dist=1)
        # print(len(matches_linker1), len(matches_linker2), len(matches_trim2))
        # Check if both linker sequences are present in read 2 sequence
        if len(matches_linker1)==1 and len(matches_linker2)==1 and len(matches_trim2)==1:
            n_reads_passed += 1
            linker1_start = matches_linker1[0].start
            linker2_start = matches_linker2[0].start
            barcode = str(read2_seq[(linker1_start-8):linker1_start] + read2_seq[(linker2_start-8):linker2_start])

            # Count barcodes
            if barcode in barcode_cnt:
                barcode_cnt[barcode] += 1
            else:
                barcode_cnt[barcode] = 1

            # process and write read 1
            new_id = barcode+"_"+read1_record.id
            new_seq = read1_record.seq
            new_qual = read1_record.letter_annotations["phred_quality"]
            adapter_matches = []
            for adapter in args.adapters:
                adapter_matches.extend(find_near_matches(adapter, str(read1_seq), max_l_dist=1))
                
            if len(adapter_matches) > 0:
                trim_index = min([x.start for x in adapter_matches])
                n_reads_trimmed += 1
                n_bases_trimmed += len(new_seq) - trim_index
                new_seq = new_seq[:trim_index]
                new_qual = new_qual[:trim_index]
            write_seqrecord_to_fastq(new_id, new_seq, new_qual, output_read1_handle)

            # process and write read 2
            trim_index = matches_trim2[0].end
            new_id = barcode+"_"+read2_record.id
            new_seq = read2_record.seq[trim_index:]
            new_qual = read2_record.letter_annotations["phred_quality"][trim_index:]
            write_seqrecord_to_fastq(new_id, new_seq, new_qual, output_read2_handle)

sorted_barcodes = sorted(barcode_cnt.items(), key=lambda x: x[1], reverse=True)
output_bc_handle = open(args.output_prefix + "_barcodes.txt", "w")
cnt_max = sorted_barcodes[0][1]
for barcode, cnt in sorted_barcodes:
    output_bc_handle.write("{}\t{}\t{}\n".format(barcode, cnt, "T" if cnt > cnt_max/100 else "F"))
output_bc_handle.close()

output_stats_handle = open(args.output_prefix + "_stats.txt", "w")
output_stats_handle.write("Reads_detected\t{}\n".format(n_reads))
output_stats_handle.write("Reads_passed\t{}\n".format(n_reads_passed))
output_stats_handle.write("Reads_trimmed\t{}\n".format(n_reads_trimmed))
output_stats_handle.write("Bases_trimmed\t{}\n".format(n_bases_trimmed))
output_stats_handle.close()

# Close the output file handles
output_read1_handle.close()
output_read2_handle.close()
