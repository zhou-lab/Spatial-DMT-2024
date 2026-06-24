#!/usr/bin/env python3
"""Extract the published SmC-seq spatial-RNA barcode and UMI cassette."""

import argparse
import gzip
import os

import regex


MOTIF1 = "ATCGGCGTACGACT"
MOTIF2 = "ATCCACGTGCTTGAGCGCGCTGCATACTTG"
MOTIF3 = "CCCATGATCGTCCGATGCAGTCGTGCCATGAGATGTGTATAAGAGACAG"
R1_READTHROUGH_MOTIFS = (
    "AGTCGTACGCCGAT",
    "CAAGTATGCAGCGCGCTCAAGCACGTGGAT",
    "CTGTCTCTTATACACATCTCATGGCACGACTGCATCGGACGATCATGGG",
    "AAAAAAAAAAAAAAAAAAAA",
)
MIN_TRANSCRIPT_LENGTH = 25

PRIMER_PATTERN = regex.compile(rf"^.(?:{MOTIF1}){{s<=1}}")
LINKER1_PATTERN = regex.compile(
    rf"^.(?:{MOTIF1}){{s<=1}}.{{8}}(?:{MOTIF2}){{e<=1}}"
)
FULL_PATTERN = regex.compile(
    rf"^.(?:{MOTIF1}){{s<=1}}"
    rf"(?P<barcode1>.{{8}})"
    rf"(?:{MOTIF2}){{e<=1}}"
    rf"(?P<barcode2>.{{8}})"
    rf"(?:{MOTIF3}){{s<=2,i<=1,d<=1}}"
    rf"(?P<umi>.{{10}}).*$"
)


def open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def read_fastq(handle):
    while True:
        title = handle.readline()
        if not title:
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not quality:
            raise ValueError("Truncated FASTQ record")
        yield title.rstrip("\n"), sequence.rstrip("\n"), plus.rstrip("\n"), quality.rstrip("\n")


def read_id(title):
    return title.split()[0].removeprefix("@").removesuffix("/1").removesuffix("/2")


def write_record(handle, title, sequence, plus, quality):
    handle.write(f"{title}\n{sequence}\n{plus}\n{quality}\n")


def trim_r1_readthrough(record):
    """Trim SmC-seq cassette read-through from the transcript read."""
    title, sequence, plus, quality = record
    cut = len(sequence)
    for motif in R1_READTHROUGH_MOTIFS:
        pos = sequence.find(motif)
        if pos != -1:
            cut = min(cut, pos)
    if cut < MIN_TRANSCRIPT_LENGTH:
        return None
    return title, sequence[:cut], plus, quality[:cut]


def write_stats(path, inputs, total_pairs, matched_pairs):
    pct = 100.0 * matched_pairs / total_pairs if total_pairs else 0.0
    with open(path, "w") as handle:
        handle.write(f"#File\t{inputs[0]}\t{inputs[1]}\n")
        handle.write(f"#Total\t{2 * total_pairs}\n")
        handle.write(f"#Matched\t{matched_pairs}\t{pct:.5f}%\n")
        handle.write("#Name\tReads\tReadsPct\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--read-seq-output", required=True)
    parser.add_argument("--barcode-output", required=True)
    parser.add_argument("--primer-stats", required=True)
    parser.add_argument("--linker1-stats", required=True)
    parser.add_argument("--linker2-stats", required=True)
    args = parser.parse_args()

    for path in (
        args.read_seq_output,
        args.barcode_output,
        args.primer_stats,
        args.linker1_stats,
        args.linker2_stats,
    ):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    total = primer_pass = linker1_pass = full_pass = transcript_pass = 0
    with open_text(args.r1, "r") as r1_handle, open_text(args.r2, "r") as r2_handle, open_text(
        args.read_seq_output, "w"
    ) as read_seq_out, open_text(args.barcode_output, "w") as barcode_out:
        r1_records = read_fastq(r1_handle)
        r2_records = read_fastq(r2_handle)
        while True:
            r1 = next(r1_records, None)
            r2 = next(r2_records, None)
            if r1 is None or r2 is None:
                if r1 is not None or r2 is not None:
                    raise ValueError("RNA R1 and R2 contain different record counts")
                break
            if read_id(r1[0]) != read_id(r2[0]):
                raise ValueError(f"RNA read IDs do not match: {r1[0]} / {r2[0]}")

            total += 1
            sequence = r2[1]
            if PRIMER_PATTERN.match(sequence):
                primer_pass += 1
            if LINKER1_PATTERN.match(sequence):
                linker1_pass += 1
            match = FULL_PATTERN.match(sequence)
            if not match:
                continue

            barcode_start = match.start("barcode1")
            barcode2_start = match.start("barcode2")
            umi_start = match.start("umi")
            processed_sequence = (
                match.group("barcode1") + match.group("barcode2") + match.group("umi")
            )
            processed_quality = (
                r2[3][barcode_start : barcode_start + 8]
                + r2[3][barcode2_start : barcode2_start + 8]
                + r2[3][umi_start : umi_start + 10]
            )
            full_pass += 1
            trimmed_r1 = trim_r1_readthrough(r1)
            if trimmed_r1 is None:
                continue
            write_record(read_seq_out, *trimmed_r1)
            write_record(barcode_out, r2[0], processed_sequence, r2[2], processed_quality)
            transcript_pass += 1

    inputs = (args.r1, args.r2)
    write_stats(args.primer_stats, inputs, total, primer_pass)
    write_stats(args.linker1_stats, inputs, primer_pass, linker1_pass)
    write_stats(args.linker2_stats, inputs, linker1_pass, full_pass)
    print(
        f"RNA pairs: total={total} primer={primer_pass} "
        f"linker1={linker1_pass} full_structure={full_pass} "
        f"transcript_length={transcript_pass}"
    )


if __name__ == "__main__":
    main()
