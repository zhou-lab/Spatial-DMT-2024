#!/usr/bin/env python3
"""Write zero-based BED intervals for CpG dyads in a reference FASTA."""

import argparse
import gzip
from pathlib import Path
import sys


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def fasta_records(path):
    name = None
    parts = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence encountered before its header")
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=Path)
    args = parser.parse_args()

    for chrom, sequence in fasta_records(args.fasta):
        start = sequence.find("CG")
        while start >= 0:
            print(chrom, start, start + 2, sep="\t")
            start = sequence.find("CG", start + 1)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
