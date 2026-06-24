#!/usr/bin/env python3
"""Count zero, unique, and multiple fuzzy-anchor matches in SF R2 reads."""

import gzip
import sys
from collections import Counter

import regex


PATTERNS = {
    "linker1": regex.compile(r"(GTGGTTGATGTTTTGTATTGGTGTATGATT){e<=2}"),
    "linker2": regex.compile(r"(ATTTATGTGTTTGAGAGGTTAGAGTATTTG){e<=2}"),
    "tn5_me": regex.compile(r"(AGATGTGTATAAGAGATAG){e<=1}"),
}


def bucket(n):
    return "0" if n == 0 else "1" if n == 1 else ">1"


independent = Counter()
first_failure = Counter()
combinations = Counter()
total_sf = 0

with gzip.open(sys.argv[1], "rt") as handle:
    while True:
        header = handle.readline()
        if not header:
            break
        sequence = handle.readline().rstrip("\n")
        handle.readline()
        handle.readline()
        if not header.split(None, 1)[0].endswith("_SF"):
            continue

        total_sf += 1
        counts = {name: len(pattern.findall(sequence)) for name, pattern in PATTERNS.items()}
        states = tuple(bucket(counts[name]) for name in PATTERNS)
        combinations[states] += 1
        for name, n in counts.items():
            independent[(name, bucket(n))] += 1

        for name in PATTERNS:
            state = bucket(counts[name])
            if state != "1":
                first_failure[(name, state)] += 1
                break
        else:
            first_failure[("all_ok", "1")] += 1

print("section\tanchor_or_pattern\tmatch_count_class\tcount\tfraction")
for name in PATTERNS:
    for state in ("0", "1", ">1"):
        n = independent[(name, state)]
        print(f"independent\t{name}\t{state}\t{n}\t{n / total_sf:.6f}")
for name in (*PATTERNS, "all_ok"):
    for state in ("0", ">1", "1"):
        n = first_failure[(name, state)]
        if n:
            print(f"first_failure\t{name}\t{state}\t{n}\t{n / total_sf:.6f}")
for states, n in combinations.most_common():
    print(f"combination\t{'/'.join(states)}\tNA\t{n}\t{n / total_sf:.6f}")
print(f"total\tSF\tNA\t{total_sf}\t1.000000")
