#!/usr/bin/env python3
"""Exhaustively compare exact-seeded matching with fuzzysearch."""

import itertools
import pathlib
import sys

TOOLS_DIR = pathlib.Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS_DIR))

from fuzzy_match import find_near_matches_exact_seeded  # noqa: E402
from fuzzysearch import find_near_matches  # noqa: E402


def strings(alphabet, length):
    return map("".join, itertools.product(alphabet, repeat=length))


def signature(matches):
    return [
        (match.start, match.end, match.dist, match.matched)
        for match in matches
    ]


def main():
    alphabet = "ACG"
    sequences = [
        sequence
        for length in range(7)
        for sequence in strings(alphabet, length)
    ]
    checked = 0
    for max_l_dist in (1, 2):
        for pattern_length in range(max_l_dist + 1, max_l_dist + 3):
            for pattern in strings(alphabet, pattern_length):
                for sequence in sequences:
                    expected = signature(
                        find_near_matches(
                            pattern, sequence, max_l_dist=max_l_dist
                        )
                    )
                    observed = signature(
                        find_near_matches_exact_seeded(
                            pattern, sequence, max_l_dist=max_l_dist
                        )
                    )
                    if observed != expected:
                        raise AssertionError(
                            f"pattern={pattern!r} sequence={sequence!r} "
                            f"k={max_l_dist}: {observed!r} != {expected!r}"
                        )
                    checked += 1
    print(f"exact-seeded matches fuzzysearch for {checked:,} cases")


if __name__ == "__main__":
    main()
