#!/usr/bin/env python3
"""Check the precomputed Hamming-1 lookup against the original search."""

import argparse
import ast
from pathlib import Path


FUNCTIONS = {"assign_barcode_hamming1", "build_hamming1_lookup"}


def load_functions(script_path):
    tree = ast.parse(script_path.read_text(), filename=str(script_path))
    body = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in FUNCTIONS
    ]
    if {node.name for node in body} != FUNCTIONS:
        raise RuntimeError(f"could not find both barcode functions in {script_path}")
    namespace = {}
    exec(compile(ast.Module(body=body, type_ignores=[]), str(script_path), "exec"), namespace)
    return namespace["assign_barcode_hamming1"], namespace["build_hamming1_lookup"]


def load_whitelist(path, column):
    whitelist = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= max(column, 3):
                whitelist[fields[column - 1]] = (int(fields[1]), int(fields[2]))
    if not whitelist:
        raise RuntimeError(f"empty whitelist from {path} column {column}")
    return whitelist


def observations(whitelist):
    seen = set(whitelist)
    for barcode in whitelist:
        for position in range(len(barcode)):
            for base in "ACGTNR":
                seen.add(barcode[:position] + base + barcode[position + 1 :])
    return sorted(seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trim_script", type=Path)
    parser.add_argument("whitelist", type=Path)
    parser.add_argument("--column", type=int, required=True)
    args = parser.parse_args()

    assign_original, build_lookup = load_functions(args.trim_script)
    whitelist = load_whitelist(args.whitelist, args.column)
    cases = observations(whitelist)
    comparisons = 0

    for policy in ("drop", "first"):
        lookup = build_lookup(whitelist, policy)
        for barcode in cases:
            if barcode in whitelist:
                expected = (barcode, whitelist[barcode])
                observed = expected
            else:
                expected = assign_original(barcode, whitelist, policy)
                observed = lookup.get(barcode, (None, None))
                if (
                    observed[0] is None
                    and not set(barcode).issubset({"A", "C", "G", "T"})
                ):
                    observed = assign_original(barcode, whitelist, policy)
            if observed != expected:
                raise AssertionError(
                    f"{policy} mismatch for {barcode}: {observed} != {expected}"
                )
            comparisons += 1

    print(
        f"precomputed lookup matches original Hamming-1 search for "
        f"{comparisons:,} cases"
    )


if __name__ == "__main__":
    main()
