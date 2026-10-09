#!/usr/bin/env python3
"""Summarize one per-pixel YAME .cg file per spatial coordinate."""

import argparse
import csv
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


COORDINATE_RE = re.compile(r"^(\d{2})(\d{2})\.cg$")


def summarize(path, yame):
    match = COORDINATE_RE.fullmatch(path.name)
    if match is None:
        raise ValueError(f"Invalid spatial CG filename: {path}")
    completed = subprocess.run(
        [yame, "summary", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = list(csv.DictReader(completed.stdout.splitlines(), delimiter="\t"))
    if len(rows) != 1 or "N_query" not in rows[0]:
        raise RuntimeError(f"Unexpected YAME output for {path}")
    return {
        "spatial_id": path.stem,
        "x": int(match.group(1)),
        "y": int(match.group(2)),
        "covered_cpgs": int(rows[0]["N_query"]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cg-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--yame", default="yame")
    parser.add_argument("--threads", type=int, default=24)
    args = parser.parse_args()

    cg_dir = Path(args.cg_dir)
    paths = sorted(
        path
        for path in cg_dir.glob("*.cg")
        if COORDINATE_RE.fullmatch(path.name)
    )
    if not paths:
        parser.error(f"No coordinate-named CG files found under {cg_dir}")
    if args.threads < 1:
        parser.error("--threads must be at least 1")

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        rows = list(executor.map(lambda path: summarize(path, args.yame), paths))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["spatial_id", "x", "y", "covered_cpgs"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output} with {len(rows):,} pixels")


if __name__ == "__main__":
    main()
