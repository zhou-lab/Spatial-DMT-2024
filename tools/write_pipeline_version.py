#!/usr/bin/env python3
"""Write the pipeline release label bundled with a sample result."""

import argparse
import csv
from pathlib import Path


FIELDS = (
    "sample",
    "pipeline_version",
    "protocol",
    "reference",
    "assignment",
)


def read_version(path):
    version = Path(path).read_text().strip()
    if not version:
        raise ValueError(f"Empty pipeline version file: {path}")
    return version


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--version-file", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument(
        "--assignment",
        default="workflow_generated",
        choices=("workflow_generated", "backfilled_compatible"),
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    row = {
        "sample": args.sample,
        "pipeline_version": read_version(args.version_file),
        "protocol": args.protocol,
        "reference": args.reference,
        "assignment": args.assignment,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(row)
    temporary.replace(output)


if __name__ == "__main__":
    main()
