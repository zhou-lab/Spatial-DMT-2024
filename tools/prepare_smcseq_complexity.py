#!/usr/bin/env python3
"""Prepare published SmC-seq pixel counts for library-complexity analysis."""

import argparse
import json
import re
from pathlib import Path

import pandas as pd


SAMPLE_NAMES = {
    "E55_5mC_10um_rep1": "SmC_E55_10um_rep1",
    "E55_5mC_10um_rep2": "SmC_E55_10um_rep2",
    "E65_5mC_10um": "SmC_E65_10um",
    "E75_5mC_10um_rep1": "SmC_E75_10um_rep1",
    "E75_5mC_10um_rep2": "SmC_E75_10um_rep2",
}


def coordinate(series: pd.Series, prefix: str) -> pd.Series:
    values = series.astype(str).str.extract(rf"^{re.escape(prefix)}(\d+)$")[0]
    if values.isna().any():
        invalid = sorted(series.loc[values.isna()].astype(str).unique())
        raise ValueError(f"Invalid {prefix} coordinates: {invalid[:5]}")
    return values.astype(int)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    pixel = pd.read_excel(
        args.source, sheet_name="Supplementary Table 2", header=2
    )
    pixel.columns = [str(column).strip() for column in pixel.columns]
    required = {"Sample", "nA", "nB", "Mapped_reads", "Unique_reads"}
    missing = required - set(pixel.columns)
    if missing:
        raise ValueError(f"Pixel table is missing columns: {sorted(missing)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for published_name, output_name in SAMPLE_NAMES.items():
        sample = pixel.loc[pixel["Sample"].eq(published_name)].copy()
        if sample.empty:
            raise ValueError(f"No rows found for {published_name}")

        mapped = pd.to_numeric(sample["Mapped_reads"], errors="raise").astype(int)
        unique = pd.to_numeric(sample["Unique_reads"], errors="raise").astype(int)
        if (mapped <= 0).any() or (unique <= 0).any() or (unique > mapped).any():
            raise ValueError(f"Invalid mapped/unique counts for {published_name}")

        output = pd.DataFrame(
            {
                "x": coordinate(sample["nA"], "A"),
                "y": coordinate(sample["nB"], "B"),
                "mapped_reads": mapped,
                "dup_reads": mapped - unique,
            }
        )
        if output.duplicated(["x", "y"]).any():
            raise ValueError(f"Duplicate coordinates for {published_name}")
        if not output[["x", "y"]].isin(range(1, 97)).all().all():
            raise ValueError(f"Coordinates outside the 96 x 96 grid for {published_name}")

        path = args.output_dir / f"{output_name}_published_counts.tsv"
        output.sort_values(["y", "x"]).to_csv(path, sep="\t", index=False)
        manifest.append(
            {
                "sample": output_name,
                "published_sample": published_name,
                "pixels": len(output),
                "mapped_reads": int(mapped.sum()),
                "unique_reads": int(unique.sum()),
                "duplication_fraction": float(1 - unique.sum() / mapped.sum()),
                "output": str(path),
            }
        )

    (args.output_dir / "published_counts_manifest.json").write_text(
        json.dumps(
            {
                "source": str(args.source),
                "sheet": "Supplementary Table 2",
                "coordinate_mapping": {"x": "numeric part of nA", "y": "numeric part of nB"},
                "count_mapping": {
                    "mapped_reads": "Mapped_reads",
                    "dup_reads": "Mapped_reads - Unique_reads",
                },
                "samples": manifest,
            },
            indent=2,
        )
        + "\n"
    )
    print(pd.DataFrame(manifest).to_string(index=False))


if __name__ == "__main__":
    main()
