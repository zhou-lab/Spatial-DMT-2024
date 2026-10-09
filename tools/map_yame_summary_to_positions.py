#!/usr/bin/env python3
"""Map a barcode-indexed yame summary onto spatial grid coordinates."""

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--positions", type=Path)
    source.add_argument(
        "--barcode-map",
        type=Path,
        help="Five-column Spatial-DMT barcode table; yame barcode is column 4",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    summary = pd.read_csv(args.summary, sep="\t", dtype={"Query": str})
    required = {"Query", "N_query"}
    if not required.issubset(summary.columns):
        raise ValueError(f"{args.summary} is missing {sorted(required - set(summary))}")

    if args.positions:
        positions = pd.read_csv(
            args.positions,
            header=None,
            names=[
                "spatial_id",
                "in_tissue",
                "array_row",
                "array_col",
                "pixel_y",
                "pixel_x",
            ],
            dtype={"spatial_id": str},
        )
        positions["x"] = positions["array_col"].astype(int) + 1
        positions["y"] = positions["array_row"].astype(int) + 1
        extra_columns = ["in_tissue", "pixel_x", "pixel_y"]
    else:
        barcode_map = pd.read_csv(args.barcode_map, sep="\t", header=None, dtype=str)
        if barcode_map.shape[1] < 4:
            raise ValueError(f"{args.barcode_map} must contain at least four columns")
        positions = pd.DataFrame(
            {
                "spatial_id": barcode_map.iloc[:, 3],
                "x": pd.to_numeric(barcode_map.iloc[:, 1], errors="raise").astype(int),
                "y": pd.to_numeric(barcode_map.iloc[:, 2], errors="raise").astype(int),
            }
        )
        extra_columns = []
    position_source = args.positions or args.barcode_map
    if positions["spatial_id"].duplicated().any():
        raise ValueError(f"{position_source} contains duplicate barcodes")

    mapped = summary.loc[:, ["Query", "N_query"]].merge(
        positions, left_on="Query", right_on="spatial_id", how="left", validate="one_to_one"
    )
    missing = mapped["x"].isna()
    if missing.any():
        examples = ", ".join(mapped.loc[missing, "Query"].head(5))
        raise ValueError(f"{missing.sum()} yame barcodes lack positions: {examples}")
    absent = positions.loc[~positions["spatial_id"].isin(summary["Query"])]
    if not absent.empty:
        raise ValueError(f"{len(absent)} position barcodes lack yame summaries")

    output = mapped.loc[:, ["spatial_id", "x", "y", *extra_columns]].copy()
    output["x"] = output["x"].astype(int)
    output["y"] = output["y"].astype(int)
    output["covered_cpgs"] = pd.to_numeric(
        mapped["N_query"], errors="raise"
    ).astype(int)
    output = output.loc[
        :, ["spatial_id", "x", "y", "covered_cpgs", *extra_columns]
    ].sort_values(["x", "y"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
