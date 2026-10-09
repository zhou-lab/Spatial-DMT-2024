#!/usr/bin/env python3
"""Extract indexed YAME per-pixel methylation summaries for selected samples."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

TOOLS_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_ROOT))
from plot_spatial_methylation import summarize_cells_bulk  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--yame", default="yame")
    parser.add_argument("--min-covered-cpgs", type=int, default=100)
    parser.add_argument("samples", nargs="+")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = ["sample", "spatial_id", "x", "y", "covered_cpgs", "cpg_methylation", "depth", "passes_threshold"]
    for sample in args.samples:
        cg_path = args.run_root / "final_output" / sample / f"{sample}.cg"
        rows = summarize_cells_bulk(str(cg_path), args.yame)
        out = args.output_dir / f"{sample}_spatial_cpg_methylation.tsv"
        with out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {
                        "sample": sample,
                        "spatial_id": row["spatial_id"],
                        "x": row["x"],
                        "y": row["y"],
                        "covered_cpgs": row["covered_cpgs"],
                        "cpg_methylation": f'{row["beta"]:.6f}' if np.isfinite(row["beta"]) else "NA",
                        "depth": f'{row["depth"]:.6f}' if np.isfinite(row["depth"]) else "NA",
                        "passes_threshold": int(row["covered_cpgs"] >= args.min_covered_cpgs),
                    }
                )
        print(f"Wrote {out} ({len(rows)} pixels)", flush=True)


if __name__ == "__main__":
    main()
