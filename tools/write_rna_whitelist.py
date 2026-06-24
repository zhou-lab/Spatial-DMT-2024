#!/usr/bin/env python3
"""Build combined STARsolo and coordinate whitelists from ordered 8-mer barcodes."""

import argparse
import os


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--star-whitelist", required=True)
    parser.add_argument("--coordinate-whitelist", required=True)
    args = parser.parse_args()

    with open(args.barcodes) as handle:
        barcodes = [line.strip() for line in handle if line.strip()]
    if len(barcodes) != len(set(barcodes)):
        raise ValueError("RNA barcode component list contains duplicates")
    if not barcodes or any(len(barcode) != 8 for barcode in barcodes):
        raise ValueError("RNA barcode components must be unique 8-mers")

    for path in (args.star_whitelist, args.coordinate_whitelist):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    with open(args.star_whitelist, "w") as star, open(
        args.coordinate_whitelist, "w"
    ) as coordinates:
        for x, barcode1 in enumerate(barcodes, 1):
            for y, barcode2 in enumerate(barcodes, 1):
                combined = barcode1 + barcode2
                star.write(f"{combined}\n")
                coordinates.write(f"{combined}\t{x}\t{y}\n")


if __name__ == "__main__":
    main()
