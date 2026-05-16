#!/usr/bin/env python3
"""Stream-split a coord-sorted, deduped BAM into per-cell BAMs keyed on CB:Z tag.

Input BAM is expected to carry CB:Z:<8-char ACGT> on every alignment (set by
biscuit align -9 from the <origid>_<CB>_<UMI> read-name convention written by
spatialmeth_trimadapters.py). The 8-char CB is decoded back to (X, Y) using
4-nt-per-axis base-4 (A=0 C=1 G=2 T=3), and the per-cell BAM is named
<out_dir>/<XXYY>.bam (zero-padded). The all-N sentinel goes to UNMATCHED.bam.

Output BAMs are coord-sorted since the input is and we make a single pass.
Reads with no CB tag go to UNTAGGED.bam (separate from UNMATCHED, which marks
"barcode was parsed but missed the whitelist").

Raises the per-process fd ulimit to accommodate up to ~9216 cells + headers.
"""
import argparse, os, resource, sys
import pysam


def raise_fd_limit(target=65535):
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    new = min(target, hard)
    if new > soft:
        resource.setrlimit(resource.RLIMIT_NOFILE, (new, hard))


_REV = {"A": 0, "C": 1, "G": 2, "T": 3}
_UNMATCHED_CB = "NNNNNNNN"


def cb_to_cell_name(cb):
    """Return XXYY (e.g. '0142') for a valid 8-char ACGT CB, 'UNMATCHED' for
    the all-N sentinel, or None if the CB is malformed."""
    if cb == _UNMATCHED_CB:
        return "UNMATCHED"
    if len(cb) != 8 or any(c not in _REV for c in cb):
        return None
    x = 1 + sum(_REV[c] * (4 ** (3 - i)) for i, c in enumerate(cb[:4]))
    y = 1 + sum(_REV[c] * (4 ** (3 - i)) for i, c in enumerate(cb[4:]))
    return "%02d%02d" % (x, y)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("in_bam")
    ap.add_argument("out_dir")
    ap.add_argument("--tag", default="CB", help="SAM tag to split on (default CB)")
    ap.add_argument("--threads", type=int, default=1,
        help="threads for BGZF compression on each output BAM (default 1)")
    args = ap.parse_args()

    raise_fd_limit()
    os.makedirs(args.out_dir, exist_ok=True)

    inb = pysam.AlignmentFile(args.in_bam, "rb")
    out_handles = {}
    n_reads = n_tagged = n_untagged = n_malformed = 0

    for read in inb:
        n_reads += 1
        try:
            cb = read.get_tag(args.tag)
            n_tagged += 1
        except KeyError:
            cell = "UNTAGGED"
            n_untagged += 1
        else:
            cell = cb_to_cell_name(cb)
            if cell is None:
                cell = "MALFORMED"
                n_malformed += 1
        h = out_handles.get(cell)
        if h is None:
            out_path = os.path.join(args.out_dir, f"{cell}.bam")
            h = pysam.AlignmentFile(out_path, "wb", template=inb, threads=args.threads)
            out_handles[cell] = h
        h.write(read)

    for h in out_handles.values():
        h.close()
    inb.close()

    sys.stderr.write(
        f"reads_total={n_reads}  tagged={n_tagged}  untagged={n_untagged}  "
        f"malformed={n_malformed}  cells={len(out_handles)}\n")


if __name__ == "__main__":
    main()
