"""8-char ACGT encoding of (X, Y) spatial coords for the CB SAM tag.

The encoding is mandatory: dupsifter -B silently collapses any non-{ACGTN+-}
CB tag to packed_barcode=0 (dupsifter.c:get_packed_barcode at v1.3.0 L779-786),
so a 4-digit numeric XXYY would over-dedup all bad-CB reads as one cell.

Layout: 4 nt per axis, base-4 (A=0 C=1 G=2 T=3, MSB-first); X then Y. Total
8 chars per CB; (1,1)->"AAAAAAAA", (96,96)->"CCTTCCTT".

The merged BAM only ever contains reads with a real whitelist coord -- reads
that fail structure or whitelist match are split out at the trim/tag stage
into the separate Unmatched fastq pair and aligned to bam_unmatched/. There
is no sentinel-CB cell in this codec.
"""

_B4 = "ACGT"
_REV = {"A": 0, "C": 1, "G": 2, "T": 3}


def coord_to_acgt8(x, y):
    """Encode (X, Y) in [1..96] as an 8-char ACGT string (4 nt per axis)."""
    def enc(n):
        return _B4[(n >> 6) & 3] + _B4[(n >> 4) & 3] + _B4[(n >> 2) & 3] + _B4[n & 3]
    return enc(x - 1) + enc(y - 1)


def acgt8_to_coord(s):
    """Decode 8-char ACGT back to (X, Y). Returns None if the CB is malformed."""
    if len(s) != 8 or any(c not in _REV for c in s):
        return None
    def dec(s4):
        n = 0
        for c in s4:
            n = n * 4 + _REV[c]
        return n + 1
    return (dec(s[:4]), dec(s[4:]))


def cb_to_cell_name(cb):
    """Return 'XXYY' (e.g. '0142') for a valid CB, or None if malformed."""
    xy = acgt8_to_coord(cb)
    if xy is None:
        return None
    return "%02d%02d" % xy
