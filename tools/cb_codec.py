"""8-char ACGT encoding of (X, Y) spatial coords for the CB SAM tag.

The encoding is mandatory: dupsifter -B silently collapses any non-{ACGTN+-}
CB tag to packed_barcode=0 (dupsifter.c:get_packed_barcode at v1.3.0 L779-786),
so a 4-digit numeric XXYY would over-dedup all bad-CB reads as one cell.

Layout: 4 nt per axis, base-4 (A=0 C=1 G=2 T=3, MSB-first); X then Y. Total
8 chars per CB; (1,1)->"AAAAAAAA", (96,96)->"CCTTCCTT". Reserves the all-N
sentinel "NNNNNNNN" for reads that passed structural parse but missed the
whitelist (dupsifter accepts 'N', all such reads group as one pseudo-cell).
"""

_B4 = "ACGT"
_REV = {"A": 0, "C": 1, "G": 2, "T": 3}
UNMATCHED_CB = "NNNNNNNN"


def coord_to_acgt8(x, y):
    """Encode (X, Y) in [1..96] as an 8-char ACGT string (4 nt per axis)."""
    def enc(n):
        return _B4[(n >> 6) & 3] + _B4[(n >> 4) & 3] + _B4[(n >> 2) & 3] + _B4[n & 3]
    return enc(x - 1) + enc(y - 1)


def acgt8_to_coord(s):
    """Decode 8-char ACGT back to (X, Y). Returns None for the all-N sentinel."""
    if s == UNMATCHED_CB:
        return None
    def dec(s4):
        n = 0
        for c in s4:
            n = n * 4 + _REV[c]
        return n + 1
    return (dec(s[:4]), dec(s[4:]))


def cb_to_cell_name(cb):
    """Return 'XXYY' (e.g. '0142') for a valid CB, 'UNMATCHED' for the all-N
    sentinel, or None if the CB is malformed (wrong length or unknown char)."""
    if cb == UNMATCHED_CB:
        return "UNMATCHED"
    if len(cb) != 8 or any(c not in _REV for c in cb):
        return None
    x, y = acgt8_to_coord(cb)
    return "%02d%02d" % (x, y)
