# Bundled reference assets

These compact YAME state masks are the exact ChromHMM masks named by the
checked-in hg38 and mm10 profiles. They are distributed with the repository so
the standard `final_output` report does not depend on an undocumented lab
filesystem.

| File | Assembly | CpG rows | SHA-256 |
|---|---|---:|---|
| `mm10/ChromHMM.20220414.cm` | UCSC mm10 | 21,867,837 | `49a6a884892d92145b6f1e2ef5120e8338f1228a9ecdc9dd5edf9fbdb0d66bd8` |
| `hg38/ChromHMM.20220303.cm` | UCSC hg38 primary chromosomes | 29,401,795 | `c38ff4081fad37d4678e97c72264e1b80fd2cfb2a5377a5392c07e37f4f3e160` |

The masks use YAME format 2 (categorical state data). Their row order matches
the primary-chromosome CpG BED derived from the pinned BISCUIT QC asset as
documented in `docs/reference_setup.md`. A mask must not be used with a CpG
universe having a different row count or order.
