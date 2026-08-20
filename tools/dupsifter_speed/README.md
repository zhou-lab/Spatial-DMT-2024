# Optimized vendored dupsifter

This directory vendors the MIT-licensed dupsifter 1.3.0 source used by the
production pipeline. The upstream source files are unchanged between the
`v1.3.0.20241113` tag and the current 1.3.1 development branch, apart from the
version string.

Local changes are intentionally limited to:

- hash every duplicate-signature field instead of a truncated first-read
  coordinate;
- use 64-bit counters for samples with more than 4.29 billion alignments;
- optionally emit uncompressed BAM with `--uncompressed-bam` when streaming
  directly into `samtools sort`; and
- return a failure when input, alignment output, or statistics I/O is
  truncated instead of allowing a pipeline to accept partial output. Invalid
  paired-template errors return through HTSlib cleanup so every previously
  completed record is flushed before the nonzero exit.

Build with `make`. Override `HTSLIB_ROOT` when the shared lab HTSlib 1.17 build
is not available.
