# HTSlib CB splitter

`split_bam_by_cb_hts` streams a BAM into per-cell BAMs using the `CB:Z` tag.
The eight-base ACGT tag is decoded to the pipeline's `XXYY` filename. Missing
and invalid tags are written to `UNTAGGED.bam` and `MALFORMED.bam`.

Build against the lab HTSlib installation:

```bash
make
```

Override `HTSLIB_ROOT` for another static HTSlib build. Direct mode opens one
output per observed cell and is intended for compatibility/testing:

```bash
split_bam_by_cb_hts -@ 8 input.bam output_dir
```

Production uses `../split_bam_by_cb_bucketed.sh`. Its first pass writes a
bounded number of deterministic buckets; bucket splits then run concurrently.
This bounds open files and HTSlib thread state while preserving record order
within every cell. The wrapper also plans local/fallback scratch capacity and
removes buckets as they are consumed.
