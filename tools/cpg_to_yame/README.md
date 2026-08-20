# cpg_to_yame

`cpg_to_yame` maps BISCUIT methylation calls directly into sparse YAME
format-3 CpG records. It avoids expanding every cell to every reference CpG,
repeated reference scans, and intermediate context BED files.

Build the shared coordinate index once:

```bash
cpg_to_yame index cpg.bed.gz cpg.idx
```

Then convert one cell directly from the pileup VCF:

```bash
cpg_to_yame pack-vcf -i cpg.idx -o cell.cg -c cell.conv.tsv cell.vcf.gz
```

`pack-vcf` reads the first VCF sample, matching `biscuit vcf2bed`'s default.
`-@ 2` enables a second BGZF reader thread for latency-sensitive large cells;
one thread generally gives better throughput when many cells run concurrently.

Convert every sample in a joint VCF in one pass:

```bash
mkdir cells
cpg_to_yame pack-vcf-batch -i cpg.idx -o cells batch.vcf.gz
```

`pack-vcf` and `pack-vcf-batch` also accept `-`, allowing BISCUIT output to be
converted as it is generated:

```bash
biscuit pileup reference.fa cell.bam |
    cpg_to_yame pack-vcf -i cpg.idx -o cell.cg -c cell.conv.tsv -
```

Batch outputs are named from the VCF sample columns. Names are restricted to
letters, digits, `.`, `_`, and `-`; duplicate or unsafe names are rejected.
Each output is staged before publication.

The `pack` subcommand remains available for an existing nine-column
`biscuit vcf2bed -e -t c` stream.

The index stores one 32-bit start coordinate per reference CpG and is intended
to be memory-mapped by concurrent workers. The output record uses YAME's
documented format-3 run encoding and is readable by the standard `yame` CLI.

Build and run the synthetic parser/encoding tests with:

```bash
make test
```
