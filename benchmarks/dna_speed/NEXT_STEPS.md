# DNA speed optimization: resume notes

This file records work intentionally left after the 2026-08-20 closeout. The
implemented and measured changes are summarized in `README.md`. Resume from
branch `perf/dna-pipeline-speed-20260818` in
`/home/suc3/Lab/Spatial-DMT-2024-speedopt`; do not modify or cancel an older
production job merely to obtain benchmark capacity.

## Highest-value next experiments

1. **Measure BISCUIT scaling at production thread counts.** After the old 23DM
   job releases its 48 CPUs, run one serial fixed-input matrix on the same node
   at 48, 64, and 69 alignment threads. Use at least 50M real unmatched pairs
   so startup, final partial chunks, and sort tails are small relative to the
   alignment. Record `sacct` CPU/RSS, BISCUIT CPU/real time, primary counts,
   flagstat/stats, and coordinate-normalized records. The 69-CPU deep estimate
   is currently an extrapolation from 20-thread fixed input and historical
   24/48-thread runs.

2. **Complete a real shallow sample, then one optimized deep sample.** The
   closeout launches a real shallow SpMET sample only far enough to rule out
   startup/integration failures. Let that sample finish and compare trim
   statistics, primary mapping, duplicate statistics, CpG dimensions, and QC
   metrics with its previous pipeline run. Then rerun 21DM or 23DM from trim
   through final QC using the live capacity helper and retain per-rule Slurm
   accounting. This is required before replacing the production branch.

3. **Distribute per-cell pileup across two nodes.** Refactor the monolithic
   pileup rule into a checkpoint plus moderate batches (roughly 10-25 cells per
   Slurm job), not 2,500 one-cell scheduler jobs. With matched/unmatched
   alignment complete, the current 144-CPU QOS could run about 138 one-thread
   cells across two nodes instead of 69 on one node. Preserve atomic VCF/CG
   outputs, restartability, deterministic cell order, one final `yame index -s`,
   and the exact deep-cell regression.

4. **Remove repeated reference-CpG scans during cell conversion.** Each cell
   still coordinate-sorts calls and asks bedtools to scan the full CpG BED.
   Extend the context helper or add a batch converter that loads/reference-maps
   CpG coordinates once and processes several cell streams. Benchmark a mix of
   sparse and dense 23DM cells, and require byte-identical `.cg` output. This is
   likely the strongest remaining algorithmic pileup opportunity.

5. **Benchmark BAM compression end to end, not at sort alone.** Compression
   levels 1/3 cut standalone sort time by 50%/39% but enlarge BAMs by 28.5%/
   17.3%; they were not adopted because sort is a tiny tail and larger BAMs can
   slow split, flagstat, and pileup on shared storage. Time all downstream
   consumers on equivalent deep BAMs. Adopt a lower level only if total wall
   time improves under representative NFS load; moderate extra storage is
   acceptable, but another storage exhaustion is not.

6. **Validate the time-first trim defaults at full concurrency.** The adopted
   512-MiB blocks and gzip level 1 are exact and fastest in the final
   two-worker deep-concurrency fixture. Run at 65 workers on at least a 50M-pair
   real stream, capture aggregate Slurm RSS, and confirm the 160-GiB request has
   useful headroom. Compare against 512-MiB/gzip-3 on the same node. Do not infer
   total RSS by multiplying `/usr/bin/time` maximum child RSS.

## Harder algorithmic options

- **Compiled trim/parser implementation.** Port structure parsing, exact-seeded
  matching, barcode lookup, and FASTQ emission to a compiled streaming tool.
  The current matcher projects to roughly 160-165 CPU-hours on 21DM, so a
  rigorously equivalent compiled path could still save hours. Preserve every
  trim-stat definition and compare decompressed outputs on exhaustive synthetic
  and real fixtures.

- **Multi-node alignment scatter.** Equal-core four-shard BISCUIT was slower
  than one monolithic process. A useful scatter design must use otherwise idle
  QOS capacity and retain global duplicate semantics. One viable but
  storage-heavy design is independent alignment shards followed by original-
  order concatenation into one global dupsifter/sort pass. Quantify the large
  intermediate peak before implementation; do not deduplicate shards
  independently because duplicate signatures can cross shard boundaries.

- **Overlap downstream work safely.** Cell splitting and pileup cannot start
  until a final coordinate-sorted/deduplicated BAM exists today. A redesign
  around deterministic coordinate partitions might overlap stages, but it must
  preserve mate-aware duplicate calls and cell-level record conservation.

## Validation gaps and guardrails

- Run the integrated pileup equivalence on hg38 as well as mm10 and on cells
  containing unusual contigs and zero CpG calls.
- Add an injected upstream BISCUIT failure test for the direct bgzip pipeline;
  per-worker `pipefail` is present, but the failure path should be exercised.
- Install compatible ASan/UBSan runtime libraries or use a container/compiler
  with matching runtimes. Valgrind passed; sanitizer binaries could not run on
  the current host because the required runtime libraries were absent.
- Recheck `tools/slurm_dna_capacity.sh` before every production submission.
  The observed QOS changed during this work; the 2026-08-20 closeout snapshot
  was 144 CPUs/1,787 GiB after memory had reported 1,340 GiB earlier, so neither
  value is a constant. Keep six CPUs and 40 GiB reserved unless there is a
  measured reason to change those margins.
- Keep benchmark outputs under a benchmark-owned scratch root. Delete large
  raw/sorted VCF and duplicate BAM fixtures after preserving summaries,
  checksums, command lines, and Slurm accounting; never clean production paths.

## Resume checklist

1. Read `benchmarks/dna_speed/README.md` and inspect the branch diff/commit.
2. Check live production jobs, QOS CPU/memory limits, node CPU counts, and free
   scratch space.
3. Rebuild all three native helpers and rerun native failure tests.
4. Use fixed inputs and change one variable per benchmark.
5. Require exact outputs where the tool contract is deterministic. Across
   BISCUIT thread counts, compare primary/aggregate metrics because secondary
   output can be nondeterministic.
6. Update both performance documents with accepted and rejected results.
