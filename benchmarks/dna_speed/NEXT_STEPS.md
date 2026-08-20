# DNA speed optimization: resume notes

This file records work intentionally left after the 2026-08-20 closeout. The
implemented and measured changes are summarized in `README.md`. The pileup
follow-up is on branch `perf/pileup-speed-20260820` in
`/home/suc3/Lab/Spatial-DMT-2024-pileupopt`; it extends
`perf/dna-pipeline-speed-20260818`. Do not modify or cancel an older production
job merely to obtain benchmark capacity.

## Completed pileup follow-up (2026-08-20)

- Replaced per-cell dense CpG/reference intersections with an 87.5-MB shared
  mmap coordinate index and direct BISCUIT-VCF to YAME format-3 conversion.
  Deep conversion fell from 94 to 13.44 seconds; shallow conversion fell from
  25.28 to 0.09 seconds. Real outputs are exact.
- Added direct multi-sample VCF conversion and deterministic size-balanced
  multi-BAM pileup. Defaults are 16 cells, a 384-MiB aggregate BAM target, two
  BISCUIT threads per batch, and one converter thread.
- Removed new intermediate VCFs. BISCUIT streams into a staging directory;
  outputs publish only after `pipefail`, dimension, and sample-set checks pass.
  Valid legacy VCFs are still consumed on restart.
- Fixed two-core direct batches reduced wall/CPU by 38%/37% for four 40-MB
  cells and 21%/22% for four 80-MB cells. A 1.03-GB joint batch was beyond the
  useful range after its 95-second conversion, which is why byte limiting is
  required.
- Completed native synthetic tests, exact 16-sample and deep-cell regressions,
  a real four-cell Snakemake run, a 1.22-second restart, legacy-VCF recovery,
  and injected upstream-failure validation. See `pileup_20260820.tsv`.

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

3. **Scatter planned pileup batches across two Slurm allocations.** The new
   runner creates restart-safe content-keyed manifests but still executes them
   inside one monolithic Snakemake rule/allocation. Refactor manifest planning
   into a checkpoint and submit moderate batch groups, not 2,500 cell jobs.
   Under a free 144-CPU QOS, two 69-CPU allocations can consume all planned
   batches. Preserve staged outputs, deterministic final cell order, and one
   aggregate `yame index -s`.

4. **Replace the full-sample pileup model with production accounting.** Run all
   2,500 cells for shallow 10D and one 23DM-class sample at 48 or 69 CPUs. Record
   per-batch GNU Parallel timing, Slurm CPU/RSS, filesystem throughput, restart
   time, and exact aggregate outputs. The current model projects historical
   23DM from 3:37 at 48 CPUs to about 32 minutes at 48 or 23 minutes at 69;
   shared NFS load is not represented. Re-test 256/384/512-MiB byte targets if
   the model misses materially.

5. **Reduce multi-sample VCF parsing overhead.** A compressed 16-sample VCF
   converts in 9.39 seconds, while a 998-MB four-deep-cell VCF takes 95 seconds.
   Profile field splitting and per-sample `CV`/`BT` extraction. Promising work
   includes parsing the fixed BISCUIT FORMAT layout without repeated scans,
   flushing completed chromosomes incrementally, and adding a compact BISCUIT
   output mode or callback that bypasses VCF text entirely. Every path must
   remain exact on C/G strand rounding and conversion summaries.

6. **Benchmark BAM compression end to end, not at sort alone.** Compression
   levels 1/3 cut standalone sort time by 50%/39% but enlarge BAMs by 28.5%/
   17.3%; they were not adopted because sort is a tiny tail and larger BAMs can
   slow split, flagstat, and pileup on shared storage. Time all downstream
   consumers on equivalent deep BAMs. Adopt a lower level only if total wall
   time improves under representative NFS load; moderate extra storage is
   acceptable, but another storage exhaustion is not.

7. **Validate the time-first trim defaults at full concurrency.** The adopted
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
- Add an automated runner fixture around the already validated injected
  BISCUIT failure and legacy-VCF recovery paths. Failed staging directories are
  deliberately retained for diagnosis; add a scoped cleanup command only
  after proving it cannot touch current manifests or completed cell outputs.
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

1. Read `benchmarks/dna_speed/README.md`, `pileup_20260820.tsv`, and inspect the
   branch diff/commit.
2. Check live production jobs, QOS CPU/memory limits, node CPU counts, and free
   scratch space.
3. Rebuild all native helpers, run `make test` in `tools/cpg_to_yame`, run the
   planner unit test, and rerun native failure tests.
4. Use fixed inputs and change one variable per benchmark.
5. Require exact outputs where the tool contract is deterministic. Across
   BISCUIT thread counts, compare primary/aggregate metrics because secondary
   output can be nondeterministic.
6. Update both performance documents with accepted and rejected results.
