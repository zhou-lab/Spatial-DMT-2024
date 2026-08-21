# DNA pipeline performance work

This directory contains fixed-input benchmarks for the deep Spatial-DMT DNA
pipeline. Work was performed in the isolated worktree
`/home/suc3/Lab/Spatial-DMT-2024-speedopt` on branch
`perf/dna-pipeline-speed-20260818`. Benchmark outputs live under
`/scr1/users/suc3/20260818_dna_speedopt`; production outputs are never used as
benchmark destinations. The active `SpMETSLE23DM` job was not modified.

## Scheduler envelope

Checked on 2026-08-18, 2026-08-19, and 2026-08-20 with `sacctmgr`, `scontrol`,
and the partition QOS:

- partition/QOS: `defq`;
- the enforced QOS is dynamic: observed CPU limits were 361, 181, and then
  144; at 2026-08-20 00:12 the authoritative value was
  `cpu=144,mem=1787G`, after memory had changed from 1,340 GiB earlier that
  evening;
- values in parentheses in `assoc_mgr` output are current use, not limits;
- compute nodes are mostly 80 CPUs, with some 70- and 94-CPU nodes; and
- checked compute-node local `/tmp` capacity is about 400 GiB.

The benchmark session intentionally stayed at or below the requested temporary
70-CPU ceiling: 48 CPUs for live `SpMETSLE23DM` plus at most 20 benchmark CPUs.
The historical 100-CPU value was not the current enforced limit, while 144 is
the latest QOS value rather than a constant to embed in submission scripts.
Run `tools/slurm_dna_capacity.sh` immediately before each production
submission; it reads the partition QOS, current running CPUs, and a reserve,
then prints paired-rule and controller CPU/memory budgets after current job
allocations. It also prints a `single_rule_cpus_*` budget for targeted restarts
that have no concurrent matched/unmatched partner.

At the current 144-CPU cap and an otherwise empty account, the helper chooses
a 138-CPU controller ledger so the matched and unmatched 69-CPU alignment
rules can run together while retaining six spare CPUs. Use its live values:

```bash
capacity=$(tools/slurm_dna_capacity.sh)
value() { awk -F '\t' -v key="$1" '$1 == key { print $2 }' <<< "$capacity"; }
rule_cpus=$(value paired_rule_cpus_now)
controller_cpus=$(value controller_cpus_now)
controller_mem=$(value controller_memory_mb)
snakemake --cores "$controller_cpus" --jobs 4 \
  --resources cpus_per_task="$controller_cpus" mem_mb="$controller_mem" \
  --config \
    DNA_TRIM_CPUS="$rule_cpus" \
    DNA_ALIGN_CPUS="$rule_cpus" DNA_ALIGN_SORT_CPUS=16 \
    DNA_PILEUP_CPUS="$rule_cpus" \
    UNMATCHED_ALIGN_CPUS="$rule_cpus" UNMATCHED_SORT_CPUS=16
```

Heavy DNA rules declare the same allocation as native Snakemake `threads` and
as `cpus_per_task`; this prevents the controller from over-submitting while
still requesting the correct Slurm allocation. Use the helper's
`paired_rule_cpus_*` value for all three DNA CPU settings and its
`controller_cpus_*` value for both controller options and its
`controller_memory_mb` value for the global `mem_mb` resource. Do not hardcode
a multi-sample budget: the QOS can change while jobs are running. Node
availability and filesystem load can still be the practical limits.

The deepest observed old unmatched job reached 268 GB RSS at 48 CPUs. Matched
and unmatched defaults are therefore 400 GB when using 69 BISCUIT threads.
Two such rules request 800 GB, below both observed user memory caps and the
1-TB memory available on each full-node target. The capacity helper reports
the live memory limit and a `controller_memory_mb` value with 40 GiB reserved;
pass that value as the global `mem_mb` resource before production submission.

## Historical bottlenecks

Ranked by deep-sample critical-path wall time in the old implementation:

1. `dna_biscuit_align_unmatched`: `SpMETSLE21DM` timed out at 120 hours after
   processing about 3.98 of 4.35 billion individual alignment records. BISCUIT
   itself reported about 42.1 summed real hours; the remainder was dominated by
   downstream dupsifter hash collisions/backpressure and the final tail.
2. `dna_biscuit_align`: the 21DM attempt ran 51.2 hours before its sort scratch
   was externally deleted during storage cleanup. BISCUIT reported about 27.4
   summed real hours; old dedup/sort/split behavior accounted for much of the
   remaining wall time.
3. `dna_trim_and_tag`: 21DM took 42.3 wall hours and 1,498.8 CPU-hours on 48
   CPUs. It read 11.5 TB and wrote 11.3 TB because both FASTQs were fully
   decompressed into raw chunks and every Python worker started six pigz
   subprocesses.
4. `dna_biscuit_pileup`: historically hours rather than days and already
   parallel across cells, so it was not on the first optimization path.

At the 2026-08-18 12:08 snapshot, the untouched old 23DM unmatched job had
processed 3.55 of 7.84 billion individual records after about 64 wall hours,
showing the same backpressure pattern.

## Ranked changes

| Rank | Change | Difficulty | Decision | Fixed-input evidence |
|---:|---|---|---|---|
| 1 | Raise heavy rules to 69 CPUs under the current 144-CPU QOS and give Snakemake an explicit resource ledger | Easy | Adopted/configurable | A 2M-pair alignment was fastest with all 20 available cores assigned to BISCUIT; comparable historical E65_20um runs gained 1.67x from 24 to 48 threads |
| 2 | Use node-local sort scratch only when input-aware capacity checks pass, otherwise fall back to the output filesystem | Easy | Adopted | Forced-spill sort: 12 s local versus 14 s NFS, exact records |
| 3 | Precompute Hamming-1 barcode assignments and use a proof-safe exact-seeded fuzzy matcher | Moderate | Adopted | 5M pairs: matcher CPU 1,710.90 to 494.55 s (3.46x reduction); output FASTQs and stats exact |
| 4 | Stream paired FASTQs directly into bounded GNU Parallel blocks; remove raw uncompressed chunk files and compressor fan-out | Moderate | Adopted | 5M pairs: 36 s at 64 MiB blocks; all tested block sizes produced identical decompressed FASTQs and stats |
| 5 | Fix dupsifter's pathological signature hash, use 64-bit counters, and stream uncompressed BAM into sort | Moderate | Adopted | 438.8M real 21DM records: 6,474 to 1,144 s (5.66x); graceful SAM/BAM prefixes exact |
| 6 | Replace the thousands-of-output Python splitter with bounded C/HTSlib bucket fan-out and local/NFS tiering | Moderate | Adopted | 6.1M records/2,500 cells: 89.77 to 15.56 s (5.8x), exact files and records |
| 7 | Replace per-cell BAI files with CSI indexes | Easy | Adopted | A sparse cell index shrank from 1.32 MB to 7.8 KB with identical BISCUIT pileup records; a 347 MB deep-cell index shrank from 4.13 to 1.09 MB |
| 8 | Scatter BISCUIT over multiple independently aligned chunks/nodes | Hard | Deferred | Four 4-core shards took 532 s versus 401 s monolithic at the same 20-core ceiling and require large materialized intermediates |
| 9 | Stream pileup directly to bgzip/tabix, convert contexts once, and merge CpG strands without repeated reference loads | Moderate | Adopted | One deep 23DM cell: VCF generation 174 to 88 s; downstream conversion 164 to 94 s; combined 338 to 182 s (1.86x), exact final `.cg` and conversion tables |
| 10 | Build BAM indexes during sort and thread flagstat | Easy | Adopted | 2M records: sort plus index 22 to 16 s; 2.16-GB BAM flagstat 49 to 26 s, exact records/statistics |
| 11 | Distribute one sample's pileups across multiple Slurm nodes | Hard | Deferred | Current one-thread-per-cell scheduling uses a single allocation efficiently; cross-node scatter needs atomic collation and full-sample validation |

## Adopted implementation

### Trim and barcode assignment

`dna_trim_and_tag` now decompresses R1/R2 once, collapses each FASTQ record to
one tab-separated line, and streams paired records to workers. It no longer
materializes raw FASTQ chunks. Output shards are concatenated in deterministic
block order and deleted immediately.

Four allocated CPUs are reserved for the two decompressors, `paste`, and GNU
Parallel dispatch; by default a 69-CPU rule therefore starts 65 Python
workers. A two-worker diagnostic consumed more than two CPU equivalents once
stream overhead was included, so worker count must not equal the allocation.
The rule also suppresses the optional per-block observed-barcode tables; no
downstream rule consumes them, and writing then deleting thousands of sorted
tables wastes CPU, metadata operations, and space.

The default block size is adaptive by combined compressed input size:

- up to 20 GiB: 64 MiB;
- up to 100 GiB: 256 MiB; and
- over 100 GiB: 512 MiB.

On 5M fixed pairs, 64/128/256/512 MiB blocks took 36/43/51/84 seconds. The
larger deep-sample block limits avoid hundreds of thousands of tiny output
files while retaining parallel work. Synchronous Python gzip was as fast as
per-output pigz (86 versus 89 seconds in the old chunk harness) and avoids up
to six hidden two-thread compressors per worker.

At the two-worker deep-sample concurrency model, 512 MiB/1 GiB/2 GiB blocks
took 398/377/373 seconds, while maximum single-worker RSS grew from about
1.1/2.1/4.2 GiB. In the final same-job gzip-1 comparison, 512 MiB and 1 GiB
took 243 and 247 seconds with exact outputs, while 1 GiB doubled maximum
single-worker RSS. The time-first production setting therefore remains
512 MiB with a 160-GiB allocation for 65 workers.

With two workers and 512 MiB blocks (the deep-sample concurrency model), gzip
levels 1/3/6 took 256/269/316 seconds and produced 335.7/294.5/251.1 MB.
Decompressed FASTQs and aggregate stats were identical. Level 1 is the
time-first default: it was 5.1% faster than level 3 with 14.0% larger output,
a moderate transient-storage tradeoff. Suppressible observed-barcode tables
were disabled in all three cases.

The exact-seeded matcher uses the Levenshtein pigeonhole property: any match
within distance `k` retains one of `k+1` disjoint pattern seeds exactly. It
returns the unique exact hit only when every seed occurrence overlaps that hit;
otherwise it calls the original fuzzysearch implementation. The committed
exhaustive check covers 157,392 synthetic cases; 250K/5M real-pair tests also
matched the original output.

The Hamming-1 lookup is generated in reverse from whitelist entries and marks
multi-hit observations ambiguous under the existing policy. Exhaustive
generated-neighbor checks matched the original per-read search for 389,200
Spatial-DMT cases and 2,045,952 SmC-seq cases, including `N`/IUPAC fallback.

### Alignment, deduplication, and sort

Historical runs provide the evidence above the temporary benchmark ceiling.
Two comparable E65_20um unmatched runs processed 165,761,902 and 165,840,910
BISCUIT records. Job 10109990 used 24 threads and summed 256,738 CPU-seconds
over 10,828 real seconds (23.7 effective cores); job 14259413 used 48 threads
and summed 301,061 CPU-seconds over 6,467 real seconds (46.6 effective cores).
The 48-thread run was 1.67x faster for 0.05% more records, with a 17% CPU-cost
increase. This supports scaling beyond 20 threads, but 69 remains an
extrapolation because no equivalent fixed input has yet been measured above
48 threads.

The vendored MIT-licensed dupsifter source changes the hash only, not signature
equality. The old hash used a truncated portion of the first coordinate and
ignored the mate, barcode, and high coordinate bits, causing extreme collision
chains on deep inputs. The replacement hashes every equality field. Counters
are 64-bit so runs above 4.29 billion records report correctly, and
`--uncompressed-bam` avoids a compress/decompress cycle before `samtools sort`.

On a real 21DM name-grouped fragment containing 54,928,574 records:

| Case | Wall | User CPU | SAM-record MD5 |
|---|---:|---:|---|
| Production dupsifter/SAM | 338 s | 810.83 s | `d1eb4d02ec4476ab0dca8ec653625fa8` |
| Full hash/SAM | 204 s | 702.32 s | same |
| Full hash/uncompressed BAM | 187 s | 675.16 s | same |

All dupsifter statistics and flagstats also matched. A larger eight-fragment
21DM validation is recorded under `dupsifter_8chunk_v1`. That name-grouped
fixture contains 438,806,599 records and intentionally ends with one malformed
singleton. The inherited production error path called `exit()` before HTSlib
could flush buffered output. The optimized error path returns through cleanup,
rejects only that final record, and exits nonzero with the same diagnostic.

| Case | Pipeline wall | User CPU | Completed records | SAM-record MD5 |
|---|---:|---:|---:|---|
| Production hash/SAM, abrupt error | 6,474 s | 10,861.50 s | 438,806,593 | buffered-prefix artifact |
| Full hash/SAM, graceful error | 1,284 s | 5,575.66 s | 438,806,598 | `f27f2f70b5a52a5eba6fe980b6479050` |
| Full hash/uncompressed BAM, graceful error | 1,144 s | 5,171.62 s | 438,806,598 | same |

Graceful SAM and BAM flagstats and sorted records matched exactly. On the
complete 54.9M-record fixture, production and both optimized transports also
matched dupsifter statistics exactly. The final end-to-end unmatched smoke run
made the new bottleneck boundary explicit: BISCUIT used 6,610.6 CPU-seconds
over 365.5 wall seconds, while dupsifter used only 5.3 CPU-seconds over the
same stream.

Sort scratch uses local storage only when free space exceeds three times the
combined compressed FASTQ input plus 20 GiB, with a 250-GiB minimum. The 3x
estimate covers the observed 2.37x matched BAM/input ratio. If local scratch
does not fit, output preflight reserves two estimated BAM footprints on the
output filesystem so sort fragments and the final BAM can coexist. These
values are configurable as `SORT_LOCAL_INPUT_FACTOR`,
`SORT_LOCAL_HEADROOM_GB`, and `SORT_LOCAL_MIN_GB`.

Trim and alignment also check output capacity before doing expensive work.
Trim reserves combined compressed input size plus configured headroom;
alignment reserves one estimated 3x BAM footprint, plus a second footprint
when sort scratch shares the output filesystem. The bucketed splitter then
checks its actual merged BAM size and combined local/fallback/output peak.

### Per-cell BAM splitting

The compiled splitter first routes records into a bounded number of
deterministic BAM buckets, then splits buckets concurrently. This avoids about
2,500 simultaneous BGZF writers/HTS thread pools. The wrapper reserves 20 GiB
headroom, uses local scratch first, can tier excess buckets to NFS, deletes each
bucket after consumption, and checks combined peak space when fallback and
final outputs share a filesystem.

A forced local/fallback smoke test and the 6.1M-record benchmark matched the
Python splitter's filename set and every SAM record exactly.

### Per-cell indexes

Per-cell BAMs now use CSI rather than BAI. BISCUIT requires an index even when
processing a whole BAM, but HTSlib can load either format. On a representative
133,865-byte cell BAM, the BAI was 1,320,312 bytes and the standard-min-shift
CSI was 7,791 bytes. BISCUIT pileup took 29.97 and 29.47 seconds respectively,
and non-header VCF records had the same MD5. A separate four-thread CSI-only
run reproduced non-header VCF MD5 `abb6b6097dca85e8a14389b3a3a6011b`.
On the largest retained 23DM cell
(346,653,614-byte BAM), CSI was 1,089,972 bytes versus 4,131,064 bytes for BAI
and both indexes produced identical `samtools idxstats`. CSI construction was
7.49 versus 6.80 seconds, a small local cost in exchange for much lower index
storage and NFS writes. The existing 2,500 23DM BAI files occupy 10.06 GB.

### Pileup transport and conversion

Independent cells scale better than extra threads within one BISCUIT pileup.
Across 12 large 23DM cells under a fixed two-CPU allocation, two one-thread
cells took 627 seconds while one two-thread cell at a time took 625 seconds.
Outputs were identical, but the two-thread layout used more memory and has a
longer full-sample tail, so one thread per cell remains the default. Four
unequal deep cells took 151/146/147 seconds with 100-kb/1-Mb/5-Mb dispatch
windows; 1 Mb was fastest, while 5 Mb raised maximum RSS by 27% over 1 Mb.

The old rule materialized plain VCF, naturally sorted it, then compressed it.
On four 23DM cells, streaming sorted VCF took 148 seconds, materializing and
sorting took 174 seconds, and direct BISCUIT-to-bgzip/tabix took 88 seconds.
Direct VCF text has reference-contig order rather than natural chromosome
order, but tabix accepted every output and canonical records matched. The CpG
path still explicitly coordinate-sorts before the sorted BED intersection.

The conversion path now calls `biscuit vcf2bed -e -t c` once, tees that context
stream to the CA/CC/CG/CT accumulator, and uses a small context-aware merger
for the CpG strand pair. This removes a second VCF decode and repeated whole
reference loads by `biscuit mergecg`. On deep cell 0501:

| Path | Wall | `.cg` | conversion table |
|---|---:|---|---|
| Sorted VCF, two vcf2bed passes, BISCUIT merge | 164 s | reference | reference |
| Sorted VCF, one pass, BISCUIT merge | 105 s | exact | exact |
| Direct VCF, one pass, BISCUIT merge | 104 s | exact | exact |
| Direct VCF, one pass, context merger | 94 s | exact | exact |

Including VCF generation, the old and adopted paths took 338 and 182 seconds,
respectively: 46.2% less wall time, or 1.86x throughput. The actual integrated
Snakefile rule completed the same retained 23DM cell in 208 seconds, produced
no raw VCF, passed tabix validation, and matched fixed legacy `.cg` and
conversion outputs byte for byte. The native merger also passed malformed
input, `/dev/full`, synthetic C/G-pair, real-input, and Valgrind checks.

### Short alignment and BAM-tail measurements

On the fixed 10M-pair real 21DM unmatched subset, assigning all available CPUs
to BISCUIT and allowing the lightweight sort stages to share the cpuset was
fastest:

| BISCUIT threads | Sort threads | Wall | Max RSS | Primary records | Primary mapped |
|---:|---:|---:|---:|---:|---:|
| 20 | 8 | 2,542 s | 51.1 GiB | 20,000,000 | 14,459,424 |
| 16 | 4 | 2,851 s | 42.1 GiB | 20,000,000 | 14,459,496 |
| 12 | 8 | 3,717 s | 35.0 GiB | 20,000,000 | 14,459,454 |

Reserving four CPUs for sort was 12.2% slower; reserving eight was 46.2%
slower. BISCUIT used 46,907 CPU-seconds over 2,514 real seconds in the
20-thread case, while dupsifter consumed only 42 CPU-seconds and the complete
sort tail was about 25 seconds. This supports keeping alignment as the CPU
priority. BAM record checksums can differ when the BISCUIT thread count
changes because its multithreaded secondary-hit output is nondeterministic.
All cases retained exactly 20 million primary records, and primary-mapped
counts differed by at most 72 reads (0.0005%); primary and aggregate metrics
are therefore the comparison contract for this matrix.

Building a BAM index during `samtools sort` reduced a fixed 2M-record sort plus
index from 22 to 16 seconds with exact records and idxstats. Threading flagstat
reduced a 2.16-GB BAM scan from 49 to 26 seconds with byte-identical output.
Lower BAM compression was deliberately not adopted: levels 1 and 3 reduced a
standalone 18-second sort to 9 and 11 seconds, but enlarged every retained BAM
by 28.5% and 17.3%. Since sort contributes only seconds after a 42-minute
alignment and the larger BAM must be split, scanned, and piled up, that setting
does not optimize total pipeline time.

## End-to-end regression

A clean fixed 5M-pair 21DM subset completed every optimized front-end rule
from source at 20 CPUs. Native `threads` intentionally serialized the three
20-CPU branches under a 20-CPU controller ledger:

| Rule | Wall | Total CPU | Max RSS |
|---|---:|---:|---:|
| Trim/tag, auto 64 MiB blocks | 45s | 7m57s | 2.74 GiB |
| Matched align/dedup/sort/split/CSI | 4m28s | 1h02m30s | 36.26 GiB |
| Unmatched align/dedup/sort/index | 6m23s | 1h50m54s | 40.13 GiB |
| Lambda align/dedup/sort/index | 19s | 33s | below accounting resolution |

All six decompressed FASTQs, trim statistics, and matched/unmatched dupsifter
statistics matched the fixed baseline checksums. All 2,500 cell BAMs passed
`samtools quickcheck`; merged and summed per-cell counts were both 5,274,032
records. The unmatched BAM contained 3,850,680 records, of which 3,728,358
were primary, and lambda contained 735,690 records. The run produced 2,500 CSI
indexes, zero BAI cell indexes, and no optional per-block barcode tables.

Native failure tests also passed: dupsifter returned nonzero for `/dev/full`
and truncated BAM input, the C splitter rejected truncated BAM input, and the
context merger rejected malformed input and failed writes.

## Completed real shallow sample

`SpMETSLE10D` completed from raw FASTQs through copied final QC on 2026-08-20
using commit `c896544`, a 20-CPU controller ledger, and 20 CPUs for each heavy
DNA rule. The controller ran for 3:24:58; the interval from the first Slurm job
starting through both alignment branches was 37:23, and unmatched diagnostics
finished at 38:58. The run processed 14,199,274 raw pairs, retained 9,168,847
pairs, produced exactly 18,337,694 primary matched-BAM records, passed
`samtools quickcheck`, and produced 2,500 cell BAMs with 2,500 CSI and zero BAI
indexes. The final QC report is present. Pileup and downstream time are excluded
from this speed comparison because this run was not the validation run for the
final pileup changes.

There is no retained pre-optimization run of the same sample. The closest valid
pre-optimization comparator is `SpMETSLB16D`, run on 2026-07-04 with 24,573,333
raw pairs and the same 2,500-cell layout. It used 48 CPUs for trim, matched
alignment, 24 CPUs for unmatched alignment, and ran matched and unmatched
branches concurrently. The table therefore reports both wall
throughput and CPU efficiency rather than treating raw wall ratios as
algorithmic gains.

| Stage | Optimized 10D | Old 16D comparator | Input-normalized result |
|---|---:|---:|---:|
| Trim/tag | 2m20s, 20 CPUs | 55m30s, 48 CPUs | 13.74x pair throughput; 3.01x pairs/CPU-second |
| Matched align/dedup/sort/split | 15m16s, 20 CPUs | 19m13s, 48 CPUs | 0.66x pair throughput; 1.22x pairs/CPU-second |
| Unmatched align/dedup/sort | 16m23s, 20 CPUs | 19m41s, 24 CPUs | 0.80x pair throughput; 0.98x pairs/CPU-second |
| Unmatched diagnostics | 1m18s, 20 CPUs | 22m16s, 2 CPUs | 11.37x pair throughput; 0.95x pairs/CPU-second |
| Trim through both alignments | 37m23s | 1h17m37s | 52% less wall; 1.20x raw-pair throughput |
| Trim through unmatched diagnostics | 38m58s | 1h38m17s | 60% less wall; 1.46x raw-pair throughput |

The real trim result confirms the optimization: despite 58% fewer allocated
CPUs, it processed raw pairs 13.74x faster than the closest old sample. The
diagnostics wall gain came almost entirely from raising allocation from two to
20 CPUs. The shallow alignment branches did not improve input-normalized wall
throughput at the tested allocations; only matched alignment showed a modest
22% CPU-efficiency gain. The deep dupsifter hash improvement is not expected to
dominate a shallow run where BISCUIT remains the main cost.

The closest old comparator shows a real 1.20x input-normalized front-end
throughput gain through both alignment branches. Extending the endpoint through
unmatched diagnostics gives 1.46x, but that larger number includes the 10x
diagnostics CPU-allocation increase. The larger old `SpMETSLE15D` comparator
gives 1.00x through alignment and 1.02x through diagnostics after raw-input
normalization, so the defensible observed whole-front-end gain is approximately
1.0-1.2x before diagnostics rather than the component-level trim gain.

The intended 144-CPU production layout would run matched and unmatched 69-CPU
branches concurrently. Those production-count gains remain projections; this
20-CPU serialized run must not be cited as validation of them. Pileup must be
measured separately after its final changes are applied.

## Deep-run estimate

The matcher result projects to about 160-165 CPU-hours for 21DM's 6.01 billion
raw pairs, before decompression, compression, and NFS I/O. At 65 trim workers,
the compute floor is about 2.5 hours; allow roughly 3-6 hours for deep-input
streaming and filesystem effects rather than the old 42.3 hours.

Scaling the clean unmatched smoke run's BISCUIT time and CPU use to 4.35
billion records gives a 30-33-hour compute projection at 69 CPUs. The measured
dupsifter/sort pipeline sustains about 384,000 records/s at 20 CPUs and is no
longer expected to throttle BISCUIT. Allowing for less-than-linear scaling,
sort completion, and NFS variation gives a practical 30-45-hour range for a
21DM-sized alignment branch. Matched and unmatched can run concurrently under
the 138-CPU ledger, so plan roughly 1.5-2 days for the DNA front end and about
2 days plus downstream pileup/QC. This remains an extrapolation because the
temporary 70-CPU session ceiling prevented a full-depth 69-CPU validation;
heavy shared-filesystem load can extend it.

### Time-reduction accounting

The following separates measured component speedups from whole-deep-sample
projections. Projected rows are planning ranges, not completed 21DM timings.

| Step | Old | Adopted/projected | Reduction | Main source |
|---|---:|---:|---:|---|
| Barcode matcher, 5M pairs | 1,710.9 CPU-s | 494.6 CPU-s | 3.46x | exact-seeded matching and Hamming lookup |
| 21DM trim/tag | 42.3 h observed | 3-6 h projected | 7.1-14.1x; 36-39 h | matcher, no 11-TB raw-chunk round trip, 65 vs 48 workers |
| Dupsifter, 438.8M records | 6,474 s | 1,144 s | 5.66x; 5,330 s | full-signature hash and uncompressed BAM transport |
| BAM split, 6.1M records/2,500 cells | 89.77 s | 15.56 s | 5.77x; 74.21 s | bounded bucket fan-out in C/HTSlib |
| Deep-cell VCF plus conversion | 338 s | 182 s | 1.86x; 46.2% | direct compression, one decode, context merger |
| Sort plus BAM index, 2M records | 22 s | 16 s | 1.38x; 27.3% | `samtools sort --write-index` |
| Flagstat, 2.16-GB BAM | 49 s | 26 s | 1.88x; 46.9% | threaded BAM decode |
| 21DM unmatched branch | about 131 h old completion estimate | 30-45 h at 69 CPUs | 2.9-4.4x; 86-101 h | removal of hash backpressure plus CPU increase |
| 21DM DNA front-end critical path | about 173 h old completion estimate | 33-51 h | 3.4-5.2x; 122-140 h | trim plus concurrent matched/unmatched branches |

The CPU contribution should not be confused with the algorithmic contribution.
The dupsifter, splitter, matcher, and pileup figures above are fixed-allocation
algorithm/data-path gains. Historical comparable unmatched runs improved 1.67x
when BISCUIT increased from 24 to 48 threads. Increasing 48 to 69 threads has
an ideal ceiling of 1.44x but is not measured at fixed input; the practical
deep estimate allows for sublinear scaling. The 20-versus-16 fixed-input test
improved only 1.12x for 1.25x more alignment threads, while 20 versus 12
improved 1.46x for 1.67x more threads. This is why the final 30-45-hour range
remains wider than an ideal CPU-ratio calculation.

## Benchmark scripts

- `benchmark_unmatched.sh`: BISCUIT/sort allocations and scratch placement.
- `benchmark_sort_scratch.sh`: forced sort spills on NFS versus local scratch.
- `benchmark_trim_compression.sh`: compressor and worker fan-out.
- `benchmark_trim_streaming.sh`: paired streaming trim.
- `benchmark_trim_blocks.sh`: streaming block-size matrix.
- `benchmark_trim_levels.sh`: gzip-level matrix with exact decompressed checks.
- `check_matcher_equivalence.py`: exhaustive exact-seeded/fuzzysearch check.
- `benchmark_dupsifter_fixed.sh`: exact old/new dupsifter comparison.
- `benchmark_dupsifter_multichunk.sh`: larger globally name-sorted real data.
- `benchmark_dupsifter_boundary_fixture.sh`: real 21DM boundary-prefix test.
- `benchmark_cb_splitter.sh`: Python, direct C, and bucketed C splitters.
- `benchmark_pileup_parallelism.sh`: fixed-CPU cells/jobs/thread matrix.
- `benchmark_pileup_transport.sh`: raw/sorted/direct VCF transport matrix.
- `check_pileup_transport_downstream.sh`: one/two-pass and merge-backend
  equivalence through packed CpG output.
- `validate_integrated_pileup_cell.sh`: actual Snakefile deep-cell regression.
- `benchmark_sort_compression.sh`: BAM compression speed/space tradeoff.
- `benchmark_sort_inline_index.sh`: separate versus inline indexing.
- `benchmark_flagstat_threads.sh`: BAM scan threading.
- `check_barcode_lookup_equivalence.py`: exhaustive Hamming-1 lookup check.
- `run_e2e_pipeline.sh`: fixed shallow end-to-end rule regression.
- `validate_e2e_pipeline.sh`: checks checksums, BAM/index integrity, and record
  conservation across all three DNA branches.
- `test_native_tool_failures.sh`: failed-output and truncated-input handling.

Every adopted data-path optimization was required to match decompressed FASTQ
records or SAM records plus rule-level statistics, as applicable.
