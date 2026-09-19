# Graph, depth matrix and BAM input optimization checks

These measurements use MiDAS with 23 coordinate-sorted BAMs, 72 threads,
seed `1789826864`, and default binning parameters. The assembly contains
1,371,496 large and 4,635,674 small contigs. Dual depth produces 46 columns.
Input BAM files total 286.915 GiB; unmapped tails can be skipped, so file size
is not the number of bytes actually scanned.

The changes use the same settings for all inputs. They add no dataset-specific
parameters, CPU model dispatch or storage-device tuning. Existing build flags
are unchanged.

## Graph and matrix changes

Comparison of `5596fed` with `98cbf63`, one sequential run of each:

| Measurement | Before | After |
| --- | ---: | ---: |
| Complete binning | 207.88 s | 129.29 s |
| Graph pair pass | 127.3865 s | 50.5233 s |
| BAM scan and boundary merge | 41.3793 s | 40.7734 s |
| Large and small depth row copies | 2.3167 s | 0.0664 s |
| Peak RSS | 9.341 GiB | 9.134 GiB |

This pair showed a 37.8% reduction in complete runtime. It is one pair, with
uncontrolled page cache, rather than a general speedup guarantee. FASTA and
BAM work overlap, and row copying is a substage; these times are not additive.

The conservative batched abundance filter preserves the original final dot
test for survivors. PMH score caching and packed heap keys preserve score and
tie order. Depth shards write exclusively owned contig interiors directly;
boundary sums are merged afterwards. A failed shard causes its BAM's partial
writes to be cleared before full-file fallback. The matrix join verifies
reference names before using source order and retains the name-based fallback.

## Additional BAM input changes

Each active BAM scan stream requests a 1 MiB HTSlib input buffer, grouping
compressed reads across BGZF blocks. HTSlib retains responsibility for virtual
offsets, decompression, CRC checks and record decoding. The extra input buffer
is bounded per active stream, not per entire input file.

NM lookup traverses scalar and string tags in an inlined loop with bounds
checks. Arrays, unknown encodings, missing tags and malformed data fall back
to HTSlib. No tag position or ordering is assumed; duplicate NM tags retain
first-match behavior.

Compared with `98cbf63`, using the same inputs and binning parameters:

| Measurement | Before | Buffered reads + NM lookup |
| --- | ---: | ---: |
| Warm-cache BAM scan, two runs each | 21.1783 / 20.7622 s | 20.1561 / 20.2370 s |
| Warm-cache BAM scan mean | 20.9703 s | 20.1966 s |
| Input-stage read system calls | 31,967,872 | about 242,108 |
| Cold-cache BAM scan, one run each | 39.5473 s | 39.5509 s |
| Cold-cache input-stage physical reads | 232.149 GiB | 234.758 GiB |
| Cold-cache input-stage user CPU time | 1,631.88 s | 1,595.51 s |
| Cold-cache input-stage kernel CPU time | 969.42 s | 849.50 s |

The warm-cache comparison ran new/old/old/new in sequence. Physical reads were
1.89–2.16 GiB per run, with most input served from cache. Mean BAM scan time
fell about 3.7%; read calls fell about 99.2%. These are short local measurements.

For the cold-cache pair, `POSIX_FADV_DONTNEED` was applied only to the 24 input
files before each run. No system cache controls or device settings were
changed. Actual physical reads were recorded to check the cache condition.
Both cold scans took about 39.55 seconds: reduced CPU and syscall overhead did
not produce a measurable cold-scan wall-time improvement in this pair.
Buffered reads can fetch a little more data beyond the last consumed block.

The diagnostic runs stopped after depth merge, before graph pair processing.
Process CPU and I/O counts in the table cover the sampled input stage,
including the concurrent FASTA work. They must not be attributed solely to
BAM scanning. Earlier user-mode profiling placed 74.2% of sampled cycles in
libdeflate; this percentage excludes kernel execution and waiting.

## Result equivalence and regression coverage

- All 13 configured CTests passed on the final source tree.
- Both versions in the graph/matrix comparison produced 2,379 bins containing
  7,120,632,336 bp. The smallest bin was 200,110 bp. Member tables were byte
  identical; bin statistics were identical after excluding output filenames.
- Membership SHA-256:
  `0983c4ae43d5b7bda97cb1e97f883daf60a23b4dc1dfae8e81baaf0f100a51c5`.
- A complete MiDAS run with buffered reads and NM lookup also reproduced this
  membership hash and the same bin statistics. It took 104.43 s with warm
  input cache; this is an equivalence run, not a cache-matched comparison with
  the earlier 129.29 s run.
- Six mutual-neighbor graph comparisons and one union-neighbor comparison
  matched the original implementation, including tied scores. The union
  comparison also checked edge output order.
- Eight BAM comparisons cover shard boundaries, compact filtering, single
  and dual channels, reordered or extra names, duplicate names, clipping,
  CIGAR variants and an unmapped tail. A ninth comparison forces resync
  failure after partial writes and verifies complete fallback.
- NM lookup matches HTSlib's returned pointer and `errno` in 20,103 cases,
  including truncated fields, unknown encodings, arrays, missing or duplicate
  NM and long strings. The same differential test passed AddressSanitizer.

The graph tests and BAM comparisons are included under `test/`. BAM integration
tests require Python 3 and samtools; CMake registers them when those tools and
the fused BAM dependency are available. These checks compare exact results;
they do not introduce a new AMBER or CheckM2 quality evaluation.

Run the configured test suite after building:

```bash
cmake --build build -j 8
ctest --test-dir build --output-on-failure -j 1
```

For a separate-version BAM comparison:

```bash
python3 test/test_bam_depth_equivalence.py \
    build/src/rabbitbin build/depth-comparison samtools /path/to/old/rabbitbin
```

Enable phase reporting with `RB_TIMING=1 RB_DEPTH_PROF=1`; use a fixed seed and
record cache state, physical input bytes and thread count for timing comparisons.
