# Performance diagnostics

[Back to the README](../README.md)

Set `RB_TIMING=1` to report wall-clock phase timestamps and
`RB_GRAPH_PROF=1` to count candidate pairs, exact early exits and heap updates.
`RB_DEPTH_PROF=1` reports depth preparation, scanning and matrix merge phases.
FASTA/sketch construction overlaps in-process BAM depth calculation; these
times must not be added when computing the total runtime.

The fused graph's portable abundance filter batches dot products over a
transposed copy of the rank vectors. It rejects a pair only outside a
conservative floating-point error bound; every survivor still uses the
original dot calculation before PMH and top-k selection. Set
`RABBIT_NO_ABD_BATCH=1` for a differential run with the original dot path.
Existing integer/VNNI filters retain precedence when enabled. PMH score
conversion is cached for each possible integer winner-match count using the
same correction and clamp as the direct calculation. These optimizations do
not change thresholds, sample handling, random seeds or neighbour tie-breaking.

For coordinate-sorted BAM input, byte-range workers write complete interior
contigs directly to their depth column and retain only boundary contigs for
merging. Integer depth sums and per-sample floating-point normalization are
unchanged. Matrix merging verifies names before reusing reference order and
copies independent rows in parallel; duplicate names retain the original
name-based handling. The graph heap reuses its packed threshold key for
comparisons, caps allocation at the configured top-k, and sorts its existing
storage for mutual-neighbour lookup without popping and copying every edge.

BAM scan streams use a bounded 1 MiB compressed-input buffer to batch reads
across BGZF blocks. HTSlib still performs seeks, decompression, CRC checks and
record decoding. NM lookup traverses ordinary scalar/string tags in one loop;
arrays, unknown types, missing tags and malformed data use HTSlib's lookup.
Neither optimization assumes a particular tag order, dataset or storage device.
Measured results, cache conditions and equivalence checks are documented in
[the optimization validation note](performance-20260919.md).
