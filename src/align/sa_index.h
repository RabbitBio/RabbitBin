#ifndef RB_SA_INDEX_H_
#define RB_SA_INDEX_H_

// Self-contained CPU read-alignment index for `rabbitbin bwa` / `rabbitbin map`.
//
// Design (lightweight, optimizable; no external aligner source):
//   * Reference is loaded as a concatenated 2-bit-ish base array (0..3, 4=N).
//   * Seeds are open minimizers (k, w) of the FORWARD reference strand stored in
//     a sorted (hash, pos) table so lookups are a binary search and the layout
//     is cache friendly (CSR-like, no per-key vectors).
//   * Strand of a read is resolved at map time by seeding both the read and its
//     reverse complement against the forward table.
// Reference base array + minimizer table are the only large structures.

#include <cstdint>
#include <string>
#include <vector>

#include "sa_seed.h"  // SaRsParams, sa_hash64, syncmer/randstrobe extraction

struct SaContig {
  std::string name;
  uint64_t offset;  // start of this contig in the concatenated base array
  uint32_t len;
};

struct SaIndex {
  int k = 19;  // k-mer length (minimizers; also strobe k when use_randstrobe)
  int w = 19;  // minimizer window (unused when use_randstrobe)

  // Seeding mode. When true, the table is built from randstrobes (strobealign-
  // style) instead of open minimizers; query-time seeding switches to match.
  // The (hash -> ref_pos) table layout, query(), diagonal voting and banded DP
  // are identical either way -- only what fills the table / what a read emits
  // changes. randstrobe parameters live in `rs` (k mirrors the field above).
  bool use_randstrobe = false;
  SaRsParams rs;

  std::vector<uint8_t> seq;       // concatenated bases: 0=A 1=C 2=G 3=T 4=N
  std::vector<SaContig> contigs;  // contig directory (sorted by offset)
  uint64_t total_len = 0;

  // Sorted minimizer table: mm_hash[i] ascending, mm_pos[i] is the ref position
  // (into `seq`) of that minimizer's k-mer start. Equal hashes are contiguous.
  std::vector<uint64_t> mm_hash;
  std::vector<uint32_t> mm_pos;

  // Prefix-bucket directory over the top `mm_pbits` bits of the hash: an exact
  // acceleration structure (no effect on results) that bounds the binary search
  // in query() to a tiny contiguous slice of mm_hash instead of probing the
  // whole multi-GB table -- the dominant cache-miss cost when mapping against a
  // large reference. mm_bucket[b] = index of the first entry whose top bits >= b
  // so all entries with prefix b live in [mm_bucket[b], mm_bucket[b+1]).
  std::vector<uint32_t> mm_bucket;  // size (1<<mm_pbits)+1
  int mm_pbits = 0;
  int mm_shift = 64;

  // Build from a FASTA file. Returns false on error.
  bool build(const std::string &fasta, int k_, int w_, int threads);

  // Map a global position in `seq` to a contig index (or -1) via binary search.
  int contig_of(uint64_t global_pos) const;

  // Look up all reference positions for a minimizer hash; appends to `out`.
  void query(uint64_t hash, std::vector<uint32_t> &out) const;

  // Batched, software-pipelined lookup (latency hiding via prefetch, the
  // minibwa technique). For each of the `n` query hashes, computes the equal
  // range [lo[i], hi[i]) into mm_pos so that mm_hash[lo[i]..hi[i]) == h[i]; the
  // caller reads positions directly from mm_pos[lo[i]..hi[i]). Issuing all the
  // bucket / hash-slice / position prefetches a phase ahead keeps many random
  // memory requests in flight at once, which hides DRAM latency on the multi-GB
  // seed table (the dominant cost of seeding on a large reference). The ranges
  // are identical to calling query() per hash, so downstream voting/results are
  // unchanged. lo/hi must each have room for `n` entries.
  void query_ranges_batch(const uint64_t *h, int n, uint32_t *lo,
                          uint32_t *hi) const;
};

// 2-bit base encoding helpers.
static inline uint8_t sa_base_code(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;  // N / ambiguous
  }
}
static inline char sa_base_char(uint8_t code) {
  static const char t[5] = {'A', 'C', 'G', 'T', 'N'};
  return t[code < 5 ? code : 4];
}

#endif  // RB_SA_INDEX_H_
