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

struct SaContig {
  std::string name;
  uint64_t offset;  // start of this contig in the concatenated base array
  uint32_t len;
};

struct SaIndex {
  int k = 15;  // k-mer length for minimizers
  int w = 10;  // minimizer window

  std::vector<uint8_t> seq;       // concatenated bases: 0=A 1=C 2=G 3=T 4=N
  std::vector<SaContig> contigs;  // contig directory (sorted by offset)
  uint64_t total_len = 0;

  // Sorted minimizer table: mm_hash[i] ascending, mm_pos[i] is the ref position
  // (into `seq`) of that minimizer's k-mer start. Equal hashes are contiguous.
  std::vector<uint64_t> mm_hash;
  std::vector<uint32_t> mm_pos;

  // Build from a FASTA file. Returns false on error.
  bool build(const std::string &fasta, int k_, int w_, int threads);

  // Map a global position in `seq` to a contig index (or -1) via binary search.
  int contig_of(uint64_t global_pos) const;

  // Look up all reference positions for a minimizer hash; appends to `out`.
  void query(uint64_t hash, std::vector<uint32_t> &out) const;
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
