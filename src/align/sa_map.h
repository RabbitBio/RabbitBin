#ifndef RB_SA_MAP_H_
#define RB_SA_MAP_H_

// Per-read mapping for the CPU aligner: minimizer seeding (both strands),
// diagonal clustering, and banded extension into a CIGAR. Output is enough to
// emit a SAM/BAM record (contig, pos, strand, CIGAR, NM, MAPQ).

#include <cstdint>
#include <string>
#include <vector>

#include "sa_index.h"

struct SaAln {
  bool mapped = false;
  int contig = -1;       // contig index in SaIndex
  int64_t pos = 0;       // 0-based local position on the contig (leftmost)
  bool rev = false;      // aligned to reverse strand
  int mapq = 0;
  int nm = 0;            // edit distance (mismatches + indel bases)
  int score = 0;        // alignment score (match=1, mismatch/gap penalised)
  std::vector<uint32_t> cigar;  // htslib-packed CIGAR ops (len<<4 | op)
  // number of distinct candidate diagonals considered (for MAPQ)
  int n_cand = 0;
  int sub_score = 0;     // second-best score (for MAPQ)
};

// Tunables for the mapper (filled with defaults; CLI may override).
struct SaOpt {
  int match = 1;
  int mismatch = 4;
  int gap_open = 6;
  int gap_ext = 1;
  int band = 31;          // banded DP half-width
  int max_seed_occ = 500; // ignore minimizers occurring more than this
  int min_score_frac = 0; // (reserved)
};

// Align a single read sequence (ACGTN string) against the index. Returns the
// best alignment (mapped=false if none found).
SaAln sa_align_read(const SaIndex &idx, const SaOpt &opt, const char *seq,
                    int len);

#endif  // RB_SA_MAP_H_
