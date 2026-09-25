#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <htslib/sam.h>

using RbFragmentCounts = std::array<uint64_t, 4>; // left/right all, left/right unique

// Non-overlapping sequence halves with the SAME per-sample edge trim as the
// ordinary depth calculation, applied separately to each half. Count reference
// bases covered by M/= /X only; deletions and reference skips are not coverage.
inline std::array<uint64_t, 2> rb_fragment_bases(const bam1_t *b,
                                               int32_t length, int edge) {
  std::array<uint64_t, 2> out{{0, 0}};
  const int32_t mid = length / 2;
  const int32_t lo[2] = {edge, mid + edge};
  const int32_t hi[2] = {mid - edge, length - edge};
  int64_t pos = b->core.pos;
  const uint32_t *cigar = bam_get_cigar(b);
  for (uint32_t k = 0; k < b->core.n_cigar; ++k) {
    const int op = bam_cigar_op(cigar[k]);
    const int64_t count = bam_cigar_oplen(cigar[k]);
    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
      for (int half = 0; half < 2; ++half) {
        const int64_t overlap = std::min<int64_t>(pos + count, hi[half]) -
                                std::max<int64_t>(pos, lo[half]);
        if (overlap > 0) out[half] += (uint64_t)overlap;
      }
    }
    if (bam_cigar_type(op) & 2) pos += count;
  }
  return out;
}
