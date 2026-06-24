#ifndef RB_SA_SEED_H_
#define RB_SA_SEED_H_

// Shared seed extraction for the CPU aligner: open syncmers + randstrobes.
//
// Both the reference index build (sa_index.cpp) and the per-read query
// (sa_map.cpp) include this header so the SAME syncmer/randstrobe logic and
// hash run on both sides -- a randstrobe generated for a read is byte-for-byte
// identical to the one the reference produced at the matching locus, which is
// what makes the (hash -> ref_pos) table lookups hit.
//
// Design notes (strobealign-derived, adapted to RabbitBin's forward-only,
// non-canonical, re-seed-the-revcomp model):
//   * Syncmers are OPEN syncmers over the FORWARD strand: a k-mer is a syncmer
//     when its minimum s-mer hash sits at the centre offset (k-s)/2 of the
//     k-mer (t = (k-s)/2 + 1, matching strobealign's t_syncmer). Non-canonical
//     (forward 2-bit k-mer only) -- the reverse strand is handled, as for the
//     old minimizers, by re-seeding the reverse-complemented read.
//   * A randstrobe links syncmer i to a second syncmer chosen from the window
//     [i+w_min, i+w_max] (genomic distance <= max_dist) by minimising
//     popcount((h_i ^ h_j) & q) (strobealign "method 3'"). The randstrobe hash
//     is h_i + h_j; the ANCHOR position stored/voted is strobe1's k-mer start,
//     so diag = ref_anchor - read_anchor works exactly like a minimizer.

#include <cstdint>
#include <utility>
#include <vector>

// 64-bit integer hash (SplitMix64 finaliser); identical to the one used by the
// old minimizer path so hashing stays consistent project-wide.
static inline uint64_t sa_hash64(uint64_t x) {
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

// Randstrobe parameters (set on the index; threaded into both build & query).
struct SaRsParams {
  int k = 20;            // strobe k-mer length
  int s = 16;            // syncmer s-mer length (k-s should be even & >0)
  int w_min = 5;         // window start (in syncmer indices) for strobe2
  int w_max = 11;        // window end   (in syncmer indices) for strobe2
  int max_dist = 80;     // max genomic distance (bp) strobe1 -> strobe2
  uint64_t q = 255;      // mask for the method-3' popcount strobe2 selector
};

// Collect open syncmers from encoded bases seq[scan_lo, hi). A syncmer is
// EMITTED (appended to `out` as {kmer_hash, kmer_start_pos}) only when its
// k-mer start position lies in [own_lo, hi); the [scan_lo, own_lo) prefix only
// warms the rolling state so a chunked parallel scan converges to exactly what
// a single sequential pass would hold (each syncmer emitted by one chunk only).
// For a single read pass scan_lo == own_lo == 0 and hi == read length.
//
// `seq` codes: 0=A 1=C 2=G 3=T, anything >3 (N) resets the rolling state.
static inline void sa_collect_syncmers(
    const uint8_t *seq, uint64_t scan_lo, uint64_t own_lo, uint64_t hi, int k,
    int s, std::vector<std::pair<uint64_t, uint32_t>> &out,
    std::vector<uint64_t> &ring_scratch) {
  if (hi < (uint64_t)k) return;
  const int M = k - s + 1;  // s-mers per k-mer window
  if (M < 1) return;
  const int t = (k - s) / 2 + 1;  // syncmer position (1-based offset+1)
  const uint64_t kmask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
  const uint64_t smask = (s < 32) ? ((1ULL << (2 * s)) - 1) : ~0ULL;

  if ((int)ring_scratch.size() < M) ring_scratch.resize(M);
  uint64_t *qs = ring_scratch.data();  // ring buffer of s-mer hashes

  uint64_t kmer = 0, smer = 0;
  int valid = 0;            // consecutive valid (non-N) bases
  int head = 0, cnt = 0;    // ring buffer head + element count
  uint64_t qs_min_val = UINT64_MAX;
  long long qs_min_pos = -1;  // start position of the current minimum s-mer

  for (uint64_t i = scan_lo; i < hi; ++i) {
    uint8_t b = seq[i];
    if (b > 3) {  // N resets everything
      valid = 0;
      head = 0;
      cnt = 0;
      kmer = 0;
      smer = 0;
      qs_min_val = UINT64_MAX;
      qs_min_pos = -1;
      continue;
    }
    kmer = ((kmer << 2) | b) & kmask;
    smer = ((smer << 2) | b) & smask;
    if (++valid < s) continue;
    // s-mer ending at i, starting at i-s+1.
    uint64_t hs = sa_hash64(smer);
    if (cnt < M) {
      qs[(head + cnt) % M] = hs;
      ++cnt;
      if (cnt < M) continue;  // window not full yet
      // First full window: brute-force minimum (rightmost on ties via <=? we
      // mirror strobealign's '<' so the leftmost wins here, matched by the
      // recompute below which scans in reverse for the rightmost).
      qs_min_val = UINT64_MAX;
      for (int j = 0; j < M; ++j) {
        uint64_t v = qs[(head + j) % M];
        if (v < qs_min_val) {
          qs_min_val = v;
          qs_min_pos = (long long)i - k + j + 1;
        }
      }
    } else {
      // Window already full: evict front (s-mer starting at i-k), append new.
      bool popped_min = (qs_min_pos == (long long)i - k);
      head = (head + 1) % M;            // pop_front
      qs[(head + (M - 1)) % M] = hs;    // push_back at tail
      if (popped_min) {
        qs_min_val = UINT64_MAX;
        qs_min_pos = (long long)i - s + 1;
        for (int j = M - 1; j >= 0; --j) {  // reverse: rightmost min on ties
          uint64_t v = qs[(head + j) % M];
          if (v < qs_min_val) {
            qs_min_val = v;
            qs_min_pos = (long long)i - k + j + 1;
          }
        }
      } else if (hs < qs_min_val) {
        qs_min_val = hs;
        qs_min_pos = (long long)i - s + 1;
      }
    }
    if (qs_min_pos == (long long)i - k + t) {  // syncmer at centre offset
      uint64_t kstart = i - k + 1;
      if (kstart >= own_lo)
        out.emplace_back(sa_hash64(kmer), (uint32_t)kstart);
    }
  }
}

// Link randstrobes over a position-sorted syncmer array `sync` for indices
// [lo, hi). Appends {randstrobe_hash, strobe1_pos} to `out`. strobe1_pos is the
// anchor used for diagonal voting (k-mer start of the first strobe).
static inline void sa_link_randstrobes(
    const std::pair<uint64_t, uint32_t> *sync, size_t nsync, size_t lo,
    size_t hi, const SaRsParams &p,
    std::vector<std::pair<uint64_t, uint32_t>> &out) {
  if (nsync == 0) return;
  if (hi > nsync) hi = nsync;
  for (size_t i = lo; i < hi; ++i) {
    uint64_t h1 = sync[i].first;
    uint32_t pos1 = sync[i].second;
    size_t w_start = i + (size_t)p.w_min;
    size_t w_end = i + (size_t)p.w_max;
    if (w_end >= nsync) w_end = nsync - 1;
    uint64_t max_pos = (uint64_t)pos1 + (uint64_t)p.max_dist;
    uint64_t best_pop = ~0ULL;
    uint64_t h2 = h1;  // default: strobe2 = strobe1 (degenerate)
    for (size_t j = w_start; j <= w_end && j < nsync; ++j) {
      if ((uint64_t)sync[j].second > max_pos) break;  // beyond genomic window
      uint64_t pc = (uint64_t)__builtin_popcountll((h1 ^ sync[j].first) & p.q);
      if (pc < best_pop) {
        best_pop = pc;
        h2 = sync[j].first;
      }
    }
    out.emplace_back(h1 + h2, pos1);
  }
}

#endif  // RB_SA_SEED_H_
