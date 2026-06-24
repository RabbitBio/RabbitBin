// sa_index.cpp - reference loading + minimizer table for the CPU aligner.

#include "sa_index.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <vector>

#include <omp.h>
#include <parallel/algorithm>

#include <zlib.h>

// sa_hash64 + syncmer/randstrobe extraction come from sa_seed.h (via sa_index.h)
// so the reference build and the per-read query share byte-identical seeding.

bool SaIndex::build(const std::string &fasta, int k_, int w_, int threads) {
  if (threads < 1) threads = 1;
  k = k_;
  w = w_;
  rs.k = k;  // strobe k mirrors the index k
  if (k < 8 || k > 31) {
    fprintf(stderr, "[Error!] aligner k must be in [8,31] (got %d)\n", k);
    return false;
  }
  if (use_randstrobe && (rs.s < 1 || rs.s >= k)) {
    fprintf(stderr, "[Error!] randstrobe s must be in [1,k-1] (got %d, k=%d)\n",
            rs.s, k);
    return false;
  }

  // Read the whole reference into memory, then parse + 2-bit encode in parallel.
  // The previous version did a single-threaded char-by-char pass (~13s on a
  // 4 Gbp reference). We split the raw bytes at LINE boundaries -- each chunk
  // begins right after a '\n', where the FASTA header state is always reset, so
  // no record/line spans a split -- run the EXACT same state machine per chunk,
  // then merge in file order. The resulting `seq` and `contigs` are byte-for-
  // byte identical to the serial parse (verified end-to-end against the BAM).
  seq.clear();
  contigs.clear();
  std::vector<char> raw;
  {
    gzFile fp = gzopen(fasta.c_str(), "rb");
    if (!fp) {
      fprintf(stderr, "[Error!] cannot open reference FASTA: %s\n",
              fasta.c_str());
      return false;
    }
    const size_t CH = 32u << 20;
    size_t used = 0;
    for (;;) {
      raw.resize(used + CH);
      int n = gzread(fp, raw.data() + used, (unsigned)CH);
      if (n <= 0) break;
      used += (size_t)n;
    }
    gzclose(fp);
    raw.resize(used);
  }

  const size_t N = raw.size();
  int nchunks = threads * 4;
  if (nchunks < 1) nchunks = 1;
  if (N > 0 && (size_t)nchunks > N) nchunks = (int)N;
  if (nchunks < 1) nchunks = 1;
  std::vector<size_t> bnd(nchunks + 1);
  bnd[0] = 0;
  bnd[nchunks] = N;
  for (int c = 1; c < nchunks; ++c) {
    size_t p = (size_t)((unsigned long long)N * (unsigned)c / (unsigned)nchunks);
    while (p < N && raw[p - 1] != '\n') ++p;  // snap forward to a line start
    bnd[c] = p;
  }

  struct ChunkRes {
    std::vector<uint8_t> enc;
    std::vector<std::pair<std::string, uint64_t>> hdrs;  // (name, local offset)
  };
  std::vector<ChunkRes> cres(nchunks);

#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (int c = 0; c < nchunks; ++c) {
    size_t s = bnd[c], e = bnd[c + 1];
    if (s >= e) continue;
    ChunkRes &cr = cres[c];
    cr.enc.reserve(e - s);
    bool in_header = false;
    std::string name;
    uint64_t cur_off = 0;
    for (size_t i = s; i < e; ++i) {
      char ch = raw[i];
      if (ch == '>') {  // header start (matches serial: fires even in-header)
        in_header = true;
        name.clear();
        cur_off = cr.enc.size();
        continue;
      }
      if (in_header) {
        if (ch == '\n' || ch == '\r') {
          in_header = false;
          size_t sp = name.find_first_of(" \t");
          if (sp != std::string::npos) name.resize(sp);
          cr.hdrs.emplace_back(name, cur_off);
        } else {
          name.push_back(ch);
        }
        continue;
      }
      if (ch == '\n' || ch == '\r' || ch == ' ' || ch == '\t') continue;
      cr.enc.push_back(sa_base_code(ch));
    }
  }

  // Prefix-sum per-chunk encoded lengths -> global base offsets.
  std::vector<uint64_t> base(nchunks + 1, 0);
  for (int c = 0; c < nchunks; ++c) base[c + 1] = base[c] + cres[c].enc.size();
  total_len = base[nchunks];

  // Assemble seq (parallel copy) and contigs (serial, in file order).
  seq.resize(total_len);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (int c = 0; c < nchunks; ++c) {
    if (!cres[c].enc.empty())
      memcpy(seq.data() + base[c], cres[c].enc.data(), cres[c].enc.size());
  }
  for (int c = 0; c < nchunks; ++c)
    for (auto &h : cres[c].hdrs)
      contigs.push_back({h.first, base[c] + h.second, 0});
  for (size_t i = 0; i < contigs.size(); ++i) {
    uint64_t end = (i + 1 < contigs.size()) ? contigs[i + 1].offset : total_len;
    contigs[i].len = (uint32_t)(end - contigs[i].offset);
  }
  std::vector<char>().swap(raw);  // free the raw bytes

  if (contigs.empty() || total_len == 0) {
    fprintf(stderr, "[Error!] reference FASTA had no sequence\n");
    return false;
  }
  fprintf(stderr,
          "[rabbitbin align] reference: %zu contigs, %llu bp (k=%d w=%d)\n",
          contigs.size(), (unsigned long long)total_len, k, w);

  std::vector<std::pair<uint64_t, uint32_t>> seeds;
  if (use_randstrobe) {
    // ----- Randstrobe seeding (strobealign-style) -----------------------
    // Two parallel passes over the forward reference:
    //   (1) extract open syncmers into a position-sorted array (chunked with a
    //       warm-up margin so each syncmer is emitted exactly once, identical to
    //       a serial scan);
    //   (2) link each syncmer to a second strobe in its window -> randstrobe
    //       {hash, strobe1_pos}. Pass 2 is a pure function of the syncmer array
    //       (embarrassingly parallel by index), so no warm-up is needed.
    // The anchor stored is strobe1's k-mer start, so diagonal voting at query
    // time is identical to the minimizer path.
    std::vector<std::pair<uint64_t, uint32_t>> sync;
    if (total_len >= (uint64_t)k) {
      const uint64_t kMargin = 4096;  // converges syncmer rolling state
      int nchunks = threads * 8;
      if (nchunks < 1) nchunks = 1;
      uint64_t chunk = (total_len + nchunks - 1) / (uint64_t)nchunks;
      if (chunk < 1) chunk = 1;
      nchunks = (int)((total_len + chunk - 1) / chunk);
      std::vector<std::vector<std::pair<uint64_t, uint32_t>>> locals(nchunks);

#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
      for (int c = 0; c < nchunks; ++c) {
        uint64_t own_lo = (uint64_t)c * chunk;
        uint64_t own_hi = own_lo + chunk;
        if (own_hi > total_len) own_hi = total_len;
        uint64_t scan_lo = (own_lo > kMargin) ? (own_lo - kMargin) : 0;
        std::vector<uint64_t> ring;
        sa_collect_syncmers(seq.data(), scan_lo, own_lo, own_hi, rs.k, rs.s,
                            locals[c], ring);
      }
      // Concatenate in chunk (= position) order -> sorted by position.
      size_t nsync = 0;
      for (auto &v : locals) nsync += v.size();
      sync.reserve(nsync);
      for (auto &v : locals) {
        sync.insert(sync.end(), v.begin(), v.end());
        std::vector<std::pair<uint64_t, uint32_t>>().swap(v);
      }
    }
    fprintf(stderr, "[rabbitbin align] %zu syncmers; linking randstrobes\n",
            sync.size());

    // Pass 2: link randstrobes (parallel by syncmer index range).
    if (!sync.empty()) {
      int nchunks = threads * 8;
      if (nchunks < 1) nchunks = 1;
      size_t chunk = (sync.size() + nchunks - 1) / (size_t)nchunks;
      if (chunk < 1) chunk = 1;
      nchunks = (int)((sync.size() + chunk - 1) / chunk);
      std::vector<std::vector<std::pair<uint64_t, uint32_t>>> locals(nchunks);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
      for (int c = 0; c < nchunks; ++c) {
        size_t lo = (size_t)c * chunk;
        size_t hi = lo + chunk;
        if (hi > sync.size()) hi = sync.size();
        locals[c].reserve(hi - lo);
        sa_link_randstrobes(sync.data(), sync.size(), lo, hi, rs, locals[c]);
      }
      size_t total_seeds = 0;
      for (auto &v : locals) total_seeds += v.size();
      seeds.reserve(total_seeds);
      for (auto &v : locals) {
        seeds.insert(seeds.end(), v.begin(), v.end());
        std::vector<std::pair<uint64_t, uint32_t>>().swap(v);
      }
    }
    std::vector<std::pair<uint64_t, uint32_t>>().swap(sync);  // free syncmers
    fprintf(stderr, "[rabbitbin align] %zu randstrobe seeds; sorting\n",
            seeds.size());
  } else {
    // ----- Open minimizer seeding (legacy) ------------------------------
    // A seed is emitted whenever the minimum-hash k-mer in the sliding window
    // of w consecutive k-mers changes. Parallelised by partitioning the base
    // axis into chunks that each emit window events for the positions they OWN,
    // while scanning a `kMargin`-base warm-up prefix so the rolling k-mer /
    // window / dedup state converges to exactly what a single sequential pass
    // would hold on entry to the owned range. Bit-identical to a serial scan.
    const uint64_t kmask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
    if (total_len >= (uint64_t)k) {
      const uint64_t kMargin = 4096;
      int nchunks = threads * 8;
      if (nchunks < 1) nchunks = 1;
      uint64_t chunk = (total_len + nchunks - 1) / (uint64_t)nchunks;
      if (chunk < 1) chunk = 1;
      nchunks = (int)((total_len + chunk - 1) / chunk);

      std::vector<std::vector<std::pair<uint64_t, uint32_t>>> locals(nchunks);

#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
      for (int c = 0; c < nchunks; ++c) {
        uint64_t own_lo = (uint64_t)c * chunk;
        uint64_t own_hi = own_lo + chunk;
        if (own_hi > total_len) own_hi = total_len;
        uint64_t scan_lo = (own_lo > kMargin) ? (own_lo - kMargin) : 0;

        auto &out = locals[c];
        out.reserve((own_hi - own_lo) / (w > 0 ? w : 1) + 16);

        std::vector<std::pair<uint64_t, uint32_t>> win(w);
        uint64_t kmer = 0;
        int valid = 0;   // consecutive valid (non-N) bases
        int wfill = 0;   // how many k-mers are in the window ring
        uint32_t last_emitted = UINT32_MAX;

        for (uint64_t i = scan_lo; i < own_hi; ++i) {
          uint8_t b = seq[i];
          if (b > 3) {  // N resets the rolling k-mer (last_emitted persists)
            valid = 0;
            wfill = 0;
            kmer = 0;
            continue;
          }
          kmer = ((kmer << 2) | b) & kmask;
          if (++valid < k) continue;
          uint32_t kpos = (uint32_t)(i - k + 1);  // k-mer start
          uint64_t h = sa_hash64(kmer);
          int slot = wfill % w;
          win[slot] = {h, kpos};
          ++wfill;
          if (wfill < w) continue;
          uint64_t best_h = UINT64_MAX;
          uint32_t best_p = 0;
          int start = wfill - w;
          for (int j = 0; j < w; ++j) {
            auto &e = win[(start + j) % w];
            if (e.first < best_h) {
              best_h = e.first;
              best_p = e.second;
            }
          }
          if (best_p != last_emitted) {
            if (i >= own_lo) out.emplace_back(best_h, best_p);
            last_emitted = best_p;
          }
        }
      }

      size_t total_seeds = 0;
      for (auto &v : locals) total_seeds += v.size();
      seeds.reserve(total_seeds);
      for (auto &v : locals) {
        seeds.insert(seeds.end(), v.begin(), v.end());
        std::vector<std::pair<uint64_t, uint32_t>>().swap(v);  // free early
      }
    }
    fprintf(stderr, "[rabbitbin align] %zu minimizer seeds; sorting\n",
            seeds.size());
  }
  omp_set_num_threads(threads);
  __gnu_parallel::sort(seeds.begin(), seeds.end());

  mm_hash.resize(seeds.size());
  mm_pos.resize(seeds.size());
  const long long ns = (long long)seeds.size();
#pragma omp parallel for num_threads(threads) schedule(static)
  for (long long i = 0; i < ns; ++i) {
    mm_hash[i] = seeds[i].first;
    mm_pos[i] = seeds[i].second;
  }

  // Build the prefix-bucket directory (top mm_pbits bits of the hash). Chosen so
  // there are ~10 entries per bucket on a full reference, turning query()'s deep
  // binary search over the multi-GB table into a single bucket lookup + a probe
  // of one or two cache lines. Bit-exact: it only narrows the search interval.
  mm_pbits = 26;
  mm_shift = 64 - mm_pbits;
  const size_t nb = ((size_t)1 << mm_pbits);
  mm_bucket.assign(nb + 1, 0);
  size_t bidx = 0;
  for (size_t i = 0; i < (size_t)ns; ++i) {
    uint64_t pre = mm_hash[i] >> mm_shift;
    while (bidx <= pre) mm_bucket[bidx++] = (uint32_t)i;
  }
  while (bidx <= nb) mm_bucket[bidx++] = (uint32_t)ns;
  return true;
}

int SaIndex::contig_of(uint64_t global_pos) const {
  // contigs sorted by offset; find last contig with offset <= global_pos
  int lo = 0, hi = (int)contigs.size() - 1, ans = -1;
  while (lo <= hi) {
    int mid = (lo + hi) / 2;
    if (contigs[mid].offset <= global_pos) {
      ans = mid;
      lo = mid + 1;
    } else {
      hi = mid - 1;
    }
  }
  if (ans < 0) return -1;
  if (global_pos >= contigs[ans].offset + contigs[ans].len) return -1;
  return ans;
}

void SaIndex::query(uint64_t hash, std::vector<uint32_t> &out) const {
  // Narrow the search to this hash's prefix bucket, then binary-search that tiny
  // contiguous slice for the equal range (identical result to a full-table
  // lower_bound, just far fewer cache-cold probes).
  if (mm_bucket.empty()) {  // no directory (empty index)
    return;
  }
  uint64_t b = hash >> mm_shift;
  size_t bl = mm_bucket[b];
  size_t bh = mm_bucket[b + 1];
  size_t lo = std::lower_bound(mm_hash.begin() + bl, mm_hash.begin() + bh,
                               hash) -
              mm_hash.begin();
  for (size_t i = lo; i < bh && mm_hash[i] == hash; ++i)
    out.push_back(mm_pos[i]);
}

void SaIndex::query_ranges_batch(const uint64_t *h, int n, uint32_t *out_lo,
                                 uint32_t *out_hi) const {
  if (n <= 0) return;
  if (mm_bucket.empty()) {  // empty index: all ranges empty
    for (int i = 0; i < n; ++i) out_lo[i] = out_hi[i] = 0;
    return;
  }
  // Per-seed bucket bounds, carried between phases (reused per thread).
  static thread_local std::vector<uint32_t> bl, bh;
  if ((int)bl.size() < n) { bl.resize(n); bh.resize(n); }

  const uint32_t *bucket = mm_bucket.data();
  const uint64_t *hash = mm_hash.data();
  const uint32_t *pos = mm_pos.data();

  // Phase 1: prefetch every seed's bucket-directory entry. All n misses are
  // issued before any is consumed, so they overlap in the memory system.
  for (int i = 0; i < n; ++i)
    __builtin_prefetch(&bucket[h[i] >> mm_shift]);

  // Phase 2: read bucket bounds (now warm) and prefetch each hash slice.
  for (int i = 0; i < n; ++i) {
    uint64_t b = h[i] >> mm_shift;
    uint32_t lo = bucket[b], hi = bucket[b + 1];
    bl[i] = lo;
    bh[i] = hi;
    if (lo < hi) __builtin_prefetch(&hash[lo]);
  }

  // Phase 3: search each (warm) slice for the equal range; prefetch positions.
  for (int i = 0; i < n; ++i) {
    uint32_t lo = bl[i], hi = bh[i];
    uint64_t key = h[i];
    // lower_bound over the tiny prefix-bucket slice (usually ~1 cache line).
    while (lo < hi) {
      uint32_t mid = lo + ((hi - lo) >> 1);
      if (hash[mid] < key) lo = mid + 1; else hi = mid;
    }
    uint32_t e = lo;
    uint32_t top = bh[i];
    while (e < top && hash[e] == key) ++e;
    out_lo[i] = lo;
    out_hi[i] = e;
    if (lo < e) __builtin_prefetch(&pos[lo]);
  }
}
