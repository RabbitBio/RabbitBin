// sa_index.cpp - reference loading + minimizer table for the CPU aligner.

#include "sa_index.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <vector>

#include <omp.h>
#include <parallel/algorithm>

#include <zlib.h>

// 64-bit integer hash (SplitMix64 finaliser) for k-mer hashing.
static inline uint64_t sa_hash64(uint64_t x) {
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

bool SaIndex::build(const std::string &fasta, int k_, int w_, int threads) {
  if (threads < 1) threads = 1;
  k = k_;
  w = w_;
  if (k < 8 || k > 31) {
    fprintf(stderr, "[Error!] aligner k must be in [8,31] (got %d)\n", k);
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

  // Collect open minimizers over the forward strand (a seed is emitted whenever
  // the minimum-hash k-mer in the sliding window of w consecutive k-mers
  // changes). This is parallelised by partitioning the base axis into chunks
  // that each emit the window events for the positions they OWN, while scanning
  // a `kMargin`-base warm-up prefix so the rolling k-mer / window / dedup state
  // converges to exactly what a single sequential pass would hold on entry to
  // the owned range. Emission is keyed on the window-evaluation position `i`
  // (each `i` belongs to exactly one chunk), and the per-cell logic is byte-for-
  // byte the sequential recurrence below, so the resulting multiset of seeds is
  // identical to the serial scan (order is irrelevant: the table is sorted).
  const uint64_t kmask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;

  std::vector<std::pair<uint64_t, uint32_t>> seeds;
  if (total_len >= (uint64_t)k) {
    // Warm-up length: enough to refill kmer (k) + window (w) several times over
    // so the dedup `last_emitted` state is identical on entry to the owned span.
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
          // Warm-up positions (i < own_lo) only advance dedup state; the owned
          // span records, so each event is emitted by exactly one chunk.
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
