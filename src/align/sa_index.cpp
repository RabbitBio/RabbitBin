// sa_index.cpp - reference loading + minimizer table for the CPU aligner.

#include "sa_index.h"

#include <algorithm>
#include <cstdio>
#include <cstring>

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
  (void)threads;
  k = k_;
  w = w_;
  if (k < 8 || k > 31) {
    fprintf(stderr, "[Error!] aligner k must be in [8,31] (got %d)\n", k);
    return false;
  }

  gzFile fp = gzopen(fasta.c_str(), "rb");
  if (!fp) {
    fprintf(stderr, "[Error!] cannot open reference FASTA: %s\n",
            fasta.c_str());
    return false;
  }

  // Stream the FASTA, concatenating bases into `seq` and recording contigs.
  seq.clear();
  seq.reserve(1u << 20);
  contigs.clear();
  const int BUF = 1 << 20;
  std::vector<char> buf(BUF);
  std::string cur_name;
  bool in_header = false;
  uint64_t cur_off = 0;
  int n;
  while ((n = gzread(fp, buf.data(), BUF)) > 0) {
    for (int i = 0; i < n; ++i) {
      char c = buf[i];
      if (c == '>') {
        // close previous contig
        if (!contigs.empty()) {
          contigs.back().len = (uint32_t)(seq.size() - contigs.back().offset);
        }
        in_header = true;
        cur_name.clear();
        cur_off = seq.size();
        continue;
      }
      if (in_header) {
        if (c == '\n' || c == '\r') {
          in_header = false;
          // name = up to first whitespace
          size_t sp = cur_name.find_first_of(" \t");
          if (sp != std::string::npos) cur_name.resize(sp);
          contigs.push_back({cur_name, cur_off, 0});
        } else {
          cur_name.push_back(c);
        }
        continue;
      }
      if (c == '\n' || c == '\r' || c == ' ' || c == '\t') continue;
      seq.push_back(sa_base_code(c));
    }
  }
  gzclose(fp);
  if (!contigs.empty())
    contigs.back().len = (uint32_t)(seq.size() - contigs.back().offset);
  total_len = seq.size();

  if (contigs.empty() || total_len == 0) {
    fprintf(stderr, "[Error!] reference FASTA had no sequence\n");
    return false;
  }
  fprintf(stderr,
          "[rabbitbin align] reference: %zu contigs, %llu bp (k=%d w=%d)\n",
          contigs.size(), (unsigned long long)total_len, k, w);

  // Collect open minimizers over the forward strand.
  // Rolling k-mer hash; a minimizer is the minimum hash in each window of w
  // consecutive k-mers. We emit a seed whenever the window minimum changes.
  std::vector<std::pair<uint64_t, uint32_t>> seeds;
  seeds.reserve(total_len / (w > 0 ? w : 1) + 16);

  const uint64_t kmask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
  // Ring buffer of the last w k-mer (hash,pos) for the window minimum.
  std::vector<std::pair<uint64_t, uint32_t>> win(w);
  uint64_t kmer = 0;
  int valid = 0;  // consecutive valid (non-N) bases
  int wfill = 0;  // how many k-mers are in the window ring
  uint32_t last_emitted = UINT32_MAX;

  for (uint64_t i = 0; i < total_len; ++i) {
    uint8_t b = seq[i];
    if (b > 3) {  // N resets the rolling k-mer
      valid = 0;
      wfill = 0;
      kmer = 0;
      continue;
    }
    kmer = ((kmer << 2) | b) & kmask;
    if (++valid < k) continue;
    uint32_t kpos = (uint32_t)(i - k + 1);  // k-mer start
    uint64_t h = sa_hash64(kmer);
    // push into window ring
    int slot = wfill % w;
    win[slot] = {h, kpos};
    ++wfill;
    if (wfill < w) continue;
    // window full: find minimum over the w entries
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
      seeds.emplace_back(best_h, best_p);
      last_emitted = best_p;
    }
  }

  fprintf(stderr, "[rabbitbin align] %zu minimizer seeds; sorting\n",
          seeds.size());
  std::sort(seeds.begin(), seeds.end());

  mm_hash.resize(seeds.size());
  mm_pos.resize(seeds.size());
  for (size_t i = 0; i < seeds.size(); ++i) {
    mm_hash[i] = seeds[i].first;
    mm_pos[i] = seeds[i].second;
  }
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
  // mm_hash is sorted; find the equal range.
  size_t lo = std::lower_bound(mm_hash.begin(), mm_hash.end(), hash) -
              mm_hash.begin();
  for (size_t i = lo; i < mm_hash.size() && mm_hash[i] == hash; ++i)
    out.push_back(mm_pos[i]);
}
