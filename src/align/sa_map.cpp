// sa_map.cpp - minimizer seeding + banded extension for the CPU aligner.

#include "sa_map.h"

#include <algorithm>
#include <cstring>
#include <unordered_map>

#include <htslib/sam.h>  // BAM_CMATCH / BAM_CINS / BAM_CDEL CIGAR op codes

static inline uint64_t sa_hash64(uint64_t x) {
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

// Encode a read into 2-bit codes (0..3, 4=N). Returns codes in `out`.
static void encode_read(const char *seq, int len, std::vector<uint8_t> &out) {
  out.resize(len);
  for (int i = 0; i < len; ++i) out[i] = sa_base_code(seq[i]);
}

static void revcomp_codes(const std::vector<uint8_t> &in,
                          std::vector<uint8_t> &out) {
  int n = (int)in.size();
  out.resize(n);
  for (int i = 0; i < n; ++i) {
    uint8_t b = in[n - 1 - i];
    out[i] = (b > 3) ? 4 : (uint8_t)(3 - b);
  }
}

// Collect minimizer (hash, read_pos) seeds for one strand of an encoded read.
static void read_minimizers(const std::vector<uint8_t> &codes, int k, int w,
                            std::vector<std::pair<uint64_t, uint32_t>> &out) {
  out.clear();
  int n = (int)codes.size();
  if (n < k) return;
  const uint64_t kmask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
  // Reused across reads (window size w is fixed for the run): avoid a per-call
  // allocation. Indexed modulo w, so any extra capacity is harmless.
  static thread_local std::vector<std::pair<uint64_t, uint32_t>> win;
  if ((int)win.size() < w) win.resize(w);
  uint64_t kmer = 0;
  int valid = 0, wfill = 0;
  uint32_t last = UINT32_MAX;
  for (int i = 0; i < n; ++i) {
    uint8_t b = codes[i];
    if (b > 3) {
      valid = 0;
      wfill = 0;
      kmer = 0;
      continue;
    }
    kmer = ((kmer << 2) | b) & kmask;
    if (++valid < k) continue;
    uint32_t kpos = (uint32_t)(i - k + 1);
    uint64_t h = sa_hash64(kmer);
    win[wfill % w] = {h, kpos};
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
    if (best_p != last) {
      out.emplace_back(best_h, best_p);
      last = best_p;
    }
  }
}

// Banded global alignment of read codes vs a reference window (ref codes),
// returning CIGAR (htslib-packed), edit distance and score. The read is aligned
// end-to-end; the reference window is [ref, ref+rlen). `band` bounds indel size.
// Linear-ish gap model (gap_open applied per gap run via affine approximation
// is omitted for simplicity: we use a fixed per-base gap penalty).
struct BandResult {
  int score;
  int nm;
  int ref_consumed;  // ref bases consumed from window start to alignment end
  int ref_start;     // offset within window where alignment begins
  std::vector<uint32_t> cigar;
};

static BandResult banded_align(const uint8_t *qry, int qn, const uint8_t *ref,
                               int rn, const SaOpt &opt) {
  // DP over (i=query 0..qn, j=ref 0..rn) restricted to |j - i| <= band.
  // We want a semi-global alignment: full query, free choice of ref start/end
  // within the window (so leading/trailing ref is not penalised). Implemented
  // as: first row free (ref start), last query row picks best ref end.
  const int band = opt.band;
  const int NEG = -1000000000;
  // Banded semi-global DP. The matrix is logically (qn+1) x (rn+1) but only a
  // diagonal band of half-width `band` is ever consulted (qn ~ read length,
  // rn ~ qn + 2*band). We therefore (a) reuse thread-local buffers across reads
  // instead of allocating per call, and (b) write ONLY the cells the recurrence
  // and traceback actually read: row 0, the per-row band [lo..hi], and the two
  // just-outside-band border cells that the recurrence reads as NEG sentinels.
  //
  // This is bit-identical to a full NEG-initialised matrix: every cell read by
  // the recurrence is either inside a computed band, on row 0 (all 0), or one of
  // the explicit NEG borders below; out-of-band cells are never read and never
  // win a max, so their stale contents (from a prior read) cannot affect H, TB,
  // the chosen best_j, or the traceback.
  const int rowstride = rn + 1;
  static thread_local std::vector<int> Hbuf;
  static thread_local std::vector<uint8_t> TBbuf;  // 0=diag 1=up(del) 2=left(ins)
  const size_t need = (size_t)(qn + 1) * rowstride;
  if (Hbuf.size() < need) Hbuf.resize(need);
  if (TBbuf.size() < need) TBbuf.resize(need);
  int *H = Hbuf.data();
  uint8_t *TB = TBbuf.data();
  auto Hidx = [&](int i, int j) { return (size_t)i * rowstride + j; };

  // Loop-invariant penalties / scores hoisted out of the inner cell loop.
  const int gap = opt.gap_open + opt.gap_ext;
  const int match = opt.match;
  const int mismatch = opt.mismatch;

  // First query row: aligning 0 query bases; any ref prefix is free (score 0).
  for (int j = 0; j <= rn; ++j) {
    H[Hidx(0, j)] = 0;
    TB[Hidx(0, j)] = 3;  // start
  }
  for (int i = 1; i <= qn; ++i) {
    int lo = std::max(0, i - band);
    int hi = std::min(rn, i + band);
    // Query base for this row is constant across the inner j loop; `qc <= 3`
    // implies the original `ref==qc` check already covers `ref <= 3`.
    const uint8_t qc = qry[i - 1];
    const bool qok = qc <= 3;
    const int *Hprev = H + (size_t)(i - 1) * rowstride;
    int *Hcur = H + (size_t)i * rowstride;
    uint8_t *TBcur = TB + (size_t)i * rowstride;
    // Cell just below the band (read as `left` at j==lo when lo>=1): sentinel.
    if (lo >= 1) Hcur[lo - 1] = NEG;
    // column 0 (no ref consumed) for this query row: must insert i query bases
    if (lo == 0) {
      Hcur[0] = -(opt.gap_open + opt.gap_ext * i);
      TBcur[0] = 2;  // ins (query base unaligned to ref)
    }
    for (int j = std::max(1, lo); j <= hi; ++j) {
      int sc = (qok && ref[j - 1] == qc) ? match : -mismatch;
      int diag = Hprev[j - 1] + sc;
      int up = Hprev[j] - gap;       // query base, no ref => I
      int left = Hcur[j - 1] - gap;  // ref base, no query => D
      int best = diag;
      uint8_t tb = 0;
      if (up > best) { best = up; tb = 2; }
      if (left > best) { best = left; tb = 1; }
      Hcur[j] = best;
      TBcur[j] = tb;
    }
    // Cell just above the band (read as `up` by row i+1 at j==hi+1): sentinel.
    if (hi + 1 <= rn) Hcur[hi + 1] = NEG;
  }
  // Find best ref end on the last query row.
  int best_j = -1, best_score = NEG;
  int lo = std::max(0, qn - band), hi = std::min(rn, qn + band);
  for (int j = lo; j <= hi; ++j) {
    if (H[Hidx(qn, j)] > best_score) {
      best_score = H[Hidx(qn, j)];
      best_j = j;
    }
  }
  BandResult R;
  R.score = best_score;
  R.nm = 0;
  R.ref_consumed = best_j;
  // Traceback from (qn, best_j) to query row 0.
  std::vector<uint32_t> ops;  // temporary, will reverse
  int i = qn, j = best_j;
  int run_len = 0, run_op = -1;  // op: 0=M 1=I 2=D
  auto push = [&](int op) {
    if (op == run_op) {
      run_len++;
    } else {
      if (run_op >= 0 && run_len > 0) {
        int hop = (run_op == 0) ? BAM_CMATCH
                                : (run_op == 1 ? BAM_CINS : BAM_CDEL);
        ops.push_back(((uint32_t)run_len << 4) | hop);
      }
      run_op = op;
      run_len = 1;
    }
  };
  while (i > 0) {
    uint8_t tb = TB[Hidx(i, j)];
    if (tb == 0) {  // diagonal: M (match/mismatch)
      if (!(qry[i - 1] <= 3 && ref[j - 1] <= 3 && qry[i - 1] == ref[j - 1]))
        R.nm++;
      push(0);
      --i; --j;
    } else if (tb == 2) {  // up: query base consumed, no ref => insertion
      R.nm++;
      push(1);
      --i;
    } else if (tb == 1) {  // left: ref consumed, no query => deletion
      R.nm++;
      push(2);
      --j;
    } else {
      break;  // start
    }
  }
  // flush last run
  if (run_op >= 0 && run_len > 0) {
    int hop = (run_op == 0) ? BAM_CMATCH : (run_op == 1 ? BAM_CINS : BAM_CDEL);
    ops.push_back(((uint32_t)run_len << 4) | hop);
  }
  R.ref_start = j;  // ref offset where alignment begins
  std::reverse(ops.begin(), ops.end());
  R.cigar = std::move(ops);
  return R;
}

// One-strand attempt: seed, pick best diagonal cluster, banded-extend.
static SaAln align_one_strand(const SaIndex &idx, const SaOpt &opt,
                              const std::vector<uint8_t> &codes, bool rev,
                              int &out_best_diag_votes) {
  SaAln a;
  // Reused scratch (cleared on use): avoids per-read allocations. The diagonal
  // vote maps below are intentionally left as fresh per-call unordered_maps:
  // their iteration order decides the winning bin on vote ties (which sets the
  // alignment window), so reusing/reserving them could change results.
  static thread_local std::vector<std::pair<uint64_t, uint32_t>> mm;
  read_minimizers(codes, idx.k, idx.w, mm);
  if (mm.empty()) return a;

  // Collect (diag, ref_pos, read_pos). diag = ref_pos - read_pos.
  std::unordered_map<int64_t, int> diag_votes;
  std::unordered_map<int64_t, uint32_t> diag_minref;
  static thread_local std::vector<uint32_t> hits;
  const int64_t binw = opt.band > 0 ? opt.band : 1;  // hoisted bucket divisor
  for (auto &s : mm) {
    hits.clear();
    idx.query(s.first, hits);
    if ((int)hits.size() > opt.max_seed_occ) continue;  // repetitive seed
    for (uint32_t rp : hits) {
      int64_t diag = (int64_t)rp - (int64_t)s.second;
      // bucket diagonals into bins of width (band) to tolerate small indels
      int64_t bin = diag / binw;
      diag_votes[bin]++;
      auto it = diag_minref.find(bin);
      if (it == diag_minref.end() || rp < it->second) diag_minref[bin] = rp;
    }
  }
  if (diag_votes.empty()) return a;

  // Pick the top-2 diagonal bins by votes.
  int64_t best_bin = 0, second_bin = 0;
  int best_v = 0, second_v = 0;
  for (auto &kv : diag_votes) {
    if (kv.second > best_v) {
      second_v = best_v; second_bin = best_bin;
      best_v = kv.second; best_bin = kv.first;
    } else if (kv.second > second_v) {
      second_v = kv.second; second_bin = kv.first;
    }
  }
  out_best_diag_votes = best_v;

  int qn = (int)codes.size();
  // Estimate the read's leftmost ref position from the representative seed.
  int64_t approx_ref = (int64_t)diag_minref[best_bin];
  // Window: start a little before approx_ref, length qn + 2*band.
  int64_t win_start = approx_ref - opt.band;
  if (win_start < 0) win_start = 0;
  int64_t win_len = qn + 2 * opt.band;
  if (win_start + win_len > (int64_t)idx.total_len)
    win_len = (int64_t)idx.total_len - win_start;
  if (win_len < qn) return a;

  // Reject if window crosses a contig boundary at its used region.
  int cidx = idx.contig_of((uint64_t)approx_ref);
  if (cidx < 0) return a;
  const SaContig &ctg = idx.contigs[cidx];
  // Clamp window to this contig.
  int64_t cstart = ctg.offset;
  int64_t cend = ctg.offset + ctg.len;
  if (win_start < cstart) win_start = cstart;
  if (win_start + win_len > cend) win_len = cend - win_start;
  if (win_len < qn - opt.band) return a;

  BandResult br = banded_align(codes.data(), qn,
                               idx.seq.data() + win_start, (int)win_len, opt);
  if (br.cigar.empty()) return a;

  // Global leftmost ref position of the alignment.
  int64_t gpos = win_start + br.ref_start;
  int rc = idx.contig_of((uint64_t)gpos);
  if (rc != cidx) return a;  // alignment drifted out of contig

  a.mapped = true;
  a.contig = cidx;
  a.pos = gpos - ctg.offset;
  a.rev = rev;
  a.nm = br.nm;
  a.score = br.score;
  a.cigar = std::move(br.cigar);
  a.n_cand = (int)diag_votes.size();
  a.sub_score = second_v;  // proxy; refined by caller comparing strands
  return a;
}

SaAln sa_align_read(const SaIndex &idx, const SaOpt &opt, const char *seq,
                    int len) {
  static thread_local std::vector<uint8_t> fwd, rev;
  encode_read(seq, len, fwd);
  revcomp_codes(fwd, rev);

  int fv = 0, rv = 0;
  SaAln af = align_one_strand(idx, opt, fwd, false, fv);
  SaAln ar = align_one_strand(idx, opt, rev, true, rv);

  SaAln best;
  int best_score = -1000000000, sub = -1000000000;
  if (af.mapped) { best = af; best_score = af.score; }
  if (ar.mapped && ar.score > best_score) {
    sub = best_score;
    best = ar;
    best_score = ar.score;
  } else if (af.mapped && ar.mapped) {
    sub = ar.score;
  }
  if (!best.mapped) return best;

  // MAPQ heuristic: scaled by margin between best and runner-up.
  int margin = best_score - (sub > -1000000000 ? sub : 0);
  int mapq = 0;
  if (best.n_cand <= 1) {
    mapq = 60;
  } else if (margin >= best_score / 2) {
    mapq = 60;
  } else if (margin > 0) {
    mapq = 30;
  } else {
    mapq = 3;
  }
  // Downweight low-identity alignments.
  if (len > 0) {
    double idfrac = 1.0 - (double)best.nm / (double)len;
    if (idfrac < 0.8) mapq = std::min(mapq, 10);
  }
  best.mapq = mapq;
  return best;
}
