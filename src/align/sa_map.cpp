// sa_map.cpp - minimizer seeding + banded extension for the CPU aligner.

#include "sa_map.h"

#include <algorithm>
#include <cstring>

#include "../phmap.h"  // parallel-hashmap: flat_hash_map (open addressing)

#include <htslib/sam.h>  // BAM_CMATCH / BAM_CINS / BAM_CDEL CIGAR op codes

// sa_hash64 + syncmer/randstrobe extraction come from sa_seed.h (via sa_index.h)
// so reference-build and read-query seeding are byte-identical.

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

// Collect randstrobe (hash, read_pos) seeds for one strand of an encoded read.
// read_pos is strobe1's k-mer start, matching the reference anchor so that
// diag = ref_pos - read_pos is identical to the minimizer path. Uses the EXACT
// same syncmer + linking code as the index build (sa_seed.h), so a read seed
// equals the reference seed at the matching locus.
static void read_randstrobes(const std::vector<uint8_t> &codes,
                             const SaRsParams &rs,
                             std::vector<std::pair<uint64_t, uint32_t>> &out) {
  out.clear();
  static thread_local std::vector<std::pair<uint64_t, uint32_t>> sync;
  static thread_local std::vector<uint64_t> ring;
  sync.clear();
  uint64_t n = (uint64_t)codes.size();
  sa_collect_syncmers(codes.data(), 0, 0, n, rs.k, rs.s, sync, ring);
  if (sync.empty()) return;
  sa_link_randstrobes(sync.data(), sync.size(), 0, sync.size(), rs, out);
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

// One strand's seed/vote result and the reference window to banded-extend. The
// expensive banded DP (which makes a cache-cold access to a random ~200 bp slice
// of the multi-GB reference) is deferred to extend_plan() so the batched driver
// can prefetch many windows before any is consumed -- hiding that DRAM latency.
struct StrandPlan {
  bool has = false;        // a window worth extending was found
  int64_t win_start = 0;   // global ref offset of the DP window
  int win_len = 0;         // window length
  int cidx = -1;           // contig index of the window
  int n_cand = 0;          // # distinct diagonal bins (for MAPQ)
  int second_v = 0;        // runner-up diagonal votes (for MAPQ)
  // Anchor of the winning diagonal: read base `anchor_q` maps to global ref
  // position `anchor_r`. Lets extend_plan try a vectorizable ungapped placement
  // (read base 0 -> ref anchor_r - anchor_q) before the scalar banded DP.
  int64_t anchor_r = 0;
  int anchor_q = 0;
  int qlen = 0;            // read length (for prefetching the ungapped span)
};

// Phase 1: vote on diagonals from already-looked-up seed ranges and compute the
// reference window to extend. No DP, no reference access -- pure work on the
// (warm) vote maps -- so the only cold access left for phase 2 is the window,
// which the caller prefetches. Voting is byte-for-byte the old per-seed path.
static void vote_and_plan(const SaIndex &idx, const SaOpt &opt, int qn,
                          const std::vector<std::pair<uint64_t, uint32_t>> &mm,
                          const uint32_t *rlo, const uint32_t *rhi,
                          StrandPlan &plan) {
  plan = StrandPlan{};
  if (mm.empty()) return;

  static thread_local phmap::flat_hash_map<int64_t, int> diag_votes;
  static thread_local phmap::flat_hash_map<int64_t, uint32_t> diag_minref;
  static thread_local phmap::flat_hash_map<int64_t, uint32_t> diag_minq;
  diag_votes.clear();
  diag_minref.clear();
  diag_minq.clear();
  const int64_t binw = opt.band > 0 ? opt.band : 1;  // hoisted bucket divisor
  const uint32_t *pos = idx.mm_pos.data();
  for (size_t si = 0; si < mm.size(); ++si) {
    uint32_t lo = rlo[si], hi = rhi[si];
    if ((int)(hi - lo) > opt.max_seed_occ) continue;  // repetitive seed
    uint32_t rpos = mm[si].second;
    for (uint32_t j = lo; j < hi; ++j) {
      uint32_t rp = pos[j];
      int64_t diag = (int64_t)rp - (int64_t)rpos;
      int64_t bin = diag / binw;
      diag_votes[bin]++;
      auto it = diag_minref.find(bin);
      if (it == diag_minref.end() || rp < it->second) {
        diag_minref[bin] = rp;
        diag_minq[bin] = rpos;  // read pos paired with this min ref pos
      }
    }
  }
  if (diag_votes.empty()) return;

  // Pick the top-2 diagonal bins by votes.
  int64_t best_bin = 0;
  int best_v = 0, second_v = 0;
  for (auto &kv : diag_votes) {
    if (kv.second > best_v) {
      second_v = best_v;
      best_v = kv.second;
      best_bin = kv.first;
    } else if (kv.second > second_v) {
      second_v = kv.second;
    }
  }

  int64_t approx_ref = (int64_t)diag_minref[best_bin];
  int64_t win_start = approx_ref - opt.band;
  if (win_start < 0) win_start = 0;
  int64_t win_len = qn + 2 * opt.band;
  if (win_start + win_len > (int64_t)idx.total_len)
    win_len = (int64_t)idx.total_len - win_start;
  if (win_len < qn) return;

  int cidx = idx.contig_of((uint64_t)approx_ref);
  if (cidx < 0) return;
  const SaContig &ctg = idx.contigs[cidx];
  int64_t cstart = ctg.offset;
  int64_t cend = ctg.offset + ctg.len;
  if (win_start < cstart) win_start = cstart;
  if (win_start + win_len > cend) win_len = cend - win_start;
  if (win_len < qn - opt.band) return;

  plan.has = true;
  plan.win_start = win_start;
  plan.win_len = (int)win_len;
  plan.cidx = cidx;
  plan.n_cand = (int)diag_votes.size();
  plan.second_v = second_v;
  plan.anchor_r = approx_ref;
  plan.anchor_q = (int)diag_minq[best_bin];
  plan.qlen = qn;
}

// Phase 2: banded-extend a planned window (its reference bytes should already be
// prefetched by the batched driver). Identical DP/result to the old fused path.
// Try a vectorizable ungapped placement: read base 0 -> ref (anchor_r-anchor_q),
// full read laid straight on the reference (all-M). For the common case of a
// read with only substitutions (no indel) at its anchored locus this IS the
// banded-DP optimum, produced ~10x cheaper and WITHOUT touching the 38 KB DP
// matrix (cutting both compute and the H-matrix bandwidth that bounds align).
// Returns true and fills `a` when accepted; false to fall back to banded DP.
// Gated on high identity so a gapped alignment cannot score higher in accepted
// cases; ambiguous reads fall through to the exact scalar path.
static bool try_ungapped(const SaIndex &idx, const SaOpt &opt,
                         const std::vector<uint8_t> &codes, bool rev,
                         const StrandPlan &plan, SaAln &a) {
  const int qn = (int)codes.size();
  const int64_t rs = plan.anchor_r - plan.anchor_q;  // read[0] -> ref rs
  if (rs < 0 || rs + qn > (int64_t)idx.total_len) return false;
  const SaContig &ctg = idx.contigs[plan.cidx];
  if (rs < (int64_t)ctg.offset ||
      rs + qn > (int64_t)(ctg.offset + ctg.len))
    return false;  // would cross a contig boundary; let banded handle it

  const uint8_t *q = codes.data();
  const uint8_t *t = idx.seq.data() + rs;
  // Vectorizable mismatch count (N in read or ref counts as a mismatch, exactly
  // as the banded DP's diagonal term scores a non-equal/ambiguous pair). Early
  // exit once we exceed the acceptance bound keeps this cheap on noisy reads.
  int nm = 0;
  for (int i = 0; i < qn && nm <= 1; ++i) nm += (q[i] != t[i]) | (q[i] > 3);
  // Accept only nm <= 1: a 1 bp gap costs gap_open+gap_ext and can recover at
  // most (match+mismatch) per fixed mismatch, so with <= 1 mismatch NO gapped
  // alignment can outscore the straight all-M placement -- this is provably the
  // banded-DP optimum, so the result is equivalent and binning is unaffected.
  // (Reads with >= 2 mismatches, where an indel might win, fall back to DP.)
  if (nm > 1) return false;

  a.mapped = true;
  a.contig = plan.cidx;
  a.pos = rs - (int64_t)ctg.offset;
  a.rev = rev;
  a.nm = nm;
  a.score = (qn - nm) * opt.match - nm * opt.mismatch;
  a.cigar.assign(1, ((uint32_t)qn << 4) | BAM_CMATCH);
  a.n_cand = plan.n_cand;
  a.sub_score = plan.second_v;
  return true;
}

static SaAln extend_plan(const SaIndex &idx, const SaOpt &opt,
                         const std::vector<uint8_t> &codes, bool rev,
                         const StrandPlan &plan) {
  SaAln a;
  if (!plan.has) return a;
  if (try_ungapped(idx, opt, codes, rev, plan, a)) return a;  // fast path
  int qn = (int)codes.size();
  BandResult br = banded_align(codes.data(), qn,
                               idx.seq.data() + plan.win_start, plan.win_len,
                               opt);
  if (br.cigar.empty()) return a;

  int64_t gpos = plan.win_start + br.ref_start;
  int rc = idx.contig_of((uint64_t)gpos);
  if (rc != plan.cidx) return a;  // alignment drifted out of contig

  const SaContig &ctg = idx.contigs[plan.cidx];
  a.mapped = true;
  a.contig = plan.cidx;
  a.pos = gpos - ctg.offset;
  a.rev = rev;
  a.nm = br.nm;
  a.score = br.score;
  a.cigar = std::move(br.cigar);
  a.n_cand = plan.n_cand;
  a.sub_score = plan.second_v;  // proxy; refined below comparing strands
  return a;
}

// Combine the two strands' alignments and assign MAPQ (identical to the old
// sa_align_read tail).
static SaAln combine_strands(SaAln af, SaAln ar, int len) {
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
  if (len > 0) {
    double idfrac = 1.0 - (double)best.nm / (double)len;
    if (idfrac < 0.8) mapq = std::min(mapq, 10);
  }
  best.mapq = mapq;
  return best;
}

// Prefetch the reference bytes of a planned window into cache (the read-side
// codes and DP buffers are already warm; only the random ref slice is cold).
static inline void prefetch_window(const SaIndex &idx, const StrandPlan &plan) {
  if (!plan.has) return;
  // Prefetch the union of the ungapped span [anchor_r-anchor_q, +qlen] (the
  // common fast path) and the banded window [win_start, +win_len] (fallback).
  int64_t rs = plan.anchor_r - plan.anchor_q;
  int64_t lo = rs < plan.win_start ? rs : plan.win_start;
  int64_t hi = rs + plan.qlen;
  int64_t wend = plan.win_start + plan.win_len;
  if (wend > hi) hi = wend;
  if (lo < 0) lo = 0;
  if (hi > (int64_t)idx.total_len) hi = (int64_t)idx.total_len;
  const uint8_t *p = idx.seq.data();
  for (int64_t o = lo; o < hi; o += 64) __builtin_prefetch(p + o);
}

// Seed + look up + vote + plan one read's window (both strands). Fills planF /
// planR and prefetches both windows. mmf/mmr/hbuf/rlo/rhi are caller-owned
// scratch (reused across reads). fcodes/rcodes receive the encoded read (needed
// later by extend_plan), so they must persist until phase 2.
static void plan_one_read(const SaIndex &idx, const SaOpt &opt, const char *seq,
                          int len, std::vector<uint8_t> &fcodes,
                          std::vector<uint8_t> &rcodes,
                          std::vector<std::pair<uint64_t, uint32_t>> &mmf,
                          std::vector<std::pair<uint64_t, uint32_t>> &mmr,
                          std::vector<uint64_t> &hbuf,
                          std::vector<uint32_t> &rlo,
                          std::vector<uint32_t> &rhi, StrandPlan &planF,
                          StrandPlan &planR) {
  encode_read(seq, len, fcodes);
  revcomp_codes(fcodes, rcodes);
  if (idx.use_randstrobe) {
    read_randstrobes(fcodes, idx.rs, mmf);
    read_randstrobes(rcodes, idx.rs, mmr);
  } else {
    read_minimizers(fcodes, idx.k, idx.w, mmf);
    read_minimizers(rcodes, idx.k, idx.w, mmr);
  }
  const int nf = (int)mmf.size(), nr = (int)mmr.size(), nt = nf + nr;
  if ((int)hbuf.size() < nt) { hbuf.resize(nt); rlo.resize(nt); rhi.resize(nt); }
  for (int i = 0; i < nf; ++i) hbuf[i] = mmf[i].first;
  for (int i = 0; i < nr; ++i) hbuf[nf + i] = mmr[i].first;
  if (nt > 0) idx.query_ranges_batch(hbuf.data(), nt, rlo.data(), rhi.data());
  vote_and_plan(idx, opt, len, mmf, rlo.data(), rhi.data(), planF);
  vote_and_plan(idx, opt, len, mmr, rlo.data() + nf, rhi.data() + nf, planR);
  prefetch_window(idx, planF);
  prefetch_window(idx, planR);
}

void sa_align_reads_batch(const SaIndex &idx, const SaOpt &opt,
                          const char *const *seqs, const int *lens, int n,
                          SaAln *out) {
  if (n <= 0) return;
  // Sub-batch size: large enough that each read's window prefetch (issued in
  // phase 1) has ~S reads of other work before phase 2 consumes it, hiding the
  // DRAM latency; small enough that the per-slot scratch stays in cache.
  const int S = 256;

  // Per-slot scratch that must persist between phase 1 and phase 2.
  static thread_local std::vector<std::vector<uint8_t>> fcodes, rcodes;
  static thread_local std::vector<StrandPlan> planF, planR;
  if ((int)fcodes.size() < S) {
    fcodes.resize(S);
    rcodes.resize(S);
    planF.resize(S);
    planR.resize(S);
  }
  // Transient per-read scratch (reused within phase 1).
  static thread_local std::vector<std::pair<uint64_t, uint32_t>> mmf, mmr;
  static thread_local std::vector<uint64_t> hbuf;
  static thread_local std::vector<uint32_t> rlo, rhi;

  for (int base = 0; base < n; base += S) {
    int s = std::min(S, n - base);
    // Phase 1: seed + vote + plan + prefetch windows for the whole sub-batch.
    for (int t = 0; t < s; ++t) {
      int r = base + t;
      plan_one_read(idx, opt, seqs[r], lens[r], fcodes[t], rcodes[t], mmf, mmr,
                    hbuf, rlo, rhi, planF[t], planR[t]);
    }
    // Phase 2: extend (windows now warm) + combine strands.
    for (int t = 0; t < s; ++t) {
      int r = base + t;
      SaAln af = extend_plan(idx, opt, fcodes[t], false, planF[t]);
      SaAln ar = extend_plan(idx, opt, rcodes[t], true, planR[t]);
      out[r] = combine_strands(std::move(af), std::move(ar), lens[r]);
    }
  }
}

SaAln sa_align_read(const SaIndex &idx, const SaOpt &opt, const char *seq,
                    int len) {
  SaAln a;
  const char *s = seq;
  sa_align_reads_batch(idx, opt, &s, &len, 1, &a);
  return a;
}
