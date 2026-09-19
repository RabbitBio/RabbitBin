// RabbitBin module: rb_graph.cpp

#include "rb_abundance_filter.h"

// ═══════════════════════════════════════════════════════════════════════════
// build_similarity_graph  –  build edge list using KmerSketch Jaccard
// ═══════════════════════════════════════════════════════════════════════════
// Reference all-pairs graph build (O(N^2) Jaccard); used when inverted index is off.
static void build_graph_allpairs(Graph &g) {
  ProgressTracker progress(nobs);
  std::vector<GraphNodeId> &from = g.from;
  std::vector<GraphNodeId> &to   = g.to;
  auto &sComp = g.sComp;

  size_t TILE = 10;
  try {
    TILE = std::max(
        (size_t)((CacheSize() * 1024.) /
                 (2 * (sketch_size / 8 + sizeof(uint64_t)) +
                  maxEdges * (2 * sizeof(size_t) + sizeof(StoredDistance)))),
        (size_t)10);
  } catch (...) {}

  verbose_message("Starting Building Similarity Graph (all-pairs reference). "
                  "TILE=%d nobs=%d maxEdges=%d\n", TILE, nobs, maxEdges);

#pragma omp parallel for schedule(dynamic, 1)                                  \
    reduction(merge_graphnode : from) reduction(merge_graphnode : to)          \
    reduction(merge_storeddist : sComp)
  for (size_t ii = 0; ii < nobs; ii += TILE) {
    std::vector<std::priority_queue<Edge, std::vector<Edge>, CompareEdge>>
        tmp_iedges(TILE);
    auto i_stop = std::min(ii + TILE, nobs);
    for (size_t jj = 0; jj < nobs; jj += TILE) {
      auto j_stop = std::min(jj + TILE, nobs);
      const size_t stride = (size_t)g_sig_nw * g_sig_np;
      for (size_t i = ii; i < i_stop; ++i) {
        auto &edges = tmp_iedges[i - ii];
        const uint64_t* sig_i = g_sig_flat.data() + i * stride;
        (void)sig_i;
        for (size_t j = jj; j < j_stop; ++j) {
          if (i == j || !is_nz(i, j)) continue;
          StoredDistance sv = (StoredDistance)graph_sim(i, j);
          if (sv > 0.0f &&
              (edges.size() < maxEdges ||
               (edges.size() == maxEdges && sv > edges.top().second))) {
            if (edges.size() == maxEdges) edges.pop();
            edges.push(std::make_pair(j, sv));
          }
        }
      }
    }
    for (size_t k = 0; k < TILE; ++k) {
      auto &edges = tmp_iedges[k];
      auto i_idx = ii + k;
      if (i_idx > i_stop) break;
      while (!edges.empty()) {
        Edge edge = edges.top();
        if (i_idx < edge.first) {
          sComp.push_back(edge.second);
          from.push_back(i_idx);
          to.push_back(edge.first);
        }
        edges.pop();
      }
    }
    if (verbose && omp_get_thread_num() == 0) { progress.track(TILE); }
  }
  verbose_message("Finished Building Similarity Graph (%d edges) "
                  "[%.1fGb / %.1fGb]                                      \n",
                  g.getEdgeCount(), getUsedPhysMem(),
                  getTotalPhysMem() / 1024 / 1024);
  g.sComp.shrink_to_fit(); g.to.shrink_to_fit(); g.from.shrink_to_fit();
}

void build_similarity_graph(Graph &g) {
  std::vector<GraphNodeId> &from = g.from;
  std::vector<GraphNodeId> &to   = g.to;
  auto &sComp = g.sComp;

  if (nobs == 0) return;

  // Default to the all-pairs b-bit popcount build: the OPH signature is so
  // compact (m/8 bytes, cache-resident) that the all-pairs SIMD popcount kernel
  // is faster than inverted-index candidate generation at these contig counts,
  // where skewed posting lists (popular (bucket,min) keys) cause a sum-of-L^2
  // traversal blow-up. Set RABBIT_GRAPH_INDEX=1 to use the inverted-index build
  // (wins when N is much larger and/or the per-pair metric is expensive).
  bool use_index = false;
  if (const char *e = rb_getenv("RABBIT_GRAPH_INDEX"))
    use_index = (e[0] == '1');
  if (!use_index) { build_graph_allpairs(g); return; }

  // The optional inverted-index path considers pairs sharing at least one
  // sketch bucket. RABBIT_MINCOMMON can request a stricter experimental
  // prefilter without changing the production default.
  const uint32_t M = g_sketches[0]->getK();
  size_t minCommon = 1;
  if (const char *e = rb_getenv("RABBIT_MINCOMMON")) {
    long v = std::atol(e);
    if (v >= 1) minCommon = (size_t)v;
  }

  // ── Use the index that was built inline during sketch construction ────────
  if (!g_inv_idx) {
    verbose_message("WARN: inverted index not ready, falling back to all-pairs\n");
    build_graph_allpairs(g); return;
  }
  const rabbit_invidx::InvertedIndex& idx = *g_inv_idx;

  verbose_message("Starting Building Similarity Graph (inverted index). "
                  "nobs=%d maxEdges=%d buckets=%d minCommon=%d "
                  "postings=%zu keys=%zu\n",
                  nobs, maxEdges, M, minCommon,
                  idx.totalPostings, idx.postIdx.size());

  // ── Posting-list traversal: count collisions per candidate, keep the
  //    top-maxEdges neighbours per node, emit i<j edges.
  //    Per-sketch keys are generated on-the-fly (getKeys) – no global skKeys.
  const uint32_t* csrPtr = idx.csrPosts.data();

  std::vector<std::vector<GraphNodeId>>    tl_from(numThreads);
  std::vector<std::vector<GraphNodeId>>    tl_to(numThreads);
  std::vector<std::vector<StoredDistance>> tl_sComp(numThreads);

  ProgressTracker progress(nobs);

#pragma omp parallel num_threads(numThreads)
  {
    const int tid = omp_get_thread_num();
    // stamp/isect: per-thread scratch indexed by candidate node id. isect is
    // the collision count, bounded by m (<= sketch size) so uint16_t suffices.
    std::vector<int>      stamp(nobs, 0);
    std::vector<uint16_t> isect(nobs, 0);
    int ep = 0;
    std::vector<size_t> cand;
    cand.reserve(4096);
    std::priority_queue<Edge, std::vector<Edge>, CompareEdge> heap;

    auto &lfrom = tl_from[tid];
    auto &lto   = tl_to[tid];
    auto &lsComp = tl_sComp[tid];
    // Per-thread scratch key buffer: keys are generated on-the-fly so we never
    // need the global N×m skKeys array.
    std::vector<uint64_t> keys;

#pragma omp for schedule(dynamic, 16)
    for (size_t i = 0; i < nobs; ++i) {
      cand.clear();
      ++ep;
      if (__builtin_expect(ep == INT_MAX, 0)) {
        std::fill(stamp.begin(), stamp.end(), 0);
        ep = 1;
      }

      // Generate this sketch's keys on-the-fly (O(m), no global storage).
      g_sketches[i]->getKeys(keys);
      const size_t ksz = keys.size();
      // Cap posting-list traversal: low-complexity 21-mer buckets can have
      // enormous lists (all N contigs), turning candidate-gen into O(N^2).
      // Skip any key whose list is larger than this fraction of N — those keys
      // carry no discriminating power for proximity anyway.
      const uint32_t plCap = (uint32_t)std::max((size_t)256, nobs / 8);
      for (size_t ki = 0; ki < ksz; ++ki) {
        auto it = idx.postIdx.find(keys[ki]);
        if (it == idx.postIdx.end()) continue;
        const uint32_t  plSz = it->second.cnt;
        if (plSz > plCap) continue;             // skip near-universal keys
        const uint32_t* pl   = csrPtr + it->second.off;
        for (uint32_t pi = 0; pi < plSz; ++pi) {
          size_t j = (size_t)pl[pi];
          if (j == i) continue;
          if (stamp[j] != ep) {
            stamp[j]  = ep;
            isect[j]  = 1;
            cand.push_back(j);
          } else {
            ++isect[j];
          }
        }
      }

      // Keep the top-maxEdges neighbours by similarity.
      // In PMH mode: score every candidate with the k=6 PMH winner-match
      // (graph_sim), not OPH Jaccard — the OPH index is used only as a
      // sparse prefilter to generate candidates, not to score them.
      while (!heap.empty()) heap.pop();
      for (size_t j : cand) {
        const size_t common = isect[j];
        if (common < minCommon) continue;
        if (!is_nz(i, j)) continue;
        StoredDistance sim = (StoredDistance)graph_sim(i, j);
        if (sim > 0.0f &&
            (heap.size() < maxEdges ||
             (heap.size() == maxEdges && sim > heap.top().second))) {
          if (heap.size() == maxEdges) heap.pop();
          heap.push(std::make_pair(j, sim));
        }
      }

      while (!heap.empty()) {
        Edge edge = heap.top();
        if (i < edge.first) {
          lsComp.push_back(edge.second);
          lfrom.push_back(i);
          lto.push_back(edge.first);
        }
        heap.pop();
      }

      if (verbose && tid == 0) {
        progress.track(numThreads);
        if (progress.isStepMarker())
          verbose_message("Building Similarity Graph %s [%.1fGb / %.1fGb]    "
                          "                \r",
                          progress.getProgress(),
                          getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
      }
    }
  }

  // ── 4. Concatenate the per-thread edge lists into the graph. ─────────────
  size_t total = 0;
  for (size_t t = 0; t < numThreads; ++t) total += tl_from[t].size();
  from.reserve(total);
  to.reserve(total);
  sComp.reserve(total);
  for (size_t t = 0; t < numThreads; ++t) {
    from.insert(from.end(), tl_from[t].begin(), tl_from[t].end());
    to.insert(to.end(),     tl_to[t].begin(),   tl_to[t].end());
    sComp.insert(sComp.end(), tl_sComp[t].begin(), tl_sComp[t].end());
  }

  verbose_message("Finished Building Similarity Graph (%d edges) "
                  "[%.1fGb / %.1fGb]                                      \n",
                  g.getEdgeCount(),
                  getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

  g.sComp.shrink_to_fit();
  g.to.shrink_to_fit();
  g.from.shrink_to_fit();
}

// ── Fused top-k graph build ─────────────────────────────────────────────
// A tiled N²/2 pass accumulates the top-maxEdges PMH neighbours for every
// contig. Candidate pairs must pass the coverage feasibility checks below;
// graph emission subsequently applies the mutual-neighbour rule.
// ── Winner-banding LSH candidate generation (RABBIT_LSH=1) ─────────────────
// The O(N²/2) all-pairs scan is compute-bound and scales out to the core count,
// so the only way to go faster is to evaluate fewer pairs.  Kept graph edges
// have a raw winner-match fraction ≳0.90 while random pairs sit at the PMH
// baseline (~0.76); banding the m winners into B bands of R rows turns that gap
// into a strong candidate filter: two contigs are candidates iff some band's R
// winners match exactly (prob s^R per band).  With R=20 a same-genome pair
// (s≈0.9) is a candidate w.h.p. while a random pair (s≈0.76) is pruned ~10×.
// Default OFF — the all-pairs path stays the reference.
static bool     rb_lsh_on() { const char *e = getenv("RABBIT_LSH"); return e && e[0]=='1'; }
static uint32_t rb_lsh_R()  { const char *e = getenv("RABBIT_LSH_R");
                              uint32_t v = e ? (uint32_t)atoi(e) : 20; return v ? v : 20; }
static size_t   rb_lsh_maxbucket() { const char *e = getenv("RABBIT_LSH_MAXBUCKET");
                              return e ? (size_t)atoll(e) : 6000; }

static inline uint64_t lsh_band_hash(const uint32_t *w, uint32_t r) {
  uint64_t h = 1469598103934665603ULL;
  for (uint32_t k = 0; k < r; ++k) { h ^= w[k]; h *= 1099511628211ULL; }
  return h;
}

#if defined(__AVX512VNNI__) && defined(__AVX512BW__)
// ── rank8 row filter ────────────────────────────────────────────────────────
// Exact Spearman screen of contig i against j ∈ [j_begin, j_stop) using the
// block-transposed int8 centred-rank layout (see gen_fused_graph).  G = number
// of 4-sample groups is a template parameter so the dpbusd chain is fully
// unrolled with the broadcast operands held in registers, and the function is
// kept out of the giant OpenMP body so the compiler does not spill its state.
// Writes j's whose exact correlation is clearly ≥ the threshold to `pass`, and
// j's inside the ±band sliver (decided by the caller with unit_dot_f) to
// `band`.  Returns the number of pruned lanes (profiling only).
struct Rank8Row {
  const int8_t *blk_base;     // r8_blk.data()
  const uint32_t *rowu;       // (r8+128) packed groups of contig i
  const float *invn;          // r8_invn.data()
  const uint16_t *nzmask;     // r8_nzmask.data(), consulted only if !nz_i
  float invn_i, lo, hi;
  bool nz_i;
};

template <int G>
static uint64_t rank8_row_filter(const Rank8Row &R, size_t j_begin,
                                 size_t j_stop, uint32_t *__restrict__ pass,
                                 size_t &npass, uint32_t *__restrict__ band,
                                 size_t &nband) {
  __m512i bc[G];
  for (int g = 0; g < G; ++g) bc[g] = _mm512_set1_epi32((int)R.rowu[g]);
  const __m512 v_lo = _mm512_set1_ps(R.lo);
  const __m512 v_hi = _mm512_set1_ps(R.hi);
  const __m512 v_ii = _mm512_set1_ps(R.invn_i);
  uint64_t pruned = 0;
  npass = 0;
  nband = 0;
  const size_t b_begin = j_begin >> 4, b_end = (j_stop + 15u) >> 4;
  for (size_t b = b_begin; b < b_end; ++b) {
    const size_t j0 = b << 4;
    __mmask16 valid = 0xFFFFu;
    if (j0 < j_begin) valid &= (__mmask16)(0xFFFFu << (j_begin - j0));
    if (j0 + 16 > j_stop) valid &= (__mmask16)((1u << (j_stop - j0)) - 1u);
    if (!R.nz_i) valid &= (__mmask16)R.nzmask[b];
    if (!valid) continue;
    const int8_t *blk = R.blk_base + b * (size_t)G * 64u;
    __m512i a0 = _mm512_setzero_si512(), a1 = a0, a2 = a0, a3 = a0;
    int g = 0;
    for (; g + 4 <= G; g += 4) {
      a0 = _mm512_dpbusd_epi32(a0, bc[g],     _mm512_loadu_si512((const void *)(blk + (size_t)g * 64u)));
      a1 = _mm512_dpbusd_epi32(a1, bc[g + 1], _mm512_loadu_si512((const void *)(blk + (size_t)(g + 1) * 64u)));
      a2 = _mm512_dpbusd_epi32(a2, bc[g + 2], _mm512_loadu_si512((const void *)(blk + (size_t)(g + 2) * 64u)));
      a3 = _mm512_dpbusd_epi32(a3, bc[g + 3], _mm512_loadu_si512((const void *)(blk + (size_t)(g + 3) * 64u)));
    }
    if (g < G)     a0 = _mm512_dpbusd_epi32(a0, bc[g],     _mm512_loadu_si512((const void *)(blk + (size_t)g * 64u)));
    if (g + 1 < G) a1 = _mm512_dpbusd_epi32(a1, bc[g + 1], _mm512_loadu_si512((const void *)(blk + (size_t)(g + 1) * 64u)));
    if (g + 2 < G) a2 = _mm512_dpbusd_epi32(a2, bc[g + 2], _mm512_loadu_si512((const void *)(blk + (size_t)(g + 2) * 64u)));
    const __m512i acc = _mm512_add_epi32(_mm512_add_epi32(a0, a1),
                                         _mm512_add_epi32(a2, a3));
    const __m512 rho = _mm512_mul_ps(
        _mm512_mul_ps(_mm512_cvtepi32_ps(acc), _mm512_loadu_ps(R.invn + j0)),
        v_ii);
    const __mmask16 ge_lo = _mm512_mask_cmp_ps_mask(valid, rho, v_lo, _CMP_GE_OQ);
    if (__builtin_expect(ge_lo == 0, 1)) {
      pruned += (uint64_t)__builtin_popcount((unsigned)valid);
      continue;
    }
    const __mmask16 ge_hi = _mm512_mask_cmp_ps_mask(ge_lo, rho, v_hi, _CMP_GE_OQ);
    pruned += (uint64_t)__builtin_popcount((unsigned)(valid & ~ge_lo));
    __mmask16 p = ge_hi, q = (__mmask16)(ge_lo & ~ge_hi);
    while (p) {
      const unsigned l = (unsigned)__builtin_ctz((unsigned)p);
      p = (__mmask16)(p & (p - 1u));
      pass[npass++] = (uint32_t)(j0 + l);
    }
    while (q) {
      const unsigned l = (unsigned)__builtin_ctz((unsigned)q);
      q = (__mmask16)(q & (q - 1u));
      band[nband++] = (uint32_t)(j0 + l);
    }
  }
  return pruned;
}

typedef uint64_t (*Rank8RowFn)(const Rank8Row &, size_t, size_t, uint32_t *,
                               size_t &, uint32_t *, size_t &);
static Rank8RowFn rank8_row_filter_for(size_t G) {
  switch (G) {
#define RB_R8_CASE(n) case n: return &rank8_row_filter<n>;
    RB_R8_CASE(1)  RB_R8_CASE(2)  RB_R8_CASE(3)  RB_R8_CASE(4)
    RB_R8_CASE(5)  RB_R8_CASE(6)  RB_R8_CASE(7)  RB_R8_CASE(8)
    RB_R8_CASE(9)  RB_R8_CASE(10) RB_R8_CASE(11) RB_R8_CASE(12)
    RB_R8_CASE(13) RB_R8_CASE(14) RB_R8_CASE(15) RB_R8_CASE(16)
    RB_R8_CASE(17) RB_R8_CASE(18) RB_R8_CASE(19) RB_R8_CASE(20)
    RB_R8_CASE(21) RB_R8_CASE(22) RB_R8_CASE(23) RB_R8_CASE(24)
    RB_R8_CASE(25) RB_R8_CASE(26) RB_R8_CASE(27) RB_R8_CASE(28)
    RB_R8_CASE(29) RB_R8_CASE(30) RB_R8_CASE(31) RB_R8_CASE(32)
#undef RB_R8_CASE
    default: return nullptr;
  }
}
#endif

static void gen_fused_graph(Graph &g) {

  // Per-contig neighbor heap (min-heap of size ≤ maxEdges; top = weakest kept).
  // Triangle pass: each unordered pair {i,j} is evaluated exactly once and used
  // to update BOTH heaps[i] and heaps[j] (graph_sim is symmetric), halving the
  // similarity computations vs the old full-N² scan.  Because every heap can now
  // be touched by multiple threads, each row gets a lock + an atomic threshold:
  // row_thresh[r] mirrors heaps[r].top() once the heap is full, so the common
  // low-similarity candidate is rejected with a single relaxed atomic read and
  // never contends for the lock.  Correctness is re-verified inside the lock.
  // Reuse the threshold's packed order in the heap itself. For graph_sim's
  // nonnegative similarities, IEEE float bits have the same order as values;
  // the complemented id retains the original smaller-id-wins tie-break.
  // Sifting therefore needs one integer comparison per level, with no float
  // comparison or second branch for equal scores. Both the kept set and the
  // pop order are unchanged. Each edge still occupies eight bytes.
  static_assert(sizeof(StoredDistance) == 4 &&
                    std::numeric_limits<StoredDistance>::is_iec559,
                "packed edge key requires an IEEE-754 32-bit similarity");
  struct Edge32 {
    uint64_t key;
    static uint64_t pack(StoredDistance sv, uint32_t id) noexcept {
      uint32_t bits;
      std::memcpy(&bits, &sv, sizeof(bits));
      return ((uint64_t)bits << 32) | (uint32_t)~id;
    }
    uint32_t id() const noexcept { return ~(uint32_t)key; }
    StoredDistance sv() const noexcept {
      const uint32_t bits = (uint32_t)(key >> 32);
      StoredDistance value;
      std::memcpy(&value, &bits, sizeof(value));
      return value;
    }
  };
  struct CompareEdge32 {
    constexpr bool operator()(Edge32 const &a, Edge32 const &b) const noexcept {
      return a.key > b.key;
    }
  };
  // Minimal binary min-heap (CompareEdge32 order) with a single-pass
  // replace_top: a kept candidate on a FULL heap used to cost pop()+push()
  // (sift-down then sift-up); overwriting the root and sifting down once does
  // the same job.  Element set and pop order are unchanged (strict total order).
  struct Heap {
    std::vector<Edge32> v;
    size_t size() const noexcept { return v.size(); }
    bool empty() const noexcept { return v.empty(); }
    const Edge32 &top() const noexcept { return v.front(); }
    void push(const Edge32 &e, size_t limit) {
      if (v.size() == v.capacity()) {
        // Grow only as far as the configured top-k. Default vector growth can
        // otherwise keep almost twice as many slots as the heap can ever use.
        const size_t next = v.empty() ? 1 :
            v.size() + std::min(v.size(), limit - v.size());
        v.reserve(next);
      }
      v.push_back(e);
      std::push_heap(v.begin(), v.end(), CompareEdge32{});
    }
    void pop() {
      std::pop_heap(v.begin(), v.end(), CompareEdge32{});
      v.pop_back();
    }
    void replace_top(const Edge32 &e) {
      static const CompareEdge32 c{};
      const size_t n = v.size();
      size_t k = 0;
      for (;;) {
        size_t l = 2 * k + 1;
        if (l >= n) break;
        size_t r = l + 1;
        // child that must move up = the one "kept-preferred" by CompareEdge32
        // (std::push_heap puts the comparator-maximal element at the root).
        size_t m = (r < n && c(v[l], v[r])) ? r : l;
        if (!c(e, v[m])) break;
        v[k] = v[m];
        k = m;
      }
      v[k] = e;
    }
  };
  std::vector<Heap> heaps(nobs);

  // Per-row sync state (threshold + spinlock) co-located on its OWN cache line.
  // The triangle tiling makes every row r's heap a target for every thread that
  // processes a tile-pair containing r, so the spinlock's acquire RMW and the
  // threshold store are hammered concurrently.  Packing them as separate 1-byte
  // / 4-byte entries in two flat arrays put 64 locks (or 16 thresholds) per
  // cache line, so an update to row r falsely invalidated the line for up to 63
  // neighbouring rows.  alignas(64) gives each row a private line: a row update
  // touches exactly one line and never false-shares with other rows.  True
  // same-row contention (correct and unavoidable) is unchanged.
  // The critical section is a single heap push/pop, so a futex-backed omp_lock
  // is overkill: a relaxed test-and-test-and-set spinlock keeps the hot path in
  // user space.  Most candidates never reach here (rejected by the threshold
  // fast-path), so contention is brief.
  // thresh holds the packed (similarity, neighbour id) key of the weakest kept
  // edge (heap top). Packing both into one atomic lets the lock-free fast path
  // reject a candidate iff it cannot beat the current top under the SAME total
  // order as CompareEdge — including the equal-similarity id tie-break — so the
  // deterministic tie-break adds zero extra spinlock traffic vs the sv-only
  // threshold. (The two-field sv/id alternative would need either an extra atomic
  // read that can race, or letting every equal-sv candidate take the lock.)
  // thresh_sv: similarity of the weakest kept edge (heap top) — a pure float for
  // the cheap common-case reject (sv strictly below threshold), identical cost to
  // the original sv-only test. thresh_key: the FULL packed (similarity, id) key of
  // that same top, consulted only for the rare near-threshold candidates so the
  // deterministic equal-similarity id tie-break is resolved lock-free too. Both
  // are written together under the lock; both are monotonic so the relaxed reads
  // never reject a true winner.
  struct alignas(64) RowSync {
    std::atomic<float>    thresh_sv;
    std::atomic<uint64_t> thresh_key;
    std::atomic<uint8_t>  spin;
  };
  std::vector<RowSync> rowsync(nobs);
  for (size_t r = 0; r < nobs; ++r) {
    rowsync[r].thresh_sv.store(std::numeric_limits<float>::lowest(),
                               std::memory_order_relaxed);
    rowsync[r].thresh_key.store(0, std::memory_order_relaxed);
    rowsync[r].spin.store(0, std::memory_order_relaxed);
  }

  // Optional hot-loop diagnostics.  Counters are strictly thread-local inside
  // the parallel region; enabling RB_GRAPH_PROF therefore avoids introducing
  // the very atomics/locks that this profile is intended to measure.
  struct alignas(64) GraphPassProfile {
    uint64_t tiles = 0;
    uint64_t geometricPairs = 0;
    uint64_t similarityPairs = 0;
    uint64_t upperBoundPruned = 0;
    uint64_t q8Pruned = 0;
    uint64_t q8FloatChecked = 0;
    uint64_t q8FloatRejected = 0;
    uint64_t partialThenFull = 0;
    uint64_t graphSimCalls = 0;
    uint64_t rowUpdateCalls = 0;
    uint64_t similarityRejects = 0;
    uint64_t keyRejects = 0;
    uint64_t lockAcquires = 0;
    uint64_t lockContentions = 0;
    uint64_t heapPushes = 0;
    uint64_t heapReplacements = 0;
    uint64_t eeStop[8] = {0, 0, 0, 0, 0, 0, 0, 0};  // early exits per checkpoint
  };
  const bool graphProfileOn = getenv("RB_GRAPH_PROF") != nullptr;
  std::vector<GraphPassProfile> graphProfile(numThreads);

  // Tile size: keep one tile-pair's *j-block* resident in a core's L2.  The loop
  // structure is `for i in i-tile: screen j-tile; process survivors`, so the data
  // reused across the whole i-loop is the j-tile's winner rows (and, on the
  // abundance-per-pair paths, its depth rows).  The i-side is streamed one row at
  // a time, so only a single block — not two — needs to stay resident.  (The old
  // heuristic charged 2×, which halved the block and left it undersized.)  The
  // production PMH path packs winners to 16 bits, so charging sizeof(uint32_t)
  // here was also too conservative; row heaps are sparse cross-graph writes, not
  // tile-local reuse, so maxEdges is excluded.  We fill ~2/3 of L2 with the
  // j-block and leave the rest for the streamed i-row, the abundance/rank8
  // filter's own j-block, rowsync lines, and set-associativity headroom (winner
  // rows are ~1 KB-strided, which stresses a set-associative L2).  Everything is
  // derived from CacheSize(); no absolute sizes are baked in.  Rounded down to a
  // SIMD-friendly multiple of 16 rows.
  size_t TILE = 10;
  try {
    const size_t winner_bytes = (g_win_bits == 16)
                                    ? sizeof(uint16_t)
                                    : sizeof(uint32_t);
    const size_t winner_row = (size_t)g_pmh_m * winner_bytes;
    const size_t depth_row = num_depth_samples > 1
                                 ? (size_t)num_depth_samples * sizeof(float)
                                 : 0;
    const size_t block_row_bytes = winner_row + depth_row;
    if (block_row_bytes > 0) {
      // The j winner block is *pinned* in L2 for the whole (tall) i-sweep while
      // thousands of transient i-rows and rank8 rows stream past it, so it must
      // fit with generous headroom or capacity/conflict misses evict it.  Fill
      // ~45% of L2 with the pinned block; the rest absorbs the streaming churn.
      const size_t usable_l2 = (size_t)(CacheSize() * 1024. * 0.45);
      TILE = usable_l2 / block_row_bytes;
      TILE = std::max((size_t)10, std::min((size_t)1024, TILE));
      // With hundreds of abundance dimensions the dot-product stream, dynamic
      // work balance, and early heap-threshold formation benefit from more tile
      // pairs than the capacity-only estimate provides.  On this path cap the
      // block at 256 rows; low/medium-dimensional data still use the larger L2
      // block. Lower-dimensional data retain the capacity-derived size.
      if (depth_row > winner_row / 2) TILE = std::min((size_t)256, TILE);
      if (TILE >= 16) TILE = (TILE / 16) * 16;
    }
  } catch (...) {}
  // TILE only affects blocking (never the result), so it is a safe speed knob.
  // The default sizes a tile-pair's j winner block for L2; larger tiles let more
  // i-rows reuse that resident block, but once it no longer fits L2 the reuse
  // spills to DRAM and the pass slows again.  RABBIT_TILE forces an override.
  if (const char *e = getenv("RABBIT_TILE")) {
    long v = atol(e);
    if (v >= 10) TILE = (size_t)v;
  }

  // Asymmetric blocking.  The j-block (winner rows) is what must stay L2-resident
  // across a tile-pair, and it is *reused once per i-row* in the tile.  Survivors
  // of the abundance screen are sparse, so with a square tile a given j winner
  // row is only revisited by a couple of i-rows before eviction — the PMH stage
  // then re-streams winner rows from DRAM and becomes bandwidth-bound.  Keeping
  // the j-block narrow (L2-resident) but making the i-block *tall* multiplies how
  // many i-rows reuse each resident winner row, cutting winner-row DRAM traffic
  // several-fold.  The i-block costs no extra L2 (its rows are streamed one at a
  // time), so its height is bounded only by load-balance granularity, not cache.
  size_t J_TILE = TILE;                 // L2-resident winner-block width
  size_t iMul = 32;                     // i-block is iMul × taller than j-block
  if (const char *e = getenv("RABBIT_ITILE_MUL")) {
    long v = atol(e);
    if (v >= 1) iMul = (size_t)v;
  }
  size_t I_TILE = J_TILE * iMul;
  if (const char *e = getenv("RABBIT_JTILE")) {
    long v = atol(e);
    if (v >= 10) J_TILE = (size_t)v;
  }
  bool i_tile_forced = false;
  if (const char *e = getenv("RABBIT_ITILE")) {
    long v = atol(e);
    if (v >= 10) { I_TILE = (size_t)v; i_tile_forced = true; }
  }
  if (I_TILE < J_TILE) I_TILE = J_TILE;
  // Load-balance floor: the auto height is capped only if it would leave too few
  // tile-pairs to keep the thread pool busy.  Because the j-dimension is finely
  // tiled (many j-tiles), this rarely binds — a tall i-block is fine.  An explicit
  // RABBIT_ITILE bypasses the cap for experimentation.
  if (!i_tile_forced) {
    const size_t nj_tiles = (nobs + J_TILE - 1) / J_TILE;
    const size_t want_tiles = (size_t)numThreads * 8;
    if (nj_tiles > 0) {
      size_t min_i_tiles = (want_tiles + nj_tiles - 1) / nj_tiles;
      if (min_i_tiles < 1) min_i_tiles = 1;
      size_t cap = (nobs + min_i_tiles - 1) / min_i_tiles;
      if (cap < J_TILE) cap = J_TILE;
      if (I_TILE > cap) I_TILE = cap;
    }
  }

  // Flatten the upper-triangle tile-pairs into a list so dynamic scheduling
  // balances load evenly.  Tiles are I_TILE (rows) × J_TILE (cols); a tile is
  // kept only if it can contain a pair with j>i, i.e. its top-right corner lies
  // above the diagonal.  Straddling tiles clamp j_begin per row to max(jj,i+1).
  std::vector<std::pair<uint32_t, uint32_t>> tilepairs;
  {
    const size_t ni_tiles = (nobs + I_TILE - 1) / I_TILE;
    const size_t nj_tiles = (nobs + J_TILE - 1) / J_TILE;
    tilepairs.reserve(ni_tiles * nj_tiles);
    for (size_t ii = 0; ii < nobs; ii += I_TILE)
      for (size_t jj = 0; jj < nobs; jj += J_TILE) {
        const size_t j_stop = std::min(jj + J_TILE, nobs);
        if (j_stop <= ii + 1) continue;  // entirely below diagonal (no j>i)
        tilepairs.emplace_back((uint32_t)ii, (uint32_t)jj);
      }
  }

  verbose_message(
      "Fused graph (triangle): nobs=%zu I_TILE=%zu J_TILE=%zu tilepairs=%zu\n",
      nobs, I_TILE, J_TILE, tilepairs.size());

  const bool rb_timing = (getenv("RB_TIMING") != nullptr);
  std::chrono::steady_clock::time_point _t_pass0;
  if (rb_timing) _t_pass0 = std::chrono::steady_clock::now();

  // ── Abundance-first exact upper-bound prune (RABBIT_ABDFIRST=1) ───────────
  // The default coverage weight is w = min(max(corr,0), Jcov), so a pair
  // with corr < min_edge_weight cannot pass the final cutoff.  Skip its PMH
  // kernel before it can consume a bounded candidate-neighbour slot.
  // This bound also holds for corr-only and product fusion, but not for
  // weighted Jaccard alone, geometric-mean fusion, or dual-channel correlation.
  // Adaptive cutoffs are unknown here and must not use the fixed-cutoff bound.
  // Abundance-first prune is ON by default for multi-sample data (it is exact
  // w.r.t. the fused-edge filter and substantially improves multi-sample
  // binning by stopping high-composition / low-abundance pairs from starving
  // the composition top-k). RABBIT_NO_ABDFIRST=1 disables this pruning.
  const bool corr_bounds_weight = !g_dual_conj &&
      (g_depth_sim == 0 ||
       (g_depth_sim == 2 && (g_depth_fuse == 1 || g_depth_fuse == 2)));
  const bool abdfirst = (getenv("RABBIT_NO_ABDFIRST") == nullptr) &&
                        g_coverage_samples >= 3 && corr_bounds_weight &&
                        rb_env_edge_cut_mode() == 0;
  const double abd_corr_min = abdfirst ? (double)min_edge_weight : -2.0;
  // Precompute per-contig unit rank vectors u_i so the abundance correlation
  // corr(i,j) = Σ_k u_i[k]·u_j[k] is a single length-S dot product — bit-equal
  // to the Pearson-on-ranks that cal_depth_corr() computes (depth_matrix is
  // already rank-transformed), but ~100× cheaper than its per-pair Welford+sqrt.
  // A small safety margin keeps the prune conservative against float rounding:
  // we only skip when the fast estimate is clearly below the exact bound, so no
  // pair that could pass the fused filter is ever dropped.
  const uint32_t ABD_S = (uint32_t)num_depth_samples;
  const float abd_corr_min_eps = (float)(abd_corr_min - 1e-4);
#if defined(__AVX512F__)
  const bool preload_short_abd = abdfirst && ABD_S >= 16 && ABD_S <= 64 &&
      (getenv("RABBIT_NO_ABD_PRELOAD") == nullptr);
#else
  const bool preload_short_abd = false;
#endif
  // The abd-first prune needs centered + L2-normalised rank vectors so that the
  // depth correlation is a plain dot product.  The edge-weight stage already
  // builds exactly this in the global g_depth_unit (identical layout and values;
  // it is built before build_similarity_graph() is called and is not modified
  // afterwards).  Reuse it to avoid a duplicate O(nobs·S) build pass and a second
  // ~nobs·S float matrix; fall back to a local build only if it is unavailable.
  std::vector<float> abd_unit;   // fallback storage only
  const float *abd_u = nullptr;
  bool abd_reused = false;
  bool use_abd_q8 = false;
  size_t abd_q8_stride = 0;
  int32_t abd_q8_prune_max = std::numeric_limits<int32_t>::min();
  std::vector<int8_t> abd_q8;
  std::vector<int32_t> abd_q8_sum;
  // ── Exact integer Spearman filter (rank8, 1×16 transposed VNNI kernel) ──
  // depth_matrix rows are average-tie Spearman ranks over S columns, so
  // 2·(rank − mean) is an exact integer in [−(S−1), S−1] and fits int8 for
  // S ≤ 128.  corr(i,j) = dot(r8_i, r8_j) · invn_i · invn_j EXACTLY (no
  // quantisation); the only rounding is the final float scale, covered by a
  // tiny fallback band around the prune threshold in which the original
  // unit_dot_f decides.  Layout: contigs in blocks of 16; block b, sample group
  // g (4 samples) is one 64-byte vector whose lane l holds the 4 int8 ranks of
  // contig 16b+l — so one vpdpbusd advances 16 pairs by 4 samples and the whole
  // 16-pair dot needs G = ceil(S/4) instructions and NO horizontal reduction.
  // Σ r8 = 0 per row, so biasing the broadcast left operand by +128 (u8) needs
  // no correction term.
  bool use_rank8 = false;
  // Fallback band: |float unit_dot_f − exact corr| is bounded by ~8·2⁻²⁴·S⁰
  // (≤ 6e-7 for S ≤ 128) and the rank8 float scale adds ≤ 3e-7; 4e-6 leaves
  // a wide margin, and the band is decided by unit_dot_f itself.
  constexpr float R8_BAND = 4e-6f;
  size_t r8_G = 0;                    // sample groups of 4
  std::vector<int8_t>  r8_blk;        // nblocks × G × 64 bytes (transposed)
  std::vector<uint32_t> r8_row_u8;    // nobs × G packed (r8+128) bytes, broadcast source
  std::vector<float>   r8_invn;       // nobs: 1/sqrt(Σ r8²), 0 for constant rows
  std::vector<uint16_t> r8_nzmask;    // nblocks: lane l set iff g_anynz[16b+l]
  if (abdfirst) {
    if (!g_depth_unit.empty() && g_depth_unit.size() == (size_t)nobs * ABD_S) {
      abd_u = g_depth_unit.data();
      abd_reused = true;
    } else {
      abd_unit.assign((size_t)nobs * ABD_S, 0.0f);
#pragma omp parallel for num_threads(numThreads) schedule(static)
      for (size_t r = 0; r < nobs; ++r) {
        double mean = 0.0;
        for (uint32_t k = 0; k < ABD_S; ++k) mean += depth_matrix(r, k);
        mean /= ABD_S;
        double ss = 0.0;
        for (uint32_t k = 0; k < ABD_S; ++k) {
          double d = (double)depth_matrix(r, k) - mean; ss += d * d;
        }
        if (ss > 0.0) {
          const double inv = 1.0 / std::sqrt(ss);
          float *u = abd_unit.data() + r * ABD_S;
          for (uint32_t k = 0; k < ABD_S; ++k)
            u[k] = (float)(((double)depth_matrix(r, k) - mean) * inv);
        }
      }
      abd_u = abd_unit.data();
    }
    verbose_message("Abundance-first prune: skip pairs with depth corr < %.12g "
                    "(coverage-only, min_edge=%.12g)%s\n",
                    abd_corr_min, (double)min_edge_weight,
                    abd_reused ? " [reused g_depth_unit]" : "");
#if defined(__AVX512VNNI__) && defined(__AVX512BW__)
    {
      const char *r8Env = getenv("RABBIT_ABD_RANK8");
      const bool r8Requested = r8Env ? atoi(r8Env) != 0 : true;
      const bool anynz_ok = depth_file.empty() || g_anynz.size() == nobs;
      if (r8Requested && ABD_S >= 3 && ABD_S <= 128 && anynz_ok) {
        r8_G = ((size_t)ABD_S + 3u) / 4u;
        const size_t nblk = (nobs + 15u) / 16u;
        r8_blk.assign(nblk * r8_G * 64u, (int8_t)0);
        r8_row_u8.assign((size_t)nobs * r8_G, 0u);
        r8_invn.assign(nobs, 0.0f);
        r8_nzmask.assign(nblk, 0u);
        int bad = 0;
#pragma omp parallel for num_threads(numThreads) schedule(static) \
    reduction(+ : bad)
        for (size_t b = 0; b < nblk; ++b) {
          int8_t *blk = r8_blk.data() + b * r8_G * 64u;
          uint16_t nzm = 0;
          for (size_t l = 0; l < 16; ++l) {
            const size_t r = b * 16u + l;
            if (r >= nobs) break;
            if (depth_file.empty() || g_anynz[r]) nzm |= (uint16_t)(1u << l);
            double mean = 0.0;
            for (uint32_t k = 0; k < ABD_S; ++k) mean += depth_matrix(r, k);
            mean /= (double)ABD_S;
            long long ss = 0, sum = 0;
            uint32_t *rowu = r8_row_u8.data() + r * r8_G;
            for (uint32_t k = 0; k < ABD_S; ++k) {
              const double v = 2.0 * ((double)depth_matrix(r, k) - mean);
              const long q = std::lround(v);
              if (std::fabs(v - (double)q) > 1e-6 || q < -127 || q > 127) ++bad;
              ss += (long long)q * q;
              sum += q;
              blk[(size_t)(k >> 2) * 64u + l * 4u + (k & 3u)] = (int8_t)q;
              rowu[k >> 2] |= (uint32_t)((uint8_t)(q + 128)) << (8u * (k & 3u));
            }
            if (sum != 0) ++bad;
            // Padded samples (k ≥ S) stay 0 in the block; the broadcast side
            // holds +128 there so the pad contributes 128·0 = 0.
            for (uint32_t k = ABD_S; k < r8_G * 4u; ++k)
              rowu[k >> 2] |= 128u << (8u * (k & 3u));
            r8_invn[r] = (ss > 0) ? (float)(1.0 / std::sqrt((double)ss)) : 0.0f;
          }
          r8_nzmask[b] = nzm;
        }
        if (bad == 0) {
          use_rank8 = true;
          verbose_message(
              "Abundance prefilter: exact int8 Spearman rank dot (1x16 VNNI, "
              "S=%u groups=%zu, %.1f MB)\n",
              ABD_S, r8_G,
              (double)(r8_blk.size() + r8_row_u8.size() * 4 +
                       r8_invn.size() * 4) / 1048576.0);
        } else {
          verbose_message("Abundance prefilter: rank8 disabled (%d rows not "
                          "integer-centered ranks); using q8\n", bad);
          std::vector<int8_t>().swap(r8_blk);
          std::vector<uint32_t>().swap(r8_row_u8);
          std::vector<float>().swap(r8_invn);
          std::vector<uint16_t>().swap(r8_nzmask);
        }
      }
    }
    // ISA-gated exact-safe int8 prefilter. Quantisation is never used as the
    // similarity value: an analytically conservative error radius only rejects
    // pairs that cannot reach the original float-dot threshold. All remaining
    // pairs still execute unit_dot_f and therefore keep the original graph
    // semantics. The dimension limit follows directly from int32 accumulator
    // range for dpbusd's worst-case (255 * 127) product, not from input data.
    const char *abdQ8Env = getenv("RABBIT_ABD_Q8");
    constexpr int Q8_SCALE = 127;
    const size_t q8MinDefaultDim = sizeof(__m512i) / sizeof(float);
    const size_t q8SafeDim =
        (size_t)std::numeric_limits<int32_t>::max() /
        (255u * (size_t)Q8_SCALE);
    // Below one 512-bit float vector, padding to a 64-byte int8 row can move
    // more data than the original dot. Keep that small-S case on float unless
    // explicitly forced; RABBIT_ABD_Q8=0 disables the prefilter everywhere.
    const bool q8Requested = abdQ8Env ? atoi(abdQ8Env) != 0
                                      : (size_t)ABD_S >= q8MinDefaultDim;
    if (!use_rank8 && q8Requested && ABD_S > 0 && (size_t)ABD_S <= q8SafeDim) {
      abd_q8_stride = ((size_t)ABD_S + 63u) & ~(size_t)63u;
      if (abd_q8_stride != 0 &&
          nobs <= std::numeric_limits<size_t>::max() / abd_q8_stride) {
        abd_q8.assign(nobs * abd_q8_stride, (int8_t)0);
        abd_q8_sum.assign(nobs, 0);
        double maxXNorm = 0.0, maxErrNorm = 0.0;
#pragma omp parallel for num_threads(numThreads) schedule(static) \
    reduction(max : maxXNorm, maxErrNorm)
        for (size_t r = 0; r < nobs; ++r) {
          const float *src = abd_u + r * (size_t)ABD_S;
          int8_t *dst = abd_q8.data() + r * abd_q8_stride;
          int32_t qsum = 0;
          double xss = 0.0, ess = 0.0;
          for (uint32_t k = 0; k < ABD_S; ++k) {
            long q = std::lround((double)src[k] * Q8_SCALE);
            q = std::max(-Q8_SCALE, std::min(Q8_SCALE, (int)q));
            dst[k] = (int8_t)q;
            qsum += (int32_t)q;
            const double x = (double)q / Q8_SCALE;
            const double e = (double)src[k] - x;
            xss += x * x;
            ess += e * e;
          }
          abd_q8_sum[r] = qsum;
          maxXNorm = std::max(maxXNorm, std::sqrt(xss));
          maxErrNorm = std::max(maxErrNorm, std::sqrt(ess));
        }

        // For u=x+e and v=y+f:
        // u.v <= x.y + ||x||||f|| + ||y||||e|| + ||e||||f||.
        // Global maxima make one pair-independent bound, leaving the hot loop
        // with only an integer dot and compare. Add a forward-error bound for
        // the original float FMA/reduction so this prefilter cannot reject a
        // pair that unit_dot_f would round upward across the threshold.
        const double eps = std::numeric_limits<float>::epsilon();
        const double dotOps = (double)(((size_t)ABD_S + 15u) / 16u + 4u);
        const double gamma =
            (dotOps * eps < 0.5) ? (dotOps * eps) / (1.0 - dotOps * eps)
                                 : 1.0;
        const double maxUnitNorm = maxXNorm + maxErrNorm;
        const double slack = 2.0 * maxXNorm * maxErrNorm +
                             maxErrNorm * maxErrNorm +
                             gamma * maxUnitNorm * maxUnitNorm + 1e-7;
        const double qCutReal =
            ((double)abd_corr_min_eps - slack) * Q8_SCALE * Q8_SCALE;
        const double qCutCeil = std::ceil(qCutReal);
        if (qCutCeil > (double)std::numeric_limits<int32_t>::min() &&
            qCutCeil <= (double)std::numeric_limits<int32_t>::max()) {
          abd_q8_prune_max = (int32_t)qCutCeil - 1;
          use_abd_q8 = true;
          verbose_message(
              "Abundance prefilter: exact-safe int8/VNNI S=%u stride=%zu "
              "error_bound=%.6g q_prune_max=%d\n",
              ABD_S, abd_q8_stride, slack, abd_q8_prune_max);
        } else {
          std::vector<int8_t>().swap(abd_q8);
          std::vector<int32_t>().swap(abd_q8_sum);
        }
      }
    }
#endif
    if (!use_abd_q8 && preload_short_abd)
      verbose_message("Abundance dot: batched AVX-512 left-row preload (S=%u)\n",
                      ABD_S);
  }

  rabbit_abundance::DotFilter abd_batch;
  const bool use_abd_batch = abdfirst && !use_rank8 && !use_abd_q8 &&
      getenv("RABBIT_NO_ABD_BATCH") == nullptr &&
      abd_batch.build(abd_u, nobs, ABD_S, abd_corr_min_eps, (int)numThreads);
  if (use_abd_batch)
    verbose_message("Abundance prefilter: conservative batched dot (S=%u, "
                    "%.1f MB); survivors use original dot\n", ABD_S,
                    (double)abd_batch.bytes() / 1048576.0);

  // ── Composition early-exit (exact) ───────────────────────────────────────
  // For the default PMH winners metric a pair can only become an edge if its
  // similarity beats at least one endpoint's current heap threshold.  After
  // matching the first EE_P winners we have an upper bound on the full match
  // fraction (= all remaining winners also match), hence — after the monotone
  // baseline correction — an upper bound on the similarity sv.  If that upper
  // bound is below BOTH endpoints' (monotonically rising) thresholds the pair
  // can never be kept, so the remaining (g_pmh_m − EE_P) winner comparisons are
  // skipped. This is exact with respect to the retained top-k heaps.
  // RABBIT_NO_COMP_EE=1 disables this optimization.
  const bool win16 = (g_win_bits == 16);
  const bool comp_ee = g_pmh_mode && !g_exact_cos_cmp && g_pmh_m >= 16 &&
                       (win16 ? !g_win16.empty() : !g_win_flat.empty()) &&
                       (getenv("RABBIT_NO_COMP_EE") == nullptr);
  // A PMH score depends only on an integer match count in [0,m]. Compute the
  // SAME correction/clamp once per count, rather than for every checkpoint
  // of every pair. This table also preserves float ties in the heap bound.
  std::vector<StoredDistance> pmh_count_scores;
  if (comp_ee) {
    pmh_count_scores.resize((size_t)g_pmh_m + 1);
    for (size_t c = 0; c < pmh_count_scores.size(); ++c)
      pmh_count_scores[c] = (StoredDistance)pmh_sv_from_count((uint32_t)c);
  }
  const uint32_t EE_P = (g_pmh_m >= 64) ? 64u : g_pmh_m;
  // Progressive re-check granularity (registers) after the EE_P prefix.
  // RABBIT_COMP_EE_STEP overrides; 0 or >= m restores the single check.
  uint32_t EE_STEP = 128u;
  if (const char *e = getenv("RABBIT_COMP_EE_STEP")) {
    long v = atol(e);
    EE_STEP = (v <= 0 || v >= (long)g_pmh_m) ? g_pmh_m : (uint32_t)v;
  }

  // Spread the large read-mostly winner array over both memory controllers so
  // the all-pairs scan is not bottlenecked on the first-touch socket's DRAM.
  if (win16) {
    if (!g_win16.empty())
      numa_interleave_buffer(g_win16.data(), g_win16.size() * sizeof(uint16_t));
  } else if (!g_win_flat.empty()) {
    numa_interleave_buffer(g_win_flat.data(),
                           g_win_flat.size() * sizeof(uint32_t));
  }

  // ── Optional LSH candidate generation ───────────────────────────────────
  // Build B band buckets: for each band, contig ids sorted by the band's
  // R-winner hash so equal-hash contigs (candidates) form contiguous ranges.
  const bool     use_lsh = rb_lsh_on() && g_pmh_mode && !g_win_flat.empty();
  const uint32_t LSH_R   = rb_lsh_R();
  const uint32_t LSH_B   = use_lsh ? (g_pmh_m / LSH_R) : 0;
  const size_t   LSH_MAXB= rb_lsh_maxbucket();
  std::vector<std::vector<uint64_t>> band_sorted_key;  // [B][nobs] ascending
  std::vector<std::vector<uint32_t>> band_order;       // [B][nobs] contig ids by key
  std::vector<std::vector<uint64_t>> band_key_of;      // [B][nobs] key per contig
  if (use_lsh) {
    band_sorted_key.assign(LSH_B, std::vector<uint64_t>(nobs));
    band_order.assign(LSH_B, std::vector<uint32_t>(nobs));
    band_key_of.assign(LSH_B, std::vector<uint64_t>(nobs));
    const uint32_t *winbase = g_win_flat.data();
#pragma omp parallel for num_threads(numThreads) schedule(dynamic, 1)
    for (uint32_t b = 0; b < LSH_B; ++b) {
      auto &ord = band_order[b]; auto &skey = band_sorted_key[b];
      auto &kof = band_key_of[b];
      const uint32_t off = b * LSH_R;
      for (size_t c = 0; c < nobs; ++c) {
        kof[c] = lsh_band_hash(winbase + c * g_pmh_m + off, LSH_R);
        ord[c] = (uint32_t)c;
      }
      std::sort(ord.begin(), ord.end(),
                [&](uint32_t a, uint32_t c) { return kof[a] < kof[c]; });
      for (size_t p = 0; p < nobs; ++p) skey[p] = kof[ord[p]];
    }
    verbose_message("LSH candidate gen: B=%u bands × R=%u rows, maxbucket=%zu\n",
                    LSH_B, LSH_R, LSH_MAXB);
  }

#pragma omp parallel num_threads(numThreads)
  {
    const int tid = omp_get_thread_num();
    GraphPassProfile &my_prof = graphProfile[tid];
#if defined(__AVX512VNNI__) && defined(__AVX512BW__)
    // rank8 per-row survivor scratch (bounded by the tile width).
    std::vector<uint32_t> r8_pass(use_rank8 ? J_TILE + 16 : 0);
    std::vector<uint32_t> r8_band(use_rank8 ? J_TILE + 16 : 0);
    std::vector<uint32_t> r8_merge;
    const Rank8RowFn r8_fn = use_rank8 ? rank8_row_filter_for(r8_G) : nullptr;
#endif
    std::vector<uint32_t> visited;       // LSH per-i dedup stamp
    uint32_t visit_gen = 0;
    if (use_lsh) visited.assign(nobs, 0u);

    // Heap keys and atomic thresholds use the same total order. The weakest
    // kept key only increases, so a stale relaxed load cannot reject a winner.
    static const CompareEdge32 cmp_edge{};
    auto store_thresh = [&](RowSync &rs, const Edge32 &top) {
      rs.thresh_sv.store(top.sv(), std::memory_order_relaxed);
      rs.thresh_key.store(top.key, std::memory_order_relaxed);
    };
    // t_seen: a threshold value for row r read earlier by the caller.  Row
    // thresholds only ever rise, so rejecting on a stale (lower) snapshot is
    // exact and saves re-touching the contended rowsync line in the common
    // case; the snapshot is refreshed only when the fast test does not reject.
    auto update_row = [&](size_t r, size_t other, StoredDistance sv,
                          float t_seen) {
      if (graphProfileOn) ++my_prof.rowUpdateCalls;
      RowSync &rs = rowsync[r];
      // Common case: similarity strictly below the kept threshold — rejected with
      // one float compare, identical cost to the original sv-only fast path.
      if (sv < t_seen || sv < rs.thresh_sv.load(std::memory_order_relaxed)) {
        if (graphProfileOn) ++my_prof.similarityRejects;
        return;
      }
      // Near/at threshold (rare): resolve the deterministic equal-similarity
      // smaller-id-wins tie-break lock-free via the packed key.
      const uint64_t ck = Edge32::pack(sv, (uint32_t)other);
      if (ck <= rs.thresh_key.load(std::memory_order_relaxed)) {
        if (graphProfileOn) ++my_prof.keyRejects;
        return;
      }
      // ── acquire spinlock (test-and-test-and-set) ──
      if (rs.spin.exchange(1, std::memory_order_acquire)) {
        if (graphProfileOn) ++my_prof.lockContentions;
        do {
        while (rs.spin.load(std::memory_order_relaxed))
#if defined(__x86_64__) || defined(__i386__)
          __builtin_ia32_pause();
#else
          ;
#endif
        } while (rs.spin.exchange(1, std::memory_order_acquire));
      }
      if (graphProfileOn) ++my_prof.lockAcquires;
      Heap &h = heaps[r];
      const Edge32 cand{ck};
      if (h.size() < (size_t)maxEdges) {
        h.push(cand, (size_t)maxEdges);
        if (graphProfileOn) ++my_prof.heapPushes;
        if (h.size() == (size_t)maxEdges) store_thresh(rs, h.top());
      } else if (cmp_edge(cand, h.top())) {  // cand kept-preferred over weakest
        h.replace_top(cand);
        if (graphProfileOn) ++my_prof.heapReplacements;
        store_thresh(rs, h.top());
      }
      rs.spin.store(0, std::memory_order_release);  // release spinlock
    };

    auto process_pair = [&](size_t i, size_t j) {
      if (graphProfileOn) ++my_prof.similarityPairs;
      StoredDistance sv;
      float t_i = std::numeric_limits<float>::lowest();
      float t_j = std::numeric_limits<float>::lowest();
      if (comp_ee) {
        // Partial winner match → exact upper bound on sv; prune if it cannot
        // beat either endpoint's heap threshold.  Result is bit-identical to
        // graph_sim(i,j) on the non-pruned path (same total match count).
        uint32_t cP, cFull;
        if (win16) {
          // Progressive exact early exit: after p registers with c matches the
          // full count is ≤ c + (m − p), so the baseline-corrected similarity
          // is bounded above; once that bound is below BOTH endpoints' current
          // heap thresholds (monotonically rising, read once) the pair can
          // never be kept and the remaining registers are skipped.  Checked
          // after EE_P, then every EE_STEP registers.
          const uint16_t *aw = g_win16.data() + i * (size_t)g_pmh_m;
          const uint16_t *bw = g_win16.data() + j * (size_t)g_pmh_m;
          const float ti = rowsync[i].thresh_sv.load(std::memory_order_relaxed);
          const float tj = rowsync[j].thresh_sv.load(std::memory_order_relaxed);
          t_i = ti; t_j = tj;
          const float tmin = ti < tj ? ti : tj;
          auto bound_below = [&](uint32_t c, uint32_t p) -> bool {
            return pmh_count_scores[c + (g_pmh_m - p)] < tmin;
          };
          uint32_t c = pmh_match_count16(aw, bw, EE_P);
          uint32_t p = EE_P;
          if (bound_below(c, p)) {
            if (graphProfileOn) { ++my_prof.upperBoundPruned; ++my_prof.eeStop[0]; }
            return;  // exact: unkeepable
          }
          uint32_t stage = 1;
          while (p < g_pmh_m) {
            const uint32_t step = std::min<uint32_t>(EE_STEP, g_pmh_m - p);
            c += pmh_match_count16(aw + p, bw + p, step);
            p += step;
            if (p < g_pmh_m && bound_below(c, p)) {
              if (graphProfileOn) {
                ++my_prof.upperBoundPruned;
                ++my_prof.eeStop[stage < 7 ? stage : 7];
              }
              return;  // exact: unkeepable
            }
            ++stage;
          }
          cP = c;
          cFull = c;
          (void)cP;
        } else {
          const uint32_t *aw = g_win_flat.data() + i * (size_t)g_pmh_m;
          const uint32_t *bw = g_win_flat.data() + j * (size_t)g_pmh_m;
          cP = pmh_match_count(aw, bw, EE_P);
          const StoredDistance up = pmh_count_scores[cP + (g_pmh_m - EE_P)];
          const float ti = rowsync[i].thresh_sv.load(std::memory_order_relaxed);
          const float tj = rowsync[j].thresh_sv.load(std::memory_order_relaxed);
          t_i = ti; t_j = tj;
          if (up < (ti < tj ? ti : tj)) {
            if (graphProfileOn) ++my_prof.upperBoundPruned;
            return;  // exact: unkeepable
          }
          cFull = cP + pmh_match_count(aw + EE_P, bw + EE_P, g_pmh_m - EE_P);
        }
        if (graphProfileOn) ++my_prof.partialThenFull;
        sv = pmh_count_scores[cFull];
      } else {
        if (graphProfileOn) ++my_prof.graphSimCalls;
        sv = (StoredDistance)graph_sim(i, j);
      }
      update_row(i, j, sv, t_i);
      update_row(j, i, sv, t_j);
    };

    if (!use_lsh) {
#pragma omp for schedule(dynamic, 1)
    for (size_t p = 0; p < tilepairs.size(); ++p) {
      const size_t ii = tilepairs[p].first;
      const size_t jj = tilepairs[p].second;
      const size_t i_stop = std::min(ii + I_TILE, nobs);
      const size_t j_stop = std::min(jj + J_TILE, nobs);
      if (graphProfileOn) ++my_prof.tiles;
      for (size_t i = ii; i < i_stop; ++i) {
        // Only j>i: clamp the column start to max(jj, i+1).  For tiles wholly
        // above the diagonal this is just jj; for straddling tiles the early
        // rows may have no work.
        const size_t j_begin = jj > i + 1 ? jj : i + 1;
        if (j_begin >= j_stop) continue;
        if (graphProfileOn)
          my_prof.geometricPairs += (uint64_t)(j_stop - j_begin);
#if defined(__AVX512VNNI__) && defined(__AVX512BW__)
        if (use_rank8) {
          // Exact int8 Spearman screen for row i over [j_begin, j_stop):
          // survivors (clearly >= threshold) go to PMH in ascending j order;
          // the +-R8_BAND sliver is decided by unit_dot_f (bit-identical).
          Rank8Row R;
          R.blk_base = r8_blk.data();
          R.rowu = r8_row_u8.data() + i * r8_G;
          R.invn = r8_invn.data();
          R.nzmask = r8_nzmask.data();
          R.invn_i = r8_invn[i];
          R.lo = abd_corr_min_eps - R8_BAND;
          R.hi = abd_corr_min_eps + R8_BAND;
          R.nz_i = depth_file.empty() || g_anynz[i];
          size_t npass = 0, nband = 0;
          const uint64_t pruned = r8_fn(R, j_begin, j_stop, r8_pass.data(), npass,
                                        r8_band.data(), nband);
          if (graphProfileOn) {
            my_prof.q8Pruned += pruned;
            my_prof.q8FloatChecked += nband;
          }
          if (nband) {
            // Merge band survivors (rare) into the pass list, keeping j order.
            size_t kept = 0;
            for (size_t q = 0; q < nband; ++q) {
              const size_t j = r8_band[q];
              const float c = unit_dot_f(abd_u + i * ABD_S, abd_u + j * ABD_S, ABD_S);
              if (c < abd_corr_min_eps) {
                if (graphProfileOn) ++my_prof.q8FloatRejected;
              } else {
                r8_band[kept++] = (uint32_t)j;
              }
            }
            if (kept) {
              std::vector<uint32_t> &m = r8_merge;
              m.resize(npass + kept);
              std::merge(r8_pass.begin(), r8_pass.begin() + npass,
                         r8_band.begin(), r8_band.begin() + kept, m.begin());
              std::copy(m.begin(), m.end(), r8_pass.begin());
              npass += kept;
            }
          }
          if (npass) {
            // Warm the lines the PMH stage will miss on: each survivor's
            // contended rowsync line and the head of its winner row.
            for (size_t q = 0; q < npass; ++q) {
              const size_t j = r8_pass[q];
              _mm_prefetch((const char *)&rowsync[j], _MM_HINT_T0);
              if (win16) {
                const char *wr = (const char *)(g_win16.data() + j * (size_t)g_pmh_m);
                _mm_prefetch(wr, _MM_HINT_T0);
                _mm_prefetch(wr + 64, _MM_HINT_T0);
              }
            }
            for (size_t q = 0; q < npass; ++q) process_pair(i, r8_pass[q]);
          }
          continue;
        }
        if (use_abd_q8) {
          const int8_t *qi = abd_q8.data() + i * abd_q8_stride;
          for (size_t j = j_begin; j < j_stop; ++j) {
            if (!is_nz(i, j)) continue;
            const int32_t qdot = unit_dot_q8_vnni(
                qi, abd_q8.data() + j * abd_q8_stride,
                abd_q8_stride, abd_q8_sum[j]);
            if (qdot <= abd_q8_prune_max) {
              if (graphProfileOn) ++my_prof.q8Pruned;
              continue;
            }
            if (graphProfileOn) ++my_prof.q8FloatChecked;
            const float c = unit_dot_f(
                abd_u + i * ABD_S, abd_u + j * ABD_S, ABD_S);
            if (c < abd_corr_min_eps) {
              if (graphProfileOn) ++my_prof.q8FloatRejected;
              continue;
            }
            process_pair(i, j);
          }
          continue;
        }
#endif
        if (use_abd_batch) {
          constexpr size_t width = rabbit_abundance::DotFilter::width;
          const float *left = abd_u + i * ABD_S;
          const size_t block_end = j_stop / width + (j_stop % width != 0);
          for (size_t b = j_begin / width; b < block_end; ++b) {
            const size_t j0 = b * width;
            uint32_t valid = ~uint32_t(0);
            if (j0 < j_begin) valid &= ~uint32_t(0) << (j_begin - j0);
            if (j_stop - j0 < width) valid &= (uint32_t(1) << (j_stop - j0)) - 1;
            uint32_t pass = abd_batch.candidates(left, b) & valid;
            while (pass) {
              const size_t j = j0 + (unsigned)__builtin_ctz(pass);
              pass &= pass - 1;
              if (!is_nz(i, j)) continue;
              // Preserve the exact previous prune decision, including the
              // floating-point boundary and the order of surviving j values.
              if (unit_dot_f(left, abd_u + j * ABD_S, ABD_S) < abd_corr_min_eps)
                continue;
              process_pair(i, j);
            }
          }
          continue;
        }
#if defined(__AVX512F__)
        if (preload_short_abd) {
          UnitDotLeft64 left;
          unit_dot_left64_init(left, abd_u + i * ABD_S, ABD_S);
          static constexpr size_t DOT_BATCH = 8;
          for (size_t jb = j_begin; jb < j_stop; jb += DOT_BATCH) {
            const size_t nb = std::min(DOT_BATCH, j_stop - jb);
            uint8_t nz[DOT_BATCH];
            float corr[DOT_BATCH];
            for (size_t q = 0; q < nb; ++q) {
              const size_t j = jb + q;
              nz[q] = (uint8_t)is_nz(i, j);
              if (nz[q])
                corr[q] = unit_dot_f_left64(left, abd_u + j * ABD_S);
            }
            // Keep PMH/heap processing in the original j order.  Only the
            // independent abundance tests are computed ahead within a batch.
            for (size_t q = 0; q < nb; ++q) {
              if (!nz[q] || corr[q] < abd_corr_min_eps) continue;
              process_pair(i, jb + q);
            }
          }
          continue;
        }
#endif
        for (size_t j = j_begin; j < j_stop; ++j) {
          if (!is_nz(i, j)) continue;
          // Abundance-first exact upper-bound prune: skip the PMH kernel for
          // pairs whose best-case fused score is already below the cutoff.
          if (abdfirst) {
            const float c = unit_dot_f(abd_u + i * ABD_S, abd_u + j * ABD_S, ABD_S);
            if (c < abd_corr_min_eps) continue;
          }
          process_pair(i, j);
        }
      }
      if (verbose && tid == 0 && (p & 0x3FFu) == 0)
        verbose_message("Fusion D %zu/%zu\r", p, tilepairs.size());
    }
    } else {
      // ── LSH path: evaluate only co-banded candidate pairs ────────────────
#pragma omp for schedule(dynamic, 16)
      for (size_t i = 0; i < nobs; ++i) {
        ++visit_gen;
        visited[i] = visit_gen;                 // exclude self
        for (uint32_t b = 0; b < LSH_B; ++b) {
          const uint64_t key   = band_key_of[b][i];
          const auto    &skey  = band_sorted_key[b];
          const auto    &ord   = band_order[b];
          const size_t lo = (size_t)(std::lower_bound(skey.begin(), skey.end(), key) - skey.begin());
          const size_t hi = (size_t)(std::upper_bound(skey.begin(), skey.end(), key) - skey.begin());
          if (hi - lo > LSH_MAXB) continue;     // skip non-discriminative (repeat) buckets
          for (size_t p = lo; p < hi; ++p) {
            const uint32_t j = ord[p];
            if (visited[j] == visit_gen) continue;  // dedup across bands
            visited[j] = visit_gen;
            if ((size_t)j <= i) continue;           // each unordered pair once (i<j)
            if (!is_nz(i, j)) continue;
            StoredDistance sv = (StoredDistance)graph_sim(i, j);
            update_row(i, j, sv, std::numeric_limits<float>::lowest());
            update_row(j, i, sv, std::numeric_limits<float>::lowest());
          }
        }
        if (verbose && tid == 0 && (i & 0xFFFu) == 0)
          verbose_message("LSH cand %zu/%zu\r", i, nobs);
      }
    }
  }

  if (graphProfileOn) {
    GraphPassProfile total{};
    uint64_t minPairs = std::numeric_limits<uint64_t>::max();
    uint64_t maxPairs = 0;
    for (const auto &p : graphProfile) {
      total.tiles += p.tiles;
      total.geometricPairs += p.geometricPairs;
      total.similarityPairs += p.similarityPairs;
      total.upperBoundPruned += p.upperBoundPruned;
      total.q8Pruned += p.q8Pruned;
      total.q8FloatChecked += p.q8FloatChecked;
      total.q8FloatRejected += p.q8FloatRejected;
      total.partialThenFull += p.partialThenFull;
      total.graphSimCalls += p.graphSimCalls;
      total.rowUpdateCalls += p.rowUpdateCalls;
      total.similarityRejects += p.similarityRejects;
      total.keyRejects += p.keyRejects;
      total.lockAcquires += p.lockAcquires;
      total.lockContentions += p.lockContentions;
      total.heapPushes += p.heapPushes;
      total.heapReplacements += p.heapReplacements;
      for (int s = 0; s < 8; ++s) total.eeStop[s] += p.eeStop[s];
      if (p.geometricPairs != 0) {
        minPairs = std::min(minPairs, p.geometricPairs);
        maxPairs = std::max(maxPairs, p.geometricPairs);
      }
    }
    if (minPairs == std::numeric_limits<uint64_t>::max()) minPairs = 0;
    uint64_t nzPairs = total.geometricPairs;
    if (g_anynz.size() == nobs) {
      const uint64_t nzRows = (uint64_t)std::count(
          g_anynz.begin(), g_anynz.end(), (uint8_t)1);
      const uint64_t zeroRows = (uint64_t)nobs - nzRows;
      nzPairs -= zeroRows * (zeroRows - 1) / 2;
    }
    fprintf(stderr,
            "[RB_GRAPH_PROF] tiles=%llu geometric_pairs=%llu nz_pairs=%llu "
            "similarity_pairs=%llu upper_pruned=%llu partial_full=%llu "
            "graph_sim=%llu q8_pruned=%llu q8_float_checked=%llu "
            "q8_float_rejected=%llu thread_pair_min=%llu thread_pair_max=%llu\n",
            (unsigned long long)total.tiles,
            (unsigned long long)total.geometricPairs,
            (unsigned long long)nzPairs,
            (unsigned long long)total.similarityPairs,
            (unsigned long long)total.upperBoundPruned,
            (unsigned long long)total.partialThenFull,
            (unsigned long long)total.graphSimCalls,
            (unsigned long long)total.q8Pruned,
            (unsigned long long)total.q8FloatChecked,
            (unsigned long long)total.q8FloatRejected,
            (unsigned long long)minPairs,
            (unsigned long long)maxPairs);
    fprintf(stderr,
            "[RB_GRAPH_PROF] row_updates=%llu sim_rejects=%llu key_rejects=%llu "
            "lock_acquires=%llu lock_contentions=%llu heap_pushes=%llu "
            "heap_replacements=%llu\n",
            (unsigned long long)total.rowUpdateCalls,
            (unsigned long long)total.similarityRejects,
            (unsigned long long)total.keyRejects,
            (unsigned long long)total.lockAcquires,
            (unsigned long long)total.lockContentions,
            (unsigned long long)total.heapPushes,
            (unsigned long long)total.heapReplacements);
    fprintf(stderr,
            "[RB_GRAPH_PROF] pmh_early_exit@checkpoint: p=%u:%llu",
            EE_P, (unsigned long long)total.eeStop[0]);
    for (int s = 1; s < 8; ++s)
      if (total.eeStop[s])
        fprintf(stderr, " p=%u:%llu",
                (unsigned)std::min<uint64_t>(EE_P + (uint64_t)s * EE_STEP, g_pmh_m),
                (unsigned long long)total.eeStop[s]);
    fprintf(stderr, " full=%llu\n", (unsigned long long)total.partialThenFull);
  }

  if (rb_timing) {
    double ms = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - _t_pass0).count();
    fprintf(stderr, "[RB_TIMING] fusedD pair-pass: %.1f ms\n", ms);
  }

  { std::vector<RowSync>().swap(rowsync); }

  // Emit positive-similarity top-k candidates; no absolute PMH cutoff.
  verbose_message("Starting Building Similarity Graph (Fusion D). "
                  "nobs=%zu maxEdges=%zu\n", nobs, maxEdges);

  std::vector<GraphNodeId>    &from  = g.from;
  std::vector<GraphNodeId>    &to    = g.to;
  std::vector<StoredDistance> &sComp  = g.sComp;

  if (!g_mutual_knn) {
    // Union k-NN (original): emit edge (i,j) if j ∈ top-k(i) OR i ∈ top-k(j)
    for (size_t i = 0; i < nobs; ++i) {
      while (!heaps[i].empty()) {
        auto e = heaps[i].top(); heaps[i].pop();
        size_t j   = e.id();
        StoredDistance sv = e.sv();
        if (sv <= 0.0f) continue;
        if (i < j) {
          from.push_back(i); to.push_back(j); sComp.push_back(sv);
        }
        // i > j edges emitted by the j iteration
      }
    }
  } else {
    // Mutual k-NN (RABBIT_MUTUAL_KNN=1): keep edge (i,j) only if both
    // j ∈ top-k(i) AND i ∈ top-k(j).
    //
    // Each heap already owns a contiguous array of all its neighbours. Filter
    // and sort that array in place: popping every edge before sorting it again
    // performed an unnecessary O(k log k) heap traversal plus a second buffer
    // allocation and copy for every row. No heap operation is needed after the
    // pair pass. Membership remains a binary search over ascending ids.
#pragma omp parallel for num_threads(numThreads) schedule(dynamic, 64)
    for (size_t i = 0; i < nobs; ++i) {
      std::vector<Edge32> &row = heaps[i].v;
      row.erase(std::remove_if(row.begin(), row.end(),
                               [](const Edge32 &e) { return !(e.sv() > 0.0f); }),
                row.end());
      std::sort(row.begin(), row.end(),
                [](const Edge32 &a, const Edge32 &b) { return a.id() < b.id(); });
    }

    // Parallel emit: keep (i,j) with i<j iff row j also contains i.
    std::vector<std::vector<GraphNodeId>>    tl_from(numThreads);
    std::vector<std::vector<GraphNodeId>>    tl_to(numThreads);
    std::vector<std::vector<StoredDistance>> tl_sv(numThreads);
#pragma omp parallel num_threads(numThreads)
    {
      const int tid = omp_get_thread_num();
      auto &lf = tl_from[tid]; auto &lt = tl_to[tid]; auto &ls = tl_sv[tid];
#pragma omp for schedule(dynamic, 64)
      for (size_t i = 0; i < nobs; ++i) {
        const std::vector<Edge32> &row = heaps[i].v;
        for (const Edge32 &nj : row) {
          const uint32_t j = nj.id();
          if ((size_t)j <= i) continue;          // emit each undirected pair once
          const std::vector<Edge32> &rj = heaps[j].v; // is i a neighbor of j?
          auto it = std::lower_bound(
              rj.begin(), rj.end(), (uint32_t)i,
              [](const Edge32 &p, uint32_t v) { return p.id() < v; });
          if (it != rj.end() && it->id() == (uint32_t)i) {
            lf.push_back(i); lt.push_back((size_t)j); ls.push_back(nj.sv());
          }
        }
      }
    }

    // (3) Concatenate per-thread edge lists (deterministic thread order).
    size_t total = 0;
    for (size_t t = 0; t < numThreads; ++t) total += tl_from[t].size();
    from.reserve(total); to.reserve(total); sComp.reserve(total);
    for (size_t t = 0; t < numThreads; ++t) {
      from.insert(from.end(), tl_from[t].begin(), tl_from[t].end());
      to.insert(to.end(),     tl_to[t].begin(),   tl_to[t].end());
      sComp.insert(sComp.end(), tl_sv[t].begin(),   tl_sv[t].end());
    }
  }

  verbose_message("Finished Building Similarity Graph (%zu edges) "
                  "[%.1fGb / %.1fGb]                         \n",
                  g.getEdgeCount(), getUsedPhysMem(),
                  getTotalPhysMem() / 1024 / 1024);
  g.sComp.shrink_to_fit(); g.to.shrink_to_fit(); g.from.shrink_to_fit();
}
