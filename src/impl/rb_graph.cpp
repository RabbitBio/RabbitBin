// RabbitBin module: rb_graph.cpp

// ═══════════════════════════════════════════════════════════════════════════
// build_similarity_graph  –  build edge list using KmerSketch Jaccard
// ═══════════════════════════════════════════════════════════════════════════
// All-pairs reference implementation (O(N^2) Jaccard), kept for A/B validation.
// Reference all-pairs graph build (O(N^2) Jaccard); used when inverted index is off.
static void build_graph_allpairs(Graph &g, Similarity cutoff) {
  ProgressTracker progress(nobs);
  std::vector<size_t> &from = g.from;
  std::vector<size_t> &to   = g.to;
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
    reduction(merge_size_t : from) reduction(merge_size_t : to)                \
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
          if (sv > cutoff &&
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

// ── PMH-winner inverted-index graph build ─────────────────────────────────
// Replaces the all-pairs PMH winner-match scan.  For each contig i:
//   1. Walk its m winners → look up posting lists → accumulate hit counts.
//   2. Hit count for candidate j  ==  PMH match numerator (= #shared winners).
//   3. raw_sim  = count / m;  corrected sim = (raw - b0)/(1 - b0)
//   4. Only candidates with corrected_sim > cutoff are inserted as edges.
// Complexity: O(N × avg_posting_size × avg_candidates), typically sub-quadratic.
// Average posting size for k=6 ≈ nobs/4096 ≈ 7.5 → very short lists.
static void gen_pmh_graph_index(Graph &g, Similarity cutoff) {
  if (!g_pmh_idx) {
    verbose_message("WARN: PMH index not ready, falling back to all-pairs\n");
    build_graph_allpairs(g, cutoff); return;
  }
  const rabbit_invidx::InvertedIndex &idx = *g_pmh_idx;
  const uint32_t *csrPtr = idx.csrPosts.data();
  const double    Md     = (double)g_pmh_m;

  // Minimum raw match count guaranteed to yield corrected_sim > cutoff.
  // corrected = (raw/m - b0)/(1-b0) > cutoff  ↔  raw > m*(cutoff*(1-b0)+b0)
  // Use a small safety margin (0.95×) to avoid false-negative edge loss
  // from estimator variance around the threshold.
  double b0 = (g_pmh_base_on && g_pmh_baseline > 0.0 && g_pmh_baseline < 1.0)
              ? g_pmh_baseline : 0.0;
  size_t minCount = (size_t)std::floor(Md * (cutoff * (1.0 - b0) + b0) * 0.95);
  if (minCount < 1) minCount = 1;

  verbose_message("Starting Building Similarity Graph (PMH winner index). "
                  "nobs=%zu cutoff=%.4f b0=%.4f minCount=%zu "
                  "postings=%zu keys=%zu\n",
                  nobs, (double)cutoff, b0, minCount,
                  idx.totalPostings, idx.postIdx.size());

  std::vector<size_t>         &from  = g.from;
  std::vector<size_t>         &to    = g.to;
  std::vector<StoredDistance> &sComp  = g.sComp;

  std::vector<std::vector<size_t>>         tl_from(numThreads);
  std::vector<std::vector<size_t>>         tl_to(numThreads);
  std::vector<std::vector<StoredDistance>> tl_sComp(numThreads);

  ProgressTracker progress(nobs);

#pragma omp parallel num_threads(numThreads)
  {
    const int tid = omp_get_thread_num();
    // Per-thread epoch-stamped scratch (stamp-based clear → O(candidates), not O(N)).
    std::vector<int>      stamp(nobs, 0);
    std::vector<uint32_t> isect(nobs, 0);  // collision count; bounded by m
    int ep = 0;
    std::vector<size_t> cand;
    cand.reserve(8192);
    std::priority_queue<Edge, std::vector<Edge>, CompareEdge> heap;

    auto &lf = tl_from[tid], &lt = tl_to[tid];
    std::vector<StoredDistance> &ls = tl_sComp[tid];

#pragma omp for schedule(dynamic, 8)
    for (size_t i = 0; i < nobs; ++i) {
      cand.clear();
      if (++ep == INT_MAX) { std::fill(stamp.begin(), stamp.end(), 0); ep = 1; }

      // Walk every non-zero winner register of contig i.
      // Look up winner-index keys using the full 64-bit g_win64_flat.
      // Key = winners64[pos] XOR (pos * GOLDEN) — must match getWinnerIndexKeys().
      static constexpr uint64_t GOLDEN = UINT64_C(0x9E3779B97F4A7C15);
      const uint64_t *wi64 = g_win64_flat.data() + i * g_pmh_m;
      for (uint32_t pos = 0; pos < g_pmh_m; ++pos) {
        if (wi64[pos] == 0ULL) continue;
        const uint64_t key = wi64[pos] ^ ((uint64_t)pos * GOLDEN);
        auto it = idx.postIdx.find(key);
        if (it == idx.postIdx.end()) continue;
        const uint32_t *pl = csrPtr + it->second.off;
        const uint32_t  sz = it->second.cnt;
        for (uint32_t pi = 0; pi < sz; ++pi) {
          size_t j = pl[pi];
          if (j == i) continue;
          if (stamp[j] != ep) { stamp[j] = ep; isect[j] = 1; cand.push_back(j); }
          else                 { ++isect[j]; }
        }
      }

      // Score candidates; only emit i<j edges above cutoff.
      while (!heap.empty()) heap.pop();
      for (size_t j : cand) {
        if (isect[j] < (uint32_t)minCount) continue;
        if (!is_nz(i, j)) continue;
        StoredDistance sim = (StoredDistance)graph_sim(i, j);
        if (sim > (Similarity)cutoff &&
            (heap.size() < (size_t)maxEdges ||
             (heap.size() == (size_t)maxEdges && sim > heap.top().second))) {
          if (heap.size() == (size_t)maxEdges) heap.pop();
          heap.push(std::make_pair(j, sim));
        }
      }
      while (!heap.empty()) {
        auto e = heap.top(); heap.pop();
        if (i < e.first) { ls.push_back(e.second); lf.push_back(i); lt.push_back(e.first); }
      }

      if (verbose && tid == 0) {
        progress.track(numThreads);
        if (progress.isStepMarker())
          verbose_message("Building Similarity Graph %s [%.1fGb / %.1fGb]\r",
                          progress.getProgress(),
                          getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
      }
    }
  }

  // Merge per-thread edge lists.
  size_t total = 0;
  for (int t = 0; t < (int)numThreads; ++t) total += tl_from[t].size();
  from.reserve(total); to.reserve(total); sComp.reserve(total);
  for (int t = 0; t < (int)numThreads; ++t) {
    from.insert(from.end(), tl_from[t].begin(), tl_from[t].end());
    to  .insert(to  .end(), tl_to[t]  .begin(), tl_to[t]  .end());
    sComp.insert(sComp.end(), tl_sComp[t].begin(), tl_sComp[t].end());
  }
  verbose_message("Finished Building Similarity Graph (%zu edges) "
                  "[%.1fGb / %.1fGb]                         \n",
                  g.getEdgeCount(),
                  getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
  g.sComp.shrink_to_fit(); g.to.shrink_to_fit(); g.from.shrink_to_fit();
}

void build_similarity_graph(Graph &g, Similarity cutoff) {
  std::vector<size_t> &from = g.from;
  std::vector<size_t> &to   = g.to;
  auto &sComp = g.sComp;

  if (nobs == 0) return;

  // NOTE: the PMH winner index (g_pmh_idx) is NOT used for graph build.
  // With k=6 there are only 4096 distinct k-mers, so posting lists average
  // ~1748 contigs — dense, not sparse. Candidate-gen work exceeds the all-pairs
  // SIMD kernel. The all-pairs AVX-512 path (below) is optimal for this regime.

  // Default to the all-pairs b-bit popcount build: the OPH signature is so
  // compact (m/8 bytes, cache-resident) that the all-pairs SIMD popcount kernel
  // is faster than inverted-index candidate generation at these contig counts,
  // where skewed posting lists (popular (bucket,min) keys) cause a sum-of-L^2
  // traversal blow-up. Set RABBIT_GRAPH_INDEX=1 to use the inverted-index build
  // (wins when N is much larger and/or the per-pair metric is expensive).
  bool use_index = false;
  if (const char *e = rb_getenv("RABBIT_GRAPH_INDEX"))
    use_index = (e[0] == '1');
  if (!use_index) { build_graph_allpairs(g, cutoff); return; }

  // Number of OPH buckets m (= sketch size). The collision count between two
  // sketches divided by m is the Jaccard estimate, so an edge with similarity
  // > cutoff requires at least ceil(cutoff * m) colliding buckets. That lower
  // bound drives the inverted index: only pairs sharing >= minCommon buckets
  // are even considered, replacing the O(N^2) all-pairs Jaccard scan.
  const uint32_t M = g_sketches[0]->getK();
  const double   Md = (double)M;
  // A pair whose true Jaccard exceeds `cutoff` collides in roughly cutoff*m
  // buckets. We relax this lower bound by a safety margin so estimator noise
  // can't drop a genuine edge below the candidate threshold; the exact Jaccard
  // (same metric the cutoff was calibrated against) is then evaluated only for
  // the surviving candidates — never for all N^2 pairs.
  size_t minCommon = (size_t)std::floor((double)cutoff * Md * 0.7);
  if (minCommon < 1) minCommon = 1;
  // Hybrid experiment: the cutoff-derived bound assumes the candidate metric IS
  // the OPH Jaccard, but in PMH mode the index is only a sparse k-mer prefilter
  // (candidate = shares >= minCommon exact OPH buckets) and the edge weight is
  // the k=6 PMH winner-match.  Allow forcing minCommon to probe that regime.
  if (const char *e = rb_getenv("RABBIT_MINCOMMON")) {
    long v = std::atol(e);
    if (v >= 1) minCommon = (size_t)v;
  }

  // ── Use the index that was built inline during sketch construction ────────
  if (!g_inv_idx) {
    verbose_message("WARN: inverted index not ready, falling back to all-pairs\n");
    build_graph_allpairs(g, cutoff); return;
  }
  const rabbit_invidx::InvertedIndex& idx = *g_inv_idx;

  verbose_message("Starting Building Similarity Graph (inverted index). "
                  "nobs=%d maxEdges=%d buckets=%d cutoff=%.4f minCommon=%d "
                  "postings=%zu keys=%zu\n",
                  nobs, maxEdges, M, (double)cutoff, minCommon,
                  idx.totalPostings, idx.postIdx.size());

  // ── Posting-list traversal: count collisions per candidate, keep the
  //    top-maxEdges neighbours per node, emit i<j edges.
  //    Per-sketch keys are generated on-the-fly (getKeys) – no global skKeys.
  const uint32_t* csrPtr = idx.csrPosts.data();

  std::vector<std::vector<size_t>>         tl_from(numThreads);
  std::vector<std::vector<size_t>>         tl_to(numThreads);
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
        if (sim > cutoff &&
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

// ═══════════════════════════════════════════════════════════════════════════
// ── calibrate_sim_cutoff_converge ──────────────────────────────────────────
// Inner: given maxsim[] for _nobs sampled contigs, binary-search for the
// similarity cutoff p such that exactly coverage fraction are "connected"
// (i.e. have at least one neighbour ≥ cutoff). Returns p×1000.
static size_t calib_converge(const std::vector<Similarity>& maxsim,
                              Distance coverage) {
  size_t _nobs = maxsim.size();
  size_t p = 999, pp = 1000;
  Distance cov = 0, pcov = 0;
  for (; p > RB_SIM_FLOOR;) {
    Distance cutoff = (Distance)p / 1000.;
    size_t counton = 0;
    for (size_t i = 0; i < _nobs; ++i)
      if (maxsim[i] >= cutoff) counton++;
    cov = (Distance)counton / _nobs;
    if (cov >= coverage) {
      if (cov - coverage > coverage - pcov) { p = pp; cov = pcov; }
      break;
    } else {
      verbose_message("Preparing Similarity Graph Building [pSim = %2.1f; "
                      "%zu / %zu (P = %2.2f%%)]               \r",
                      p / 10., counton, _nobs, cov * 100.);
    }
    pp = p; pcov = cov;
    if (p > 990)       p -= rand() % 3 + 1;
    else if (p > 900)  p -= rand() % 3 + 3;
    else               p -= rand() % 3 + 9;
  }
  return p;
}

// pmh_calib_maxsim – compute maxsim[] for a sample of contigs using the PMH
// winner index.  For each sampled contig, walks its winner registers, accumulates
// hit counts per candidate, and records the maximum corrected similarity.
// O(sample × avg_posting_size) instead of O(sample × N).
static void pmh_calib_maxsim(const std::vector<size_t>& sample,
                              std::vector<Similarity>& maxsim) {
  const rabbit_invidx::InvertedIndex &idx = *g_pmh_idx;
  const uint32_t *csrPtr = idx.csrPosts.data();
  const double    Md     = (double)g_pmh_m;
  const double    b0     = (g_pmh_base_on && g_pmh_baseline > 0.0 && g_pmh_baseline < 1.0)
                           ? g_pmh_baseline : 0.0;
  const size_t    S      = sample.size();

#pragma omp parallel num_threads(numThreads)
  {
    std::vector<int>      stamp(nobs, 0);
    std::vector<uint32_t> isect(nobs, 0);
    int ep = 0;
    std::vector<size_t> cand;
    cand.reserve(8192);

#pragma omp for schedule(dynamic, 8)
    for (size_t si = 0; si < S; ++si) {
      const size_t i = sample[si];
      cand.clear();
      if (++ep == INT_MAX) { std::fill(stamp.begin(), stamp.end(), 0); ep = 1; }

      static constexpr uint64_t GOLDEN_C = UINT64_C(0x9E3779B97F4A7C15);
      const uint64_t *wi64 = g_win64_flat.data() + i * g_pmh_m;
      for (uint32_t pos = 0; pos < g_pmh_m; ++pos) {
        if (wi64[pos] == 0ULL) continue;
        const uint64_t key = wi64[pos] ^ ((uint64_t)pos * GOLDEN_C);
        auto it = idx.postIdx.find(key);
        if (it == idx.postIdx.end()) continue;
        const uint32_t *pl = csrPtr + it->second.off;
        const uint32_t  sz = it->second.cnt;
        for (uint32_t pi = 0; pi < sz; ++pi) {
          size_t j = pl[pi];
          if (j == i) continue;
          if (stamp[j] != ep) { stamp[j] = ep; isect[j] = 1; cand.push_back(j); }
          else                 { ++isect[j]; }
        }
      }

      // maxsim = highest corrected similarity among all candidates.
      uint32_t best = 0;
      for (size_t j : cand) if (isect[j] > best) best = isect[j];
      double raw = (double)best / Md;
      double sim = (b0 < 1.0) ? (raw - b0) / (1.0 - b0) : 0.0;
      if (sim < 0.0) sim = 0.0;
      if (sim > 1.0 - 1e-6) sim = 1.0 - 1e-6;
      maxsim[si] = (Similarity)sim;
    }
  }
}

// ── Fusion D: fused calibration + graph build ─────────────────────────────
// Replaces the two sequential O(N²) passes:
//   • calibrate_sim_cutoff_fused : 25,000 × N = 771M pair-sims
//   • build_graph_allpairs       : N²/2       = 475M pair-sims
//
// A single tiled N²/2 pass simultaneously:
//   1. Accumulates per-contig top-maxEdges neighbor heaps (uncutoff).
//   2. Accumulates maxsim[] for 10 × 2500 calibration samples.
// Then calibrates (→ simCutoff) and emits edges from stored heaps.
//
// Pair-sim count: 475M vs 1,246M → ~2.6× reduction in the two heaviest phases.
// Memory: N × maxEdges heaps ≈ 30841 × 200 × 12B ≈ 74 MB (negligible).
//
// Returns the calibrated simCutoff×10 and fills Graph g with edges.
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

static Distance gen_fused_calib_graph(Graph &g, Distance coverage) {
  static constexpr size_t NROUNDS = 10;
  static constexpr size_t SAMP    = 2500;

  // Draw NROUNDS independent 2500-element samples (permuted prefix).
  std::vector<std::vector<size_t>> idxs(NROUNDS, std::vector<size_t>(nobs));
  for (size_t r = 0; r < NROUNDS; ++r) {
    std::iota(idxs[r].begin(), idxs[r].end(), 0);
    random_unique(idxs[r].begin(), idxs[r].end(), SAMP);
  }

  // Flat sample membership: sample_round[i] = round index if contig i is in
  // that round's sample, else -1.  Used for O(1) lookup in the inner loop.
  // A contig may appear in multiple rounds; we record the lowest-index round.
  std::vector<int8_t> sample_round(nobs, -1);   // -1 = not sampled
  std::vector<size_t> sample_local(nobs, 0);    // local index within sample_round
  for (int r = (int)NROUNDS - 1; r >= 0; --r)  // lowest round wins
    for (size_t li = 0; li < SAMP; ++li) {
      const size_t ci = idxs[r][li];
      sample_round[ci] = (int8_t)r;
      sample_local[ci] = li;
    }

  // maxsim[r][li]: max similarity seen by sample contig li in round r.
  std::vector<std::vector<float>> maxsim(NROUNDS, std::vector<float>(SAMP, 0.f));

  // Per-contig neighbor heap (min-heap of size ≤ maxEdges; top = weakest kept).
  // Triangle pass: each unordered pair {i,j} is evaluated exactly once and used
  // to update BOTH heaps[i] and heaps[j] (graph_sim is symmetric), halving the
  // similarity computations vs the old full-N² scan.  Because every heap can now
  // be touched by multiple threads, each row gets a lock + an atomic threshold:
  // row_thresh[r] mirrors heaps[r].top() once the heap is full, so the common
  // low-similarity candidate is rejected with a single relaxed atomic read and
  // never contends for the lock.  Correctness is re-verified inside the lock.
  using Heap = std::priority_queue<Edge, std::vector<Edge>, CompareEdge>;
  std::vector<Heap> heaps(nobs);

  std::vector<std::atomic<float>> row_thresh(nobs);
  for (size_t r = 0; r < nobs; ++r)
    row_thresh[r].store(std::numeric_limits<float>::lowest(),
                        std::memory_order_relaxed);
  // Per-row spinlock (1 byte) instead of omp_lock_t. The critical section is a
  // single heap push/pop, so a futex-backed omp_lock is overkill: its acquire
  // path is a syscall under contention and its cache line bounces expensively
  // across sockets. A relaxed test-and-test-and-set spinlock keeps the hot path
  // entirely in user space and only touches one cache line per row. Most
  // candidates never reach here (rejected by the row_thresh atomic fast-path),
  // so contention is brief.
  std::vector<std::atomic<uint8_t>> row_spin(nobs);
  for (size_t r = 0; r < nobs; ++r)
    row_spin[r].store(0, std::memory_order_relaxed);

  // Per-thread maxsim accumulators — both endpoints of every pair feed these;
  // merged into the global maxsim[][] after the pass (max is order-independent,
  // so calibration is bit-identical regardless of pair processing order).
  std::vector<std::vector<std::vector<float>>> thread_jmax(
      numThreads,
      std::vector<std::vector<float>>(NROUNDS, std::vector<float>(SAMP, 0.f)));

  // Tile size: same cache heuristic as before.
  size_t TILE = 10;
  try {
    TILE = std::max(
        (size_t)((CacheSize() * 1024.) /
                 (2 * (g_pmh_m * sizeof(uint32_t)) +
                  maxEdges * (sizeof(size_t) + sizeof(StoredDistance)))),
        (size_t)10);
  } catch (...) {}

  // Flatten the upper-triangle tile-pairs (jj ≥ ii) into a list so dynamic
  // scheduling balances load evenly (each tile-pair is ~equal work).
  std::vector<std::pair<uint32_t, uint32_t>> tilepairs;
  {
    const size_t ntiles = (nobs + TILE - 1) / TILE;
    tilepairs.reserve(ntiles * (ntiles + 1) / 2);
    for (size_t ii = 0; ii < nobs; ii += TILE)
      for (size_t jj = ii; jj < nobs; jj += TILE)
        tilepairs.emplace_back((uint32_t)ii, (uint32_t)jj);
  }

  verbose_message("Fusion D (triangle): single-pass calib+graph nobs=%zu "
                  "SAMP=%zu×%zu TILE=%zu tilepairs=%zu\n",
                  nobs, NROUNDS, SAMP, TILE, tilepairs.size());

  const bool rb_timing = (getenv("RB_TIMING") != nullptr);
  std::chrono::steady_clock::time_point _t_pass0;
  if (rb_timing) _t_pass0 = std::chrono::steady_clock::now();

  // ── Abundance-first exact upper-bound prune (RABBIT_ABDFIRST=1) ───────────
  // The final edge weight is w = g_w_comp·sComp + (1−g_w_comp)·corr, kept only
  // if w ≥ min_edge_weight.  sComp ≤ 1, so the best-case score for a pair is
  // g_w_comp·1 + (1−g_w_comp)·corr.  If even that is below the cutoff the pair
  // can NEVER become an edge — so we can skip the (≈100× more expensive) PMH
  // composition kernel after only an O(num_samples) abundance correlation.
  // This is exact w.r.t. the fused-edge filter: no pruned pair could survive.
  // corr threshold: corr < (min_edge_weight − g_w_comp)/(1 − g_w_comp).
  // Abundance-first prune is ON by default for multi-sample data (it is exact
  // w.r.t. the fused-edge filter and substantially improves multi-sample
  // binning by stopping high-composition / low-abundance pairs from starving
  // the composition top-k).  RABBIT_NO_ABDFIRST=1 restores the old all-pairs
  // composition pass for comparison.
  const bool abdfirst = (getenv("RABBIT_NO_ABDFIRST") == nullptr) &&
                        num_depth_samples > 1 && g_w_comp < 1.0 &&
                        g_w_comp < min_edge_weight;
  const double abd_corr_min = abdfirst
      ? (min_edge_weight - g_w_comp) / (1.0 - g_w_comp) : -2.0;
  // Precompute per-contig unit rank vectors u_i so the abundance correlation
  // corr(i,j) = Σ_k u_i[k]·u_j[k] is a single length-S dot product — bit-equal
  // to the Pearson-on-ranks that cal_depth_corr() computes (depth_matrix is
  // already rank-transformed), but ~100× cheaper than its per-pair Welford+sqrt.
  // A small safety margin keeps the prune conservative against float rounding:
  // we only skip when the fast estimate is clearly below the exact bound, so no
  // pair that could pass the fused filter is ever dropped.
  const uint32_t ABD_S = (uint32_t)num_depth_samples;
  const float abd_corr_min_eps = (float)(abd_corr_min - 1e-4);
  std::vector<float> abd_unit;
  if (abdfirst) {
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
    verbose_message("Abundance-first prune: skip pairs with depth corr < %.4f "
                    "(g_w_comp=%.2f, min_edge=%.2f)\n",
                    abd_corr_min, g_w_comp, (double)min_edge_weight);
  }
  const float *abd_u = abd_unit.empty() ? nullptr : abd_unit.data();

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
    auto &my_jmax = thread_jmax[tid];
    std::vector<uint32_t> visited;       // LSH per-i dedup stamp
    uint32_t visit_gen = 0;
    if (use_lsh) visited.assign(nobs, 0u);

    // Thread-safe update of row r's neighbor heap with candidate (other, sv).
    auto update_row = [&](size_t r, size_t other, StoredDistance sv) {
      if (sv <= row_thresh[r].load(std::memory_order_relaxed)) return;
      // ── acquire spinlock (test-and-test-and-set) ──
      while (row_spin[r].exchange(1, std::memory_order_acquire)) {
        while (row_spin[r].load(std::memory_order_relaxed))
#if defined(__x86_64__) || defined(__i386__)
          __builtin_ia32_pause();
#else
          ;
#endif
      }
      Heap &h = heaps[r];
      if (h.size() < (size_t)maxEdges) {
        h.push({other, sv});
        if (h.size() == (size_t)maxEdges)
          row_thresh[r].store(h.top().second, std::memory_order_relaxed);
      } else if (sv > h.top().second) {
        h.pop();
        h.push({other, sv});
        row_thresh[r].store(h.top().second, std::memory_order_relaxed);
      }
      row_spin[r].store(0, std::memory_order_release);  // release spinlock
    };

    if (!use_lsh) {
#pragma omp for schedule(dynamic, 1)
    for (size_t p = 0; p < tilepairs.size(); ++p) {
      const size_t ii = tilepairs[p].first;
      const size_t jj = tilepairs[p].second;
      const size_t i_stop = std::min(ii + TILE, nobs);
      const size_t j_stop = std::min(jj + TILE, nobs);
      const bool diag = (ii == jj);
      for (size_t i = ii; i < i_stop; ++i) {
        const int    i_round = sample_round[i];
        const size_t i_local = sample_local[i];
        // Diagonal tile: only j>i; off-diagonal: full j range.
        for (size_t j = (diag ? i + 1 : jj); j < j_stop; ++j) {
          if (!is_nz(i, j)) continue;
          // Abundance-first exact upper-bound prune: skip the PMH kernel for
          // pairs whose best-case fused score is already below the cutoff.
          if (abdfirst) {
            const float *ui = abd_u + i * ABD_S;
            const float *uj = abd_u + j * ABD_S;
            float c = 0.0f;
            for (uint32_t k = 0; k < ABD_S; ++k) c += ui[k] * uj[k];
            if (c < abd_corr_min_eps) continue;
          }
          StoredDistance sv = (StoredDistance)graph_sim(i, j);
          update_row(i, j, sv);
          update_row(j, i, sv);
          if (i_round >= 0) {
            float &m = my_jmax[i_round][i_local];
            if ((float)sv > m) m = (float)sv;
          }
          const int j_round = sample_round[j];
          if (j_round >= 0) {
            float &m = my_jmax[j_round][sample_local[j]];
            if ((float)sv > m) m = (float)sv;
          }
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
        const int    i_round = sample_round[i];
        const size_t i_local = sample_local[i];
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
            update_row(i, j, sv);
            update_row(j, i, sv);
            if (i_round >= 0) {
              float &m = my_jmax[i_round][i_local];
              if ((float)sv > m) m = (float)sv;
            }
            const int j_round = sample_round[j];
            if (j_round >= 0) {
              float &m = my_jmax[j_round][sample_local[j]];
              if ((float)sv > m) m = (float)sv;
            }
          }
        }
        if (verbose && tid == 0 && (i & 0xFFFu) == 0)
          verbose_message("LSH cand %zu/%zu\r", i, nobs);
      }
    }
  }

  if (rb_timing) {
    double ms = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - _t_pass0).count();
    fprintf(stderr, "[RB_TIMING] fusedD pair-pass: %.1f ms\n", ms);
  }

  { std::vector<std::atomic<float>>().swap(row_thresh); }
  { std::vector<std::atomic<uint8_t>>().swap(row_spin); }

  // Merge per-thread maxsim into the global array.
  for (size_t r = 0; r < NROUNDS; ++r)
    for (size_t li = 0; li < SAMP; ++li)
      for (int tid = 0; tid < (int)numThreads; ++tid)
        if (thread_jmax[tid][r][li] > maxsim[r][li])
          maxsim[r][li] = thread_jmax[tid][r][li];
  thread_jmax.clear();

  // ── Calibrate: converge each round independently. ─────────────────────────
  Distance sum_p = 0;
  for (size_t r = 0; r < NROUNDS; ++r) {
    std::vector<Similarity> sub(SAMP);
    for (size_t li = 0; li < SAMP; ++li) sub[li] = (Similarity)maxsim[r][li];
    Distance _minp = (Distance)calib_converge(sub, coverage);
    if (_minp < (Distance)(RB_SIM_FLOOR + 1)) _minp = (Distance)RB_SIM_FLOOR;
    sum_p += _minp;
    if (r == 1 && sum_p / 2 < (Distance)(RB_SIM_FLOOR + 1))
      { sum_p = (Distance)RB_SIM_FLOOR * (Distance)NROUNDS; break; }
  }
  Distance simCutoff = sum_p / (Distance)NROUNDS;

  verbose_message("Finished Preparing Similarity Graph Building [pSim = %2.2f] "
                  "[%.1fGb / %.1fGb]                                      \n",
                  simCutoff / 10., getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

  // ── Emit edges: drain heaps, apply cutoff, emit only i<j. ────────────────
  const Similarity cutoff = (Similarity)(simCutoff / 1000.);
  verbose_message("Starting Building Similarity Graph (Fusion D). "
                  "nobs=%zu cutoff=%.4f\n", nobs, (double)cutoff);

  std::vector<size_t>         &from  = g.from;
  std::vector<size_t>         &to    = g.to;
  std::vector<StoredDistance> &sComp  = g.sComp;

  if (!g_mutual_knn) {
    // Union k-NN (original): emit edge (i,j) if j ∈ top-k(i) OR i ∈ top-k(j)
    for (size_t i = 0; i < nobs; ++i) {
      while (!heaps[i].empty()) {
        auto e = heaps[i].top(); heaps[i].pop();
        size_t j   = e.first;
        StoredDistance sv = e.second;
        if (sv <= cutoff) continue;
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
    // Parallel + compact: instead of one ~16M-entry hash map (serial build +
    // serial probe, plus a fat 24 B/entry footprint and cross-NUMA reads on the
    // thread-scattered heaps), we drain each row's heap into a per-row sorted
    // vector of (neighborId, sim). Membership "is i a neighbor of j?" becomes a
    // branch-light binary search in contiguous 8 B pairs (≈128 MB vs 200-400 MB
    // for the map). Both the drain and the mutual-emit passes parallelise over
    // rows, and the edge order is now deterministic (sorted by i then j).
    typedef std::pair<uint32_t, StoredDistance> Nbr;   // (neighborId, sim)
    std::vector<std::vector<Nbr>> nbr(nobs);

    // (1) Parallel drain: heaps[i] → sorted nbr[i] (ascending neighbor id).
#pragma omp parallel for num_threads(numThreads) schedule(dynamic, 64)
    for (size_t i = 0; i < nobs; ++i) {
      Heap &h = heaps[i];
      std::vector<Nbr> &row = nbr[i];
      row.reserve(h.size());
      while (!h.empty()) {
        const Edge &e = h.top();
        if (e.second > cutoff)
          row.emplace_back((uint32_t)e.first, e.second);
        h.pop();
      }
      Heap().swap(h);   // release this row's heap buffer early
      std::sort(row.begin(), row.end(),
                [](const Nbr &a, const Nbr &b) { return a.first < b.first; });
    }

    // (2) Parallel emit: for each i, keep (i,j) with i<j iff i ∈ nbr[j] too.
    std::vector<std::vector<size_t>>         tl_from(numThreads);
    std::vector<std::vector<size_t>>         tl_to(numThreads);
    std::vector<std::vector<StoredDistance>> tl_sv(numThreads);
#pragma omp parallel num_threads(numThreads)
    {
      const int tid = omp_get_thread_num();
      auto &lf = tl_from[tid]; auto &lt = tl_to[tid]; auto &ls = tl_sv[tid];
#pragma omp for schedule(dynamic, 64)
      for (size_t i = 0; i < nobs; ++i) {
        const std::vector<Nbr> &row = nbr[i];
        for (const Nbr &nj : row) {
          const uint32_t j = nj.first;
          if ((size_t)j <= i) continue;          // emit each undirected pair once
          const std::vector<Nbr> &rj = nbr[j];   // is i a neighbor of j?
          auto it = std::lower_bound(
              rj.begin(), rj.end(), (uint32_t)i,
              [](const Nbr &p, uint32_t v) { return p.first < v; });
          if (it != rj.end() && it->first == (uint32_t)i) {
            lf.push_back(i); lt.push_back((size_t)j); ls.push_back(nj.second);
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
  return simCutoff;
}

// calibrate_sim_cutoff  –  sample-based auto-calibration of similarity cutoff
// ═══════════════════════════════════════════════════════════════════════════
// `full=true` uses all nobs contigs as both sample and reference (medium path).
// `full=false` uses a 2500-sample subset (high path, called 10 times).
size_t calibrate_sim_cutoff(Distance coverage, bool full) {
  size_t _nobs = full ? nobs : std::min(nobs, (size_t)2500);

  std::vector<size_t> idx(nobs);
  std::iota(idx.begin(), idx.end(), 0);
  random_unique(idx.begin(), idx.end(), _nobs);
  idx.resize(_nobs);

  std::vector<Similarity> maxsim(_nobs, (Similarity)0);

#pragma omp parallel for schedule(dynamic, 4)
  for (size_t i = 0; i < _nobs; ++i) {
    const size_t contig_i = idx[i];
    Similarity mx = (Similarity)0;
    for (size_t j = 0; j < nobs; ++j) {
      if (j == i) continue;
      Similarity s = (Similarity)graph_sim(contig_i, idx[j]);
      if (s > mx) mx = s;
    }
    maxsim[i] = mx;
  }
  return calib_converge(maxsim, coverage);
}

// calibrate_sim_cutoff_fused ── Fusion C ───────────────────────────────
// Replaces the 10-serial-call loop used when nobs > 25000.
// Draws 10 independent 2500-sample subsets at once, computes all 25000
// maxsim values in a SINGLE parallel OMP pass, then independently converges
// each sample to find its cutoff.  Benefits vs 10 serial calls:
//   • 1 OMP region / barrier instead of 10   (less scheduling overhead)
//   • 10× better load balancing (25000 tasks vs 2500 tasks per region)
//   • g_sig_flat access pattern more linear (less cache thrashing)
// Returns the same average-p that the 10-call loop would produce.
static Distance calibrate_sim_cutoff_fused(Distance coverage) {
  static constexpr size_t NROUNDS = 10;
  static constexpr size_t SAMP    = 2500;  // per round

  // Draw NROUNDS independent random samples (each a permuted prefix of idx)
  std::vector<std::vector<size_t>> idxs(NROUNDS, std::vector<size_t>(nobs));
  for (size_t r = 0; r < NROUNDS; ++r) {
    std::iota(idxs[r].begin(), idxs[r].end(), 0);
    random_unique(idxs[r].begin(), idxs[r].end(), SAMP);
  }

  const size_t total = NROUNDS * SAMP;
  std::vector<Similarity> maxsim(total, (Similarity)0);

  // Single parallel pass over all NROUNDS*SAMP sampled contigs.
  // Task i = (round r, local index li): r = i/SAMP, li = i%SAMP.
#pragma omp parallel for schedule(dynamic, 8)
  for (size_t i = 0; i < total; ++i) {
    const size_t r  = i / SAMP;
    const size_t li = i % SAMP;
    const size_t contig_i = idxs[r][li];
    Similarity mx = (Similarity)0;
    for (size_t j = 0; j < nobs; ++j) {
      if (j == li) continue;
      Similarity s = (Similarity)graph_sim(contig_i, idxs[r][j]);
      if (s > mx) mx = s;
    }
    maxsim[i] = mx;
  }

  // Converge each round independently and accumulate.
  Distance sum_p = 0;
  for (size_t r = 0; r < NROUNDS; ++r) {
    const std::vector<Similarity> sub(maxsim.begin() + r * SAMP,
                                       maxsim.begin() + (r + 1) * SAMP);
    Distance _minp = (Distance)calib_converge(sub, coverage);
    if (_minp < (Distance)(RB_SIM_FLOOR + 1)) _minp = (Distance)RB_SIM_FLOOR;
    sum_p += _minp;
    if (r == 1 && sum_p / 2 < (Distance)(RB_SIM_FLOOR + 1)) {
      return (Distance)RB_SIM_FLOOR;
    }
  }
  return sum_p / (Distance)NROUNDS;
}

// ═══════════════════════════════════════════════════════════════════════════
// Similarity calibration helpers
// ═══════════════════════════════════════════════════════════════════════════

