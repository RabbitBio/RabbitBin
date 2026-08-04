// RabbitBin module: rb_split.cpp

static void marker_guided_split(BinMap &cls) {
  if (marker_seed_file.empty()) return;

  // contig name -> unified index (large: i; small: j + nobs)
  std::unordered_map<std::string, size_t> name2idx;
  name2idx.reserve((nobs + nobs1) * 2);
  auto add_name = [&](const std::string &raw, size_t idx) {
    std::string nm = raw;
    size_t sp = nm.find_first_of(" \t");
    if (sp != std::string::npos) nm.resize(sp);
    name2idx[nm] = idx;
  };
  for (size_t i = 0; i < nobs; ++i)  add_name(contig_names[i], i);
  for (size_t j = 0; j < nobs1; ++j) add_name(small_contig_names[j], j + nobs);

  std::ifstream fin(marker_seed_file);
  if (!fin) {
    cerr << "[Warn] cannot open --marker-seed file: " << marker_seed_file
         << " (skipping marker-guided split)\n";
    return;
  }
  std::unordered_map<size_t, std::vector<int>> contig_markers;
  std::string line;
  int marker_id = 0;
  size_t hits = 0;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::string tok;
    bool first = true, any = false;
    while (std::getline(ss, tok, '\t')) {
      if (first) { first = false; continue; }     // marker name column
      if (tok.empty()) continue;
      auto it = name2idx.find(tok);
      if (it != name2idx.end()) {
        contig_markers[it->second].push_back(marker_id);
        ++hits; any = true;
      }
    }
    if (any) ++marker_id;
  }
  if (marker_id == 0 || num_depth_samples < 1) {
    verbose_message("Marker-guided split: no usable markers (%d) — skipped\n",
                    marker_id);
    return;
  }

  // Both depth_matrix and small_depth_matrix are rank-transformed in place during the pipeline,
  // so read the raw-mean snapshots taken before those transforms.
  auto depth_at = [&](size_t c, size_t i) -> double {
    if (c < nobs) {
      if (!g_large_means.empty() && c * (size_t)num_depth_samples + i < g_large_means.size())
        return (double)g_large_means[c * num_depth_samples + i];
      return (double)depth_matrix(c, i);
    }
    size_t s = c - nobs;
    if (!g_small_means.empty() && s * (size_t)num_depth_samples + i < g_small_means.size())
      return (double)g_small_means[s * num_depth_samples + i];
    return (double)small_depth_matrix(s, i);
  };

  auto bin_bp = [&](const ContigVector &contigs) -> size_t {
    size_t bp = 0;
    for (size_t c : contigs)
      bp += (c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs];
    return bp;
  };

  std::mt19937 rng((unsigned)(seed ? seed : 1ULL));
  BinMap out;
  int next = 0;
  size_t n_split = 0, n_kept = 0, n_drop = 0;
  for (auto &kv : cls) {
    const ContigVector &contigs = kv.second;
    // Size gate: only refine bins that would pass the min_bin_bp output filter.
    if (bin_bp(contigs) < min_bin_bp) { ++n_drop; continue; }
    // marker multiplicity within this bin
    std::unordered_map<int, int> mc;
    int mult = 1;
    for (size_t c : contigs) {
      auto it = contig_markers.find(c);
      if (it == contig_markers.end()) continue;
      for (int m : it->second) { int v = ++mc[m]; if (v > mult) mult = v; }
    }
    if (mult <= 1 || (int)contigs.size() < splitMinContigs) {
      out[next++] = contigs; ++n_kept; continue;
    }
    int k = std::min((int)contigs.size(), std::min(mult, splitMaxK));
    const size_t n = contigs.size();
    std::vector<float> X(n * num_depth_samples);
    for (size_t r = 0; r < n; ++r)
      for (size_t i = 0; i < num_depth_samples; ++i)
        X[r * num_depth_samples + i] = (float)std::log(depth_at(contigs[r], i) + 1.0);
    std::vector<int> labels = rb_kmeans(X.data(), n, num_depth_samples, k, rng);
    std::vector<ContigVector> sub(k);
    for (size_t r = 0; r < n; ++r) sub[labels[r]].push_back(contigs[r]);
    for (int t = 0; t < k; ++t)
      if (!sub[t].empty()) out[next++] = std::move(sub[t]);
    ++n_split;
  }
  cls.swap(out);
  verbose_message("Marker-guided split: %d markers (%zu hits), %zu kept, "
                  "%zu split, %zu dropped(<min_bin_bp) -> %d bins\n",
                  marker_id, hits, n_kept, n_split, n_drop, next);
  // Split products may be < min_bin_bp; they are legitimate genome bins, so
  // disable the size filter in output_bins (the base bins were already gated).
  min_bin_bp = 0;
}

// ═══════════════════════════════════════════════════════════════════════════
// Marker-FREE bin splitting (Phase 2; --split-bins) — abundance multimodality
// ═══════════════════════════════════════════════════════════════════════════
// RB_KM_PROF accumulators: which part of split_bin the time actually goes to.
static std::atomic<double> g_ms_feat{0}, g_ms_kmeans{0}, g_ms_sil{0}, g_ms_coh{0};
static inline void rb_prof_add(std::atomic<double> &a,
                               std::chrono::steady_clock::time_point t0) {
  const double v = std::chrono::duration<double, std::milli>(
                       std::chrono::steady_clock::now() - t0).count();
  double cur = a.load(std::memory_order_relaxed);
  while (!a.compare_exchange_weak(cur, cur + v, std::memory_order_relaxed)) {}
}

// Mean silhouette of a labeling over feature rows X (Euclidean).  O(n^2 d); for
// large bins we sample up to g_sil_sample_cap rows so this stays cheap.
static double fkmv_silhouette(const float *X, size_t n, size_t d,
                              const std::vector<int> &labels, int k,
                              std::mt19937 &rng) {
  if (n < 3 || k < 2) return -1.0;
  // sample indices if too large
  std::vector<size_t> idx;
  if (n > g_sil_sample_cap) {
    idx.resize(n);
    for (size_t i = 0; i < n; ++i) idx[i] = i;
    std::shuffle(idx.begin(), idx.end(), rng);
    idx.resize(g_sil_sample_cap);
  } else {
    idx.resize(n);
    for (size_t i = 0; i < n; ++i) idx[i] = i;
  }
  const size_t m = idx.size();
  auto dist = [&](size_t a, size_t b) -> double {
    return std::sqrt(rb_sqdist(X + a * d, X + b * d, d));
  };
  // per-cluster sizes within the sample
  std::vector<int> csz(k, 0);
  for (size_t a : idx) csz[labels[a]]++;
  double sil_sum = 0;
  size_t cnt = 0;
  std::vector<double> sumd(k);   // hoisted: was one allocation per sampled row
  for (size_t ii = 0; ii < m; ++ii) {
    size_t a = idx[ii];
    int la = labels[a];
    if (csz[la] <= 1) continue;  // singleton cluster: s=0, skip
    std::fill(sumd.begin(), sumd.end(), 0.0);
    for (size_t jj = 0; jj < m; ++jj) {
      if (jj == ii) continue;
      sumd[labels[idx[jj]]] += dist(a, idx[jj]);
    }
    double ai = sumd[la] / (double)(csz[la] - 1);
    double bi = std::numeric_limits<double>::infinity();
    for (int c = 0; c < k; ++c) {
      if (c == la || csz[c] == 0) continue;
      double mc = sumd[c] / (double)csz[c];
      if (mc < bi) bi = mc;
    }
    if (!std::isfinite(bi)) continue;
    double denom = std::max(ai, bi);
    if (denom > 0) { sil_sum += (bi - ai) / denom; ++cnt; }
  }
  return cnt ? sil_sum / (double)cnt : -1.0;
}

// Re-split internally multi-modal bins using per-sample log-abundance KMeans,
// choosing k by silhouette.  No markers / gene prediction needed.
static void abundance_guided_split_current(BinMap &cls) {
  if (num_depth_samples < 1) {
    verbose_message("Abundance split: needs >=1 abundance sample — skipped\n");
    return;
  }
  auto depth_at = [&](size_t c, size_t i) -> double {
    if (c < nobs) {
      if (!g_large_means.empty() && c * (size_t)num_depth_samples + i < g_large_means.size())
        return (double)g_large_means[c * num_depth_samples + i];
      return (double)depth_matrix(c, i);
    }
    size_t s = c - nobs;
    if (!g_small_means.empty() && s * (size_t)num_depth_samples + i < g_small_means.size())
      return (double)g_small_means[s * num_depth_samples + i];
    return (double)small_depth_matrix(s, i);
  };
  auto bin_bp = [&](const ContigVector &contigs) -> size_t {
    size_t bp = 0;
    for (size_t c : contigs)
      bp += (c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs];
    return bp;
  };

  // Mean intra-bin depth correlation over (up to 64) large-contig members,
  // reusing the precomputed unit-rank vectors via depth_corr_fast.  Small contigs
  // (index >= nobs) have no unit vector, so coherence is estimated on the large
  // members that anchor the bin.  Returns -2 when it cannot be estimated.
  auto bin_coherence = [&](const ContigVector &cs) -> double {
    constexpr size_t CAP = 64;
    std::vector<size_t> L;
    L.reserve(cs.size());
    for (size_t c : cs) if (c < nobs) L.push_back(c);
    if (L.size() < 2) return -2.0;
    if (L.size() > CAP) {                            // deterministic stride sample
      std::vector<size_t> s; s.reserve(CAP);
      const double step = (double)L.size() / (double)CAP;
      for (size_t t = 0; t < CAP; ++t) s.push_back(L[(size_t)(t * step)]);
      L.swap(s);
    }
    double sum = 0.0; size_t cnt = 0;
    for (size_t a = 0; a < L.size(); ++a)
      for (size_t b = a + 1; b < L.size(); ++b) {
        double c = depth_corr_fast(L[a], L[b]);
        if (std::isfinite(c)) { sum += c; ++cnt; }
      }
    return cnt ? sum / (double)cnt : -2.0;
  };

  // KMeans+silhouette split decision for one bin.  Returns true and fills `sub`
  // with >=2 sub-clusters when the bin is multi-modal (best silhouette >= thr);
  // false otherwise.  Identical maths to the historical split so the main path
  // stays byte-for-byte unchanged.
  auto split_bin = [&](const ContigVector &contigs, size_t bi,
                       std::vector<ContigVector> &sub) -> bool {
    ContigVector canonical;
    const ContigVector *work = &contigs;
    if (g_stable_split_kmeans) {
      canonical = contigs;
      std::sort(canonical.begin(), canonical.end());
      work = &canonical;
    }
    const ContigVector &items = *work;
    const size_t n = items.size();
    if ((int)n < splitMinContigs) return false;
    const size_t sd = (size_t)num_depth_samples;
    auto tk0 = std::chrono::steady_clock::now();
    std::vector<float> X(n * sd);
    for (size_t r = 0; r < n; ++r)
      for (size_t i = 0; i < sd; ++i)
        X[r * sd + i] = (float)std::log(depth_at(items[r], i) + 1.0);
    if (g_km_prof) rb_prof_add(g_ms_feat, tk0);
    uint64_t rng_seed = (seed ? (uint64_t)seed : 1ULL) + bi * 2654435761ULL;
    unsigned rng_seed32 = (unsigned)rng_seed; // preserve the established path
    if (g_stable_split_kmeans) {
      // FNV-1a over canonical contig indices: the same biological bin gets the
      // same RNG stream even when hash-map iteration changes its transient ID.
      uint64_t h = 1469598103934665603ULL;
      for (size_t c : items) {
        uint64_t x = (uint64_t)c;
        for (int b = 0; b < 8; ++b) {
          h ^= (unsigned char)(x & 0xffU);
          h *= 1099511628211ULL;
          x >>= 8;
        }
      }
      rng_seed = (seed ? (uint64_t)seed : 1ULL) ^ h;
      rng_seed32 = (unsigned)(rng_seed ^ (rng_seed >> 32));
    }
    std::mt19937 rng(rng_seed32);
    int best_k = 1; double best_sil = -1.0; std::vector<int> best_labels;
    int kmax = std::min((int)n - 1, splitMaxK);
    for (int k = 2; k <= kmax; ++k) {
      tk0 = std::chrono::steady_clock::now();
      std::vector<int> lab = rb_kmeans(X.data(), n, sd, k, rng);
      if (g_km_prof) rb_prof_add(g_ms_kmeans, tk0);
      int seen = 0; { std::vector<char> u(k, 0); for (int l : lab) if (!u[l]) { u[l] = 1; ++seen; } }
      if (seen < 2) continue;
      if (g_split_reject_small) {
        std::vector<size_t> child_n((size_t)k, 0), child_bp((size_t)k, 0);
        for (size_t r = 0; r < n; ++r) {
          const size_t lab_r = (size_t)lab[r];
          ++child_n[lab_r];
          const size_t c = items[r];
          child_bp[lab_r] += (c < nobs) ? seq_lens[c]
                                           : small_seq_lens[c - nobs];
        }
        bool valid = true;
        for (int j = 0; j < k; ++j) {
          if (child_n[(size_t)j] < splitMinSubContigs ||
              (splitMinSubBp > 0 && child_bp[(size_t)j] < splitMinSubBp)) {
            valid = false;
            break;
          }
        }
        if (!valid) continue;
      }
      tk0 = std::chrono::steady_clock::now();
      double s = fkmv_silhouette(X.data(), n, sd, lab, k, rng);
      if (g_km_prof) rb_prof_add(g_ms_sil, tk0);
      if (s > best_sil) { best_sil = s; best_k = k; best_labels = std::move(lab); }
    }
    if (best_k <= 1 || best_sil < g_split_sil) return false;
    sub.assign(best_k, ContigVector());
    for (size_t r = 0; r < n; ++r) sub[best_labels[r]].push_back(items[r]);
    return true;
  };

  // ── Parallel over bins ──────────────────────────────────────────────────
  // Each bin's KMeans+silhouette sweep is independent, so distribute bins
  // across threads (this was a serial Amdahl bottleneck). Determinism is kept
  // independent of thread count by seeding each bin's RNG from (seed, binIdx)
  // instead of one shared, sequentially-advanced generator.
  std::vector<const ContigVector*> binList;
  binList.reserve(cls.size());
  for (auto &kv : cls) binList.push_back(&kv.second);
  const size_t nbins = binList.size();

  // Coherence machinery usable only with multi-sample abundance + unit vectors;
  // otherwise the legacy absolute-bp discard runs and behaviour is unchanged.
  const bool coh_gate = g_split_keep_coherent && num_depth_samples > 1 &&
                        !g_depth_unit.empty();

  // kind: 0 = captured sub-floor candidate (post-processed below),
  //       1 = kept whole, 2 = split into emit[]
  struct BinResult { std::vector<ContigVector> emit; int kind = 1; };
  std::vector<BinResult> results(nbins);

  // ── Largest-bin-first scheduling order ────────────────────────────────────
  // split_bin's cost is dominated by uncapped KMeans over the bin's OWN contig
  // count. Under schedule(dynamic,1), a straggler megabin can start after smaller
  // bins are drained and leave other workers idle. Visiting
  // bins largest-first lets that bin start at t=0 and run concurrently with the
  // small-bin sweep, hiding most of its cost. `order` only changes ITERATION
  // ORDER: `bi` (used for both `results[bi]` and the per-bin RNG seed below) is
  // untouched, so every bin's KMeans/silhouette output is byte-identical to the
  // bi-ascending schedule — this cannot change accuracy, only wall time.
  std::vector<size_t> order(nbins);
  for (size_t i = 0; i < nbins; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return binList[a]->size() > binList[b]->size();
  });

  // RB_TIMING diagnostic: per-bin wall time, so the sum (total work) can be
  // compared against the loop's wall time (critical path) to tell a genuinely
  // saturated loop apart from one straggler bin holding 63 threads idle.
  const bool prof_bins = getenv("RB_SPLIT_PROF") != nullptr;
  std::vector<double> bin_ms(prof_bins ? nbins : 0, 0.0);

#pragma omp parallel for schedule(dynamic, 1) num_threads(numThreads)
  for (size_t oi = 0; oi < nbins; ++oi) {
    const size_t bi = order[oi];
    const ContigVector &contigs = *binList[bi];
    BinResult &res = results[bi];
    const auto t_bin0 = prof_bins ? std::chrono::steady_clock::now()
                                  : std::chrono::steady_clock::time_point{};
    struct BinTimer {
      const bool on; const std::chrono::steady_clock::time_point t0;
      double *slot;
      ~BinTimer() {
        if (on) *slot = std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - t0).count();
      }
    } bin_timer{prof_bins, t_bin0, prof_bins ? &bin_ms[bi] : nullptr};

    // Confidently sized bins follow the regular split path.
    if (bin_bp(contigs) < min_bin_bp) {
      // Historically dropped (contigs lost).  With the coherence gate, capture
      // them for an independent purification pass instead of discarding.
      if (coh_gate) res.kind = 0;                    // post-processed below
      else          res.kind = 1, res.emit.clear();  // legacy: emit nothing → drop
      continue;
    }
    std::vector<ContigVector> sub;
    bool do_split = split_bin(contigs, bi, sub);
    // ── Conservative split guard for large, internally-coherent bins ─────────
    // A >= g_split_guard_bp bin is usually one near-complete genome.  Veto the
    // proposed split unless it IMPROVES coherence (min child coherence exceeds
    // the parent's by g_split_guard_margin) — i.e. it separated two co-binned
    // genomes rather than cutting one genome in half.  Needs the same multi-
    // sample coherence machinery as the retention gate; otherwise behaviour is
    // unchanged.  Cost: a few bin_coherence calls on the rare large split
    // candidates only (O(64^2 * S) each), so no measurable runtime impact.
    if (do_split && g_split_guard && coh_gate &&
        bin_bp(contigs) >= g_split_guard_bp) {
      auto tc0 = std::chrono::steady_clock::now();
      double pc = bin_coherence(contigs);
      if (pc > -1.5) {                                // parent coherence usable
        double min_child = 2.0;
        for (auto &s : sub) {
          if (s.empty()) continue;
          double cc = bin_coherence(s);              // -2 (unscoreable) → veto
          if (cc < min_child) min_child = cc;
        }
        if (!(min_child >= pc + g_split_guard_margin)) do_split = false;
      }
      if (g_km_prof) rb_prof_add(g_ms_coh, tc0);
    }
    if (do_split) {
      for (auto &s : sub) {
        if (s.empty()) continue;
        if (bin_bp(s) < min_bin_bp) continue;        // legacy sub-product drop
        res.emit.push_back(std::move(s));
      }
      if (!res.emit.empty()) res.kind = 2;
      else { res.kind = 1; res.emit.push_back(contigs); }  // nothing survived
    } else {
      res.kind = 1; res.emit.push_back(contigs);
    }
  }

  if (g_km_prof) {
    fprintf(stderr,
            "[RB_KM_PROF] split_bin CPU-ms: features=%.0f kmeans=%.0f "
            "silhouette=%.0f coherence-guard=%.0f\n",
            g_ms_feat.load(), g_ms_kmeans.load(), g_ms_sil.load(),
            g_ms_coh.load());
    const uint64_t rs = g_km_restarts.load(), it = g_km_iters.load();
    fprintf(stderr,
            "[RB_KM_PROF] %lu kmeans restarts, %lu Lloyd iters (avg %.1f), "
            "%lu hit the 200 cap (%.1f%%)\n",
            (unsigned long)rs, (unsigned long)it, rs ? (double)it / rs : 0.0,
            (unsigned long)g_km_capped.load(),
            rs ? 100.0 * g_km_capped.load() / rs : 0.0);
  }
  if (prof_bins) {
    std::vector<size_t> byt(nbins);
    for (size_t i = 0; i < nbins; ++i) byt[i] = i;
    std::sort(byt.begin(), byt.end(),
              [&](size_t a, size_t b) { return bin_ms[a] > bin_ms[b]; });
    double tot = 0.0;
    for (double v : bin_ms) tot += v;
    fprintf(stderr, "[RB_SPLIT_PROF] %zu bins, total work %.1f ms\n", nbins, tot);
    for (size_t i = 0; i < std::min<size_t>(8, nbins); ++i)
      fprintf(stderr, "[RB_SPLIT_PROF]   #%zu contigs=%zu  %.1f ms\n", byt[i],
              binList[byt[i]]->size(), bin_ms[byt[i]]);
  }

  // ── Self-calibrated coherence bar over confident (kept/split) bins ────────
  // bar = chosen percentile of the intra-bin coherence of the bins that passed
  // the size floor.  Data-derived (no per-dataset constant); identical formula
  // across datasets.  Used only to admit purified sub-floor pieces below.
  double bar = 1.0;
  if (coh_gate) {
    std::vector<double> conf;
    for (size_t bi = 0; bi < nbins; ++bi) {
      if (results[bi].kind != 1 && results[bi].kind != 2) continue;
      for (auto &v : results[bi].emit) {
        double c = bin_coherence(v);
        if (c > -1.5) conf.push_back(c);
      }
    }
    if (!conf.empty()) {
      std::sort(conf.begin(), conf.end());
      size_t idx = (size_t)((g_split_coh_pct / 100.0) * (conf.size() - 1) + 0.5);
      if (idx >= conf.size()) idx = conf.size() - 1;
      bar = conf[idx];
    }
  }

  // ── Purify + retain captured sub-floor bins (additive; main path untouched) ─
  // Each captured <min_bin_bp bin is independently KMeans-split (separating any
  // co-binned small genomes), then each resulting piece (or the whole bin) is
  // kept iff it is internally at least as coherent as the bar.  Because these
  // bins were going to be discarded, every retained piece is additive.
  std::vector<size_t> captured;
  for (size_t bi = 0; bi < nbins; ++bi)
    if (results[bi].kind == 0) captured.push_back(bi);

  std::vector<std::vector<ContigVector>> recov(captured.size());
  if (coh_gate) {
    // Same largest-first scheduling fix as the main path above: `co` only
    // reorders which `ci` a thread grabs next, `ci` itself (and thus the
    // `captured[ci]` seed and `recov[ci]` output slot) is unchanged, so this
    // cannot alter any bin's split/coherence result — perf only.
    std::vector<size_t> co(captured.size());
    for (size_t i = 0; i < captured.size(); ++i) co[i] = i;
    std::sort(co.begin(), co.end(), [&](size_t a, size_t b) {
      return binList[captured[a]]->size() > binList[captured[b]]->size();
    });
#pragma omp parallel for schedule(dynamic, 1) num_threads(numThreads)
    for (size_t oi = 0; oi < captured.size(); ++oi) {
      const size_t ci = co[oi];
      const ContigVector &contigs = *binList[captured[ci]];
      std::vector<ContigVector> pieces, sub;
      if (split_bin(contigs, captured[ci], sub)) {
        for (auto &s : sub) if (!s.empty()) pieces.push_back(std::move(s));
      } else {
        pieces.push_back(contigs);
      }
      for (auto &p : pieces)
        if (bin_coherence(p) >= bar) recov[ci].push_back(std::move(p));
    }
  }

  // ── Sub-floor fragment consolidation (additive; gated; deterministic) ─────
  // Flatten the retained >=min_bin_bp cores and the recovered sub-floor pieces
  // into stable, index-ordered lists, then (optionally) attach each fragment to
  // the best-matching core by mean cross-bin depth correlation.  A fragment is
  // merged iff that correlation >= bar + delta (same self-calibrated bar as the
  // retention gate) — i.e. it is as depth-coherent with the core as the core is
  // with itself, which only same-genome tails achieve.  Unmerged fragments are
  // emitted as their own bins exactly as before, so the unfiltered behaviour and
  // every dataset where cores are rare stay unchanged.
  std::vector<ContigVector> cores;          // retained >=min_bin_bp bins
  size_t n_split = 0, n_kept = 0, n_drop = 0, n_recovered = 0, n_merged = 0;
  for (size_t bi = 0; bi < nbins; ++bi) {
    BinResult &res = results[bi];
    if (res.kind == 0) continue;                     // captured: handled below
    if (res.emit.empty()) { ++n_drop; continue; }    // legacy sub-floor drop
    if (res.kind == 2) ++n_split; else ++n_kept;
    for (auto &v : res.emit) cores.push_back(std::move(v));
  }
  std::vector<ContigVector> frags;                   // recovered sub-floor pieces
  for (size_t ci = 0; ci < captured.size(); ++ci) {
    if (recov[ci].empty()) { ++n_drop; continue; }
    for (auto &v : recov[ci]) frags.push_back(std::move(v));
  }
  n_recovered = frags.size();

  const bool do_consolidate = g_split_consolidate && coh_gate && !cores.empty();
  std::vector<long> target(frags.size(), -1);   // frag -> core index, or -1 (keep)
  if (do_consolidate) {
    // Anchor = up to CAP large members (c < nobs) per core / fragment, used to
    // estimate the mean cross correlation cheaply via depth_corr_fast.
    constexpr size_t CAP_C = 24, CAP_F = 8;
    auto anchors = [&](const ContigVector &cs, size_t cap) {
      std::vector<size_t> L;
      for (size_t c : cs) if (c < nobs) L.push_back(c);
      if (L.size() > cap) {
        std::vector<size_t> s; s.reserve(cap);
        const double step = (double)L.size() / (double)cap;
        for (size_t t = 0; t < cap; ++t) s.push_back(L[(size_t)(t * step)]);
        L.swap(s);
      }
      return L;
    };
    std::vector<std::vector<size_t>> coreAnchor(cores.size());
    std::vector<double> coreThr(cores.size());
    for (size_t k = 0; k < cores.size(); ++k) {
      coreAnchor[k] = anchors(cores[k], CAP_C);
      // Per-core adaptive gate: a fragment must be at least as depth-coherent
      // with this core as the core is INTERNALLY (own mean pairwise corr), so a
      // tight single-genome core only absorbs equally tight same-genome tails.
      double cc = bin_coherence(cores[k]);
      coreThr[k] = ((cc > -1.5) ? cc : bar) + g_split_consolidate_delta;
    }
    // For each fragment, find the core with the highest mean cross correlation
    // (read-only over cores → parallel-safe; the actual append is sequential).
#pragma omp parallel for schedule(dynamic, 16) num_threads(numThreads)
    for (size_t f = 0; f < frags.size(); ++f) {
      std::vector<size_t> fa = anchors(frags[f], CAP_F);
      if (fa.size() < 1) continue;                   // no large anchor → keep separate
      double best = -2.0; long bestk = -1;
      for (size_t k = 0; k < cores.size(); ++k) {
        const auto &ca = coreAnchor[k];
        if (ca.size() < 2) continue;
        double sum = 0.0; size_t cnt = 0;
        for (size_t a : fa)
          for (size_t b : ca) {
            double c = depth_corr_fast(a, b);
            if (std::isfinite(c)) { sum += c; ++cnt; }
          }
        if (!cnt) continue;
        double m = sum / (double)cnt;
        // gate against THIS core's own coherence (margin baked into coreThr)
        if (m >= coreThr[k] && m > best) { best = m; bestk = (long)k; }
      }
      if (bestk >= 0) target[f] = bestk;
    }
    // Sequential, deterministic append.
    for (size_t f = 0; f < frags.size(); ++f) {
      if (target[f] >= 0) {
        ContigVector &dst = cores[(size_t)target[f]];
        dst.insert(dst.end(), frags[f].begin(), frags[f].end());
        ++n_merged;
      }
    }
  }

  BinMap out;
  int next = 0;
  for (auto &v : cores) out[next++] = std::move(v);
  for (size_t f = 0; f < frags.size(); ++f) {
    if (target[f] >= 0) continue;                      // already merged into a core
    out[next++] = std::move(frags[f]);
  }
  cls.swap(out);
  verbose_message("Abundance split (sil>=%.2f, bar=%.3f): %zu kept, %zu split, "
                  "%zu recovered(purified<floor), %zu merged->core, %zu dropped -> %d bins\n",
                  g_split_sil, bar, n_kept, n_split, n_recovered, n_merged, n_drop, next);
  // Retained pieces may be < min_bin_bp; the floor was already applied here.
  min_bin_bp = 0;
}

// ═══════════════════════════════════════════════════════════════════════════
// Simplified marker-free splitter. Coverage features are log1p transformed. For
// every K, K-means uses several deterministic initialisations; only partitions
// whose every child clears the configured size guards compete by silhouette.
// Rejecting an invalid proposal leaves the parent intact, so an outlier child is
// never silently discarded.
static void abundance_guided_split_simple(BinMap &cls) {
  if (num_depth_samples < 1) {
    verbose_message("Abundance split: needs >=1 abundance sample — skipped\n");
    return;
  }
  if (cls.empty()) return;

  const size_t S = (size_t)num_depth_samples;
  auto depth_at = [&](size_t c, size_t k) -> double {
    if (c < nobs) {
      const size_t o = c * S + k;
      if (o < g_large_means.size()) return (double)g_large_means[o];
      return (double)depth_matrix(c, k);
    }
    const size_t s = c - nobs;
    const size_t o = s * S + k;
    if (o < g_small_means.size()) return (double)g_small_means[o];
    return (double)small_depth_matrix(s, k);
  };
  auto child_bp = [&](const ContigVector &cs) -> size_t {
    size_t bp = 0;
    for (size_t c : cs)
      bp += (c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs];
    return bp;
  };

  struct InputBin { int id; const ContigVector *members; };
  std::vector<InputBin> bins;
  bins.reserve(cls.size());
  int next_id = 0;
  for (auto &kv : cls) {
    bins.push_back({kv.first, &kv.second});
    next_id = std::max(next_id, kv.first + 1);
  }
  std::sort(bins.begin(), bins.end(),
            [](const InputBin &a, const InputBin &b) { return a.id < b.id; });

  struct SplitResult {
    std::vector<ContigVector> children;
    size_t rejected_small = 0;
  };
  std::vector<SplitResult> results(bins.size());

#pragma omp parallel for schedule(dynamic, 1) num_threads(numThreads)
  for (size_t bi = 0; bi < bins.size(); ++bi) {
    const ContigVector &members = *bins[bi].members;
    const size_t n = members.size();
    if ((int)n < splitMinContigs) continue;

    int kmax = std::min(splitMaxK, (int)n - 1);
    if (splitMinSubContigs > 1)
      kmax = std::min(kmax, (int)(n / splitMinSubContigs));
    if (kmax < 2) continue;

    std::vector<float> X(n * S);
    for (size_t r = 0; r < n; ++r)
      for (size_t k = 0; k < S; ++k)
        X[r * S + k] = (float)std::log1p(std::max(0.0, depth_at(members[r], k)));

    const uint64_t base_seed = seed ? (uint64_t)seed : 1ULL;
    std::mt19937 rng((uint32_t)(base_seed ^
        ((uint64_t)(uint32_t)bins[bi].id + 0x9e3779b9ULL) * 2654435761ULL));
    double best_sil = -1.0;
    std::vector<ContigVector> best_children;

    for (int k = 2; k <= kmax; ++k) {
      std::vector<int> labels = rb_kmeans(X.data(), n, S, k, rng);
      std::vector<ContigVector> children((size_t)k);
      for (size_t r = 0; r < n; ++r) {
        const int lab = labels[r];
        if (lab >= 0 && lab < k) children[(size_t)lab].push_back(members[r]);
      }

      bool valid = true;
      for (const auto &child : children) {
        if (child.size() < splitMinSubContigs ||
            (splitMinSubBp > 0 && child_bp(child) < splitMinSubBp)) {
          valid = false;
          break;
        }
      }
      if (!valid) {
        ++results[bi].rejected_small;
        continue;
      }

      const double sil = fkmv_silhouette(X.data(), n, S, labels, k, rng);
      if (sil > best_sil) {
        best_sil = sil;
        best_children = std::move(children);
      }
    }

    if (best_sil >= g_split_sil)
      results[bi].children = std::move(best_children);
  }

  BinMap out;
  out.reserve(cls.size() * 2);
  size_t n_split = 0, rejected_small = 0;
  for (size_t bi = 0; bi < bins.size(); ++bi) {
    SplitResult &res = results[bi];
    rejected_small += res.rejected_small;
    if (res.children.size() < 2) {
      out[bins[bi].id] = *bins[bi].members;
      continue;
    }
    ++n_split;
    out[bins[bi].id] = std::move(res.children[0]);
    for (size_t i = 1; i < res.children.size(); ++i)
      out[next_id++] = std::move(res.children[i]);
  }
  cls.swap(out);
  verbose_message("Abundance split (log1p, K=2..%d, restarts=%d, sil>=%.2f, "
                  "child>=%zu contigs, child>=%zu bp): %zu split, %zu kept, "
                  "%zu invalid-K proposals -> %zu bins\n",
                  splitMaxK, splitKmeansRestarts, g_split_sil,
                  splitMinSubContigs, splitMinSubBp,
                  n_split, bins.size() - n_split, rejected_small, cls.size());
}

// A frozen core stores exactly the post-split members used to derive its
// recruitment threshold and centroid.  Recruits are never folded back into it.
struct FrozenRecruitCore {
  int id = -1;
  double threshold = 1.0;
  std::vector<double> centroid;
};

static int frozen_core_unique_match(const float *unit, size_t S,
                                    const std::vector<FrozenRecruitCore> &cores,
                                    bool &ambiguous) {
  ambiguous = false;
  double norm2 = 0.0;
  for (size_t k = 0; k < S; ++k) norm2 += (double)unit[k] * unit[k];
  if (!(norm2 > 0.0)) return -1; // Spearman is undefined for a constant profile

  int match = -1;
  for (const FrozenRecruitCore &core : cores) {
    double score = 0.0;
    for (size_t k = 0; k < S; ++k) score += (double)unit[k] * core.centroid[k];
    if (score >= core.threshold) {
      if (match >= 0) {
        ambiguous = true;
        return -1;                 // exactly one passing bin is required
      }
      match = core.id;
    }
  }
  return match;
}

// Recruit in two batches: unassigned long contigs first, then retained short
// contigs.  Both batches compare independently with the same immutable cores.
// For unit-rank vectors u, the formulas below are exactly the stated mean
// Spearman scores, evaluated in O(nS) instead of materialising every core pair:
//
//   tau_b = (||sum_i u_i||^2 - sum_i ||u_i||^2) / (n(n-1))
//   r_cb   = u_c . (sum_i u_i / n)
static void recruit_to_frozen_cores(BinMap &cls, size_t core_floor) {
  const size_t S = (size_t)num_depth_samples;
  if (no_recruit || S < 2 || g_depth_unit.empty() || cls.empty()) return;

  std::vector<int> ids;
  ids.reserve(cls.size());
  for (const auto &kv : cls) ids.push_back(kv.first);
  std::sort(ids.begin(), ids.end());

  std::vector<FrozenRecruitCore> cores;
  cores.reserve(ids.size());
  for (int id : ids) {
    const ContigVector &bin = cls.find(id)->second;
    size_t span = 0;
    for (size_t c : bin)
      span += (c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs];
    // A core must already satisfy the requested output floor.  Otherwise a
    // tiny low-confidence cluster can cross the floor solely by recruitment,
    // creating a new bin whose apparent completeness comes from the recruits.
    if (span < core_floor) continue;
    std::vector<double> sum(S, 0.0);
    double self_sum = 0.0;
    size_t valid = 0;
    for (size_t c : bin) {
      if (c >= nobs) continue;
      const float *u = g_depth_unit.data() + c * S;
      double norm2 = 0.0;
      for (size_t k = 0; k < S; ++k) norm2 += (double)u[k] * u[k];
      if (!(norm2 > 0.0)) continue;
      for (size_t k = 0; k < S; ++k) sum[k] += u[k];
      self_sum += norm2;
      ++valid;
    }
    if (valid < std::max<size_t>(2, minCS)) continue;

    double sum_norm2 = 0.0;
    for (double v : sum) sum_norm2 += v * v;
    double tau = (sum_norm2 - self_sum) /
                 ((double)valid * (double)(valid - 1));
    tau = std::max(-1.0, std::min(1.0, tau));
    for (double &v : sum) v /= (double)valid;
    FrozenRecruitCore core;
    core.id = id;
    core.threshold = tau;
    core.centroid = std::move(sum);
    cores.push_back(std::move(core));
  }
  if (cores.empty()) {
    verbose_message("Frozen-core recruit: no bins have >=%zu valid coverage profiles\n",
                    minCS);
    return;
  }

  // This mask and the cores are frozen before either batch.  Committing the long
  // batch therefore cannot affect any later score or threshold.
  std::vector<char> initially_binned(nobs, 0);
  for (const auto &kv : cls)
    for (size_t c : kv.second) if (c < nobs) initially_binned[c] = 1;

  std::vector<int> long_assignment(nobs, -1);
  size_t long_assigned = 0, long_ambiguous = 0, long_unmatched = 0;
#pragma omp parallel for schedule(dynamic, 128) num_threads(numThreads) \
    reduction(+:long_assigned,long_ambiguous,long_unmatched)
  for (size_t c = 0; c < nobs; ++c) {
    if (initially_binned[c]) continue;
    bool ambiguous = false;
    const int id = frozen_core_unique_match(g_depth_unit.data() + c * S, S,
                                            cores, ambiguous);
    long_assignment[c] = id;
    if (id >= 0) ++long_assigned;
    else if (ambiguous) ++long_ambiguous;
    else ++long_unmatched;
  }
  for (size_t c = 0; c < nobs; ++c)
    if (long_assignment[c] >= 0) cls[long_assignment[c]].push_back(c);
  verbose_message("Frozen-core recruit (mean Spearman, unique-bin): long %zu assigned, "
                  "%zu ambiguous, %zu unmatched against %zu cores\n",
                  long_assigned, long_ambiguous, long_unmatched, cores.size());

  std::vector<int> short_assignment(nobs1, -1);
  size_t short_assigned = 0, short_ambiguous = 0, short_unmatched = 0;
  unsigned long long short_bp = 0;
#pragma omp parallel num_threads(numThreads) \
    reduction(+:short_assigned,short_ambiguous,short_unmatched)
  {
    std::vector<double> raw(S), ranked(S);
    std::vector<float> unit(S);
#pragma omp for schedule(dynamic, 128)
    for (size_t s = 0; s < nobs1; ++s) {
      for (size_t k = 0; k < S; ++k) raw[k] = (double)small_depth_matrix(s, k);
      rank(raw, ranked);
      double mean = 0.0;
      for (double v : ranked) mean += v;
      mean /= (double)S;
      double ss = 0.0;
      for (double v : ranked) { const double d = v - mean; ss += d * d; }
      if (!(ss > 0.0)) {
        ++short_unmatched;
        continue;
      }
      const double inv = 1.0 / std::sqrt(ss);
      for (size_t k = 0; k < S; ++k)
        unit[k] = (float)((ranked[k] - mean) * inv);
      bool ambiguous = false;
      const int id = frozen_core_unique_match(unit.data(), S, cores, ambiguous);
      short_assignment[s] = id;
      if (id >= 0) ++short_assigned;
      else if (ambiguous) ++short_ambiguous;
      else ++short_unmatched;
    }
  }

  for (size_t s = 0; s < nobs1; ++s) {
    if (short_assignment[s] < 0) continue;
    cls[short_assignment[s]].push_back(s + nobs);
    short_bp += small_seq_lens[s];
  }
  verbose_message("Frozen-core recruit (same immutable cores): short %zu assigned "
                  "(%llu bp), %zu ambiguous, %zu unmatched\n",
                  short_assigned, short_bp, short_ambiguous, short_unmatched);
}

// Shared abundance statistics — one model for merge / recruit
// ═══════════════════════════════════════════════════════════════════════════
// Every post-clustering stage answers the same question: are these two sets of
// contigs samples of ONE genome's coverage profile?  A fixed correlation cutoff
// answers it badly, because the correlation a set attains depends on how many
// contigs it holds and how long they are — coverage estimated over 1 kb is far
// noisier than over 100 kb, so a true member is penalised for being short, and a
// large set is penalised for averaging over more noise.  Model the noise instead.
//
// For contig i of length L_i in sample k let y_ik = log1p(depth_ik).  Within one
// genome
//        y_ik ~ N( mu_k , sigma_k^2(L_i) ),   sigma_k^2(L) = a_k + b_k / L
// a length-independent term (strain heterogeneity, GC / mappability bias) plus a
// sampling term decaying as 1/L, because a coverage mean over L positions has
// variance proportional to 1/L.  (a_k, b_k) are estimated from the within-bin
// scatter of the current partition, so the noise scale is calibrated to the data
// rather than assumed.
//
// With inverse-variance weights w_ik = 1/sigma_k^2(L_i) a set A has weighted mean
// mu_Ak = (sum w_ik y_ik)/W_Ak and Var(mu_Ak) = 1/W_Ak, W_Ak = sum w_ik.  Under
// the hypothesis that A and B come from the same genome
//        T2(A,B) = sum_k (mu_Ak - mu_Bk)^2 / (1/W_Ak + 1/W_Bk)   ~   chi2_S
// A single contig tested against a bin is the same expression with A = {i}, so
// merge and statistical recruitment both reduce to one statistic against one
// chi-square critical value.  The bar follows from the sample count S and a significance
// level; short contigs receive exactly the wider tolerance the model implies.

// Raw per-sample mean depth of a contig (large or small).  depth_matrix and
// small_depth_matrix are rank-transformed in place during the pipeline, so read
// the raw-mean snapshots taken before those transforms.
static inline double rb_depth_at(size_t c, size_t k) {
  const size_t S = (size_t)num_depth_samples;
  if (c < nobs) {
    if (!g_large_means.empty() && c * S + k < g_large_means.size())
      return (double)g_large_means[c * S + k];
    return (double)depth_matrix(c, k);
  }
  const size_t s = c - nobs;
  if (!g_small_means.empty() && s * S + k < g_small_means.size())
    return (double)g_small_means[s * S + k];
  return (double)small_depth_matrix(s, k);
}

static inline double rb_contig_len(size_t c) {
  return (double)((c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs]);
}

// sigma_k^2(L) = a[k] + b[k]/L, fitted once per run.
struct RbAbdVar {
  std::vector<double> a, b;
  bool ok = false;
};
static RbAbdVar g_abd_var;

static inline double rb_abd_w(size_t c, size_t k) {
  const double v = g_abd_var.a[k] + g_abd_var.b[k] / std::max(rb_contig_len(c), 1.0);
  return v > 1e-12 ? 1.0 / v : 0.0;
}
static inline double rb_abd_y(size_t c, size_t k) {
  return std::log1p(std::max(rb_depth_at(c, k), 0.0));
}

// Estimate (a_k, b_k) from the within-bin scatter of log-depth against 1/L.
// Bins of the current partition act as (imperfect) samples of single genomes;
// a single trimming pass removes the contribution of members that sit far
// outside the fitted band, so residual cross-genome contamination in the input
// partition does not inflate the noise scale.
static bool rb_fit_abd_var(const BinMap &cls) {
  const size_t S = (size_t)num_depth_samples;
  g_abd_var.ok = false;
  if (S < 1) return false;

  std::vector<const ContigVector *> use;
  for (auto &kv : cls) {
    size_t n = 0;
    for (size_t c : kv.second) if (c < nobs && ++n >= 8) break;
    if (n >= 8) use.push_back(&kv.second);
  }
  if (use.empty()) return false;

  constexpr size_t CAP = 128;         // members sampled per bin
  g_abd_var.a.assign(S, 0.0);
  g_abd_var.b.assign(S, 0.0);

  // Two passes: pass 0 fits on all sampled residuals, pass 1 refits after
  // trimming residuals above 3 sigma of the pass-0 band.
  for (int pass = 0; pass < 2; ++pass) {
    const int nthr = (int)std::max<size_t>(1, numThreads);
    // Per-thread accumulators of the 2x2 normal equations, per sample.
    std::vector<std::vector<double>> acc(nthr, std::vector<double>(S * 5, 0.0));
#pragma omp parallel num_threads(numThreads)
    {
      int tid = omp_get_thread_num();
      if (tid >= nthr) tid = nthr - 1;
      double *A = acc[tid].data();
      std::vector<size_t> mem;
      std::vector<double> mu(S);
#pragma omp for schedule(dynamic, 16)
      for (size_t u = 0; u < use.size(); ++u) {
        mem.clear();
        for (size_t c : *use[u]) if (c < nobs) mem.push_back(c);
        if (mem.size() < 8) continue;
        if (mem.size() > CAP) {
          const double step = (double)mem.size() / (double)CAP;
          std::vector<size_t> s(CAP);
          for (size_t t = 0; t < CAP; ++t) s[t] = mem[(size_t)(t * step)];
          mem.swap(s);
        }
        const double n = (double)mem.size();
        // Unbiased within-bin scatter: E[(y-mean)^2] = sigma^2 (n-1)/n.
        const double corr = n / (n - 1.0);
        for (size_t k = 0; k < S; ++k) {
          double s = 0.0;
          for (size_t c : mem) s += rb_abd_y(c, k);
          mu[k] = s / n;
        }
        for (size_t c : mem) {
          const double x = 1.0 / std::max(rb_contig_len(c), 1.0);
          for (size_t k = 0; k < S; ++k) {
            const double d = rb_abd_y(c, k) - mu[k];
            const double d2 = d * d * corr;
            if (pass == 1) {
              const double band = g_abd_var.a[k] + g_abd_var.b[k] * x;
              if (band > 0.0 && d2 > 9.0 * band) continue;   // >3 sigma → trim
            }
            double *Ak = A + k * 5;
            Ak[0] += 1.0; Ak[1] += x; Ak[2] += x * x; Ak[3] += d2; Ak[4] += d2 * x;
          }
        }
      }
    }
    for (size_t k = 0; k < S; ++k) {
      double n = 0, Sx = 0, Sxx = 0, Sd = 0, Sdx = 0;
      for (int t = 0; t < nthr; ++t) {
        const double *Ak = acc[t].data() + k * 5;
        n += Ak[0]; Sx += Ak[1]; Sxx += Ak[2]; Sd += Ak[3]; Sdx += Ak[4];
      }
      if (n < 4.0) { g_abd_var.a[k] = 1e-4; g_abd_var.b[k] = 0.0; continue; }
      const double det = n * Sxx - Sx * Sx;
      double a = 0.0, b = 0.0;
      if (std::fabs(det) > 1e-30) {
        a = (Sd * Sxx - Sdx * Sx) / det;
        b = (n * Sdx - Sx * Sd) / det;
      }
      if (!(a > 0.0) || !std::isfinite(a) || !std::isfinite(b)) { a = Sd / n; b = 0.0; }
      if (b < 0.0) b = 0.0;
      g_abd_var.a[k] = std::max(a, 1e-6);
      g_abd_var.b[k] = b;
    }
  }
  g_abd_var.ok = true;
  return true;
}

// Sufficient statistics of a contig set: W[k] = sum of weights, Z[k] = sum w*y.
static inline void rb_suff_add(const ContigVector &cs, double *W, double *Z, size_t S) {
  for (size_t c : cs) {
    const double L = std::max(rb_contig_len(c), 1.0);
    for (size_t k = 0; k < S; ++k) {
      const double v = g_abd_var.a[k] + g_abd_var.b[k] / L;
      if (!(v > 1e-12)) continue;
      const double w = 1.0 / v;
      W[k] += w;
      Z[k] += w * rb_abd_y(c, k);
    }
  }
}

static inline double rb_t2(const double *WA, const double *ZA,
                           const double *WB, const double *ZB, size_t S) {
  double t2 = 0.0;
  for (size_t k = 0; k < S; ++k) {
    if (!(WA[k] > 0.0) || !(WB[k] > 0.0)) continue;
    const double d = ZA[k] / WA[k] - ZB[k] / WB[k];
    t2 += d * d / (1.0 / WA[k] + 1.0 / WB[k]);
  }
  return t2;
}

// chi2_S critical value at family-wise level alpha over n_tests comparisons
// (Bonferroni).  alpha is a significance level, not a per-dataset constant: the
// resulting bar adapts to the sample count and to how many comparisons the
// stage actually performs.
static double rb_chi2_crit(size_t S, double alpha, size_t n_tests) {
  if (S < 1) return 0.0;
  double q = alpha / (double)std::max<size_t>(n_tests, 1);
  if (!(q > 0.0)) q = 1e-300;
  if (q >= 1.0) q = std::nextafter(1.0, 0.0);
  boost::math::chi_squared_distribution<double> d((double)S);
  return boost::math::quantile(boost::math::complement(d, q));
}

// Family-wise significance level shared by statistical merge / recruitment.
static double g_abd_alpha = 0.05;

// ═══════════════════════════════════════════════════════════════════════════
// Paired-end bin merge — rejoin fragments of ONE genome
// ═══════════════════════════════════════════════════════════════════════════
// LPA / abundance split can scatter a single genome across several pure bins
// (coverage sub-clusters, or separate LPA communities).  Each fragment is then
// individually < min_bin_bp, so the academic ">=200kb bin" protocol discards the
// genome even though the UNION of its fragments is a high-quality MAG.  This pass
// agglomerates such bins using paired-end read linkage (mates aligned to two
// different contigs) as a physical same-genome signal orthogonal to
// composition and abundance.
// A cheap depth-centroid correlation is an additional abundance sanity check.
// Safety: a merge requires PE evidence; never fuses two independent >=floor
// cores; with fragment-only mode it cannot touch any >=floor bin at all.
static void consolidate_bins(BinMap &cls, size_t floor) {
  if (!g_bin_merge) return;
  const size_t S = (size_t)num_depth_samples;
  // The abundance test needs raw per-sample depths; depth_matrix is rank-
  // transformed in place, so the raw-mean snapshot must be present.
  if (S < 1 || g_large_means.empty()) return;

  std::vector<int> ids; std::vector<ContigVector*> bp;
  ids.reserve(cls.size()); bp.reserve(cls.size());
  for (auto &kv : cls) { ids.push_back(kv.first); bp.push_back(&kv.second); }
  const size_t B = bp.size();
  if (B < 2) return;

  // GAIN-ONLY guard: only sub-floor fragments may participate in a merge. A bin
  // already above the floor is not modified, so size-filtered recovery is
  // monotone: merging can only create new above-floor bins by
  // uniting a genome's scattered fragments, never break an existing good bin.
  auto bin_bp = [&](size_t b) -> size_t {
    size_t s = 0;
    for (size_t c : *bp[b]) s += (c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs];
    return s;
  };
  // Fragment-only mode (RABBIT_BIN_MERGE_FRAGONLY=1) forbids touching any
  // >=floor bin (strictly monotone, can only build new bins from fragments).
  // When fragment-only mode is disabled, a fragment may also be merged into an
  // above-floor core. Paired-end evidence protects against cross-genome merges.
  // A merge still requires at least one participant to be a
  // fragment, so two independent >=floor cores are never fused.
  const bool frag_only = [] { const char *e = getenv("RABBIT_BIN_MERGE_FRAGONLY");
                              return e && e[0] != '0'; }();
  std::vector<char> is_frag(B, 0);
  for (size_t b = 0; b < B; ++b) is_frag[b] = (bin_bp(b) < floor) ? 1 : 0;

  // ── Paired-end linkage support per candidate bin-pair ─────────────────────
  // Convert the captured links (compact depth-row space) to contig indices, then
  // aggregate read-pair counts between distinct bins.  A merge later REQUIRES
  // pe_support >= g_bin_merge_pe_min so only physically-linked (≈same-genome)
  // bins are joined — the high-precision confirmation that TNF+depth alone lack.
  std::vector<int> binpos(nobs, -1);
  for (size_t b = 0; b < B; ++b)
    for (size_t c : *bp[b]) if (c < nobs) binpos[c] = (int)b;
  phmap::flat_hash_map<uint64_t, uint32_t> pe_pair;   // pack(min,max bin pos) -> count
  if (!g_pe_links_compact.empty() && !g_pe_names.empty()) {
    // compact row -> contig index (binning space), via name.
    phmap::flat_hash_map<std::string, int> name2idx;
    name2idx.reserve(nobs * 2);
    for (size_t i = 0; i < nobs; ++i) name2idx[contig_names[i]] = (int)i;
    auto compact2idx = [&](int32_t cr) -> int {
      if (cr < 0 || (size_t)cr >= g_pe_names.size()) return -1;
      std::string nm = g_pe_names[cr];
      size_t sp = nm.find_first_of(" \t"); if (sp != std::string::npos) nm.resize(sp);
      auto it = name2idx.find(nm);
      return it == name2idx.end() ? -1 : it->second;
    };
    for (auto &lk : g_pe_links_compact) {
      int ia = compact2idx(std::get<0>(lk));
      int ib = compact2idx(std::get<1>(lk));
      if (ia < 0 || ib < 0) continue;
      int ba = binpos[ia], bb = binpos[ib];
      if (ba < 0 || bb < 0 || ba == bb) continue;
      uint32_t lo = (uint32_t)std::min(ba, bb), hi = (uint32_t)std::max(ba, bb);
      pe_pair[((uint64_t)lo << 32) | hi] += std::get<2>(lk);
    }
  }
  uint32_t pe_min = 1;
  if (const char *e = getenv("RABBIT_BIN_MERGE_PE_MIN")) pe_min = (uint32_t)atoi(e);
  if (pe_pair.empty()) return;   // no physical linkage evidence → merge nothing

  // ── Abundance model + per-bin sufficient statistics ───────────────────────
  if (!rb_fit_abd_var(cls)) return;
  std::vector<double> Wb((size_t)B * S, 0.0), Zb((size_t)B * S, 0.0);
#pragma omp parallel for schedule(dynamic, 16) num_threads(numThreads)
  for (size_t b = 0; b < B; ++b)
    rb_suff_add(*bp[b], Wb.data() + (size_t)b * S, Zb.data() + (size_t)b * S, S);

  // The candidate set IS the set of paired-end-linked bin pairs: PE linkage is
  // the physical same-genome evidence a merge requires, so enumerating the
  // observed links replaces the quadratic all-pairs scan.
  double alpha = g_abd_alpha;
  if (const char *e = getenv("RABBIT_BIN_MERGE_ALPHA")) alpha = atof(e);
  const double crit = rb_chi2_crit(S, alpha, pe_pair.size());

  struct MPair { double t2; uint32_t a, b; };
  std::vector<MPair> pairs;
  pairs.reserve(pe_pair.size());
  for (auto &kv : pe_pair) {
    if (kv.second < pe_min) continue;
    const uint32_t a  = (uint32_t)(kv.first >> 32);
    const uint32_t b2 = (uint32_t)(kv.first & 0xffffffffu);
    if ((size_t)a >= B || (size_t)b2 >= B) continue;
    // At least one participant must be a fragment (never fuse two independent
    // cores); fragment-only mode additionally forbids touching any core.
    if (frag_only) { if (!is_frag[a] || !is_frag[b2]) continue; }
    else           { if (!is_frag[a] && !is_frag[b2]) continue; }
    const double t2 = rb_t2(Wb.data() + (size_t)a  * S, Zb.data() + (size_t)a  * S,
                            Wb.data() + (size_t)b2 * S, Zb.data() + (size_t)b2 * S, S);
    if (t2 > crit) continue;               // coverage profiles differ → distinct
    pairs.push_back({t2, a, b2});
  }
  if (pairs.empty()) return;
  // Strongest evidence first, so a genome's own fragments unite before a
  // marginally compatible neighbour gets its chance.
  std::sort(pairs.begin(), pairs.end(),
            [](const MPair &x, const MPair &y) { return x.t2 < y.t2; });

  std::vector<size_t> uf(B);
  std::vector<char> comp_core(B);            // does this component contain a >=floor core?
  for (size_t i = 0; i < B; ++i) { uf[i] = i; comp_core[i] = is_frag[i] ? 0 : 1; }
  auto find = [&](size_t x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; };
  size_t merges = 0, n_drift = 0;
  for (auto &p : pairs) {
    const size_t ra = find(p.a), rb = find(p.b);
    if (ra == rb) continue;
    // Never fuse two components that each already contain a >=floor core: that
    // would join two independent (potentially HQ) MAGs.  Fragments may still flow
    // into a single core, or unite into a new core.
    if (comp_core[ra] && comp_core[rb]) continue;
    // Re-test the CURRENT components rather than the original bins.  Pairwise
    // admissibility is not transitive: without this a core accumulates drift one
    // individually-plausible fragment at a time and ends up spanning two genomes.
    double *WA = Wb.data() + ra * S, *ZA = Zb.data() + ra * S;
    double *WB = Wb.data() + rb * S, *ZB = Zb.data() + rb * S;
    if (rb_t2(WA, ZA, WB, ZB, S) > crit) { ++n_drift; continue; }
    const size_t root = std::min(ra, rb), other = std::max(ra, rb);
    double *Wr = Wb.data() + root * S, *Zr = Zb.data() + root * S;
    const double *Wo = Wb.data() + other * S, *Zo = Zb.data() + other * S;
    for (size_t k = 0; k < S; ++k) { Wr[k] += Wo[k]; Zr[k] += Zo[k]; }
    uf[other] = root;
    comp_core[root] = comp_core[ra] || comp_core[rb];
    ++merges;
  }
  if (!merges) return;

  BinMap out;
  for (size_t b = 0; b < B; ++b) {
    size_t r = find(b);
    ContigVector &dst = out[ids[r]];
    dst.insert(dst.end(), bp[b]->begin(), bp[b]->end());
  }
  const size_t before = cls.size();
  cls.swap(out);
  verbose_message("Bin merge (PE pe>=%u, chi2_%zu crit=%.1f alpha=%.3g over %zu "
                  "candidates): %zu -> %zu bins (%zu merges, %zu rejected on "
                  "accumulated drift)\n", pe_min, S, crit, alpha, pe_pair.size(),
                  before, cls.size(), merges, n_drift);
}

// ═══════════════════════════════════════════════════════════════════════════
// Composition+depth recruit — attach UNBINNED contigs to the core they belong to
// ═══════════════════════════════════════════════════════════════════════════
// For each unbinned large contig, attach
// it to the >=floor core it matches on BOTH tetranucleotide composition (cosine)
// AND abundance (depth-corr to the core centroid) — a CONJUNCTIVE gate (rejects a
// contaminant that matches only one signal) plus a margin over the 2nd-best core
// (unambiguous).  Additive but purity-safe: a contig only joins a core when it
// looks like that core's own genome on both independent signals.
static void recruit_unbinned_to_cores(BinMap &cls, size_t floor) {
  if (!g_bin_recruit) return;
  const size_t S = (size_t)num_depth_samples;
  if (S < 1 || g_depth_unit.empty() || g_merge_tnf.empty()) return;

  std::vector<int> ids; std::vector<ContigVector*> bp;
  for (auto &kv : cls) { ids.push_back(kv.first); bp.push_back(&kv.second); }
  const size_t B = bp.size();
  auto bin_bp = [&](const ContigVector &cs) {
    size_t s = 0; for (size_t c : cs) s += (c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs];
    return s;
  };

  // Cores = >=floor bins.  Mark binned contigs.
  std::vector<size_t> coreIdx;
  std::vector<char> binned(nobs, 0);
  for (size_t b = 0; b < B; ++b) {
    if (bin_bp(*bp[b]) >= floor) coreIdx.push_back(b);
    for (size_t c : *bp[b]) if (c < nobs) binned[c] = 1;
  }
  if (coreIdx.empty()) return;
  const size_t NC = coreIdx.size();

  double thr = 0.88, margin = 0.02;   // conjunctive min(tnf,depth) bar + unambiguity margin
  int iters = 1;                      // >1 over-recruits as centroids drift (purity loss)
  if (const char *e = getenv("RABBIT_BIN_RECRUIT_THR"))    thr = atof(e);
  if (const char *e = getenv("RABBIT_BIN_RECRUIT_MARGIN")) margin = atof(e);
  if (const char *e = getenv("RABBIT_BIN_RECRUIT_ITERS"))  iters = atoi(e);

  size_t n_rec_total = 0;
  for (int it = 0; it < iters; ++it) {
    // Recompute depth + TNF centroids from the CURRENT (growing) core membership.
    std::vector<std::vector<float>> dcen(NC), tcen(NC);
#pragma omp parallel for schedule(dynamic, 8) num_threads(numThreads)
    for (size_t ci = 0; ci < NC; ++ci) {
      const ContigVector &cs = *bp[coreIdx[ci]];
      std::vector<size_t> L; for (size_t c : cs) if (c < nobs) L.push_back(c);
      if (L.empty()) continue;
      constexpr size_t CAP = 64;
      std::vector<size_t> Ls = L;
      if (Ls.size() > CAP) { std::vector<size_t> s; double st=(double)Ls.size()/CAP;
        for (size_t t=0;t<CAP;++t) s.push_back(Ls[(size_t)(t*st)]); Ls.swap(s); }
      std::vector<float> d(S, 0.f);
      for (size_t c : Ls) { const float *u = g_depth_unit.data() + c * S;
                            for (size_t k = 0; k < S; ++k) d[k] += u[k]; }
      double nd = 0; for (size_t k = 0; k < S; ++k) { d[k] /= (float)Ls.size(); nd += (double)d[k]*d[k]; }
      if (nd > 1e-12) { float iv=(float)(1.0/std::sqrt(nd)); for (size_t k=0;k<S;++k) d[k]*=iv; dcen[ci]=std::move(d); }
      std::vector<float> t(256, 0.f); double tw = 0;
      for (size_t c : cs) if (c < nobs) { const float *v = g_merge_tnf.data() + c*256;
        double w = (double)seq_lens[c]; for (int k=0;k<256;++k) t[k]+=(float)(w*v[k]); tw += w; }
      if (tw > 0) { double nt=0; for (int k=0;k<256;++k) nt+=(double)t[k]*t[k];
        if (nt>1e-12){ float iv=(float)(1.0/std::sqrt(nt)); for(int k=0;k<256;++k) t[k]*=iv; tcen[ci]=std::move(t);} }
    }
    // For each UNBINNED large contig, find best/2nd-best core by conjunctive score.
    std::vector<int> assign(nobs, -1);
#pragma omp parallel for schedule(dynamic, 256) num_threads(numThreads)
    for (size_t c = 0; c < nobs; ++c) {
      if (binned[c]) continue;
      const float *du = g_depth_unit.data() + c * S;
      const float *tu = g_merge_tnf.data() + c * 256;
      double best = -2, second = -2; long bestk = -1;
      for (size_t ci = 0; ci < NC; ++ci) {
        if (dcen[ci].empty() || tcen[ci].empty()) continue;
        const float *dc = dcen[ci].data(); double dco = 0;
        for (size_t k = 0; k < S; ++k) dco += (double)du[k]*dc[k];
        if (dco < thr) continue;                          // depth gate first (cheap reject)
        const float *tc = tcen[ci].data(); double tco = 0;
        for (int k = 0; k < 256; ++k) tco += (double)tu[k]*tc[k];
        double sc = dco < tco ? dco : tco;                // conjunctive = min
        if (sc > best) { second = best; best = sc; bestk = (long)ci; }
        else if (sc > second) second = sc;
      }
      if (bestk >= 0 && best >= thr && (second < 0 || best - second >= margin))
        assign[c] = (int)coreIdx[bestk];
    }
    size_t n_rec = 0;
    for (size_t c = 0; c < nobs; ++c)
      if (assign[c] >= 0) { bp[assign[c]]->push_back(c); binned[c] = 1; ++n_rec; }
    n_rec_total += n_rec;
    if (!n_rec) break;                                    // converged
  }
  if (n_rec_total)
    verbose_message("Composition+depth recruit (min(tnf,depth)>=%.2f margin>=%.2f): "
                    "%zu unbinned contigs into cores\n", thr, margin, n_rec_total);
}

// ═══════════════════════════════════════════════════════════════════════════
