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
    std::vector<std::vector<float>> X(n, std::vector<float>(num_depth_samples));
    for (size_t r = 0; r < n; ++r)
      for (size_t i = 0; i < num_depth_samples; ++i)
        X[r][i] = (float)std::log(depth_at(contigs[r], i) + 1.0);
    std::vector<int> labels = rb_kmeans(X, k, rng);
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
// Mean silhouette of a labeling over feature rows X (Euclidean).  O(n^2 d); for
// large bins we sample up to g_sil_sample_cap rows so this stays cheap.
static double fkmv_silhouette(const std::vector<std::vector<float>> &X,
                              const std::vector<int> &labels, int k,
                              std::mt19937 &rng) {
  const size_t n = X.size();
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
  const size_t d = X[0].size();
  auto dist = [&](size_t a, size_t b) -> double {
    double s = 0;
    for (size_t i = 0; i < d; ++i) { double t = X[a][i] - X[b][i]; s += t * t; }
    return std::sqrt(s);
  };
  // per-cluster sizes within the sample
  std::vector<int> csz(k, 0);
  for (size_t a : idx) csz[labels[a]]++;
  double sil_sum = 0;
  size_t cnt = 0;
  for (size_t ii = 0; ii < m; ++ii) {
    size_t a = idx[ii];
    int la = labels[a];
    if (csz[la] <= 1) continue;  // singleton cluster: s=0, skip
    std::vector<double> sumd(k, 0.0);
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
static void abundance_guided_split(BinMap &cls) {
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

  // ── Parallel over bins ──────────────────────────────────────────────────
  // Each bin's KMeans+silhouette sweep is independent, so distribute bins
  // across threads (this was a serial Amdahl bottleneck). Determinism is kept
  // independent of thread count by seeding each bin's RNG from (seed, binIdx)
  // instead of one shared, sequentially-advanced generator.
  std::vector<const ContigVector*> binList;
  binList.reserve(cls.size());
  for (auto &kv : cls) binList.push_back(&kv.second);
  const size_t nbins = binList.size();

  // kind: 0 = dropped (<min_bin_bp), 1 = kept whole, 2 = split into emit[]
  struct BinResult { std::vector<ContigVector> emit; int kind = 1; };
  std::vector<BinResult> results(nbins);

#pragma omp parallel for schedule(dynamic, 1) num_threads(numThreads)
  for (size_t bi = 0; bi < nbins; ++bi) {
    const ContigVector &contigs = *binList[bi];
    BinResult &res = results[bi];

    if (bin_bp(contigs) < min_bin_bp) { res.kind = 0; continue; }
    const size_t n = contigs.size();
    if ((int)n < splitMinContigs) { res.kind = 1; res.emit.push_back(contigs); continue; }

    std::vector<std::vector<float>> X(n, std::vector<float>(num_depth_samples));
    for (size_t r = 0; r < n; ++r)
      for (size_t i = 0; i < (size_t)num_depth_samples; ++i)
        X[r][i] = (float)std::log(depth_at(contigs[r], i) + 1.0);

    // Per-bin RNG: thread-count-independent, reproducible.
    std::mt19937 rng((unsigned)((seed ? (uint64_t)seed : 1ULL)
                                + bi * 2654435761ULL));
    int best_k = 1; double best_sil = -1.0; std::vector<int> best_labels;
    int kmax = std::min((int)n - 1, splitMaxK);
    for (int k = 2; k <= kmax; ++k) {
      std::vector<int> lab = rb_kmeans(X, k, rng);
      int seen = 0; { std::vector<char> u(k, 0); for (int l : lab) if (!u[l]) { u[l] = 1; ++seen; } }
      if (seen < 2) continue;
      double s = fkmv_silhouette(X, lab, k, rng);
      if (s > best_sil) { best_sil = s; best_k = k; best_labels = std::move(lab); }
    }

    if (best_k > 1 && best_sil >= g_split_sil) {
      std::vector<ContigVector> sub(best_k);
      for (size_t r = 0; r < n; ++r) sub[best_labels[r]].push_back(contigs[r]);
      for (int t = 0; t < best_k; ++t) {
        if (sub[t].empty()) continue;
        if (bin_bp(sub[t]) < min_bin_bp) continue;  // drop tiny sub-products
        res.emit.push_back(std::move(sub[t]));
      }
      if (!res.emit.empty()) res.kind = 2;
      else { res.kind = 1; res.emit.push_back(contigs); }  // nothing survived
    } else {
      res.kind = 1; res.emit.push_back(contigs);
    }
  }

  // ── Sequential merge (deterministic bin order + counts) ──────────────────
  BinMap out;
  int next = 0;
  size_t n_split = 0, n_kept = 0, n_drop = 0;
  for (size_t bi = 0; bi < nbins; ++bi) {
    BinResult &res = results[bi];
    if (res.kind == 0) { ++n_drop; continue; }
    if (res.kind == 2) ++n_split; else ++n_kept;
    for (auto &v : res.emit) out[next++] = std::move(v);
  }
  cls.swap(out);
  verbose_message("Abundance split (sil>=%.2f): %zu kept, %zu split, "
                  "%zu dropped(<min_bin_bp) -> %d bins\n",
                  g_split_sil, n_kept, n_split, n_drop, next);
  // Sub-products may be < min_bin_bp; base bins were already gated above.
  min_bin_bp = 0;
}

// ═══════════════════════════════════════════════════════════════════════════
