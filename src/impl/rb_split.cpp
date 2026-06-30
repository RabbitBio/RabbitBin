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
    const size_t n = contigs.size();
    if ((int)n < splitMinContigs) return false;
    std::vector<std::vector<float>> X(n, std::vector<float>(num_depth_samples));
    for (size_t r = 0; r < n; ++r)
      for (size_t i = 0; i < (size_t)num_depth_samples; ++i)
        X[r][i] = (float)std::log(depth_at(contigs[r], i) + 1.0);
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
    if (best_k <= 1 || best_sil < g_split_sil) return false;
    sub.assign(best_k, ContigVector());
    for (size_t r = 0; r < n; ++r) sub[best_labels[r]].push_back(contigs[r]);
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

#pragma omp parallel for schedule(dynamic, 1) num_threads(numThreads)
  for (size_t bi = 0; bi < nbins; ++bi) {
    const ContigVector &contigs = *binList[bi];
    BinResult &res = results[bi];

    // ── MAIN PATH — IDENTICAL to the historical behaviour for >=min_bin_bp ──
    // bins.  This guarantees confidently-sized bins (the only ones that matter
    // on CAMI1/2) are byte-for-byte unchanged, so recovery there cannot drop.
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
  // bins were going to be discarded, every retained piece is a pure addition —
  // it cannot lower recovery on datasets where such bins are rare (CAMI1/2).
  std::vector<size_t> captured;
  for (size_t bi = 0; bi < nbins; ++bi)
    if (results[bi].kind == 0) captured.push_back(bi);

  std::vector<std::vector<ContigVector>> recov(captured.size());
  if (coh_gate) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(numThreads)
    for (size_t ci = 0; ci < captured.size(); ++ci) {
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
// Paired-end bin merge — rejoin fragments of ONE genome
// ═══════════════════════════════════════════════════════════════════════════
// LPA / abundance split can scatter a single genome across several pure bins
// (coverage sub-clusters, or separate LPA communities).  Each fragment is then
// individually < min_bin_bp, so the academic ">=200kb bin" protocol discards the
// genome even though the UNION of its fragments is a high-quality MAG.  This pass
// agglomerates such bins using PAIRED-END read linkage (mates aligned to two
// different contigs) as the decisive same-genome signal — a physical, ~98%-
// intra-genome adjacency from the same BAM, orthogonal to composition/abundance.
// A cheap depth-centroid correlation is an additional abundance sanity check.
// Safety: a merge requires PE evidence; never fuses two independent >=floor
// cores; with fragment-only mode it cannot touch any >=floor bin at all.
static void consolidate_bins(BinMap &cls, size_t floor) {
  if (!g_bin_merge) return;
  const size_t S = (size_t)num_depth_samples;
  if (S < 1 || g_depth_unit.empty()) return;  // need multi-sample depth

  std::vector<int> ids; std::vector<ContigVector*> bp;
  ids.reserve(cls.size()); bp.reserve(cls.size());
  for (auto &kv : cls) { ids.push_back(kv.first); bp.push_back(&kv.second); }
  const size_t B = bp.size();
  if (B < 2) return;

  // GAIN-ONLY guard: only sub-floor fragments may participate in a merge.  A bin
  // that is already >= floor (a potential HQ MAG) is never touched, so the >=200kb
  // recovery count is monotone — merging can only CREATE new >=floor bins by
  // uniting a genome's scattered fragments, never break an existing good bin.
  auto bin_bp = [&](size_t b) -> size_t {
    size_t s = 0;
    for (size_t c : *bp[b]) s += (c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs];
    return s;
  };
  // Fragment-only mode (RABBIT_BIN_MERGE_FRAGONLY=1) forbids touching any
  // >=floor bin (strictly monotone, can only build new bins from fragments).
  // Default OFF: with the paired-end gate (~98% precise) a fragment may also be
  // merged INTO a >=floor core to COMPLETE it, which is where the real >=200kb HQ
  // gain is; PE prevents the cross-genome contamination that sank the TNF+depth-
  // only attempt.  A merge still requires at least one participant to be a
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
  const bool pe_required = !pe_pair.empty() || !g_pe_links_compact.empty();

  constexpr size_t CAP = 48;   // depth-centroid sample of large members
  // Per-bin depth centroid (normalised mean unit-rank row, dim S) — a cheap
  // abundance-consistency check that reuses the precomputed unit vectors.
  std::vector<std::vector<float>> cent(B);

  auto sampleLarge = [&](size_t b, size_t cap) {
    std::vector<size_t> L;
    for (size_t c : *bp[b]) if (c < nobs) L.push_back(c);
    if (L.size() > cap) {
      std::vector<size_t> s; s.reserve(cap);
      const double step = (double)L.size() / (double)cap;
      for (size_t t = 0; t < cap; ++t) s.push_back(L[(size_t)(t * step)]);
      L.swap(s);
    }
    return L;
  };

#pragma omp parallel for schedule(dynamic, 16) num_threads(numThreads)
  for (size_t b = 0; b < B; ++b) {
    std::vector<size_t> L = sampleLarge(b, CAP);
    if (!L.empty()) {
      std::vector<float> m(S, 0.f);
      for (size_t c : L) { const float *u = g_depth_unit.data() + c * S;
                           for (size_t k = 0; k < S; ++k) m[k] += u[k]; }
      double nrm = 0;
      for (size_t k = 0; k < S; ++k) { m[k] /= (float)L.size(); nrm += (double)m[k] * m[k]; }
      if (nrm > 1e-12) { float inv = (float)(1.0 / std::sqrt(nrm));
                         for (size_t k = 0; k < S; ++k) m[k] *= inv; cent[b] = std::move(m); }
    }
  }

  // Cross-bin gate: the PHYSICAL paired-end linkage (pe_pair) is the decisive
  // ~98%-precise same-genome signal; the depth-centroid correlation is only a
  // cheap abundance sanity check (a covarying but distinct genome is still
  // rejected by the absence of read-pair linkage).  No composition/TNF term:
  // empirically it never changed the >=200kb HQ count once PE was required.
  double dbar = 0.50;
  if (const char *e = getenv("RABBIT_BIN_MERGE_DCORR")) dbar = atof(e);

  struct MPair { uint32_t a, b; };
  int nthr = (int)std::max<size_t>(1, numThreads);
  std::vector<std::vector<MPair>> tl(nthr);
#pragma omp parallel num_threads(numThreads)
  {
    int tid = omp_get_thread_num();
    if (tid >= nthr) tid = nthr - 1;
#pragma omp for schedule(dynamic, 8)
    for (size_t a = 0; a < B; ++a) {
      if (cent[a].empty()) continue;
      const float *ca = cent[a].data();
      for (size_t b2 = a + 1; b2 < B; ++b2) {
        if (cent[b2].empty()) continue;
        // At least one must be a fragment (never fuse two independent cores);
        // fragment-only mode additionally forbids touching any core.
        if (frag_only) { if (!is_frag[a] || !is_frag[b2]) continue; }
        else           { if (!is_frag[a] && !is_frag[b2]) continue; }
        const float *cb = cent[b2].data();
        double dc = 0; for (size_t k = 0; k < S; ++k) dc += (double)ca[k] * cb[k];
        if (dc < dbar) continue;
        // Decisive gate: physical paired-end support confirms same organism.
        if (pe_required) {
          uint32_t lo = (uint32_t)std::min(a, b2), hi = (uint32_t)std::max(a, b2);
          auto it = pe_pair.find(((uint64_t)lo << 32) | hi);
          if (it == pe_pair.end() || it->second < pe_min) continue;
        } else continue;  // no PE evidence available → do not merge (safe)
        tl[tid].push_back({(uint32_t)a, (uint32_t)b2});
      }
    }
  }
  std::vector<MPair> pairs;
  for (auto &v : tl) pairs.insert(pairs.end(), v.begin(), v.end());
  if (pairs.empty()) return;
  std::sort(pairs.begin(), pairs.end(),
            [](const MPair &x, const MPair &y) { return x.a != y.a ? x.a < y.a : x.b < y.b; });

  std::vector<size_t> uf(B);
  std::vector<char> comp_core(B);            // does this component contain a >=floor core?
  for (size_t i = 0; i < B; ++i) { uf[i] = i; comp_core[i] = is_frag[i] ? 0 : 1; }
  auto find = [&](size_t x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; };
  size_t merges = 0;
  for (auto &p : pairs) {
    size_t ra = find(p.a), rb = find(p.b);
    if (ra == rb) continue;
    // Never fuse two components that each already contain a >=floor core: that
    // would join two independent (potentially HQ) MAGs.  Fragments may still flow
    // into a single core, or unite into a new core.
    if (comp_core[ra] && comp_core[rb]) continue;
    size_t root = std::min(ra, rb), other = std::max(ra, rb);
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
  verbose_message("Bin merge (paired-end, dcorr>=%.2f pe>=%u): %zu -> %zu bins "
                  "(%zu merges)\n", dbar, pe_min, before, cls.size(), merges);
}

// ═══════════════════════════════════════════════════════════════════════════
// Marker-free decontamination — shed foreign contigs from a bin (subtraction)
// ═══════════════════════════════════════════════════════════════════════════
// A contig that drifted into a bin by composition but belongs to another genome
// has a DIFFERENT abundance profile (low depth-correlation to the bin centroid)
// and NO paired-end linkage to the rest of the bin.  Removing such contigs raises
// purity (pushing complete-but-impure bins past 95%) without reducing the main
// genome's completeness.  PE is a PROTECT signal: any contig physically linked to
// a bin-mate is kept regardless of its depth correlation, so genuinely co-binned
// low-coverage contigs are never shed.  Subtraction-only; never enlarges a bin.
static void decontaminate_bins(BinMap &cls, size_t floor) {
  if (!g_bin_decontam) return;
  const size_t S = (size_t)num_depth_samples;
  if (S < 1 || g_depth_unit.empty() || g_pe_links_compact.empty()) return;

  // Per-contig PE neighbours (binning index space), via compact-row -> name -> idx.
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
  phmap::flat_hash_map<int, std::vector<int>> peadj;
  for (auto &lk : g_pe_links_compact) {
    int ia = compact2idx(std::get<0>(lk)), ib = compact2idx(std::get<1>(lk));
    if (ia < 0 || ib < 0 || ia == ib) continue;
    peadj[ia].push_back(ib);
    peadj[ib].push_back(ia);
  }

  // Only decontaminate a bin that is OTHERWISE tight: its members' mean depth-
  // correlation to the centroid is >= mu_min.  In such a bin a member sitting
  // below abs_thr is clearly foreign (the rest agree strongly).  A uniformly
  // low-coverage bin (mean < mu_min, e.g. CAMI2 plant's fragmented genomes) is
  // skipped entirely, so legitimate noisy contigs are never gutted.  This
  // tightness gate (not a per-dataset threshold) is what keeps the rule general.
  double abs_thr = 0.80;   // a member below this corr is a removal candidate
  double mu_min  = 0.82;   // only clean bins whose mean member corr is >= this
  if (const char *e = getenv("RABBIT_BIN_DECONTAM_DCORR"))  abs_thr = atof(e);
  if (const char *e = getenv("RABBIT_BIN_DECONTAM_MUMIN"))  mu_min  = atof(e);

  auto bin_bp = [&](const ContigVector &cs) -> size_t {
    size_t s = 0;
    for (size_t c : cs) s += (c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs];
    return s;
  };

  // Gather bins into an indexable list for parallel processing.
  std::vector<int> ids; std::vector<ContigVector*> bp;
  for (auto &kv : cls) { ids.push_back(kv.first); bp.push_back(&kv.second); }
  const size_t B = bp.size();
  std::vector<std::vector<size_t>> shed(B);   // removed contigs per bin

#pragma omp parallel for schedule(dynamic, 8) num_threads(numThreads)
  for (size_t b = 0; b < B; ++b) {
    ContigVector &cs = *bp[b];
    if (bin_bp(cs) < floor) continue;          // only the >=200kb bins matter here
    // Depth centroid over (capped) large members.
    std::vector<size_t> L;
    for (size_t c : cs) if (c < nobs) L.push_back(c);
    if (L.size() < 4) continue;                // too small to judge an outlier
    constexpr size_t CAP = 64;
    std::vector<size_t> Ls = L;
    if (Ls.size() > CAP) {
      std::vector<size_t> s; s.reserve(CAP);
      const double step = (double)Ls.size() / (double)CAP;
      for (size_t t = 0; t < CAP; ++t) s.push_back(Ls[(size_t)(t * step)]);
      Ls.swap(s);
    }
    std::vector<float> cent(S, 0.f);
    for (size_t c : Ls) { const float *u = g_depth_unit.data() + c * S;
                          for (size_t k = 0; k < S; ++k) cent[k] += u[k]; }
    double nrm = 0;
    for (size_t k = 0; k < S; ++k) { cent[k] /= (float)Ls.size(); nrm += (double)cent[k] * cent[k]; }
    if (nrm <= 1e-12) continue;
    float inv = (float)(1.0 / std::sqrt(nrm));
    for (size_t k = 0; k < S; ++k) cent[k] *= inv;
    // Member set for PE-membership test.
    phmap::flat_hash_set<int> memset;
    memset.reserve(cs.size() * 2);
    for (size_t c : cs) if (c < nobs) memset.insert((int)c);
    // Per-bin corr to centroid over ALL large members; mean = bin tightness.
    std::vector<double> corr(L.size());
    double sum = 0;
    for (size_t i = 0; i < L.size(); ++i) {
      const float *u = g_depth_unit.data() + L[i] * S;
      double dc = 0; for (size_t k = 0; k < S; ++k) dc += (double)u[k] * cent[k];
      corr[i] = dc; sum += dc;
    }
    const double mu = sum / (double)L.size();
    if (mu < mu_min) continue;                 // loose/low-coverage bin → don't touch
    for (size_t i = 0; i < L.size(); ++i) {
      if (corr[i] >= abs_thr) continue;        // abundance-consistent → keep
      const int c = (int)L[i];
      bool pe_to_bin = false;
      auto it = peadj.find(c);
      if (it != peadj.end())
        for (int nb : it->second) if (nb != c && memset.count(nb)) { pe_to_bin = true; break; }
      if (pe_to_bin) continue;                 // physically linked → keep (protect)
      shed[b].push_back((size_t)c);            // foreign by abundance AND no linkage
    }
  }

  // Apply removals: pull shed contigs out of their bins, emit each as its own bin
  // (kept in the map so nothing is lost; almost always < floor and thus filtered
  // out of the >=200kb evaluation, exactly as a contaminant should be).
  size_t n_shed = 0, n_bins_touched = 0;
  int next_id = 0;
  for (int id : ids) if (id >= next_id) next_id = id + 1;
  for (size_t b = 0; b < B; ++b) {
    if (shed[b].empty()) continue;
    ++n_bins_touched;
    phmap::flat_hash_set<size_t> rm(shed[b].begin(), shed[b].end());
    ContigVector kept;
    kept.reserve(bp[b]->size());
    for (size_t c : *bp[b]) if (!rm.count(c)) kept.push_back(c);
    *bp[b] = std::move(kept);
    for (size_t c : shed[b]) { cls[next_id++] = ContigVector{c}; ++n_shed; }
  }
  if (n_shed)
    verbose_message("Decontaminate (tight bins mu>=%.2f, shed corr<%.2f & no PE): "
                    "%zu contigs from %zu bins\n", mu_min, abs_thr, n_shed, n_bins_touched);
}

// ═══════════════════════════════════════════════════════════════════════════
// Composition+depth recruit — attach UNBINNED contigs to the core they belong to
// ═══════════════════════════════════════════════════════════════════════════
// ~38% of large contigs end up unbinned, and the missing mass of near-complete
// MAGs is mostly such unbinned contigs.  For each unbinned large contig, attach
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
