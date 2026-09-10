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

  // Large-contig depth is rank-transformed for graph scoring; short-contig
  // depth remains raw until the final recruitment pass.
  auto depth_at = [&](size_t c, size_t i) -> double {
    if (c < nobs) {
      if (!g_large_means.empty() && c * (size_t)num_depth_samples + i < g_large_means.size())
        return (double)g_large_means[c * num_depth_samples + i];
      return (double)depth_matrix(c, i);
    }
    size_t s = c - nobs;
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
  // Split products may be < min_bin_bp; output_bins() applies the final strict
  // size filter.
}

// ═══════════════════════════════════════════════════════════════════════════
// Marker-FREE bin splitting (Phase 2; --split-bins) — abundance multimodality
// ═══════════════════════════════════════════════════════════════════════════
// RB_KM_PROF accumulators: which part of split_bin the time actually goes to.
static std::atomic<double> g_ms_feat{0}, g_ms_kmeans{0}, g_ms_sil{0};
static inline void rb_prof_add(std::atomic<double> &a,
                               std::chrono::steady_clock::time_point t0) {
  const double v = std::chrono::duration<double, std::milli>(
                       std::chrono::steady_clock::now() - t0).count();
  double cur = a.load(std::memory_order_relaxed);
  while (!a.compare_exchange_weak(cur, cur + v, std::memory_order_relaxed)) {}
}

// Mean silhouette of a labeling over feature rows X (Euclidean).  O(n^2 d); for
// large bins we sample up to g_sil_sample_cap rows so this stays cheap.
#include "rb_silhouette.h"

static double fkmv_silhouette(const float *X, size_t n, size_t d,
                              const std::vector<int> &labels, int k,
                              std::mt19937 &rng) {
  auto dist = [&](size_t a, size_t b) -> double {
    return std::sqrt(rb_sqdist(X + a * d, X + b * d, d));
  };
  return rb_mean_silhouette(n, labels, k, g_sil_sample_cap, rng, dist);
}

// Re-split internally multi-modal bins using per-sample log-abundance KMeans,
// choosing k by silhouette.  No markers / gene prediction needed.

// Re-split output-sized, internally multi-modal bins using per-sample
// log-abundance K-means, choosing K by mean silhouette.
static void abundance_guided_split_current(BinMap &cls) {
  if (num_depth_samples < 1) {
    verbose_message("Abundance split: needs >=1 abundance sample — skipped\n");
    return;
  }

  auto depth_at = [&](size_t c, size_t i) -> double {
    if (c < nobs) {
      const size_t off = c * (size_t)num_depth_samples + i;
      if (off < g_large_means.size()) return (double)g_large_means[off];
      return (double)depth_matrix(c, i);
    }
    const size_t s = c - nobs;
    return (double)small_depth_matrix(s, i);
  };
  auto bin_bp = [&](const ContigVector &contigs) -> size_t {
    size_t bp = 0;
    for (size_t c : contigs)
      bp += (c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs];
    return bp;
  };

  const bool audit_splits = getenv("RB_SPLIT_AUDIT") != nullptr;
  struct SplitAudit {
    bool attempted = false, passed = false;
    int best_k = 1;
    double best_sil = -1.0;
  };

  auto split_bin = [&](const ContigVector &items, size_t bi,
                       std::vector<ContigVector> &sub,
                       SplitAudit *audit) -> bool {
    const size_t n = items.size();
    if ((int)n < splitMinContigs) return false;
    if (audit) audit->attempted = true;

    const size_t sd = (size_t)num_depth_samples;
    auto tk0 = std::chrono::steady_clock::now();
    std::vector<float> X(n * sd);
    for (size_t r = 0; r < n; ++r)
      for (size_t i = 0; i < sd; ++i)
        X[r * sd + i] = (float)std::log(depth_at(items[r], i) + 1.0);
    if (g_km_prof) rb_prof_add(g_ms_feat, tk0);

    const uint64_t rng_seed =
        (seed ? (uint64_t)seed : 1ULL) + bi * 2654435761ULL;
    std::mt19937 rng((unsigned)rng_seed);
    int best_k = 1;
    double best_sil = -1.0;
    std::vector<int> best_labels;
    const int kmax = std::min((int)n - 1, splitMaxK);
    for (int k = 2; k <= kmax; ++k) {
      tk0 = std::chrono::steady_clock::now();
      std::vector<int> lab = rb_kmeans(X.data(), n, sd, k, rng);
      if (g_km_prof) rb_prof_add(g_ms_kmeans, tk0);
      std::vector<char> seen((size_t)k, 0);
      int represented = 0;
      for (int label : lab)
        if (!seen[(size_t)label]) {
          seen[(size_t)label] = 1;
          ++represented;
        }
      if (represented < 2) continue;

      tk0 = std::chrono::steady_clock::now();
      const double sil = fkmv_silhouette(X.data(), n, sd, lab, k, rng);
      if (g_km_prof) rb_prof_add(g_ms_sil, tk0);
      if (sil > best_sil) {
        best_sil = sil;
        best_k = k;
        best_labels = std::move(lab);
      }
    }

    const bool passed = best_k > 1 && best_sil >= g_split_sil;
    if (audit) {
      audit->best_k = best_k;
      audit->best_sil = best_sil;
      audit->passed = passed;
    }
    if (!passed) return false;

    sub.assign((size_t)best_k, ContigVector());
    for (size_t r = 0; r < n; ++r)
      sub[(size_t)best_labels[r]].push_back(items[r]);
    return true;
  };

  std::vector<const ContigVector *> bins;
  bins.reserve(cls.size());
  for (auto &kv : cls) bins.push_back(&kv.second);
  const size_t nbins = bins.size();

  std::vector<size_t> order(nbins);
  for (size_t i = 0; i < nbins; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return bins[a]->size() > bins[b]->size();
  });

  struct BinResult {
    std::vector<ContigVector> emit;
    size_t dropped_children = 0;
    bool split = false;
  };
  std::vector<BinResult> results(nbins);
  std::vector<SplitAudit> audits(audit_splits ? nbins : 0);
  const bool profile_bins = getenv("RB_SPLIT_PROF") != nullptr;
  std::vector<double> bin_ms(profile_bins ? nbins : 0, 0.0);

#pragma omp parallel for schedule(dynamic, 1) num_threads(numThreads)
  for (size_t oi = 0; oi < nbins; ++oi) {
    const size_t bi = order[oi];
    const ContigVector &contigs = *bins[bi];
    const auto begin = profile_bins ? std::chrono::steady_clock::now()
                                    : std::chrono::steady_clock::time_point{};

    if (bin_bp(contigs) >= min_bin_bp) {
      std::vector<ContigVector> children;
      const bool accepted = split_bin(
          contigs, bi, children, audit_splits ? &audits[bi] : nullptr);
      if (accepted) {
        for (auto &child : children) {
          if (child.empty()) continue;
          if (bin_bp(child) >= min_bin_bp)
            results[bi].emit.push_back(std::move(child));
          else
            ++results[bi].dropped_children;
        }
        if (!results[bi].emit.empty())
          results[bi].split = true;
        else
          results[bi].emit.push_back(contigs);
      } else {
        results[bi].emit.push_back(contigs);
      }
    }

    if (profile_bins)
      bin_ms[bi] = std::chrono::duration<double, std::milli>(
                       std::chrono::steady_clock::now() - begin)
                       .count();
  }

  if (g_km_prof) {
    fprintf(stderr,
            "[RB_KM_PROF] split_bin CPU-ms: features=%.0f kmeans=%.0f "
            "silhouette=%.0f\n",
            g_ms_feat.load(), g_ms_kmeans.load(), g_ms_sil.load());
    const uint64_t restarts = g_km_restarts.load();
    const uint64_t iters = g_km_iters.load();
    fprintf(stderr,
            "[RB_KM_PROF] %lu kmeans restarts, %lu Lloyd iters (avg %.1f), "
            "%lu hit the 200 cap (%.1f%%)\n",
            (unsigned long)restarts, (unsigned long)iters,
            restarts ? (double)iters / restarts : 0.0,
            (unsigned long)g_km_capped.load(),
            restarts ? 100.0 * g_km_capped.load() / restarts : 0.0);
  }
  if (profile_bins) {
    std::vector<size_t> by_time(nbins);
    for (size_t i = 0; i < nbins; ++i) by_time[i] = i;
    std::sort(by_time.begin(), by_time.end(),
              [&](size_t x, size_t y) { return bin_ms[x] > bin_ms[y]; });
    double total_ms = 0.0;
    for (double ms : bin_ms) total_ms += ms;
    fprintf(stderr, "[RB_SPLIT_PROF] %zu bins, total work %.1f ms\n",
            nbins, total_ms);
    for (size_t i = 0; i < std::min<size_t>(8, nbins); ++i)
      fprintf(stderr, "[RB_SPLIT_PROF]   #%zu contigs=%zu  %.1f ms\n",
              by_time[i], bins[by_time[i]]->size(), bin_ms[by_time[i]]);
  }

  if (audit_splits) {
    std::ofstream out(outFile + ".split_audit.tsv");
    if (!out) throw std::runtime_error("Cannot write split audit");
    out.precision(17);
    out << "parent_index\tparent_contigs\tparent_bp\tattempted\tbest_k"
           "\tbest_silhouette\tthreshold\tthreshold_pass\tretained_pieces"
           "\tretained_contigs\tretained_bp\n";
    for (size_t bi = 0; bi < nbins; ++bi) {
      size_t retained_contigs = 0, retained_bp = 0;
      for (const auto &piece : results[bi].emit) {
        retained_contigs += piece.size();
        retained_bp += bin_bp(piece);
      }
      const SplitAudit &audit = audits[bi];
      out << bi << '\t' << bins[bi]->size() << '\t' << bin_bp(*bins[bi])
          << '\t' << audit.attempted << '\t' << audit.best_k << '\t'
          << audit.best_sil << '\t' << g_split_sil << '\t' << audit.passed
          << '\t' << results[bi].emit.size() << '\t' << retained_contigs
          << '\t' << retained_bp << '\n';
    }
  }

  BinMap refined;
  int next = 0;
  size_t kept = 0, split = 0, dropped = 0, dropped_children = 0;
  for (size_t bi = 0; bi < nbins; ++bi) {
    if (results[bi].emit.empty()) {
      ++dropped;
      continue;
    }
    if (results[bi].split) ++split;
    else ++kept;
    dropped_children += results[bi].dropped_children;
    for (auto &piece : results[bi].emit)
      refined[next++] = std::move(piece);
  }
  cls.swap(refined);
  verbose_message(
      "Abundance split (sil>=%.2f): %zu kept, %zu split, %zu dropped below "
      "minimum size (%zu split children) -> %d bins\n",
      g_split_sil, kept, split, dropped, dropped_children, next);
}

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
// statistical recruitment reduces to one statistic against one
// chi-square critical value.  The bar follows from the sample count S and a significance
// level; short contigs receive exactly the wider tolerance the model implies.

// Selectively attach all unassigned large and short contigs to immutable
// output-sized cores after abundance splitting.
// Recruitment uses coverage only, matching the coverage-driven graph weighting.
// A leave-one-out classification of the core members supplies positive and
// negative predictions, from which the run learns one ROC/Youden confidence
// boundary.  No reference labels or dataset-specific cutoffs are used.
static void recruit_unbinned_to_cores(BinMap &cls, size_t floor) {
  if (!g_bin_recruit) return;
  const size_t samples = (size_t)num_depth_samples;
  if (samples < 2 || g_depth_unit.size() < nobs * samples) return;
  const size_t candidate_count = nobs + nobs1;

  std::vector<ContigVector *> bins;
  bins.reserve(cls.size());
  for (auto &entry : cls) bins.push_back(&entry.second);

  auto bin_bp = [](const ContigVector &contigs) {
    size_t total = 0;
    for (size_t c : contigs)
      total += c < nobs ? seq_lens[c] : small_seq_lens[c - nobs];
    return total;
  };

  std::vector<size_t> core_bins;
  std::vector<char> binned(candidate_count, 0);
  for (size_t b = 0; b < bins.size(); ++b) {
    if (bin_bp(*bins[b]) >= floor) core_bins.push_back(b);
    for (size_t c : *bins[b])
      if (c < candidate_count) binned[c] = 1;
  }
  if (core_bins.size() < 2) return;
  const size_t core_count = core_bins.size();

  struct CoreFeatures {
    size_t dimension = 0;
    std::vector<double> sums;
    std::vector<double> norms;
  };

  auto row_norms = [](const std::vector<float> &features, size_t rows,
                      size_t dimension) {
    std::vector<float> norms(rows, 0.0f);
#pragma omp parallel for schedule(static) num_threads(numThreads)
    for (size_t c = 0; c < rows; ++c) {
      const float *row = features.data() + c * dimension;
      double sum_squares = 0.0;
      for (size_t k = 0; k < dimension; ++k)
        sum_squares += (double)row[k] * row[k];
      if (sum_squares > 1e-30) norms[c] = (float)std::sqrt(sum_squares);
    }
    return norms;
  };

  auto build_core_features = [&](const std::vector<float> &features,
                                 const std::vector<float> &norms,
                                 size_t dimension) {
    CoreFeatures result;
    result.dimension = dimension;
    result.sums.assign(core_count * dimension, 0.0);
    result.norms.assign(core_count, 0.0);
#pragma omp parallel for schedule(dynamic, 8) num_threads(numThreads)
    for (size_t ci = 0; ci < core_count; ++ci) {
      double *sum = result.sums.data() + ci * dimension;
      for (size_t c : *bins[core_bins[ci]]) {
        if (c >= nobs || !(norms[c] > 0.0f)) continue;
        const float *row = features.data() + c * dimension;
        const double scale = 1.0 / (double)norms[c];
        for (size_t k = 0; k < dimension; ++k)
          sum[k] += scale * (double)row[k];
      }
      double sum_squares = 0.0;
      for (size_t k = 0; k < dimension; ++k)
        sum_squares += sum[k] * sum[k];
      if (sum_squares > 1e-30) result.norms[ci] = std::sqrt(sum_squares);
    }
    return result;
  };

  const std::vector<float> depth_norms =
      row_norms(g_depth_unit, nobs, samples);
  const CoreFeatures depth_cores =
      build_core_features(g_depth_unit, depth_norms, samples);

  // Short contigs were excluded from graph construction. Build their rank-based
  // coverage vectors here so long and short candidates use the same score and
  // the same learned confidence boundary.
  std::vector<float> small_unit(nobs1 * samples, 0.0f);
#pragma omp parallel num_threads(numThreads)
  {
    std::vector<StoredDistance> row(samples);
#pragma omp for schedule(static)
    for (size_t s = 0; s < nobs1; ++s) {
      for (size_t k = 0; k < samples; ++k)
        row[k] = small_depth_matrix(s, k);
      rank(row, row);
      double mean = 0.0;
      for (size_t k = 0; k < samples; ++k) mean += row[k];
      mean /= (double)samples;
      double sum_squares = 0.0;
      for (size_t k = 0; k < samples; ++k) {
        const double d = (double)row[k] - mean;
        sum_squares += d * d;
      }
      if (!(sum_squares > 1e-30)) continue;
      const double inv = 1.0 / std::sqrt(sum_squares);
      float *out = small_unit.data() + s * samples;
      for (size_t k = 0; k < samples; ++k)
        out[k] = (float)(((double)row[k] - mean) * inv);
    }
  }
  const std::vector<float> small_norms =
      row_norms(small_unit, nobs1, samples);

  auto feature_cosine = [](size_t c, size_t ci,
                           const std::vector<float> &features,
                           const std::vector<float> &norms,
                           const CoreFeatures &cores,
                           bool leave_self_out) {
    if (!(norms[c] > 0.0f))
      return -std::numeric_limits<double>::infinity();
    const float *row = features.data() + c * cores.dimension;
    const double *sum = cores.sums.data() + ci * cores.dimension;
    double dot = 0.0;
    double sum_squares = leave_self_out ? 0.0 : cores.norms[ci] * cores.norms[ci];
    for (size_t k = 0; k < cores.dimension; ++k) {
      const double core_value = leave_self_out
          ? sum[k] - (double)row[k] / (double)norms[c]
          : sum[k];
      dot += (double)row[k] / (double)norms[c] * core_value;
      if (leave_self_out) sum_squares += core_value * core_value;
    }
    if (!(sum_squares > 1e-30))
      return -std::numeric_limits<double>::infinity();
    return std::max(-1.0, std::min(1.0, dot / std::sqrt(sum_squares)));
  };

  auto depth_score = [&](size_t c, size_t ci, bool leave_self_out) {
    if (c >= nobs)
      return feature_cosine(c - nobs, ci, small_unit, small_norms,
                            depth_cores, false);
    return feature_cosine(c, ci, g_depth_unit, depth_norms, depth_cores,
                          leave_self_out);
  };
  auto log_residual_ratio = [](double best, double second) {
    if (!std::isfinite(best) || !std::isfinite(second))
      return -std::numeric_limits<double>::infinity();
    const double winner_distance = std::max(0.0, 1.0 - best);
    const double runner_distance = std::max(0.0, 1.0 - second);
    if (winner_distance == 0.0)
      return runner_distance > 0.0
                 ? std::numeric_limits<double>::infinity() : 0.0;
    if (runner_distance == 0.0) return 0.0;
    return std::log(runner_distance / winner_distance);
  };

  struct Evidence { double confidence; bool positive; };
  std::vector<std::vector<Evidence>> evidence(core_count);
#pragma omp parallel for schedule(dynamic, 4) num_threads(numThreads)
  for (size_t ci = 0; ci < core_count; ++ci) {
    for (size_t c : *bins[core_bins[ci]]) {
      if (c >= nobs || !(depth_norms[c] > 0.0f))
        continue;
      const double own = depth_score(c, ci, true);
      double best_other = -std::numeric_limits<double>::infinity();
      double second_other = -std::numeric_limits<double>::infinity();
      for (size_t cj = 0; cj < core_count; ++cj) {
        if (cj == ci) continue;
        const double score = depth_score(c, cj, false);
        if (score > best_other) {
          second_other = best_other;
          best_other = score;
        } else if (score > second_other) {
          second_other = score;
        }
      }
      if (!std::isfinite(second_other)) continue;
      if (own > best_other) {
        evidence[ci].push_back({log_residual_ratio(own, best_other), true});
      } else if (best_other > own) {
        evidence[ci].push_back(
            {log_residual_ratio(best_other, std::max(own, second_other)), false});
      }
    }
  }

  std::vector<double> positives, negatives;
  for (const auto &core_evidence : evidence)
    for (const Evidence &item : core_evidence)
      (item.positive ? positives : negatives).push_back(item.confidence);
  if (positives.empty() || negatives.empty()) return;

  struct RocPoint { double value; bool positive; };
  std::vector<RocPoint> points;
  points.reserve(positives.size() + negatives.size());
  for (double value : positives) points.push_back({value, true});
  for (double value : negatives) points.push_back({value, false});
  std::sort(points.begin(), points.end(), [](const RocPoint &a, const RocPoint &b) {
    return a.value > b.value;
  });
  size_t true_positive = 0, false_positive = 0;
  double boundary = std::numeric_limits<double>::infinity();
  double boundary_tpr = 0.0, boundary_fpr = 0.0;
  double best_youden = -std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < points.size();) {
    size_t next = i;
    while (next < points.size() && points[next].value == points[i].value) {
      if (points[next].positive) ++true_positive; else ++false_positive;
      ++next;
    }
    const double tpr = (double)true_positive / positives.size();
    const double fpr = (double)false_positive / negatives.size();
    if (tpr - fpr > best_youden) {
      best_youden = tpr - fpr;
      boundary = points[i].value;
      boundary_tpr = tpr;
      boundary_fpr = fpr;
    }
    i = next;
  }

  struct Choice {
    int core = -1;
    double best = -std::numeric_limits<double>::infinity();
    double second = -std::numeric_limits<double>::infinity();
  };
  auto choose_core = [&](size_t c) {
    Choice choice;
    for (size_t ci = 0; ci < core_count; ++ci) {
      const double score = depth_score(c, ci, false);
      if (score > choice.best) {
        choice.second = choice.best;
        choice.best = score;
        choice.core = (int)ci;
      } else if (score > choice.second) {
        choice.second = score;
      }
    }
    return choice;
  };

  std::vector<int> assignments(candidate_count, -1);
  size_t candidates = 0, rejected = 0;
#pragma omp parallel for schedule(dynamic, 256) num_threads(numThreads) \
    reduction(+:candidates,rejected)
  for (size_t c = 0; c < candidate_count; ++c) {
    const bool has_profile = c < nobs ? depth_norms[c] > 0.0f
                                      : small_norms[c - nobs] > 0.0f;
    if (binned[c] || !has_profile)
      continue;
    ++candidates;
    const Choice choice = choose_core(c);
    const double confidence = log_residual_ratio(choice.best, choice.second);
    if (choice.core < 0 || !(confidence >= boundary)) {
      ++rejected;
      continue;
    }
    assignments[c] = (int)core_bins[(size_t)choice.core];
  }

  size_t recruited = 0, recruited_large = 0, recruited_small = 0;
  for (size_t c = 0; c < candidate_count; ++c) {
    if (assignments[c] < 0) continue;
    bins[(size_t)assignments[c]]->push_back(c);
    ++recruited;
    if (c < nobs) ++recruited_large; else ++recruited_small;
  }
  verbose_message(
      "Post-split coverage recruit (leave-one-out ROC/Youden "
      "log residual ratio>=%.4g [TPR=%.3g,FPR=%.3g]): %zu/%zu unbinned "
      "contigs recruited (%zu large, %zu short; %zu rejected; "
      "calibration=%zu)\n",
      boundary, boundary_tpr, boundary_fpr, recruited, candidates,
      recruited_large, recruited_small, rejected,
      positives.size() + negatives.size());
}
