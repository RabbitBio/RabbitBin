// RabbitBin module: rb_certify.cpp
//
// Stability-aware MAG certification (feature #6).  Re-runs the cheap binning tail
// (edge-score -> incidence -> label propagation) under parameter/seed/sample
// perturbations on the ALREADY-BUILT similarity graph, then scores how stably each
// contig stays with its baseline bin-mates.  Leverages RabbitBin's graph reuse:
// the expensive parse/sketch/graph build happens once; perturbations cost only the
// tail.  #included into the rabbitbin.cpp translation unit (shares its globals).

// Re-runs the tail under a fixed grid of parameter/seed perturbations (and, when
// there are >=4 depth samples, a sample jackknife) and fills g_cert_* with each
// large contig's per-bin assignment stability vs the baseline partition.  Clobbers
// g.edgeScore / g.incs (caller must no longer need them — invoked after the
// baseline bins are collected, just before `g` is destroyed).
static void certify_stability(Graph &g, const std::vector<size_t> &baseline,
                              bool has_depth) {
  const size_t E = g.getEdgeCount();
  const size_t N = nobs;
  static constexpr StoredDistance SSCR_MAX = 1.0f - 1e-6f;
  g_cert_prob.assign(N, -1.0f);
  g_cert_alt.assign(N, SIZE_MAX);
  g_cert_alt_prob.assign(N, 0.0f);
  if (N == 0 || E == 0) return;

  // Baseline "binned" mask from the production incidence still resident in g.
  std::vector<char> baseBinned(N, 0);
  for (size_t i = 0; i < N; ++i) baseBinned[i] = !g.incs[i].empty();

  // Mirror the production edge-score / incidence / LP exactly (so a zero-delta
  // config reproduces the baseline partition).
  auto compute_es = [&]() {
    g.edgeScore.assign(E, 0.0f);
    if (has_depth) {
#pragma omp parallel for schedule(dynamic, 1)
      for (size_t e = 0; e < E; ++e) {
        size_t i = g.from[e], j = g.to[e];
        if (edge_is_gfa(g.sComp[e])) { g.edgeScore[e] = (StoredDistance)g_gfa_weight; continue; }
        if (num_depth_samples <= 1) {
          double w = (double)g.sComp[e];
          if (g_edge_power != 1.0) w = std::pow(w, g_edge_power);
          g.edgeScore[e] = (StoredDistance)w;
        } else {
          bool depth_ok;
          double dterm = depth_edge_term(i, j, depth_ok);
          if (!depth_ok) { g.edgeScore[e] = 0.0f; continue; }
          double w = g_w_comp * (double)g.sComp[e] + (1.0 - g_w_comp) * dterm;
          if (!std::isfinite(w) || w < (double)min_edge_weight) w = 0.0;
          if (g_edge_power != 1.0 && w > 0.0) w = std::pow(w, g_edge_power);
          g.edgeScore[e] = (StoredDistance)w;
        }
      }
    } else {
      g.edgeScore = g.sComp;
      for (auto &s : g.edgeScore) {
        if (edge_is_gfa(s)) { s = (StoredDistance)g_gfa_weight; continue; }
        if (s < min_edge_weight) s = 0.0f;
      }
    }
  };
  auto rebuild_incs = [&]() {
    g.incs.assign(N, {});
    for (size_t e = 0; e < E; ++e)
      if (g.edgeScore[e] > 0) {
        if (g.edgeScore[e] > SSCR_MAX) g.edgeScore[e] = SSCR_MAX;
        g.incs[g.from[e]].push_back(e);
        g.incs[g.to[e]].push_back(e);
      }
#pragma omp parallel for schedule(dynamic, 256)
    for (size_t v = 0; v < N; ++v) {
      auto &inc = g.incs[v];
      std::sort(inc.begin(), inc.end(), [&](size_t e1, size_t e2) {
        return g.getOtherNode(e1, v) < g.getOtherNode(e2, v);
      });
    }
  };
  auto run_lp = [&](unsigned sd) {
    std::vector<size_t> mem;
    std::vector<size_t> order(N);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), std::default_random_engine(sd));
    cluster_by_propagation(g, mem, order);
    return mem;
  };

  // Perturbation grid around the baseline (composition weight, edge power,
  // min-edge-score, seed).  ~10 runs; clamps keep params in-range.
  struct Cfg { double wc, ep, mew; unsigned sd; };
  const double bwc = g_w_comp, bep = g_edge_power, bmew = (double)min_edge_weight;
  const unsigned bsd = (unsigned)seed;
  auto clampd = [](double x, double lo, double hi) { return x < lo ? lo : (x > hi ? hi : x); };
  std::vector<Cfg> cfgs = {
    {clampd(bwc - 0.10, 0.0, 1.0), bep, bmew, bsd},
    {clampd(bwc + 0.10, 0.0, 1.0), bep, bmew, bsd},
    {bwc, bep, std::max(0.01, bmew - 0.03), bsd},
    {bwc, bep, bmew + 0.03, bsd},
    {bwc, (bep == 1.0 ? 1.5 : 1.0), bmew, bsd},
    {bwc, bep, bmew, bsd * 2u + 1u},
    {bwc, bep, bmew, bsd * 7u + 13u},
    {clampd(bwc - 0.10, 0.0, 1.0), bep, bmew, bsd * 3u + 5u},
    {clampd(bwc + 0.10, 0.0, 1.0), bep, bmew, bsd * 5u + 9u},
    {bwc, bep, std::max(0.01, bmew - 0.03), bsd * 11u + 1u},
  };

  std::vector<std::vector<size_t>> runs;
  runs.reserve(cfgs.size() + 6);
  for (const Cfg &cf : cfgs) {
    g_w_comp = cf.wc; g_edge_power = cf.ep; min_edge_weight = (Similarity)cf.mew;
    compute_es();
    rebuild_incs();
    runs.push_back(run_lp(cf.sd));
  }
  // Restore baseline params before the jackknife (which keeps them fixed).
  g_w_comp = bwc; g_edge_power = bep; min_edge_weight = (Similarity)bmew;

  // ── Sample jackknife ──────────────────────────────────────────────────────
  // Drop one depth sample at a time, rebuild the depth-derived edge weights, and
  // re-bin.  Surfaces contigs whose assignment hinges on a single sample.  Needs
  // the raw per-sample depths (g_depth_raw, populated for the fuse/wjac metric);
  // rebuilds the rank-based corr view (g_depth_unit) + wjac norms over the kept
  // samples, then restores the originals.
  if (has_depth && num_depth_samples >= 4 && !g_depth_raw.empty()) {
    const size_t S = num_depth_samples;
    const size_t nJk = std::min<size_t>(S, 6);
    const std::vector<float>  sav_unit = g_depth_unit;
    const std::vector<float>  sav_raw  = g_depth_raw;
    const std::vector<double> sav_cn   = g_depth_colnorm;
    std::vector<double> rin, rout;
    for (size_t drop = 0; drop < nJk; ++drop) {
      const size_t nS = S - 1;
      std::vector<float> raw2((size_t)nobs * nS, 0.0f);
      std::vector<float> unit2((size_t)nobs * nS, 0.0f);
      rin.resize(nS);
      for (size_t r = 0; r < nobs; ++r) {
        size_t o = 0;
        for (size_t k = 0; k < S; ++k)
          if (k != drop) { float v = sav_raw[r * S + k]; raw2[r * nS + o] = v; rin[o] = v; ++o; }
        rank(rin.data(), (unsigned)nS, rout);   // per-row Spearman ranks
        double mean = 0.0; for (double x : rout) mean += x; mean /= (double)nS;
        double ss = 0.0; for (double x : rout) { double d = x - mean; ss += d * d; }
        if (ss > 0.0) {
          double iv = 1.0 / std::sqrt(ss);
          for (size_t k = 0; k < nS; ++k) unit2[r * nS + k] = (float)((rout[k] - mean) * iv);
        }
      }
      std::vector<double> cn2;
      if (g_depth_wjac_norm) {
        cn2.assign(nS, 1.0);
        for (size_t k = 0; k < nS; ++k) {
          double s = 0.0; for (size_t r = 0; r < nobs; ++r) s += raw2[r * nS + k];
          double m = s / std::max<size_t>(nobs, 1);
          cn2[k] = (m > 1e-12) ? 1.0 / m : 0.0;
        }
      }
      g_depth_unit = std::move(unit2);
      g_depth_raw  = std::move(raw2);
      g_depth_colnorm = std::move(cn2);
      num_depth_samples = nS;
      compute_es();
      rebuild_incs();
      runs.push_back(run_lp(bsd));
      g_depth_unit = sav_unit;     // restore for the next drop
      g_depth_raw  = sav_raw;
      g_depth_colnorm = sav_cn;
      num_depth_samples = S;
    }
  }

  g_cert_nruns = (int)runs.size();
  if (runs.empty()) return;

  // Baseline bins (binned contigs grouped by baseline LP label).
  std::unordered_map<size_t, std::vector<size_t>> binMembers;
  for (size_t i = 0; i < N; ++i)
    if (baseBinned[i]) binMembers[baseline[i]].push_back(i);

  std::vector<int> agree(N, 0), seen(N, 0);
  std::vector<std::unordered_map<size_t, int>> altTally(N);
  for (const auto &mem : runs) {
    // Per baseline bin: the run-cluster its members predominantly fall into.
    std::unordered_map<size_t, size_t> plurality;   // baseLabel -> runLabel
    std::unordered_map<size_t, size_t> owner;        // runLabel -> baseLabel
    for (auto &kv : binMembers) {
      std::unordered_map<size_t, int> cnt;
      int best = -1; size_t bestLab = 0;
      for (size_t i : kv.second) {
        int c = ++cnt[mem[i]];
        if (c > best) { best = c; bestLab = mem[i]; }
      }
      plurality[kv.first] = bestLab;
      // Whichever baseline bin contributes the most members to a run-cluster
      // "owns" it (for naming a contig's alternative destination).
      auto it = owner.find(bestLab);
      if (it == owner.end() || (int)binMembers[it->second].size() < (int)kv.second.size())
        owner[bestLab] = kv.first;
    }
    for (auto &kv : binMembers) {
      const size_t pl = plurality[kv.first];
      for (size_t i : kv.second) {
        seen[i]++;
        if (mem[i] == pl) agree[i]++;
        else {
          auto it = owner.find(mem[i]);
          if (it != owner.end() && it->second != kv.first)
            altTally[i][it->second]++;
        }
      }
    }
  }

  for (size_t i = 0; i < N; ++i) {
    if (!baseBinned[i] || seen[i] == 0) continue;
    g_cert_prob[i] = (float)agree[i] / (float)seen[i];
    int bestc = 0; size_t bestlab = SIZE_MAX;
    for (auto &kv : altTally[i])
      if (kv.second > bestc) { bestc = kv.second; bestlab = kv.first; }
    if (bestlab != SIZE_MAX) {
      g_cert_alt[i] = bestlab;
      g_cert_alt_prob[i] = (float)bestc / (float)seen[i];
    }
  }
}

// ── Certification output ─────────────────────────────────────────────────────
// Writes <prefix>.stable_members.tsv, .bin_stability.tsv, .ambiguous_contigs.tsv
// and (with --cert-core-fasta and sequences resident) a <prefix>.core_bins/
// directory holding each certified bin's stable-core contigs.  BinNum matches the
// final output numbering via label2binid (raw LP label -> sequential BinNum) and
// clsMap (contig -> rawLabel+1).  Called from output_bins once numbering is final.
static void rb_write_certification(
    const std::unordered_map<size_t, size_t> &label2binid,
    const std::vector<size_t> &clsMap) {
  const std::string sm = std::string(outFile) + ".stable_members.tsv";
  const std::string bs = std::string(outFile) + ".bin_stability.tsv";
  const std::string ac = std::string(outFile) + ".ambiguous_contigs.tsv";
  verbose_message("Writing MAG stability certification to: %s, %s, %s\n",
                  bs.c_str(), sm.c_str(), ac.c_str());
  std::ofstream osm(sm.c_str()), obs(bs.c_str()), oac(ac.c_str());
  osm << "ContigName\tBinNum\tLength\tAssignmentProb\tCore\tAmbiguous"
         "\tAlternativeBin\tAltProb\n";
  oac << "ContigName\tBinNum\tLength\tAssignmentProb\tAlternativeBin\tAltProb\n";
  obs << "BinNum\tTotalBp\tCoreBp\tCoreFraction\tMeanAssignmentProb\tCertified\n";

  struct Acc { size_t total = 0, core = 0; double probLw = 0.0; };
  std::map<size_t, Acc> perBin;
  std::map<size_t, std::vector<size_t>> coreContigs;   // binid -> core large idx

  for (size_t i = 0; i < nobs; ++i) {
    if (clsMap[i] == 0) continue;                       // unbinned large contig
    auto bit = label2binid.find(clsMap[i] - 1);
    if (bit == label2binid.end()) continue;
    const size_t binid = bit->second;
    const float p = (i < g_cert_prob.size()) ? g_cert_prob[i] : -1.0f;
    if (p < 0.0f) continue;                             // no stability score
    const size_t len = seq_lens[i];
    const bool core = p >= (float)g_cert_core_thr;
    const bool amb = p < (float)g_cert_amb_thr;
    long altBin = -1; float altP = 0.0f;
    if (i < g_cert_alt.size() && g_cert_alt[i] != SIZE_MAX) {
      auto ait = label2binid.find(g_cert_alt[i]);
      if (ait != label2binid.end()) altBin = (long)ait->second;
      altP = g_cert_alt_prob[i];
    }
    osm << contig_names[i] << "\t" << binid << "\t" << len << "\t"
        << std::fixed << std::setprecision(3) << p << "\t"
        << (core ? "Y" : "N") << "\t" << (amb ? "Y" : "N") << "\t"
        << altBin << "\t" << altP << "\n";
    osm.unsetf(std::ios::fixed);
    if (amb) {
      oac << contig_names[i] << "\t" << binid << "\t" << len << "\t"
          << std::fixed << std::setprecision(3) << p << "\t" << altBin
          << "\t" << altP << "\n";
      oac.unsetf(std::ios::fixed);
    }
    Acc &a = perBin[binid];
    a.total += len; a.probLw += (double)p * len;
    if (core) { a.core += len; coreContigs[binid].push_back(i); }
  }

  // Per-bin certification verdict + optional stable-core FASTA export.
  const bool want_fasta = g_cert_core_fasta && !seqs.empty();
  std::string coredir;
  if (want_fasta) {
    coredir = std::string(outFile) + ".core_bins";
    boost::system::error_code ec;
    boost::filesystem::create_directory(coredir, ec);
  }
  for (auto &kv : perBin) {
    const Acc &a = kv.second;
    double cf = a.total ? (double)a.core / (double)a.total : 0.0;
    double mp = a.total ? a.probLw / (double)a.total : 0.0;
    bool cert = (cf >= g_cert_frac_thr) && (mp >= g_cert_prob_thr);
    obs << kv.first << "\t" << a.total << "\t" << a.core << "\t"
        << std::fixed << std::setprecision(3) << cf << "\t" << mp
        << "\t" << (cert ? "Y" : "N") << "\n";
    obs.unsetf(std::ios::fixed);
    if (want_fasta && cert) {
      auto cit = coreContigs.find(kv.first);
      if (cit == coreContigs.end()) continue;
      std::string fp = coredir + "/bin." + std::to_string(kv.first) + ".core.fa";
      std::ofstream of(fp.c_str());
      for (size_t idx : cit->second)
        if (idx < seqs.size() && !seqs[idx].empty())
          printFasta(of, contig_names[idx], seqs[idx]);
    }
  }
  osm.flush(); obs.flush(); oac.flush();
}
