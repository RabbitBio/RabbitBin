// RabbitBin module: rb_strain.cpp
//
// Strain resolution (feature #5) post-processing.  The same-pass SNV extraction
// lives in rabbit_depth.cpp; here we (a) project the depth-pass SNV grid onto the
// final large-contig rows and (b) estimate, per bin, the number of strains and
// the dominant minor strain's per-sample abundance from per-site allele-frequency
// linkage.  #included into the rabbitbin.cpp translation unit (shares its globals).

// ── Strain / SNV (feature #5): materialize per-large-contig row arrays ───────
// Maps the depth-pass SNV grid (keyed by BAM tid / reference name) onto the
// final large-contig row order (contig_names[idx] aligns with depth_matrix(idx,*)
// and g_node_confidence[idx]).  Runs once, after the depth merge has compacted
// contig_names, so the resulting g_snv_* arrays survive every downstream stage.
static void materialize_snv_rows() {
  const size_t S = num_depth_samples;
  g_snv_pi.assign((size_t)nobs * S, 0.0f);
  g_snv_nsnv.assign((size_t)nobs * S, 0u);
  g_snv_cov.assign((size_t)nobs * S, 0u);
  if (!g_snv_result.enabled || g_snv_result.n_targets == 0 ||
      g_snv_result.num_bams == 0)
    return;
  const int32_t nt = g_snv_result.n_targets;
  const int nb = g_snv_result.num_bams;
  const bool haveSites = g_snv_result.resolve && !g_snv_result.sites.empty();
  if (haveSites)
    g_snv_sites.assign((size_t)nobs * S, std::vector<SnvSite>{});
  std::unordered_map<std::string, int32_t> name2tid;
  name2tid.reserve((size_t)nt * 2);
  for (int32_t t = 0; t < nt; ++t)
    name2tid.emplace(g_snv_result.names[t], t);
  const int useB = std::min<int>(nb, (int)S);  // sample order = BAM order
  for (size_t idx = 0; idx < nobs; ++idx) {
    auto it = name2tid.find(contig_names[idx]);
    if (it == name2tid.end())
      continue;
    const int32_t tid = it->second;
    for (int b = 0; b < useB; ++b) {
      const SnvContigStat &st =
          g_snv_result.stats[(size_t)b * nt + tid];
      const size_t o = idx * S + (size_t)b;
      g_snv_cov[o] = st.covered;
      g_snv_nsnv[o] = st.nsnv;
      g_snv_pi[o] = st.covered ? (float)(st.pi_sum / (double)st.covered) : 0.0f;
      if (haveSites)
        g_snv_sites[o] = std::move(g_snv_result.sites[(size_t)b * nt + tid]);
    }
  }
}

// ── Strain resolution (feature #5 iteration 2): per-bin strain count + abundance
// Estimate how many strains a bin contains and the dominant minor strain's
// per-sample relative abundance, from the per-site allele-frequency linkage.
// Sites whose per-sample alt-allele frequency trajectories point the same way
// belong to one strain pair; counting coherent trajectory clusters (cosine on
// mean-centred, sign-folded trajectories) yields strain count = 1 + #clusters.
// Returns 1 when polymorphism is insufficient.  minorAbund (length S) receives
// the largest cluster's mean per-sample minor-allele frequency (~the minor
// strain's relative abundance); empty/zero when EstStrains == 1.
static int estimate_bin_strains(const ContigVector &cluster,
                                std::vector<double> &minorAbund) {
  const size_t S = num_depth_samples;
  minorAbund.assign(S, 0.0);
  if (S < 2 || g_snv_sites.empty())
    return 1;   // linkage needs >= 2 samples and per-site data (BAM run)

  // Gather raw per-(physical site) entries: one SnvSite per reporting sample.
  struct Raw { std::vector<int> samp; std::vector<SnvSite> st; };
  std::unordered_map<uint64_t, Raw> raw;
  for (size_t ci = 0; ci < cluster.size(); ++ci) {
    const size_t idx = cluster[ci];
    if (idx >= nobs) continue;                    // small recruit: no SNV row
    for (size_t b = 0; b < S; ++b) {
      const std::vector<SnvSite> &lst = g_snv_sites[idx * S + b];
      for (const SnvSite &s : lst) {
        Raw &r = raw[((uint64_t)idx << 32) | s.pos];
        r.samp.push_back((int)b);
        r.st.push_back(s);
      }
    }
  }
  if (raw.empty()) return 1;

  // Build oriented, mean-centred, unit trajectories for sites with enough cover.
  const int minReport = std::max<int>(2, (int)((S + 1) / 2));   // >= half samples
  std::vector<std::vector<double>> dirs;     // unit, mean-centred (for clustering)
  std::vector<std::vector<double>> rawf;     // oriented alt-freq (for abundance)
  dirs.reserve(raw.size());
  rawf.reserve(raw.size());
  for (auto &kv : raw) {
    Raw &r = kv.second;
    if ((int)r.samp.size() < minReport) continue;
    // Canonical alt = alt nucleotide of the deepest-covered reporting sample.
    int deep = 0;
    for (size_t i = 1; i < r.st.size(); ++i)
      if (r.st[i].total > r.st[deep].total) deep = (int)i;
    const uint8_t canon = r.st[deep].alt;
    std::vector<double> f(S, -1.0);            // -1 = missing
    for (size_t i = 0; i < r.samp.size(); ++i) {
      const SnvSite &s = r.st[i];
      double fr = s.total ? (double)s.altCount / (double)s.total : 0.0;
      double v;
      if (s.alt == canon)      v = fr;          // freq of canonical alt
      else if (s.maj == canon) v = 1.0 - fr;    // canonical allele is the major
      else                     v = 0.0;         // a third allele: alt ~ absent
      f[r.samp[i]] = v;
    }
    // Impute missing samples with the reporting-sample mean (neutral direction).
    double sum = 0.0; int nrep = 0;
    for (size_t b = 0; b < S; ++b) if (f[b] >= 0.0) { sum += f[b]; nrep++; }
    if (nrep == 0) continue;
    double mean = sum / nrep;
    for (size_t b = 0; b < S; ++b) if (f[b] < 0.0) f[b] = mean;
    // Mean-centre + L2-normalise -> trajectory direction.
    std::vector<double> d(S);
    double dm = 0.0; for (double x : f) dm += x; dm /= S;
    double nn = 0.0; for (size_t b = 0; b < S; ++b) { d[b] = f[b] - dm; nn += d[b]*d[b]; }
    if (nn < 1e-9) continue;                     // flat: no across-sample signal
    double inv = 1.0 / std::sqrt(nn);
    for (size_t b = 0; b < S; ++b) d[b] *= inv;
    dirs.push_back(std::move(d));
    rawf.push_back(std::move(f));
  }
  const int nUsable = (int)dirs.size();
  if (nUsable < g_strain_min_sites) return 1;

  // Greedy sign-folded cosine clustering of trajectory directions.  Two distinct
  // strain trajectories over S samples are near-orthogonal, while binomial noise
  // around ONE direction stays well-correlated, so a moderate threshold both
  // merges within-strain noise and separates genuine strains.
  const double COS_THR = 0.75;
  auto absdot = [&](const std::vector<double> &x, const std::vector<double> &y) {
    double d = 0.0; for (size_t b = 0; b < S; ++b) d += x[b] * y[b];
    return d;   // signed; caller takes |.| as needed
  };
  std::vector<std::vector<double>> cent;        // cluster centroids (unit)
  std::vector<std::vector<int>> members;        // site index, sign-encoded
  for (int i = 0; i < nUsable; ++i) {
    int best = -1; double bestc = COS_THR; int bestSign = 1;
    for (size_t c = 0; c < cent.size(); ++c) {
      double dot = absdot(dirs[i], cent[c]);
      if (std::fabs(dot) >= bestc) { bestc = std::fabs(dot); best = (int)c; bestSign = (dot < 0) ? -1 : 1; }
    }
    if (best < 0) { cent.push_back(dirs[i]); members.push_back({i}); }
    else {
      members[best].push_back(bestSign > 0 ? i : -(i + 1));
      std::vector<double> &cc = cent[best];
      double nn = 0.0;
      for (size_t b = 0; b < S; ++b) { cc[b] += bestSign * dirs[i][b]; nn += cc[b] * cc[b]; }
      if (nn > 1e-12) { double inv = 1.0/std::sqrt(nn); for (double &x : cc) x *= inv; }
    }
  }
  // Merge clusters whose centroids ended up near-collinear (greedy drift can
  // spawn duplicates of the same direction).
  for (size_t a = 0; a < cent.size(); ++a) {
    if (members[a].empty()) continue;
    for (size_t c = a + 1; c < cent.size(); ++c) {
      if (members[c].empty()) continue;
      double dot = absdot(cent[a], cent[c]);
      if (std::fabs(dot) >= COS_THR) {
        const int sgn = (dot < 0) ? -1 : 1;
        for (int m : members[c]) {
          const int i = (m >= 0) ? m : (-m - 1);
          const int s = ((m >= 0) ? 1 : -1) * sgn;
          members[a].push_back(s > 0 ? i : -(i + 1));
        }
        members[c].clear();
      }
    }
  }

  // Keep clusters with enough support; each coherent minor direction is a strain.
  const int minClusterSites =
      std::max(g_strain_min_sites, (int)std::lround(0.25 * nUsable));
  int kept = 0, bestCluster = -1; size_t bestSize = 0;
  for (size_t c = 0; c < members.size(); ++c) {
    if ((int)members[c].size() < minClusterSites) continue;
    kept++;
    if (members[c].size() > bestSize) { bestSize = members[c].size(); bestCluster = (int)c; }
  }
  if (kept == 0 || bestCluster < 0) return 1;

  // Dominant minor strain's per-sample abundance = mean oriented alt-freq over
  // the largest cluster's sites (sign-aligned to the centroid).
  std::vector<double> acc(S, 0.0); int na = 0;
  for (int m : members[bestCluster]) {
    const int i = (m >= 0) ? m : (-m - 1);
    const int sign = (m >= 0) ? 1 : -1;
    for (size_t b = 0; b < S; ++b)
      acc[b] += (sign > 0) ? rawf[i][b] : (1.0 - rawf[i][b]);
    na++;
  }
  if (na > 0) for (size_t b = 0; b < S; ++b) acc[b] /= na;
  // Report the minor (rarer) side so the value is a minor-strain fraction.
  double am = 0.0; for (double x : acc) am += x; am /= S;
  if (am > 0.5) for (size_t b = 0; b < S; ++b) acc[b] = 1.0 - acc[b];
  minorAbund = std::move(acc);
  return std::min(1 + kept, g_strain_max_k);
}
