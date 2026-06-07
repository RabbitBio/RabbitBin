// RabbitBin module: rb_abundance.cpp

Distance cal_abd_dist2(Normal &p1, Normal &p2) {
  Distance k1, k2, tmp, d = 0;
  Distance m1 = p1.mean(), m2 = p2.mean();
  Distance v1 = p1.standard_deviation(); v1 *= v1;
  Distance v2 = p2.standard_deviation(); v2 *= v2;

  auto v1minusv2 = v1 - v2;
  if (FABS(v1minusv2) < 1e-4) {
    k1 = k2 = (m1 + m2) / 2;
  } else {
    auto m1minusm2 = m1 - m2;
    auto m1timesv2 = m1 * v2;
    auto m2timesv1 = m2 * v1;
    tmp = SQRT(v1 * v2 * (m1minusm2 * m1minusm2 -
                          2 * v1minusv2 * LOG(SQRT(v2 / v1))));
    k1 = (tmp - m1timesv2 + m2timesv1) / (v1minusv2);
    k2 = (tmp + m1timesv2 - m2timesv1) / (-v1minusv2);
  }
  if (k1 > k2) { tmp = k1; k1 = k2; k2 = tmp; }
  if (v1 > v2) std::swap(p1, p2);

  if (k1 == k2)
    d = FABS(boost::math::cdf(p1, k1) - boost::math::cdf(p2, k1));
  else
    d = FABS(boost::math::cdf(p1, k2) - boost::math::cdf(p1, k1) +
             boost::math::cdf(p2, k1) - boost::math::cdf(p2, k2));
  return d;
}

double cal_abd_dist(size_t r1, size_t r2, size_t i, bool &nz) {
  double d = 0;
  StoredDistance m1 = ABD(r1, i);
  StoredDistance m2 = ABD(r2, i);
  if (m1 > minCV || m2 > minCV) {
    nz = true;
    m1 = std::max(m1, (StoredDistance)1e-6);
    m2 = std::max(m2, (StoredDistance)1e-6);
    if (m1 != m2) {
      StoredDistance v1 = ABD_VAR(r1, i) < 1 ? 1 : ABD_VAR(r1, i);
      StoredDistance v2 = ABD_VAR(r2, i) < 1 ? 1 : ABD_VAR(r2, i);
      Normal p1(m1, SQRT(v1)), p2(m2, SQRT(v2));
      d = cal_abd_dist2(p1, p2);
    }
  }
  return std::min(std::max(d, (double)1e-6), (double)(1. - 1e-6));
}

size_t ncols(std::ifstream &is, int skip) {
  size_t nc = 0;
  std::string firstLine;
  while (skip-- >= 0) std::getline(is, firstLine);
  std::stringstream ss(firstLine);
  std::string col;
  while (std::getline(ss, col, tab_delim)) ++nc;
  return nc;
}

size_t ncols(const char *f, int skip) {
  std::ifstream is(f);
  if (!is.is_open()) {
    cerr << "[Error!] can't open input file " << f << "\n";
    return 0;
  }
  return ncols(is, skip);
}

std::istream &safeGetline(std::istream &is, std::string &t) {
  static int max_line = 1024;
  t.clear();
  t.reserve(max_line);
  std::istream::sentry se(is, true);
  std::streambuf *sb = is.rdbuf();
  for (;;) {
    int c = sb->sbumpc();
    switch (c) {
    case '\n': return is;
    case '\r':
      if (sb->sgetc() == '\n') sb->sbumpc();
      return is;
    case EOF:
      if (t.empty()) is.setstate(std::ios::eofbit);
      return is;
    default: t += (char)c;
    }
  }
  if (max_line < (int)t.size()) max_line = t.size() * 1.05 + 32;
  return is;
}

// Build g_anynz[i] = (any abundance sample of contig i exceeds minCV).
// O(nobs · nABD), run once; lets is_nz() avoid re-scanning per pair.
static void build_anynz_cache() {
  if (abdFile.empty()) return;
  g_anynz.assign(nobs, 0);
#pragma omp parallel for schedule(static) num_threads(numThreads)
  for (size_t r = 0; r < nobs; ++r) {
    uint8_t any = 0;
    for (size_t i = 0; i < nABD; ++i)
      if (ABD(r, i) > minCV) { any = 1; break; }
    g_anynz[r] = any;
  }
}

bool is_nz(size_t r1, size_t r2) {
  if (abdFile.empty()) return true;
  // Fast path: precomputed per-contig flags (built before the graph loop).
  if (g_anynz.size() == nobs) return g_anynz[r1] || g_anynz[r2];
  for (size_t i = 0; i < nABD; ++i) {
    if (ABD(r1, i) > minCV || ABD(r2, i) > minCV) return true;
  }
  return false;
}

template <typename D>
double _cal_abd_corr(size_t r1, size_t r2, const Matrix &ABD1,
                     const Matrix &ABD2) {
  assert(nABD > 1);
  double sum_xsq = 0.0, sum_ysq = 0.0, sum_cross = 0.0;
  D ratio, delta_x, delta_y;
  double mean_x = 0.0, mean_y = 0.0;
  D r = 0.0;

  for (int i = 0; i < (int)nABD; ++i) {
    D m1 = ABD1(r1, i);
    D m2 = ABD2(r2, i);
    if (i == 0) { mean_x = m1; mean_y = m2; continue; }
    D i_plus1 = i + 1;
    ratio = i / i_plus1;
    delta_x = m1 - mean_x;
    delta_y = m2 - mean_y;
    sum_xsq   += delta_x * delta_x * ratio;
    sum_ysq   += delta_y * delta_y * ratio;
    sum_cross += delta_x * delta_y * ratio;
    mean_x += delta_x / i_plus1;
    mean_y += delta_y / i_plus1;
  }
  if (sum_xsq <= 0.0 || sum_ysq <= 0.0) return 0.0;
  r = sum_cross / (sqrt(sum_xsq) * sqrt(sum_ysq));
  if (!std::isfinite(r)) return 0.0;
  return r;
}

double cal_abd_corr(size_t r1, size_t r2, bool second_is_small,
                    bool first_is_centroid) {
  return _cal_abd_corr<CALC_TYPE>(r1, r2,
                                  first_is_centroid ? ABD_centroids : ABD,
                                  second_is_small   ? small_ABD     : ABD);
}

void rescue_singletons(BinMap &cls) {
  verbose_message("There are %d bins already\n", cls.size());
  std::unordered_set<size_t> large_unbinned;
  for (auto i = 0; i < (int)nobs; i++)
    if (seq_lens[i] >= minClsSize) large_unbinned.insert(i);
  for (auto it = cls.begin(); it != cls.end(); ++it) {
    for (auto it2 = cls[it->first].begin(); it2 != cls[it->first].end(); ++it2)
      if (*it2 < nobs && seq_lens[*it2] >= minClsSize)
        large_unbinned.erase(*it2);
  }
  if (verbose && large_unbinned.size() > 0)
    verbose_message("Rescued %d large contig(s) into singleton bin(s)\n",
                    large_unbinned.size());
  for (auto id : large_unbinned) {
    assert(cls.find(id) == cls.end());
    cls[id].push_back(id);
  }
}


static std::string rb_bins_info_path() {
  return g_metabat_compat ? outFile + ".BinInfo.txt" : outFile + ".bins.tsv";
}
static std::string rb_members_path() {
  return g_metabat_compat ? outFile + ".BinMembers.txt" : outFile + ".members.tsv";
}
static std::string rb_matrix_path() {
  return g_metabat_compat ? outFile + ".MemberMatrix.txt" : outFile + ".members.matrix.tsv";
}
static std::string rb_bin_output_path(size_t bin_id) {
  if (g_metabat_compat) {
    std::string p = outFile + "." + boost::lexical_cast<std::string>(bin_id);
    if (!onlyLabel) p += ".fa";
    return p;
  }
  char buf[512];
  std::snprintf(buf, sizeof(buf), "%s_bin_%03zu%s", outFile.c_str(), bin_id, onlyLabel ? "" : ".fa");
  return buf;
}
static std::string rb_unbinned_path() {
  std::string p = outFile + ".unbinned";
  if (!onlyLabel) p += ".fa";
  return p;
}

