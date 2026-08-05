// RabbitBin module: rb_cluster.cpp

// Closed-form CDF of a chi-squared distribution with 2·m degrees of freedom,
// evaluated at x (here m = neighbour count, x = -2·Σlog(1-sim)).  For integer
// shape the regularized lower incomplete gamma reduces to a finite Poisson sum:
//   cdf(2m, x) = 1 − e^{−u}·Σ_{i=0}^{m−1} u^i/i!,   u = x/2.
// Accumulating the Poisson terms p_i = e^{−u}·u^i/i! (p_0=e^{−u}, p_i=p_{i-1}·u/i)
// keeps every term in [0,1] so there is no overflow for large u, and replaces
// boost's general incomplete-gamma evaluation (called tens of millions of times
// in label propagation) with a handful of multiplies.  Mathematically identical
// to boost::math::cdf(chi_squared(2m), x) up to floating-point ULPs.
static inline double chi2_2dof_sf(int m, double x) {   // upper tail Q
  const double u = 0.5 * x;
  double p = std::exp(-u);   // p_0
  double Q = p;
  for (int i = 1; i < m; ++i) { p *= u / (double)i; Q += p; }
  return Q;
}
static inline double chi2_2dof_cdf(int m, double x) {
  double v = 1.0 - chi2_2dof_sf(m, x);
  return v < 0.0 ? 0.0 : v;
}

// ── Label-score rule (ablation switch) ───────────────────────────────────────
// Propagation must reduce the edges joining v to one candidate label into a
// single comparable number.  The default reduction is a Fisher-style combination
// of the per-edge quantities log(1-w), read through a chi-squared CDF -- which
// treats w as if 1-w were a calibrated, independent p-value.  It is not: w is a
// Spearman correlation already thresholded at --min-edge-score, and incident
// edges are plainly dependent.  RABBIT_LPA_SCORE swaps the reduction so the
// aggregation rule can be ablated against an otherwise identical graph:
//   fisher (default)  chi2_{2m} CDF of -2*Sum log(1-w)   (saturates at 1.0)
//   tail              same ordering via the upper tail Q, which does not
//                     saturate: 1-Q is indistinguishable from 1.0 once Q drops
//                     below ~1e-16, while Q itself stays comparable to ~1e-308
//   logsum            Sum of -log(1-w)                   (same evidence, no CDF)
//   sum               Sum of w                            (standard weighted LPA)
//   mean              Sum of w / m
//   max               max w
enum class LpaScore { Fisher, Tail, LogSum, Sum, Mean, Max };
static LpaScore rb_lpa_score_init() {
  const char *e = getenv("RABBIT_LPA_SCORE");
  if (!e || !*e) return LpaScore::Fisher;
  if (std::strcmp(e, "sum") == 0)    return LpaScore::Sum;
  if (std::strcmp(e, "mean") == 0)   return LpaScore::Mean;
  if (std::strcmp(e, "max") == 0)    return LpaScore::Max;
  if (std::strcmp(e, "logsum") == 0) return LpaScore::LogSum;
  if (std::strcmp(e, "tail") == 0)   return LpaScore::Tail;
  return LpaScore::Fisher;
}
static const LpaScore g_lpa_score = rb_lpa_score_init();
// sum/mean/max reduce the raw weight; fisher/logsum reduce log(1-w).
static const bool g_lpa_raw_w = g_lpa_score == LpaScore::Sum ||
                                g_lpa_score == LpaScore::Mean ||
                                g_lpa_score == LpaScore::Max;
static inline void lpa_accumulate(double &dst, double x) {
  if (g_lpa_score == LpaScore::Max) { if (x > dst) dst = x; }
  else dst += x;
}
static inline double lpa_combine(int m, double acc) {
  switch (g_lpa_score) {
    case LpaScore::Fisher: return chi2_2dof_cdf(m, -2.0 * acc);
    // Maximising -Q is the same ordering as maximising the CDF 1-Q, but stays
    // discriminative in the range where 1-Q has already rounded to 1.0.
    case LpaScore::Tail:   return -chi2_2dof_sf(m, -2.0 * acc);
    case LpaScore::LogSum: return -acc;                      // acc <= 0
    case LpaScore::Mean:   return m > 0 ? acc / (double)m : 0.0;
    default:               return acc;                       // Sum, Max
  }
}

// Retain the current label when its score is tied with the best candidate up
// to floating-point noise.  Besides making ties deterministic, this prevents
// score-equivalent labels from generating spurious updates and lets the main
// propagation loop stop on a full sweep with zero actual label changes.
static constexpr double LPA_SCORE_EPS = 1e-12;
// RB_LPA_PROF: how often does the chi-squared CDF comparison actually lose
// information?  Counts decisions whose winning score has already rounded to 1.0,
// and among those, decisions where two or more candidate labels are tied there --
// those retain the current label when available, otherwise incidence-order
// first-max breaks the tie.
static const bool g_lpa_prof = getenv("RB_LPA_PROF") != nullptr;
static std::array<std::atomic<uint64_t>, 3> g_lpa_dec;  // total, saturated, tied
static const char *lpa_score_name() {
  switch (g_lpa_score) {
    case LpaScore::Tail:   return "tail(-chi2 sf)";
    case LpaScore::LogSum: return "logsum(-log(1-w))";
    case LpaScore::Sum:    return "sum(w)";
    case LpaScore::Mean:   return "mean(w)";
    case LpaScore::Max:    return "max(w)";
    default:               return "fisher(chi2 cdf)";
  }
}

// ── Alternative clustering: weighted multilevel modularity (Louvain) ─────────
// Edge weight w_e = -log(1 - edgeScore[e]) — the same per-edge "evidence" the
// Fisher label-propagation sums.  Modularity with resolution gamma: higher gamma
// yields finer communities, which matters for metagenome binning (many small,
// pure genomes) where standard modularity's resolution limit would over-merge.
// Selected via RABBIT_CLUSTER=louvain ; gamma via RABBIT_CLUSTER_GAMMA (def 1.0).
namespace rb_louvain {
struct AdjGraph {
  std::vector<std::vector<std::pair<int, double>>> adj; // (neighbor, weight)
  std::vector<double> selfloop;                         // self-loop weight
  std::vector<double> k;                                // weighted degree (incl 2*selfloop)
  double m2 = 0.0;                                      // 2m = Σ k
};

// One Louvain level: local moving on G. Fills comm (contiguous ids) and ncomm;
// returns true if any node moved.
static bool local_move(const AdjGraph &G, double gamma, std::mt19937 &rng,
                       std::vector<int> &comm, int &ncomm) {
  const int N = (int)G.adj.size();
  comm.resize(N);
  std::iota(comm.begin(), comm.end(), 0);
  std::vector<double> ktot(G.k);
  std::vector<int> order(N);
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), rng);

  const double inv2m = (G.m2 > 0.0) ? 1.0 / G.m2 : 0.0;
  std::vector<double> wToComm(N, 0.0);
  std::vector<int> touched;
  touched.reserve(64);

  bool any = false, improved = true;
  int passes = 0;
  while (improved && passes++ < 100) {
    improved = false;
    for (int idx = 0; idx < N; ++idx) {
      const int i = order[idx];
      const int ci = comm[i];
      touched.clear();
      for (const auto &pr : G.adj[i]) {
        int c = comm[pr.first];
        if (wToComm[c] == 0.0) touched.push_back(c);
        wToComm[c] += pr.second;
      }
      const double ki = G.k[i];
      ktot[ci] -= ki; // tentatively remove i from its community
      int bestc = ci;
      double bestgain = wToComm[ci] - gamma * ki * ktot[ci] * inv2m;
      for (int c : touched) {
        double gain = wToComm[c] - gamma * ki * ktot[c] * inv2m;
        if (gain > bestgain + 1e-12) { bestgain = gain; bestc = c; }
      }
      ktot[bestc] += ki;
      if (bestc != ci) { comm[i] = bestc; improved = true; any = true; }
      for (int c : touched) wToComm[c] = 0.0;
    }
  }
  std::vector<int> remap(N, -1);
  ncomm = 0;
  for (int i = 0; i < N; ++i) {
    if (remap[comm[i]] < 0) remap[comm[i]] = ncomm++;
    comm[i] = remap[comm[i]];
  }
  return any;
}

// Aggregate G into a super-graph (one node per community).
static AdjGraph aggregate(const AdjGraph &G, const std::vector<int> &comm, int ncomm) {
  AdjGraph H;
  H.adj.assign(ncomm, {});
  H.selfloop.assign(ncomm, 0.0);
  std::vector<phmap::flat_hash_map<int, double>> acc(ncomm);
  for (int i = 0; i < (int)G.adj.size(); ++i) {
    int ci = comm[i];
    H.selfloop[ci] += G.selfloop[i];
    for (const auto &pr : G.adj[i]) {
      int cj = comm[pr.first];
      if (ci == cj) H.selfloop[ci] += pr.second * 0.5; // each intra edge seen twice
      else if (ci < cj) acc[ci][cj] += pr.second;
    }
  }
  for (int c = 0; c < ncomm; ++c)
    for (const auto &kv : acc[c]) {
      H.adj[c].push_back({kv.first, kv.second});
      H.adj[kv.first].push_back({c, kv.second});
    }
  H.k.assign(ncomm, 0.0);
  H.m2 = 0.0;
  for (int c = 0; c < ncomm; ++c) {
    double kc = 2.0 * H.selfloop[c];
    for (const auto &pr : H.adj[c]) kc += pr.second;
    H.k[c] = kc;
    H.m2 += kc;
  }
  return H;
}
} // namespace rb_louvain

int cluster_by_louvain(Graph &g, std::vector<size_t> &membership, unsigned seed,
                       double gamma) {
  using namespace rb_louvain;
  const size_t n = g.getNodeCount();
  const size_t E = g.getEdgeCount();
  membership.resize(n);
  std::iota(membership.begin(), membership.end(), 0);
  if (n == 0 || E == 0) return 0;

  AdjGraph G;
  G.adj.assign(n, {});
  G.selfloop.assign(n, 0.0);
  for (size_t e = 0; e < E; ++e) {
    double s = g.edgeScore[e];
    if (s <= 0.0) continue;
    if (s > 1.0 - 1e-6) s = 1.0 - 1e-6;
    double w = -std::log(1.0 - s);
    size_t a = g.from[e], b = g.to[e];
    if (a == b) { G.selfloop[a] += w; continue; }
    G.adj[a].push_back({(int)b, w});
    G.adj[b].push_back({(int)a, w});
  }
  G.k.assign(n, 0.0);
  G.m2 = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double ki = 2.0 * G.selfloop[i];
    for (const auto &pr : G.adj[i]) ki += pr.second;
    G.k[i] = ki;
    G.m2 += ki;
  }

  std::mt19937 rng(seed);
  std::vector<int> cur(n);
  std::iota(cur.begin(), cur.end(), 0);
  for (int level = 0; level < 20; ++level) {
    std::vector<int> comm;
    int ncomm = 0;
    bool moved = local_move(G, gamma, rng, comm, ncomm);
    for (size_t i = 0; i < n; ++i) cur[i] = comm[cur[i]];
    if (!moved || ncomm == (int)G.adj.size()) break;
    G = aggregate(G, comm, ncomm);
  }
  for (size_t i = 0; i < n; ++i) membership[i] = (size_t)cur[i];
  return 0;
}

// chi-squared (2m dof) CDF with a fast path for large m: the exact Poisson sum
// is O(m), which is wasteful for the aggregated super-edges of multilevel
// merging (m = #edges between two communities, can be large).  For m>256 use the
// Wilson–Hilferty normal approximation (well within ULPs of the exact value at
// that scale, and only used for the merge tie-break).
static inline double chi2_2dof_cdf_fast(int m, double x) {
  if (m <= 256) return chi2_2dof_cdf(m, x);
  const double k = 2.0 * (double)m;
  const double t = std::cbrt(x / k);
  const double mu = 1.0 - 2.0 / (9.0 * k);
  const double sd = std::sqrt(2.0 / (9.0 * k));
  const double z = (t - mu) / sd;
  return 0.5 * std::erfc(-z / std::sqrt(2.0));
}

// Multilevel Fisher-LPA, merge phase.  After the (proven) level-0 LPA converges,
// treat each community as a super-node and CONSERVATIVELY merge super-nodes whose
// inter-community edges are collectively very significant — i.e. fragments of the
// same genome that LPA left split.  Unlike modularity (which over-merges via the
// null model), the criterion is the SAME Fisher statistic the LPA uses, gated by
// a high confidence threshold tau and a minimum supporting-edge count so only
// near-certain same-genome merges happen (no collapse).  Modifies `membership`.
static void fisher_super_merge(Graph &g, std::vector<size_t> &membership,
                               const std::vector<StoredDistance> &logsscr,
                               double tau, int min_count) {
  const size_t n = g.getNodeCount();
  const size_t E = g.getEdgeCount();
  // Relabel level-0 communities to a contiguous [0, C).
  phmap::flat_hash_map<size_t, int> remap;
  std::vector<int> comm0(n);
  int C = 0;
  for (size_t i = 0; i < n; ++i) {
    auto it = remap.find(membership[i]);
    if (it == remap.end()) { comm0[i] = C; remap.emplace(membership[i], C++); }
    else comm0[i] = it->second;
  }
  if (C <= 1) return;
  // Aggregate inter-community super-edges: (Σlog, count), additive.
  std::vector<phmap::flat_hash_map<int, std::pair<double, int>>> sadj(C);
  for (size_t e = 0; e < E; ++e) {
    if (g.edgeScore[e] <= 0) continue;
    int a = comm0[g.from[e]], b = comm0[g.to[e]];
    if (a == b) continue;
    double sl = (double)logsscr[e];
    auto &ab = sadj[a][b]; ab.first += sl; ab.second += 1;
    auto &ba = sadj[b][a]; ba.first += sl; ba.second += 1;
  }
  // Gauss-Seidel Fisher local-move on the super-graph, gated by tau/min_count so
  // a super-node only joins a neighbour when that merge is near-certain.
  std::vector<int> comm(C);
  std::iota(comm.begin(), comm.end(), 0);
  std::vector<double> sscore(C, 0.0);
  std::vector<int> scount(C, 0);
  std::vector<int> touched; touched.reserve(64);
  std::unordered_map<int, std::unordered_set<int>> visited;
  std::unordered_set<int> blacklist;
  std::vector<char> active(C, 1);
  while (true) {
    size_t changed = 0;
    for (int v = 0; v < C; ++v) {
      if (!active[v]) continue;
      active[v] = 0;
      touched.clear();
      for (const auto &kv : sadj[v]) {
        int k = comm[kv.first];
        if (k == comm[v]) continue;          // self community: no merge target
        if (scount[k] == 0) touched.push_back(k);
        sscore[k] += kv.second.first;
        scount[k] += kv.second.second;
      }
      if (touched.empty()) continue;
      double best_val = -std::numeric_limits<double>::infinity();
      int best_k = -1;
      for (int k : touched) {
        double val = chi2_2dof_cdf_fast(scount[k], -2.0 * sscore[k]);
        // Conservative gate: only consider near-certain, well-supported merges.
        if (val >= tau && scount[k] >= min_count && val > best_val) {
          best_val = val; best_k = k;
        }
        sscore[k] = 0.0; scount[k] = 0;
      }
      if (best_k < 0) continue;              // nothing confident enough → stay
      int prev = comm[v];
      if (prev != best_k && blacklist.find(v) == blacklist.end()) {
        comm[v] = best_k;
        ++changed;
        for (const auto &kv : sadj[v]) active[kv.first] = 1;
        auto &seen = visited[v];
        if (seen.empty()) seen.insert(prev);
        if (!seen.insert(best_k).second) blacklist.insert(v);
      }
    }
    if (changed == 0) break;
  }
  for (size_t i = 0; i < n; ++i) membership[i] = (size_t)comm[comm0[i]];
}

// Canonical CSR view used by the label-propagation hot loop.  Graph::incs stores
// an edge id for every incidence, which makes each visit chase that id through
// from[]/to[] (plus a branch in getOtherNode) and then through logsscr[].  Flatten
// the already neighbour-sorted incidence lists once into an 8-byte, contiguous
// {neighbour, weight} stream.  Iterating each row in the same order preserves the
// exact floating-point accumulation and first-max tie-break of the legacy path.
// RABBIT_LPA_LEGACY=1 retains the old traversal for direct A/B verification.
struct LpaArc {
  uint32_t neighbor;
  StoredDistance weight;
};
static_assert(sizeof(LpaArc) == 8, "LPA CSR arc must stay compact");

struct LpaCsr {
  std::vector<size_t> offsets;
  std::vector<LpaArc> arcs;
};

static LpaCsr build_lpa_csr(Graph &g,
                            const std::vector<StoredDistance> &edge_weight) {
  const size_t n = g.getNodeCount();
  LpaCsr csr;
  csr.offsets.resize(n + 1, 0);
  for (size_t v = 0; v < n; ++v)
    csr.offsets[v + 1] = csr.offsets[v] + g.incs[v].size();
  csr.arcs.resize(csr.offsets[n]);

#pragma omp parallel for schedule(static)
  for (size_t v = 0; v < n; ++v) {
    size_t out = csr.offsets[v];
    const std::vector<GraphEdgeId> &inc = g.incs[v];
    for (size_t j = 0; j < inc.size(); ++j) {
      const size_t edge = inc[j];
      csr.arcs[out++] = LpaArc{(uint32_t)g.getOtherNode(edge, v),
                               edge_weight[edge]};
    }
  }
  return csr;
}

int cluster_by_propagation(Graph &g, std::vector<size_t> &membership,
                      std::vector<size_t> &node_order) {
  size_t no_of_nodes = g.getNodeCount();
  size_t no_of_edges = g.getEdgeCount();

  if (no_of_nodes == 0 || no_of_edges == 0) {
    cerr << "There were " << no_of_nodes << " nodes and " << no_of_edges
         << " edges -- skipping cluster_by_propagation" << endl;
    return 0;
  }

  // Multilevel Fisher-LPA: run the standard level-0 LPA below (identity init,
  // byte-identical to the baseline), then merge communities in a super-graph.
  bool do_mlpa = false;
  double mlpa_tau = 1.0 - 1e-9;
  int mlpa_min_count = 3;

  // Algorithm dispatch (LPA stays the default).  RABBIT_CLUSTER=louvain selects
  // weighted multilevel modularity; the rng seed is derived from node_order so a
  // fixed --seed stays deterministic.
  if (const char *alg = getenv("RABBIT_CLUSTER")) {
    double gamma = 8.0;
    if (const char *gg = getenv("RABBIT_CLUSTER_GAMMA")) {
      double v = atof(gg);
      if (v > 0) gamma = v;
    }
    unsigned lseed = 0x9e3779b9u ^ (unsigned)node_order.size();
    if (!node_order.empty())
      lseed ^= (unsigned)(node_order[0] * 2654435761u);
    if (std::strcmp(alg, "louvain") == 0) {
      // Pure modularity Louvain.
      return cluster_by_louvain(g, membership, lseed, gamma);
    }
    if (std::strcmp(alg, "fuse") == 0) {
      // Fusion: Louvain global structure as the LPA initial partition, then the
      // Fisher label-propagation below refines it (membership is left sized to
      // n, so the LPA loop uses it as the seed instead of the identity).
      cluster_by_louvain(g, membership, lseed, gamma);
      // fall through to LPA refinement
    }
    if (std::strcmp(alg, "mlpa") == 0) {
      // Multilevel Fisher-LPA: standard LPA (below) + conservative super-merge.
      do_mlpa = true;
      if (const char *t = getenv("RABBIT_MLPA_TAU")) {
        double v = atof(t);
        if (v > 0 && v < 1) mlpa_tau = v;
      }
      if (const char *c = getenv("RABBIT_MLPA_MINCOUNT")) {
        int v = atoi(c);
        if (v > 0) mlpa_min_count = v;
      }
    }
  }

  if (g.edgeScore.size() != no_of_edges) {
    cerr << "edgeScore != no_of_edges" << endl;
    exit(1);
  }

  if (membership.size() != no_of_nodes) {
    membership.resize(no_of_nodes);
    std::iota(membership.begin(), membership.end(), 0);
  }

  /* Do some initial checks */
  if (*std::min_element(g.edgeScore.begin(), g.edgeScore.end()) < 0) {
    cerr << "edgeScore must be non-negative" << endl;
    exit(1);
  }

  // Precompute the per-edge quantity the propagation accumulates: LOG(1 - score)
  // for the fisher/logsum rules, the raw score for sum/mean/max.  It is invariant
  // across the (many) label-propagation rounds and across every incidence visit
  // of the edge, so hoisting it out of the hot loop removes ~10²·|E| redundant
  // log() calls.
  std::vector<StoredDistance> logsscr(no_of_edges);
  for (size_t e = 0; e < no_of_edges; ++e)
    logsscr[e] = g_lpa_raw_w ? g.edgeScore[e]
                             : (StoredDistance)LOG(1. - g.edgeScore[e]);
  if (g_lpa_score != LpaScore::Fisher)
    cerr << "[LPA] label score = " << lpa_score_name() << endl;

  const bool use_csr = no_of_nodes <= (size_t)std::numeric_limits<uint32_t>::max() &&
                       getenv("RABBIT_LPA_LEGACY") == nullptr;
  std::chrono::steady_clock::time_point csr_t0;
  if (getenv("RB_LPA_CSR_PROF")) csr_t0 = std::chrono::steady_clock::now();
  LpaCsr lpa_csr;
  if (use_csr) lpa_csr = build_lpa_csr(g, logsscr);
  if (getenv("RB_LPA_CSR_PROF")) {
    const double ms = std::chrono::duration<double, std::milli>(
                          std::chrono::steady_clock::now() - csr_t0).count();
    fprintf(stderr, "[RB_LPA_CSR] %s arcs=%zu bytes=%zu build=%.1f ms\n",
            use_csr ? "on" : "off", lpa_csr.arcs.size(),
            lpa_csr.arcs.size() * sizeof(LpaArc) +
                lpa_csr.offsets.size() * sizeof(size_t),
            ms);
  }

  const bool lpa_trace = getenv("RB_LPA_TRACE") != nullptr;

  // ── Optional parallel (Jacobi) label propagation ─────────────────────────
  // RABBIT_PAR_LPA=1 switches from the default order-dependent (Gauss-Seidel)
  // sweep to a synchronous one: every active node recomputes its best cluster
  // from the membership SNAPSHOT taken at the start of the round (read-only, so
  // the per-node work parallelises trivially), then all changes are committed
  // together in node_order.  Because neighbour updates are no longer visible
  // mid-round, the converged labelling differs from the sequential baseline,
  // so this path is opt-in. The
  // commit is serial and in node_order, so a given thread count + flag is fully
  // deterministic. Both paths retain the current label on score ties and stop
  // after a complete round with zero actual label changes.
  if (getenv("RABBIT_PAR_LPA") != nullptr) {
    const int nthreads = omp_get_max_threads();
    std::vector<std::vector<double>>  tl_nscore(nthreads, std::vector<double>(no_of_nodes, 0.0));
    std::vector<std::vector<int32_t>> tl_ncount(nthreads, std::vector<int32_t>(no_of_nodes, 0));

    std::vector<size_t> prop_k(no_of_nodes);     // proposed best cluster per node
    std::vector<char>   prop_set(no_of_nodes, 0); // 1 if node was evaluated this round
    std::vector<char>   active(no_of_nodes, 1), next_active(no_of_nodes, 0);
    std::unordered_map<size_t, std::unordered_set<size_t>> visited;
    std::unordered_set<size_t> blacklist;

    size_t round = 0;
    while (true) {
      ++round;

#pragma omp parallel
      {
        const int tid = omp_get_thread_num();
        std::vector<double>  &nscore = tl_nscore[tid];
        std::vector<int32_t> &ncount = tl_ncount[tid];
        std::vector<size_t>   touched;
        touched.reserve(256);
#pragma omp for schedule(dynamic, 256)
        for (size_t i = 0; i < node_order.size(); i++) {
          size_t v1 = node_order[i];
          if (!active[v1]) { prop_set[v1] = 0; continue; }
          touched.clear();
          if (use_csr) {
            const LpaArc *arc = lpa_csr.arcs.data() + lpa_csr.offsets[v1];
            const LpaArc *end = lpa_csr.arcs.data() + lpa_csr.offsets[v1 + 1];
            for (; arc != end; ++arc) {
              size_t k = membership[arc->neighbor];
              if (ncount[k] == 0) touched.push_back(k);
              lpa_accumulate(nscore[k], arc->weight);
              ncount[k]++;
            }
          } else {
            const std::vector<GraphEdgeId> &ineis = g.incs[v1];
            for (size_t j = 0; j < ineis.size(); ++j) {
              size_t edgeID = ineis[j];
              size_t k = membership[g.getOtherNode(edgeID, v1)];
              if (ncount[k] == 0) touched.push_back(k);
              lpa_accumulate(nscore[k], logsscr[edgeID]);
              ncount[k]++;
            }
          }
          if (touched.empty()) { prop_set[v1] = 0; continue; }
          const size_t current_k = membership[v1];
          bool current_is_candidate = false;
          double current_val = -std::numeric_limits<double>::infinity();
          double best_val = -std::numeric_limits<double>::infinity();
          size_t best_k = touched[0];
          for (size_t k : touched) {
            double val = lpa_combine((int)ncount[k], nscore[k]);
            if (k == current_k) {
              current_is_candidate = true;
              current_val = val;
            }
            if (val > best_val) { best_val = val; best_k = k; }
            nscore[k] = 0.0; ncount[k] = 0;
          }
          if (current_is_candidate &&
              current_val >= best_val - LPA_SCORE_EPS)
            best_k = current_k;
          prop_k[v1] = best_k;
          prop_set[v1] = 1;
        }
      }

      // Serial commit in node_order preserves deterministic re-activation and
      // applies the same repeated-label cycle guard as the sequential path.
      std::fill(active.begin(), active.end(), 0);
      size_t changed = 0;
      size_t frozen = 0;
      for (size_t i = 0; i < node_order.size(); i++) {
        size_t v1 = node_order[i];
        if (!prop_set[v1]) continue;
        size_t best_k = prop_k[v1];
        size_t kPrev = membership[v1];
        if (kPrev != best_k && blacklist.find(v1) == blacklist.end()) {
          membership[v1] = best_k;
          ++changed;
          if (use_csr) {
            const LpaArc *arc = lpa_csr.arcs.data() + lpa_csr.offsets[v1];
            const LpaArc *end = lpa_csr.arcs.data() + lpa_csr.offsets[v1 + 1];
            for (; arc != end; ++arc) next_active[arc->neighbor] = 1;
          } else {
            const std::vector<GraphEdgeId> &ineis = g.incs[v1];
            for (size_t j = 0; j < ineis.size(); ++j)
              next_active[g.getOtherNode(ineis[j], v1)] = 1;
          }
          auto &seen = visited[v1];
          if (seen.empty()) seen.insert(kPrev);
          if (!seen.insert(best_k).second) {
            blacklist.insert(v1);
            ++frozen;
          }
        }
      }
      if (lpa_trace)
        fprintf(stderr, "[LPA] round=%zu changed=%zu frozen=%zu\n",
                round, changed, frozen);
      if (changed == 0) break;
      active.swap(next_active);
      std::fill(next_active.begin(), next_active.end(), 0);
    }
    return 0;
  }

  // Reusable dense per-cluster accumulators (cluster IDs ∈ [0, no_of_nodes)).
  // Replaces the per-node unordered_map<cluster,(score,count)> — no hashing and
  // no per-node allocation.  Only the `touched` clusters are reset each node.
  std::vector<double>  nscore(no_of_nodes, 0.0);
  std::vector<int32_t> ncount(no_of_nodes, 0);
  std::vector<size_t>  touched;
  touched.reserve(256);

  // Active-set acceleration: a node's best cluster can only change when one of
  // its neighbours changed label since the node was last processed; otherwise
  // reprocessing it is a no-op (kPrev == best → no state change).  We therefore
  // process only "active" nodes and re-activate a node's neighbours whenever it
  // changes.  This is bit-identical to scanning every node each round (skipped
  // nodes are provably no-ops) but avoids the wasted scans once most nodes have
  // stabilised — which dominates the later rounds.
  std::vector<char> active(no_of_nodes, 1);
  std::unordered_map<size_t, std::unordered_set<size_t>> visited;
  std::unordered_set<size_t> blacklist;

  size_t round = 0;
  while (true) {
    ++round;
    size_t changed = 0;
    size_t frozen = 0;

    /* In the prescribed order, loop over the vertices and reassign labels */
    for (size_t i = 0; i < node_order.size();
         i++) { // we reconsider all nodes regardless of its previous status,
                // but is it better?
      size_t v1 = node_order[i];
      if (!active[v1]) continue;   // unchanged neighbourhood → guaranteed no-op
      active[v1] = 0;              // consume; re-activated if a neighbour changes

      // Accumulate summed log-p-value and neighbour count per neighbouring
      // cluster into the reusable dense arrays (only `touched` ones are live).
      touched.clear();
      if (use_csr) {
        const LpaArc *arc = lpa_csr.arcs.data() + lpa_csr.offsets[v1];
        const LpaArc *end = lpa_csr.arcs.data() + lpa_csr.offsets[v1 + 1];
        for (; arc != end; ++arc) {
          size_t k = membership[arc->neighbor];
          if (ncount[k] == 0) touched.push_back(k);
          lpa_accumulate(nscore[k], arc->weight);
          ncount[k]++;
        }
      } else {
        const std::vector<GraphEdgeId> &ineis = g.incs[v1];
        for (size_t j = 0; j < ineis.size(); ++j) {
          size_t edgeID = ineis[j];
          size_t k = membership[g.getOtherNode(edgeID, v1)];
          if (ncount[k] == 0) touched.push_back(k);
          lpa_accumulate(nscore[k], logsscr[edgeID]);
          ncount[k]++;
        }
      }

      if (!touched.empty()) {
        // Reduce each candidate label's incident edges to one score (rule set by
        // RABBIT_LPA_SCORE) and keep the best. If the current label is tied for
        // best within LPA_SCORE_EPS, retain it and avoid a meaningless update.
        const size_t current_k = membership[v1];
        bool current_is_candidate = false;
        double current_val = -std::numeric_limits<double>::infinity();
        double best_val = -std::numeric_limits<double>::infinity();
        size_t best_k = touched[0];
        int n_at_one = 0;
        for (size_t k : touched) {
          double val = lpa_combine((int)ncount[k], nscore[k]);
          if (k == current_k) {
            current_is_candidate = true;
            current_val = val;
          }
          if (g_lpa_prof && val >= 1.0) ++n_at_one;
          if (val > best_val) { best_val = val; best_k = k; }
          nscore[k] = 0.0; ncount[k] = 0;   // reset for the next node
        }
        if (current_is_candidate &&
            current_val >= best_val - LPA_SCORE_EPS)
          best_k = current_k;
        if (g_lpa_prof) {
          g_lpa_dec[0].fetch_add(1, std::memory_order_relaxed);
          if (n_at_one >= 1) g_lpa_dec[1].fetch_add(1, std::memory_order_relaxed);
          if (n_at_one >= 2) g_lpa_dec[2].fetch_add(1, std::memory_order_relaxed);
        }

        size_t kPrev = membership[v1];
        if (kPrev != best_k &&
            blacklist.find(v1) == blacklist.end()) {
          membership[v1] = best_k;
          ++changed;

          // v1 changed → its neighbours must reconsider (re-activate them).
          if (use_csr) {
            const LpaArc *arc = lpa_csr.arcs.data() + lpa_csr.offsets[v1];
            const LpaArc *end = lpa_csr.arcs.data() + lpa_csr.offsets[v1 + 1];
            for (; arc != end; ++arc) active[arc->neighbor] = 1;
          } else {
            const std::vector<GraphEdgeId> &ineis = g.incs[v1];
            for (size_t j = 0; j < ineis.size(); ++j)
              active[g.getOtherNode(ineis[j], v1)] = 1;
          }

          // Label retention removes score ties, but strict Fisher preferences
          // can still form a cycle. Record every assigned label and freeze only
          // a node that actually returns to one, guaranteeing that a later full
          // sweep can reach changed == 0 without a patience heuristic.
          auto &seen = visited[v1];
          if (seen.empty()) seen.insert(kPrev);
          if (!seen.insert(best_k).second) {
            blacklist.insert(v1);
            ++frozen;
          }
        }
      }
    }
    if (lpa_trace)
      fprintf(stderr, "[LPA] round=%zu changed=%zu frozen=%zu\n",
              round, changed, frozen);
    if (changed == 0) break;
  }

  // Multilevel merge phase (opt-in): conservatively merge same-genome fragments
  // the level-0 LPA left split, using the same Fisher statistic + a high gate.
  // Its tau gate is calibrated on the chi-squared CDF, so it needs the log(1-w)
  // accumuland even when the propagation above was ablated onto another rule.
  if (do_mlpa) {
    if (g_lpa_raw_w)
      for (size_t e = 0; e < no_of_edges; ++e)
        logsscr[e] = (StoredDistance)LOG(1. - g.edgeScore[e]);
    fisher_super_merge(g, membership, logsscr, mlpa_tau, mlpa_min_count);
  }

  if (g_lpa_prof) {
    const uint64_t tot = g_lpa_dec[0].load(), sat = g_lpa_dec[1].load(),
                   tie = g_lpa_dec[2].load();
    if (tot)
      fprintf(stderr,
              "[RB_LPA_PROF] %llu label decisions: winner already 1.0 in %llu "
              "(%.1f%%), >=2 labels tied at 1.0 in %llu (%.1f%%)\n",
              (unsigned long long)tot, (unsigned long long)sat,
              100.0 * sat / tot, (unsigned long long)tie, 100.0 * tie / tot);
  }

  return 0;
}

// ── Per-node assignment confidence + soft (second-choice) assignment ─────────
// Feature #8: after label propagation converges, quantify how strongly each
// node belongs to its assigned bin.  For node v we combine the per-neighbour-
// cluster Fisher p-values (exactly as the propagation step does) and report:
//   confidence[v]   = chi2cdf(assigned cluster) − chi2cdf(best OTHER cluster),
//                     clamped to [0,1].  1.0 = unambiguous; ~0 = borderline.
//   second_choice[v]= the competing cluster id (best alternative bin).
//   second_score[v] = that competitor's chi2cdf score.
// This is the soft / uncertainty information almost no binner exposes, computed
// in a single O(|E|) pass over the converged graph (negligible cost).  Nodes
// with no positive-weight edge get confidence 0 and second_choice = SIZE_MAX.
void compute_node_confidence(Graph &g,
                             const std::vector<size_t> &membership,
                             std::vector<float> &confidence,
                             std::vector<size_t> &second_choice,
                             std::vector<float> &second_score) {
  const size_t nnodes = g.getNodeCount();
  const size_t nedges = g.getEdgeCount();
  confidence.assign(nnodes, 0.0f);
  second_choice.assign(nnodes, SIZE_MAX);
  second_score.assign(nnodes, 0.0f);
  if (nnodes == 0 || nedges == 0) return;

  std::vector<StoredDistance> logsscr(nedges);
  for (size_t e = 0; e < nedges; ++e)
    logsscr[e] = (StoredDistance)LOG(1. - g.edgeScore[e]);

  std::vector<double>  nscore(nnodes, 0.0);
  std::vector<int32_t> ncount(nnodes, 0);
  std::vector<size_t>  touched;
  touched.reserve(256);

  for (size_t v1 = 0; v1 < nnodes; ++v1) {
    std::vector<GraphEdgeId> &ineis = g.incs[v1];
    if (ineis.empty()) continue;
    touched.clear();
    for (size_t j = 0; j < ineis.size(); ++j) {
      size_t edgeID = ineis[j];
      size_t k = membership[g.getOtherNode(edgeID, v1)];
      if (ncount[k] == 0) touched.push_back(k);
      nscore[k] += logsscr[edgeID];
      ncount[k]++;
    }
    const size_t kAssigned = membership[v1];
    double assigned_val = 0.0;
    double best_other = -std::numeric_limits<double>::infinity();
    size_t best_other_k = SIZE_MAX;
    for (size_t k : touched) {
      double val = chi2_2dof_cdf((int)ncount[k], -2.0 * nscore[k]);
      if (k == kAssigned) {
        assigned_val = val;
      } else if (val > best_other) {
        best_other = val; best_other_k = k;
      }
      nscore[k] = 0.0; ncount[k] = 0;   // reset for next node
    }
    double bo = (best_other_k == SIZE_MAX) ? 0.0 : best_other;
    double margin = assigned_val - bo;
    if (margin < 0.0) margin = 0.0;
    if (margin > 1.0) margin = 1.0;
    confidence[v1]   = (float)margin;
    second_choice[v1] = best_other_k;
    second_score[v1]  = (best_other_k == SIZE_MAX) ? 0.0f : (float)best_other;
  }
}

StoredDistance get_element(Matrix const &m, int i, int j) { return m(i, j); }
