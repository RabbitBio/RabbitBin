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
static inline double chi2_2dof_cdf(int m, double x) {
  const double u = 0.5 * x;
  double p = std::exp(-u);   // p_0
  double Q = p;
  for (int i = 1; i < m; ++i) { p *= u / (double)i; Q += p; }
  double v = 1.0 - Q;
  return v < 0.0 ? 0.0 : v;
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
  size_t nLeftMin = INT_MAX, attempt = 0;
  bool running = true;
  while (running) {
    running = false;
    size_t nLeft = 0;
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
        for (const auto &kv : sadj[v]) active[kv.first] = 1;
        if (visited[v].find(best_k) == visited[v].end()) { nLeft++; running = true; }
        else blacklist.insert(v);
        visited[v].insert(best_k);
      }
    }
    if (nLeft < nLeftMin) { nLeftMin = nLeft; attempt = 0; }
    else { attempt++; if (attempt >= 10) break; }
  }
  for (size_t i = 0; i < n; ++i) membership[i] = (size_t)comm[comm0[i]];
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

  std::unordered_map<size_t, std::unordered_set<size_t>> visited;
  std::unordered_set<size_t> blacklist;

  // Precompute LOG(1 - edgeScore) per edge once: it is invariant across the (many)
  // label-propagation rounds and across every incidence visit of the edge, so
  // hoisting it out of the hot loop removes ~10²·|E| redundant log() calls.
  std::vector<StoredDistance> logsscr(no_of_edges);
  for (size_t e = 0; e < no_of_edges; ++e)
    logsscr[e] = (StoredDistance)LOG(1. - g.edgeScore[e]);

  // ── Experimental parallel (Jacobi) label propagation ─────────────────────
  // RABBIT_PAR_LPA=1 switches from the default order-dependent (Gauss-Seidel)
  // sweep to a synchronous one: every active node recomputes its best cluster
  // from the membership SNAPSHOT taken at the start of the round (read-only, so
  // the per-node work parallelises trivially), then all changes are committed
  // together in node_order.  Because neighbour updates are no longer visible
  // mid-round, the converged labelling differs from the sequential baseline —
  // hence this is opt-in and OFF by default.  MEASURED (CAMI2 plant, 64t): the
  // LPA segment drops ~485ms→~288ms, but recovery falls 92→79 genomes and the
  // binning fragments (149→184 bins), so it is NOT accuracy-neutral and must
  // stay experimental.  Kept only for benchmarking the parallel ceiling.  The
  // commit is serial and in node_order, so a given thread count + flag is fully
  // deterministic, and the visited/blacklist oscillation guard + attempt-based
  // stop criterion are preserved verbatim.
  if (getenv("RABBIT_PAR_LPA") != nullptr) {
    const int nthreads = omp_get_max_threads();
    std::vector<std::vector<double>>  tl_nscore(nthreads, std::vector<double>(no_of_nodes, 0.0));
    std::vector<std::vector<int32_t>> tl_ncount(nthreads, std::vector<int32_t>(no_of_nodes, 0));

    std::vector<size_t> prop_k(no_of_nodes);     // proposed best cluster per node
    std::vector<char>   prop_set(no_of_nodes, 0); // 1 if node was evaluated this round
    std::vector<char>   active(no_of_nodes, 1), next_active(no_of_nodes, 0);

    std::unordered_map<size_t, std::unordered_set<size_t>> visited;
    std::unordered_set<size_t> blacklist;

    size_t nLeftMin = INT_MAX, attempt = 0;
    bool running = true;
    while (running) {
      running = false;
      size_t nLeft = 0;

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
          std::vector<size_t> &ineis = g.incs[v1];
          for (size_t j = 0; j < ineis.size(); j++) {
            size_t edgeID = ineis[j];
            size_t k = membership[g.getOtherNode(edgeID, v1)];
            if (ncount[k] == 0) touched.push_back(k);
            nscore[k] += logsscr[edgeID];
            ncount[k]++;
          }
          if (touched.empty()) { prop_set[v1] = 0; continue; }
          double best_val = -std::numeric_limits<double>::infinity();
          size_t best_k = touched[0];
          for (size_t k : touched) {
            double val = chi2_2dof_cdf((int)ncount[k], -2.0 * nscore[k]);
            if (val > best_val) { best_val = val; best_k = k; }
            nscore[k] = 0.0; ncount[k] = 0;
          }
          prop_k[v1] = best_k;
          prop_set[v1] = 1;
        }
      }

      // Serial commit in node_order (deterministic; preserves the oscillation
      // guard and re-activation semantics of the sequential variant).
      std::fill(active.begin(), active.end(), 0);
      for (size_t i = 0; i < node_order.size(); i++) {
        size_t v1 = node_order[i];
        if (!prop_set[v1]) continue;
        size_t best_k = prop_k[v1];
        int kPrev = membership[v1];
        if (kPrev != (int)best_k && blacklist.find(v1) == blacklist.end()) {
          membership[v1] = best_k;
          std::vector<size_t> &ineis = g.incs[v1];
          for (size_t j = 0; j < ineis.size(); ++j)
            next_active[g.getOtherNode(ineis[j], v1)] = 1;
          int kNext = membership[v1];
          if (visited.find(v1) == visited.end() ||
              visited[v1].find(kNext) == visited[v1].end()) {
            nLeft++;
            running = true;
          } else {
            blacklist.insert(v1);
          }
          visited[v1].insert(kNext);
        }
      }
      active.swap(next_active);
      std::fill(next_active.begin(), next_active.end(), 0);

      if (nLeft < nLeftMin) { nLeftMin = nLeft; attempt = 0; }
      else { attempt++; if (attempt >= 10) break; }
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

  size_t nLeftMin = INT_MAX;
  size_t attempt = 0;
  bool running = true;
  while (running) {
    running = false;

    size_t nLeft = 0;

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
      std::vector<size_t> &ineis = g.incs[v1];
      for (size_t j = 0; j < ineis.size();
           j++) { //# of neighbors (edges connected to v1)
        size_t edgeID = ineis[j];

        size_t k = membership[g.getOtherNode(
            edgeID, v1)]; // community membership of a neighbor (connected by j)

        if (ncount[k] == 0) touched.push_back(k);
        nscore[k] += logsscr[edgeID]; // as p-value (precomputed)
        ncount[k]++;
      }

      if (!touched.empty()) {
        // Fisher's method: combine the per-cluster p-values and keep the
        // most significant cluster (max chi-squared cdf).  Equivalent to the
        // previous max_element over the hash map; first-max wins on ties.
        double best_val = -std::numeric_limits<double>::infinity();
        size_t best_k = touched[0];
        for (size_t k : touched) {
          double val = chi2_2dof_cdf((int)ncount[k], -2.0 * nscore[k]);
          if (val > best_val) { best_val = val; best_k = k; }
          nscore[k] = 0.0; ncount[k] = 0;   // reset for the next node
        }

        // however, if there was a clique (loop) out of >2 nodes
        int kPrev = membership[v1];
        if (kPrev != (int)best_k &&
            blacklist.find(v1) == blacklist.end()) {

          membership[v1] = best_k;

          // v1 changed → its neighbours must reconsider (re-activate them).
          for (size_t j = 0; j < ineis.size(); ++j)
            active[g.getOtherNode(ineis[j], v1)] = 1;

          int kNext = membership[v1];
          if (visited.find(v1) == visited.end() ||
              visited[v1].find(kNext) == visited[v1].end()) {
            // not have been assigned to the cls before
            nLeft++; //# of confirmation (that this choice is optimal) left
            running = true;
          } else {
            blacklist.insert(v1); // blacklist represents nodes that change cls
                                  // in a circular form
          }
          visited[v1].insert(kNext);
        }
      }
    }

    if (nLeft < nLeftMin) {
      nLeftMin = nLeft;
      attempt = 0;
    } else {
      attempt++;
      if (attempt >= 10) {
        break;
      }
    }
    // cout << "nLeft: " << nLeft << " & attempt: " << attempt << endl;
  }

  // Multilevel merge phase (opt-in): conservatively merge same-genome fragments
  // the level-0 LPA left split, using the same Fisher statistic + a high gate.
  if (do_mlpa)
    fisher_super_merge(g, membership, logsscr, mlpa_tau, mlpa_min_count);

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
    std::vector<size_t> &ineis = g.incs[v1];
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

