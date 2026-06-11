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

int cluster_by_propagation(Graph &g, std::vector<size_t> &membership,
                      std::vector<size_t> &node_order) {
  size_t no_of_nodes = g.getNodeCount();
  size_t no_of_edges = g.getEdgeCount();

  if (no_of_nodes == 0 || no_of_edges == 0) {
    cerr << "There were " << no_of_nodes << " nodes and " << no_of_edges
         << " edges -- skipping cluster_by_propagation" << endl;
    return 0;
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

  return 0;
}

StoredDistance get_element(Matrix const &m, int i, int j) { return m(i, j); }

