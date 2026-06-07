// RabbitBin module: rb_cluster.cpp

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
          boost::math::chi_squared chi_sqr_dist(2 * ncount[k]);
          double val = boost::math::cdf(chi_sqr_dist, -2.0 * nscore[k]);
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

