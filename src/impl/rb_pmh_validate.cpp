// RabbitBin module: rb_pmh_validate.cpp
//
// Gold-aware, representation-level validation of the PMH neighbourhoods.
// This code is reached only through --validate-pmh-gold and returns from the
// bin command before graph construction.  It deliberately ignores abundance:
// the experiment asks whether the compressed sequence representation itself
// retains same-genome local neighbours.

struct RbPmhValidationCandidate {
  uint32_t id;
  float similarity;
};

// "Better" order: larger similarity first, then smaller contig id.  Supplying
// this relation to priority_queue leaves the weakest retained candidate at top.
struct RbPmhValidationBetter {
  bool operator()(const RbPmhValidationCandidate &a,
                  const RbPmhValidationCandidate &b) const noexcept {
    if (a.similarity != b.similarity) return a.similarity > b.similarity;
    return a.id < b.id;
  }
};

static int rb_validate_pmh_neighbourhoods(const std::string &gold_path,
                                          const std::string &report_path,
                                          size_t requested_queries,
                                          size_t requested_top) {
  if (!g_pmh_mode || (g_win16.empty() && g_win_flat.empty())) {
    cerr << "[Error!] --validate-pmh-gold requires the PMH winner representation\n";
    return 1;
  }
  if (nobs < 2) {
    cerr << "[Error!] PMH validation needs at least two large contigs\n";
    return 1;
  }

  std::ifstream gold(gold_path);
  if (!gold) {
    cerr << "[Error!] cannot open PMH validation gold: " << gold_path << "\n";
    return 1;
  }

  // Keys borrow the stable storage in contig_names; no second copy of all names.
  std::unordered_map<std::string_view, size_t> name_to_index;
  name_to_index.reserve(nobs * 2);
  for (size_t i = 0; i < nobs; ++i) name_to_index.emplace(contig_names[i], i);

  const uint32_t NO_LABEL = std::numeric_limits<uint32_t>::max();
  std::vector<uint32_t> label(nobs, NO_LABEL);
  std::unordered_map<std::string, uint32_t> genome_to_id;
  genome_to_id.reserve(4096);

  std::string line;
  size_t labelled = 0;
  while (std::getline(gold, line)) {
    if (line.empty() || line[0] == '@') continue;
    const size_t t1 = line.find('\t');
    if (t1 == std::string::npos) continue;
    const size_t t2 = line.find('\t', t1 + 1);
    const size_t bin_end = (t2 == std::string::npos) ? line.size() : t2;
    if (bin_end <= t1 + 1) continue;
    const std::string_view name(line.data(), t1);
    auto ni = name_to_index.find(name);
    if (ni == name_to_index.end()) continue;
    const std::string bin(line.data() + t1 + 1, bin_end - t1 - 1);
    auto gi = genome_to_id.find(bin);
    if (gi == genome_to_id.end()) {
      const uint32_t id = (uint32_t)genome_to_id.size();
      gi = genome_to_id.emplace(bin, id).first;
    }
    if (label[ni->second] == NO_LABEL) ++labelled;
    label[ni->second] = gi->second;
  }
  name_to_index.clear();

  std::vector<size_t> genome_size(genome_to_id.size(), 0);
  for (uint32_t x : label) if (x != NO_LABEL) ++genome_size[x];

  std::vector<uint32_t> eligible_queries;
  eligible_queries.reserve(labelled);
  for (size_t i = 0; i < nobs; ++i)
    if (label[i] != NO_LABEL && genome_size[label[i]] > 1)
      eligible_queries.push_back((uint32_t)i);
  if (eligible_queries.empty()) {
    cerr << "[Error!] PMH validation gold has no genome represented by >=2 "
            "large contigs\n";
    return 1;
  }

  // The same --seed and FASTA order select the same query set for every sketch
  // size, enabling paired m=100/250/500/1000 comparisons.
  std::mt19937_64 rng(seed);
  std::shuffle(eligible_queries.begin(), eligible_queries.end(), rng);
  const size_t nq = std::min(requested_queries, eligible_queries.size());
  eligible_queries.resize(nq);
  const size_t top = std::min(requested_top, nobs - 1);

  std::vector<size_t> ranks;
  for (size_t n : {size_t(50), size_t(100), size_t(200), size_t(400)})
    if (n <= top) ranks.push_back(n);
  if (ranks.empty() || ranks.back() != top) ranks.push_back(top);
  const size_t nr = ranks.size();

  std::vector<uint32_t> same_at(nq * nr, 0);
  std::vector<uint32_t> first_same_rank(nq, 0);
  std::vector<uint32_t> possible(nq, 0);
  const auto started = std::chrono::steady_clock::now();

  if (!g_win16.empty())
    numa_interleave_buffer(g_win16.data(), g_win16.size() * sizeof(uint16_t));
  else
    numa_interleave_buffer(g_win_flat.data(),
                           g_win_flat.size() * sizeof(uint32_t));

#pragma omp parallel for schedule(dynamic, 1) num_threads(numThreads)
  for (size_t qi = 0; qi < nq; ++qi) {
    const uint32_t q = eligible_queries[qi];
    const uint32_t qlabel = label[q];
    possible[qi] = (uint32_t)(genome_size[qlabel] - 1);
    std::priority_queue<RbPmhValidationCandidate,
                        std::vector<RbPmhValidationCandidate>,
                        RbPmhValidationBetter> heap;
    const RbPmhValidationBetter better{};
    for (size_t j = 0; j < nobs; ++j) {
      if (j == q) continue;
      const RbPmhValidationCandidate candidate{
          (uint32_t)j, (float)graph_sim(q, j)};
      if (heap.size() < top) {
        heap.push(candidate);
      } else if (better(candidate, heap.top())) {
        heap.pop();
        heap.push(candidate);
      }
    }

    std::vector<RbPmhValidationCandidate> row;
    row.reserve(heap.size());
    while (!heap.empty()) {
      row.push_back(heap.top());
      heap.pop();
    }
    std::sort(row.begin(), row.end(), better);

    uint32_t same = 0;
    size_t ri = 0;
    for (size_t p = 0; p < row.size(); ++p) {
      if (label[row[p].id] == qlabel) {
        ++same;
        if (first_same_rank[qi] == 0) first_same_rank[qi] = (uint32_t)(p + 1);
      }
      while (ri < nr && p + 1 == ranks[ri]) {
        same_at[qi * nr + ri] = same;
        ++ri;
      }
    }
  }

  const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  std::ofstream report(report_path);
  if (!report) {
    cerr << "[Error!] cannot write PMH validation report: " << report_path << "\n";
    return 1;
  }
  report << "sketch_m\tpmh_k\twinner_bits\tnobs\tlabelled_contigs\tgenomes"
            "\tqueries\ttop_n\tdirected_candidates\tsame_genome_hits"
            "\tprecision_at_n\trecall_all_same\trecall_capped"
            "\tmean_query_recall_capped\tquery_hit_rate"
            "\tquery_hit_rate_k5\tquery_hit_rate_k10"
            "\tmrr\tcomparison_seconds\n";
  report << std::setprecision(10);

  for (size_t ri = 0; ri < nr; ++ri) {
    uint64_t hits = 0, candidates = 0, denom_all = 0, denom_cap = 0;
    size_t query_hits = 0, query_hits_k5 = 0, query_hits_k10 = 0;
    double macro_capped = 0.0, mrr = 0.0;
    for (size_t qi = 0; qi < nq; ++qi) {
      const uint64_t h = same_at[qi * nr + ri];
      hits += h;
      candidates += ranks[ri];
      denom_all += possible[qi];
      const uint64_t cap = std::min<uint64_t>(ranks[ri], possible[qi]);
      denom_cap += cap;
      if (cap) macro_capped += (double)h / (double)cap;
      if (h) ++query_hits;
      if (h >= 5) ++query_hits_k5;
      if (h >= 10) ++query_hits_k10;
      const uint32_t fr = first_same_rank[qi];
      if (fr && fr <= ranks[ri]) mrr += 1.0 / (double)fr;
    }
    report << g_pmh_m << '\t' << g_pmh_k << '\t' << g_win_bits << '\t'
           << nobs << '\t' << labelled << '\t' << genome_to_id.size() << '\t'
           << nq << '\t' << ranks[ri] << '\t' << candidates << '\t' << hits
           << '\t' << (candidates ? (double)hits / candidates : 0.0)
           << '\t' << (denom_all ? (double)hits / denom_all : 0.0)
           << '\t' << (denom_cap ? (double)hits / denom_cap : 0.0)
           << '\t' << (nq ? macro_capped / (double)nq : 0.0)
           << '\t' << (nq ? (double)query_hits / (double)nq : 0.0)
           << '\t' << (nq ? (double)query_hits_k5 / (double)nq : 0.0)
           << '\t' << (nq ? (double)query_hits_k10 / (double)nq : 0.0)
           << '\t' << (nq ? mrr / (double)nq : 0.0)
           << '\t' << elapsed << '\n';
  }
  report.close();

  verbose_message(
      "PMH validation: %zu seed-controlled queries, top-%zu, %zu labelled "
      "large contigs / %zu genomes, %.3fs -> %s\n",
      nq, top, labelled, genome_to_id.size(), elapsed, report_path.c_str());
  return 0;
}
