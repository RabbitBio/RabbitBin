// RabbitBin module: rb_graph_audit.cpp
//
// Gold-aware audit of the graph that the production binning path actually
// builds.  Unlike --validate-pmh-gold (directed, sequence-only top-N lists),
// this diagnostic measures the mutual candidate graph after its production
// candidate gates and the retained graph after abundance edge scoring.
// Ground-truth labels are read only after edge scores have been computed and
// never affect graph construction, filtering, label propagation, or output.

struct RbGraphAuditGold {
  std::vector<uint32_t> label;
  std::vector<size_t> genome_size;
  size_t labelled = 0;
  size_t genomes = 0;
};

static bool rb_load_graph_audit_gold(const std::string &gold_path,
                                     RbGraphAuditGold &out) {
  std::ifstream gold(gold_path);
  if (!gold) {
    cerr << "[Error!] cannot open graph-audit gold: " << gold_path << "\n";
    return false;
  }

  std::unordered_map<std::string_view, size_t> name_to_index;
  name_to_index.reserve(nobs * 2);
  for (size_t i = 0; i < nobs; ++i) name_to_index.emplace(contig_names[i], i);

  const uint32_t NO_LABEL = std::numeric_limits<uint32_t>::max();
  out.label.assign(nobs, NO_LABEL);
  std::unordered_map<std::string, uint32_t> genome_to_id;
  genome_to_id.reserve(4096);

  std::string line;
  while (std::getline(gold, line)) {
    if (line.empty() || line[0] == '@' || line[0] == '#') continue;
    const size_t t1 = line.find('\t');
    if (t1 == std::string::npos) continue;
    const size_t t2 = line.find('\t', t1 + 1);
    const size_t bin_end = (t2 == std::string::npos) ? line.size() : t2;
    if (bin_end <= t1 + 1) continue;
    const std::string_view name(line.data(), t1);
    const auto ni = name_to_index.find(name);
    if (ni == name_to_index.end()) continue;

    const std::string bin(line.data() + t1 + 1, bin_end - t1 - 1);
    auto gi = genome_to_id.find(bin);
    if (gi == genome_to_id.end()) {
      const uint32_t id = (uint32_t)genome_to_id.size();
      gi = genome_to_id.emplace(bin, id).first;
    }
    if (out.label[ni->second] == NO_LABEL) ++out.labelled;
    out.label[ni->second] = gi->second;
  }

  out.genomes = genome_to_id.size();
  out.genome_size.assign(out.genomes, 0);
  for (uint32_t x : out.label)
    if (x != NO_LABEL) ++out.genome_size[x];
  return out.labelled > 0;
}

static int rb_audit_production_graph(const Graph &g,
                                     const std::string &gold_path,
                                     const std::string &report_path) {
  RbGraphAuditGold gold;
  if (!rb_load_graph_audit_gold(gold_path, gold)) {
    if (gold.label.empty()) return 1;
    cerr << "[Error!] graph-audit gold matched no large contigs\n";
    return 1;
  }

  const uint32_t NO_LABEL = std::numeric_limits<uint32_t>::max();
  std::vector<uint16_t> candidate_true_degree(nobs, 0);
  std::vector<uint16_t> retained_true_degree(nobs, 0);

  uint64_t candidate_edges = 0, candidate_true = 0, candidate_false = 0;
  uint64_t candidate_unknown = 0, retained_edges = 0, retained_true = 0;
  uint64_t retained_false = 0, retained_unknown = 0, gfa_edges = 0;

  auto increment_degree = [](uint16_t &x) {
    if (x != std::numeric_limits<uint16_t>::max()) ++x;
  };

  const size_t ne = g.from.size();
  for (size_t e = 0; e < ne; ++e) {
    if (edge_is_gfa(g.sComp[e])) {
      ++gfa_edges;
      continue;
    }
    ++candidate_edges;
    const size_t i = g.from[e], j = g.to[e];
    const bool known = i < gold.label.size() && j < gold.label.size() &&
                       gold.label[i] != NO_LABEL && gold.label[j] != NO_LABEL;
    const bool same = known && gold.label[i] == gold.label[j];
    if (!known) {
      ++candidate_unknown;
    } else if (same) {
      ++candidate_true;
      increment_degree(candidate_true_degree[i]);
      increment_degree(candidate_true_degree[j]);
    } else {
      ++candidate_false;
    }

    if (e >= g.edgeScore.size() || g.edgeScore[e] <= 0.0f) continue;
    ++retained_edges;
    if (!known) {
      ++retained_unknown;
    } else if (same) {
      ++retained_true;
      increment_degree(retained_true_degree[i]);
      increment_degree(retained_true_degree[j]);
    } else {
      ++retained_false;
    }
  }

  uint64_t eligible = 0;
  uint64_t candidate_c1 = 0, candidate_c5 = 0, candidate_c10 = 0;
  uint64_t retained_c1 = 0, retained_c5 = 0, retained_c10 = 0;
  for (size_t i = 0; i < nobs; ++i) {
    const uint32_t lab = gold.label[i];
    if (lab == NO_LABEL || gold.genome_size[lab] < 2) continue;
    ++eligible;
    const uint16_t cd = candidate_true_degree[i];
    const uint16_t rd = retained_true_degree[i];
    candidate_c1 += cd >= 1; candidate_c5 += cd >= 5; candidate_c10 += cd >= 10;
    retained_c1 += rd >= 1; retained_c5 += rd >= 5; retained_c10 += rd >= 10;
  }

  const auto ratio = [](uint64_t a, uint64_t b) {
    return b ? (double)a / (double)b : 0.0;
  };
  const uint64_t candidate_known = candidate_true + candidate_false;
  const uint64_t retained_known = retained_true + retained_false;

  std::ofstream report(report_path);
  if (!report) {
    cerr << "[Error!] cannot write graph-audit report: " << report_path << "\n";
    return 1;
  }
  report << "max_edges\tnobs\tlabelled_contigs\tgenomes\teligible_contigs"
            "\tcandidate_edges\tcandidate_true\tcandidate_false"
            "\tcandidate_unknown\tcandidate_precision\tcandidate_mean_degree"
            "\tcandidate_true_c1\tcandidate_true_c5\tcandidate_true_c10"
            "\tretained_edges\tretained_true\tretained_false"
            "\tretained_unknown\tretained_precision\tretained_mean_degree"
            "\tretained_true_c1\tretained_true_c5\tretained_true_c10"
            "\toverall_edge_survival\ttrue_edge_survival"
            "\tfalse_edge_survival\tgfa_edges\n";
  report << std::setprecision(10)
         << maxEdges << '\t' << nobs << '\t' << gold.labelled << '\t'
         << gold.genomes << '\t' << eligible << '\t'
         << candidate_edges << '\t' << candidate_true << '\t'
         << candidate_false << '\t' << candidate_unknown << '\t'
         << ratio(candidate_true, candidate_known) << '\t'
         << ratio(2 * candidate_edges, nobs) << '\t'
         << ratio(candidate_c1, eligible) << '\t'
         << ratio(candidate_c5, eligible) << '\t'
         << ratio(candidate_c10, eligible) << '\t'
         << retained_edges << '\t' << retained_true << '\t'
         << retained_false << '\t' << retained_unknown << '\t'
         << ratio(retained_true, retained_known) << '\t'
         << ratio(2 * retained_edges, nobs) << '\t'
         << ratio(retained_c1, eligible) << '\t'
         << ratio(retained_c5, eligible) << '\t'
         << ratio(retained_c10, eligible) << '\t'
         << ratio(retained_edges, candidate_edges) << '\t'
         << ratio(retained_true, candidate_true) << '\t'
         << ratio(retained_false, candidate_false) << '\t'
         << gfa_edges << '\n';
  report.close();

  verbose_message(
      "Graph audit: candidate=%llu (precision=%.4f), retained=%llu "
      "(precision=%.4f), true-neighbour C1 %.4f -> %.4f; %s\n",
      (unsigned long long)candidate_edges,
      ratio(candidate_true, candidate_known),
      (unsigned long long)retained_edges,
      ratio(retained_true, retained_known),
      ratio(candidate_c1, eligible), ratio(retained_c1, eligible),
      report_path.c_str());
  return 0;
}
