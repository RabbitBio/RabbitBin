// RabbitBin module: rb_qc.cpp
//
// `rabbitbin qc` — fast, reference-free bin quality estimation from single-copy
// marker genes (SCGs).  This is the in-house, multithreaded analogue of the
// completeness/contamination step that CheckM / MaxBin / MetaDecoder perform.
//
// Design (matches RabbitBin's "expensive once, cheap many" philosophy):
//   * The EXPENSIVE part — ORF prediction (Prodigal) + HMM search (hmmsearch)
//     against a single-copy marker set — is run ONCE per assembly by an external
//     helper (scripts/rabbitbin_markers.sh), producing a contig->marker map.
//   * `rabbitbin qc` then scores ANY binning of that assembly in milliseconds:
//     it never touches sequences, only the (binning, marker-map) pair.
//
// Per-bin metrics (CheckM-style, copy-number based over a marker set of size G):
//   present(bin)        = # distinct markers found on >=1 contig of the bin
//   completeness  (%)   = 100 * present / G
//   contamination (%)   = 100 * sum_m max(copies_m - 1, 0) / G
//   (copies_m = # contigs in the bin carrying marker m)
// MIMAG-style tiers:
//   HQ (high-quality)   : completeness > 90 AND contamination < 5
//   MQ (medium-quality) : completeness >= 50 AND contamination < 10  (and not HQ)
//
// Marker-map file (one contig per line, tab-separated):
//   contigName <TAB> markerA <TAB> markerB ...        (a contig may repeat)
//   # marker_set_size=107                             (optional header comment)
// Lines beginning with '#' are comments; `# marker_set_size=N` sets G.
//
// Binning input mirrors `rabbitbin amber`: --binning (CAMI bioboxes or 2-col
// SEQUENCEID<TAB>BINID) or --members (rabbitbin members.tsv: BinNum<TAB>Name).
//
// Reuses rb_amber:: helpers (map_file / split_tabs / Interner); this file is
// #included after rb_amber.cpp.

// ── In-pipeline marker loading + scoring (for `bin --auto` / `--ensemble`) ───
// g_contig_marker_ids / g_marker_set_size are declared (and defined) in
// rabbitbin.h; this file (part of the same TU) only implements the functions.

bool rb_load_markers_for_contigs(const std::string &path) {
  using namespace rb_amber;
  g_contig_marker_ids.assign(nobs, {});
  g_marker_set_size = 0;

  // Map whitespace-stripped large-contig name -> index.
  phmap::flat_hash_map<std::string, size_t> name2idx;
  name2idx.reserve(nobs * 2);
  for (size_t i = 0; i < nobs; ++i) {
    std::string nm = contig_names[i];
    size_t sp = nm.find_first_of(" \t");
    if (sp != std::string::npos) nm.resize(sp);
    name2idx.emplace(std::move(nm), i);
  }

  Mapped km;
  if (!map_file(path, km)) {
    cerr << "[Error!] cannot open --markers file: " << path << "\n";
    return false;
  }
  Interner markers;
  long set_size_hdr = 0;
  size_t hits = 0, matched = 0;
  const char *p = km.data, *end = km.data + km.len;
  const char *fb[64], *fe[64];
  while (p < end) {
    const char *nl = (const char *)memchr(p, '\n', (size_t)(end - p));
    const char *le = nl ? nl : end;
    if (le > p) {
      if (p[0] == '#') {
        static const char *kw = "marker_set_size";
        const size_t kl = strlen(kw);
        for (const char *q = p; q + kl < le; ++q) {
          if (strncmp(q, kw, kl) == 0) {
            const char *eq = (const char *)memchr(q + kl, '=', (size_t)(le - (q + kl)));
            if (eq) set_size_hdr = (long)parse_uint(eq + 1, le);
            break;
          }
        }
      } else {
        int n = split_tabs(p, le, fb, fe, 64);
        if (n >= 2) {
          std::string contig(fb[0], (size_t)(fe[0] - fb[0]));
          auto it = name2idx.find(contig);
          if (it != name2idx.end()) {
            auto &vec = g_contig_marker_ids[it->second];
            for (int c = 1; c < n; ++c) {
              if (fe[c] <= fb[c]) continue;
              vec.push_back(markers.intern(fb[c], fe[c]));
              ++hits;
            }
            ++matched;
          } else {
            // still count distinct markers for fallback G
            for (int c = 1; c < n; ++c)
              if (fe[c] > fb[c]) markers.intern(fb[c], fe[c]);
          }
        }
      }
    }
    if (!nl) break;
    p = nl + 1;
  }
  g_marker_set_size = set_size_hdr > 0 ? set_size_hdr : (long)markers.names.size();
  verbose_message("Loaded markers: %zu hits on %zu/%zu large contigs, "
                  "G=%ld%s\n", hits, matched, nobs, g_marker_set_size,
                  set_size_hdr > 0 ? " (header)" : " (distinct-seen)");
  return g_marker_set_size > 0 && hits > 0;
}

void rb_qc_score_membership(const std::vector<size_t> &mem,
                            const std::vector<char> &binned,
                            double comp_thr, double cont_thr,
                            long &near_complete, double &sum_comp) {
  near_complete = 0;
  sum_comp = 0.0;
  if (g_marker_set_size <= 0 || g_contig_marker_ids.empty()) return;
  // cluster id -> (marker id -> copies)
  phmap::flat_hash_map<size_t, phmap::flat_hash_map<int, int>> binmk;
  for (size_t i = 0; i < nobs && i < mem.size(); ++i) {
    if (!binned[i]) continue;
    const auto &mk = g_contig_marker_ids[i];
    if (mk.empty()) continue;
    auto &cnt = binmk[mem[i]];
    for (int m : mk) cnt[m] += 1;
  }
  const double invG = 100.0 / (double)g_marker_set_size;
  for (auto &kv : binmk) {
    long pres = 0, extra = 0;
    for (auto &mp : kv.second) {
      if (mp.second >= 1) ++pres;
      if (mp.second > 1)  extra += (mp.second - 1);
    }
    double comp = (double)pres * invG;
    double cont = (double)extra * invG;
    if (comp > comp_thr && cont < cont_thr) {
      ++near_complete;
      sum_comp += comp;
    }
  }
}

void rb_qc_ensemble(const std::vector<std::vector<size_t>> &mems_all,
                    const Graph &g, std::vector<size_t> &out_membership) {
  const size_t M = mems_all.size();
  if (M == 0 || g_marker_set_size <= 0) return;
  const double invG = 100.0 / (double)g_marker_set_size;

  // SCG quality of a contig set → (completeness, contamination).
  auto quality = [&](const std::vector<size_t> &contigs,
                     double &comp, double &cont) {
    phmap::flat_hash_map<int, int> cnt;
    for (size_t i : contigs)
      if (i < g_contig_marker_ids.size())
        for (int m : g_contig_marker_ids[i]) cnt[m] += 1;
    long pres = 0, extra = 0;
    for (auto &mp : cnt) { if (mp.second >= 1) ++pres; if (mp.second > 1) extra += mp.second - 1; }
    comp = (double)pres * invG;
    cont = (double)extra * invG;
  };

  struct Cand { double score; double comp; size_t bp; std::vector<size_t> contigs; };
  std::vector<Cand> cands;

  // Generate candidate bins from every config. Crucially we run the abundance
  // split on each config's clusters FIRST, so candidates are the same pure,
  // strain-resolved bins the production pipeline would emit — this is where most
  // HQ bins come from. We then pool post-split bins across all configs and pick
  // the best non-overlapping set. (abundance_guided_split sets min_bin_bp=0, so
  // save/restore the user's floor around each call.)
  const size_t saved_min_bin_bp = min_bin_bp;
  const bool saved_verbose = verbose;
  verbose = false;                                   // silence per-config split logs
  for (size_t c = 0; c < M; ++c) {
    const auto &mem = mems_all[c];
    BinMap cls_c;
    for (size_t i = 0; i < nobs && i < mem.size(); ++i)
      cls_c[mem[i]].push_back(i);
    min_bin_bp = saved_min_bin_bp;
    if (g_split_abundance || !g_no_split_abundance)
      abundance_guided_split(cls_c);
    for (auto &kv : cls_c) {
      if (kv.second.size() < 2) continue;          // skip singletons
      double comp, cont;
      quality(kv.second, comp, cont);
      if (comp < 30.0) continue;                    // ignore very incomplete bins
      size_t bp = 0;
      for (size_t i : kv.second) bp += (i < seq_lens.size() ? seq_lens[i] : 0);
      // DAS Tool-style SCG score: completeness penalised by contamination.
      double score = comp - 5.0 * cont;
      cands.push_back({score, comp, bp, std::move(kv.second)});
    }
  }
  min_bin_bp = saved_min_bin_bp;
  verbose = saved_verbose;

  std::sort(cands.begin(), cands.end(), [](const Cand &a, const Cand &b) {
    if (a.score != b.score) return a.score > b.score;
    if (a.comp  != b.comp)  return a.comp  > b.comp;
    return a.bp > b.bp;
  });

  // Greedy non-overlapping selection (best-first). A candidate is accepted only
  // if a majority of its contigs are still free, to avoid fragmented leftovers.
  out_membership.assign(nobs, 0);
  std::vector<char> used(nobs, 0);
  size_t next_id = 0, accepted = 0;
  for (auto &cd : cands) {
    size_t free_n = 0;
    for (size_t i : cd.contigs) if (!used[i]) ++free_n;
    if (free_n == 0 || free_n * 2 < cd.contigs.size()) continue;  // <50% free
    size_t id = next_id++;
    for (size_t i : cd.contigs) {
      if (used[i]) continue;
      used[i] = 1;
      out_membership[i] = id;
    }
    ++accepted;
  }
  // Leftover contigs → unique singleton ids (dropped/recruited downstream).
  for (size_t i = 0; i < nobs; ++i)
    if (!used[i]) out_membership[i] = next_id++;

  verbose_message("Ensemble consensus: %zu candidate bins from %zu configs -> "
                  "%zu accepted non-overlapping bins\n",
                  cands.size(), M, accepted);
}

// ── Output-time marker/taxonomy loading + per-bin metrics (features #1/#7/#8) ─
// Build a name -> unified-contig-index map (large: i; small: j+nobs), exactly as
// the membership/output stages index contigs.
static void rb_build_name2idx(phmap::flat_hash_map<std::string, size_t> &m) {
  m.reserve((nobs + nobs1) * 2);
  auto add = [&](const std::string &raw, size_t idx) {
    std::string nm = raw;
    size_t sp = nm.find_first_of(" \t");
    if (sp != std::string::npos) nm.resize(sp);
    m.emplace(std::move(nm), idx);
  };
  for (size_t i = 0; i < nobs; ++i)  add(contig_names[i], i);
  for (size_t j = 0; j < nobs1; ++j) add(small_contig_names[j], j + nobs);
}

bool rb_load_markers_unified() {
  using namespace rb_amber;
  if (g_markers_file.empty()) return false;
  g_allcontig_markers.assign(nobs + nobs1, {});
  phmap::flat_hash_map<std::string, size_t> name2idx;
  rb_build_name2idx(name2idx);

  Mapped km;
  if (!map_file(g_markers_file, km)) {
    cerr << "[Error!] cannot open --markers file: " << g_markers_file << "\n";
    return false;
  }
  Interner markers;
  long hdr = 0; size_t hits = 0, matched = 0;
  const char *p = km.data, *end = km.data + km.len;
  const char *fb[64], *fe[64];
  while (p < end) {
    const char *nl = (const char *)memchr(p, '\n', (size_t)(end - p));
    const char *le = nl ? nl : end;
    if (le > p) {
      if (p[0] == '#') {
        static const char *kw = "marker_set_size"; const size_t kl = strlen(kw);
        for (const char *q = p; q + kl < le; ++q)
          if (strncmp(q, kw, kl) == 0) {
            const char *eq = (const char *)memchr(q + kl, '=', (size_t)(le - (q + kl)));
            if (eq) hdr = (long)parse_uint(eq + 1, le);
            break;
          }
      } else {
        int n = split_tabs(p, le, fb, fe, 64);
        if (n >= 2) {
          std::string contig(fb[0], (size_t)(fe[0] - fb[0]));
          auto it = name2idx.find(contig);
          if (it != name2idx.end()) {
            auto &vec = g_allcontig_markers[it->second];
            for (int c = 1; c < n; ++c)
              if (fe[c] > fb[c]) { vec.push_back(markers.intern(fb[c], fe[c])); ++hits; }
            ++matched;
          } else {
            for (int c = 1; c < n; ++c)
              if (fe[c] > fb[c]) markers.intern(fb[c], fe[c]);
          }
        }
      }
    }
    if (!nl) break;
    p = nl + 1;
  }
  g_marker_set_size = hdr > 0 ? hdr : (long)markers.names.size();
  verbose_message("QC markers: %zu hits on %zu contigs, G=%ld%s\n",
                  hits, matched, g_marker_set_size, hdr > 0 ? " (header)" : "");
  return g_marker_set_size > 0 && hits > 0;
}

bool rb_load_taxonomy_unified() {
  using namespace rb_amber;
  if (g_taxonomy_file.empty()) return false;
  g_contig_taxon.assign(nobs + nobs1, -1);
  g_taxon_names.clear();
  phmap::flat_hash_map<std::string, size_t> name2idx;
  rb_build_name2idx(name2idx);
  phmap::flat_hash_map<std::string, int> taxid;

  Mapped tm;
  if (!map_file(g_taxonomy_file, tm)) {
    cerr << "[Error!] cannot open --taxonomy file: " << g_taxonomy_file << "\n";
    return false;
  }
  size_t matched = 0;
  const char *p = tm.data, *end = tm.data + tm.len;
  const char *fb[4], *fe[4];
  while (p < end) {
    const char *nl = (const char *)memchr(p, '\n', (size_t)(end - p));
    const char *le = nl ? nl : end;
    if (le > p && p[0] != '#') {
      int n = split_tabs(p, le, fb, fe, 2);
      if (n >= 2 && fe[1] > fb[1]) {
        std::string contig(fb[0], (size_t)(fe[0] - fb[0]));
        auto it = name2idx.find(contig);
        if (it != name2idx.end()) {
          std::string lin(fb[1], (size_t)(fe[1] - fb[1]));
          auto r = taxid.try_emplace(lin, (int)g_taxon_names.size());
          if (r.second) g_taxon_names.push_back(lin);
          g_contig_taxon[it->second] = r.first->second;
          ++matched;
        }
      }
    }
    if (!nl) break;
    p = nl + 1;
  }
  verbose_message("Taxonomy: %zu contigs labelled, %zu distinct lineages\n",
                  matched, g_taxon_names.size());
  return matched > 0;
}

void rb_bin_qc(const ContigVector &contigs, double &comp, double &cont) {
  comp = cont = 0.0;
  if (g_marker_set_size <= 0 || g_allcontig_markers.empty()) return;
  phmap::flat_hash_map<int, int> cnt;
  for (size_t c : contigs)
    if (c < g_allcontig_markers.size())
      for (int m : g_allcontig_markers[c]) cnt[m] += 1;
  long pres = 0, extra = 0;
  for (auto &mp : cnt) { if (mp.second >= 1) ++pres; if (mp.second > 1) extra += mp.second - 1; }
  const double invG = 100.0 / (double)g_marker_set_size;
  comp = (double)pres * invG;
  cont = (double)extra * invG;
}

int rb_bin_taxon(const ContigVector &contigs) {
  if (g_contig_taxon.empty()) return -1;
  phmap::flat_hash_map<int, size_t> votes;   // taxon id -> summed bp
  for (size_t c : contigs) {
    if (c >= g_contig_taxon.size()) continue;
    int t = g_contig_taxon[c];
    if (t < 0) continue;
    size_t len = (c < nobs) ? seq_lens[c] : small_seq_lens[c - nobs];
    votes[t] += len;
  }
  int best = -1; size_t best_bp = 0;
  for (auto &kv : votes)
    if (kv.second > best_bp || (kv.second == best_bp && (best < 0 || kv.first < best))) {
      best_bp = kv.second; best = kv.first;
    }
  return best;
}

void rb_purify_bins(BinMap &cls) {
  if (g_marker_set_size <= 0 || g_allcontig_markers.empty() || num_depth_samples < 1) {
    verbose_message("Purify: needs --markers and depth — skipped\n");
    return;
  }
  auto depth_at = [&](size_t c, size_t i) -> double {
    if (c < nobs) {
      if (!g_large_means.empty() && c * (size_t)num_depth_samples + i < g_large_means.size())
        return (double)g_large_means[c * num_depth_samples + i];
      return (double)depth_matrix(c, i);
    }
    size_t s = c - nobs;
    if (!g_small_means.empty() && s * (size_t)num_depth_samples + i < g_small_means.size())
      return (double)g_small_means[s * num_depth_samples + i];
    return (double)small_depth_matrix(s, i);
  };
  size_t n_removed = 0, n_bins_touched = 0;
  for (auto &kv : cls) {
    ContigVector &cs = kv.second;
    if (cs.size() < 3) continue;
    // duplicated single-copy markers in this bin?
    phmap::flat_hash_map<int, int> cnt;
    for (size_t c : cs)
      if (c < g_allcontig_markers.size())
        for (int m : g_allcontig_markers[c]) cnt[m] += 1;
    phmap::flat_hash_set<int> dup;
    for (auto &mp : cnt) if (mp.second > 1) dup.insert(mp.first);
    if (dup.empty()) continue;

    // Per-sample median log-depth centroid + MAD over the bin.
    const size_t S = num_depth_samples;
    std::vector<std::vector<double>> ld(cs.size(), std::vector<double>(S));
    for (size_t r = 0; r < cs.size(); ++r)
      for (size_t i = 0; i < S; ++i) ld[r][i] = std::log(depth_at(cs[r], i) + 1.0);
    std::vector<double> med(S), mad(S);
    std::vector<double> col(cs.size());
    for (size_t i = 0; i < S; ++i) {
      for (size_t r = 0; r < cs.size(); ++r) col[r] = ld[r][i];
      std::nth_element(col.begin(), col.begin() + col.size() / 2, col.end());
      med[i] = col[col.size() / 2];
      for (size_t r = 0; r < cs.size(); ++r) col[r] = std::fabs(ld[r][i] - med[i]);
      std::nth_element(col.begin(), col.begin() + col.size() / 2, col.end());
      mad[i] = col[col.size() / 2] * 1.4826 + 1e-9;
    }
    // A contig is a contaminant if it carries a duplicated marker AND is a
    // robust depth outlier (max |z| over samples > 3).
    ContigVector keep;
    keep.reserve(cs.size());
    bool touched = false;
    for (size_t r = 0; r < cs.size(); ++r) {
      bool has_dup = false;
      size_t c = cs[r];
      if (c < g_allcontig_markers.size())
        for (int m : g_allcontig_markers[c]) if (dup.count(m)) { has_dup = true; break; }
      double zmax = 0.0;
      for (size_t i = 0; i < S; ++i) {
        double z = std::fabs(ld[r][i] - med[i]) / mad[i];
        if (z > zmax) zmax = z;
      }
      if (has_dup && zmax > 3.0) { ++n_removed; touched = true; continue; }
      keep.push_back(c);
    }
    if (touched) { cs.swap(keep); ++n_bins_touched; }
  }
  verbose_message("Purify: removed %zu contaminant contig(s) from %zu bin(s)\n",
                  n_removed, n_bins_touched);
}

static int rb_cmd_qc(int ac, char *av[]) {
  using namespace rb_amber;

  std::string query_path, members_path, markers_path, out_path;
  size_t threads = 0;
  long set_size_opt = 0;        // 0 = auto (header comment, else distinct seen)
  double hq_comp = 90.0, hq_cont = 5.0;   // MIMAG high-quality thresholds
  double mq_comp = 50.0, mq_cont = 10.0;  // MIMAG medium-quality thresholds
  bool quiet_qc = false;

  po::options_description desc("rabbitbin qc options", 100, 50);
  desc.add_options()
      ("help,h", "Show help")
      ("binning,i", po::value<std::string>(&query_path), "Binning to score (CAMI bioboxes or 2-col SEQUENCEID<TAB>BINID)")
      ("members", po::value<std::string>(&members_path), "Binning as rabbitbin members.tsv (BinNum<TAB>SequenceName)")
      ("markers,k", po::value<std::string>(&markers_path), "Contig->marker map (contig<TAB>marker[...]; from rabbitbin_markers.sh) [required]")
      ("marker-set-size", po::value<long>(&set_size_opt)->default_value(0), "Total markers in the set G (0=auto: header comment or distinct seen)")
      ("output,o", po::value<std::string>(&out_path), "Write per-bin metrics TSV here (dir or file; optional)")
      ("hq-completeness", po::value<double>(&hq_comp)->default_value(90.0), "HQ completeness threshold (%)")
      ("hq-contamination", po::value<double>(&hq_cont)->default_value(5.0), "HQ contamination threshold (%)")
      ("mq-completeness", po::value<double>(&mq_comp)->default_value(50.0), "MQ completeness threshold (%)")
      ("mq-contamination", po::value<double>(&mq_cont)->default_value(10.0), "MQ contamination threshold (%)")
      ("threads,t", po::value<size_t>(&threads)->default_value(0), "Threads (0=all online CPUs)")
      ("quiet,q", po::value<bool>(&quiet_qc)->zero_tokens(), "Only print the summary");

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(desc).run(), vm);
  po::notify(vm);

  if (query_path.empty() && !members_path.empty()) query_path = members_path;
  const bool query_is_members = !members_path.empty();

  if (vm.count("help") || query_path.empty() || markers_path.empty()) {
    cerr << "\nrabbitbin qc: reference-free bin quality (SCG completeness/contamination)\n\n"
         << "Usage: rabbitbin qc --members out.members.tsv --markers contigs.markers.tsv\n"
         << "       rabbitbin qc -i prediction.binning -k contigs.markers.tsv -o qc/\n\n"
         << "Generate the marker map once per assembly:\n"
         << "       scripts/rabbitbin_markers.sh contigs.fa contigs.markers.tsv\n\n"
         << desc << "\n";
    return vm.count("help") ? 0 : 1;
  }

  long onln = sysconf(_SC_NPROCESSORS_ONLN);
  size_t hw = (onln > 0) ? (size_t)onln : 1;
  if (threads == 0) threads = hw; else threads = std::min(threads, hw);
  omp_set_num_threads((int)threads);

  auto t_start = std::chrono::steady_clock::now();
  auto secs = [&](std::chrono::steady_clock::time_point a) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  };

  // ── 1. Parse marker map → contig -> [markerId...], intern markers ──────────
  Mapped km;
  if (!map_file(markers_path, km)) { cerr << "[Error!] cannot open markers: " << markers_path << "\n"; return 1; }

  Interner markers;
  phmap::flat_hash_map<std::string, std::vector<int>> contig_markers;
  contig_markers.reserve(1u << 20);
  long set_size_hdr = 0;
  size_t marker_hits = 0;
  {
    const char *p = km.data, *end = km.data + km.len;
    const char *fb[64], *fe[64];
    while (p < end) {
      const char *nl = (const char *)memchr(p, '\n', (size_t)(end - p));
      const char *le = nl ? nl : end;
      if (le > p) {
        if (p[0] == '#') {
          // optional "# marker_set_size=N": scan the line for the '=' that
          // follows the keyword (no memmem; keep portable under _XOPEN_SOURCE).
          static const char *kw = "marker_set_size";
          const size_t kl = strlen(kw);
          for (const char *q = p; q + kl < le; ++q) {
            if (strncmp(q, kw, kl) == 0) {
              const char *eq = (const char *)memchr(q + kl, '=', (size_t)(le - (q + kl)));
              if (eq) set_size_hdr = (long)parse_uint(eq + 1, le);
              break;
            }
          }
        } else {
          int n = split_tabs(p, le, fb, fe, 64);
          if (n >= 2) {
            std::string contig(fb[0], (size_t)(fe[0] - fb[0]));
            auto &vec = contig_markers[contig];
            for (int c = 1; c < n; ++c) {
              if (fe[c] <= fb[c]) continue;          // skip empty fields
              vec.push_back(markers.intern(fb[c], fe[c]));
              ++marker_hits;
            }
          }
        }
      }
      if (!nl) break;
      p = nl + 1;
    }
  }

  long G = set_size_opt > 0 ? set_size_opt
         : (set_size_hdr > 0 ? set_size_hdr : (long)markers.names.size());
  if (G <= 0) { cerr << "[Error!] empty marker set\n"; return 1; }
  if (!quiet_qc) {
    fprintf(stderr, "[rabbitbin qc] markers: %zu contigs carry %zu hits, "
                    "%zu distinct; marker-set size G=%ld%s in %.2fs\n",
            contig_markers.size(), marker_hits, markers.names.size(), G,
            (set_size_opt > 0 ? " (--marker-set-size)"
                              : set_size_hdr > 0 ? " (header)" : " (distinct-seen)"),
            secs(t_start));
    if (set_size_opt == 0 && set_size_hdr == 0)
      fprintf(stderr, "[rabbitbin qc] WARNING: G inferred from distinct markers seen; "
                      "pass --marker-set-size for accurate completeness.\n");
  }

  // ── 2. Parse binning → seq2bin (intern bins) ──────────────────────────────
  auto t_q = std::chrono::steady_clock::now();
  Mapped qm;
  if (!map_file(query_path, qm)) { cerr << "[Error!] cannot open binning: " << query_path << "\n"; return 1; }

  Interner bins;
  phmap::flat_hash_map<std::string, int> seq2bin;
  seq2bin.reserve(1u << 20);
  {
    int seq_col = 0, bin_col = 1;
    if (query_is_members) { bin_col = 0; seq_col = 1; }
    const char *qp = qm.data, *qend = qm.data + qm.len;
    const char *fb2[8], *fe2[8];
    bool first_line = true;
    while (qp < qend) {
      const char *nl = (const char *)memchr(qp, '\n', (size_t)(qend - qp));
      const char *le = nl ? nl : qend;
      if (le > qp) {
        if (qp[0] == '@') {
          int sc, bc, lc;
          if (parse_cami_header(qp, le, sc, bc, lc)) { seq_col = sc; bin_col = bc; }
        } else if (first_line && (le - qp >= 6) && strncmp(qp, "BinNum", 6) == 0) {
          bin_col = 0; seq_col = 1; first_line = false;
        } else {
          int maxc = std::max(seq_col, bin_col) + 1;
          if (maxc > 8) maxc = 8;
          int n = split_tabs(qp, le, fb2, fe2, maxc);
          if (n > seq_col && n > bin_col) {
            int b = bins.intern(fb2[bin_col], fe2[bin_col]);
            std::string sk(fb2[seq_col], (size_t)(fe2[seq_col] - fb2[seq_col]));
            seq2bin.try_emplace(std::move(sk), b);   // first wins on multi-bin
          }
          first_line = false;
        }
      }
      if (!nl) break;
      qp = nl + 1;
    }
  }
  const int nbins = (int)bins.names.size();

  // Flatten contigs that have markers AND are binned: (bin, marker-list ptr).
  std::vector<int> ent_bin;
  std::vector<const std::vector<int> *> ent_mk;
  ent_bin.reserve(contig_markers.size());
  ent_mk.reserve(contig_markers.size());
  std::vector<long> bin_ncontig(nbins, 0);
  for (auto &kv : contig_markers) {
    auto it = seq2bin.find(kv.first);
    if (it == seq2bin.end()) continue;            // marker contig not in binning
    ent_bin.push_back(it->second);
    ent_mk.push_back(&kv.second);
  }
  // total contigs per bin (all binned contigs, for reporting)
  for (auto &kv : seq2bin) bin_ncontig[kv.second] += 1;
  const size_t NE = ent_bin.size();

  // ── 3. Per-bin marker copy counts (lock-free: partition bins across threads) ─
  std::vector<phmap::flat_hash_map<int, int>> bin_mk(nbins);
#pragma omp parallel
  {
    const int T = omp_get_num_threads();
    const int tid = omp_get_thread_num();
    for (size_t i = 0; i < NE; ++i) {
      int b = ent_bin[i];
      if ((b % T) != tid) continue;
      auto &cnt = bin_mk[b];
      for (int m : *ent_mk[i]) cnt[m] += 1;
    }
  }

  // ── 4. Completeness / contamination per bin ───────────────────────────────
  std::vector<double> comp(nbins, 0.0), cont(nbins, 0.0);
  std::vector<long>   present(nbins, 0);
  const double invG = 100.0 / (double)G;
#pragma omp parallel for schedule(dynamic, 256)
  for (int b = 0; b < nbins; ++b) {
    if (bin_mk[b].empty()) continue;
    long pres = 0; long extra = 0;
    for (auto &mp : bin_mk[b]) {
      if (mp.second >= 1) ++pres;
      if (mp.second > 1)  extra += (mp.second - 1);
    }
    present[b] = pres;
    comp[b] = (double)pres * invG;
    cont[b] = (double)extra * invG;
  }

  // ── 5. Tally MIMAG tiers ──────────────────────────────────────────────────
  long hq = 0, mq = 0, with_marker = 0;
  double sum_comp = 0.0, sum_cont = 0.0;
  for (int b = 0; b < nbins; ++b) {
    if (bin_mk[b].empty()) continue;
    ++with_marker;
    sum_comp += comp[b];
    sum_cont += cont[b];
    bool is_hq = comp[b] > hq_comp && cont[b] < hq_cont;
    bool is_mq = comp[b] >= mq_comp && cont[b] < mq_cont;
    if (is_hq) ++hq;
    else if (is_mq) ++mq;
  }

  // ── 6. Optional per-bin TSV ───────────────────────────────────────────────
  if (!out_path.empty()) {
    std::string fpath = out_path;
    struct stat st;
    bool is_dir = (stat(out_path.c_str(), &st) == 0 && S_ISDIR(st.st_mode));
    if (is_dir || (!out_path.empty() && out_path.back() == '/')) {
      if (!out_path.empty() && out_path.back() != '/') fpath = out_path + "/";
      else fpath = out_path;
      boost::filesystem::create_directories(fpath);
      fpath += "qc_per_bin.tsv";
    }
    FILE *f = fopen(fpath.c_str(), "w");
    if (!f) { cerr << "[Error!] cannot write metrics: " << fpath << "\n"; return 1; }
    fputs("BINID\tCompleteness\tContamination\tmarkers_present\tmarkers_total\tn_contigs\ttier\n", f);
    for (int b = 0; b < nbins; ++b) {
      const char *tier = "-";
      if (!bin_mk[b].empty()) {
        if (comp[b] > hq_comp && cont[b] < hq_cont) tier = "HQ";
        else if (comp[b] >= mq_comp && cont[b] < mq_cont) tier = "MQ";
        else tier = "LQ";
      }
      fprintf(f, "%s\t%.2f\t%.2f\t%ld\t%ld\t%ld\t%s\n",
              bins.names[b].c_str(), comp[b], cont[b], present[b], G,
              bin_ncontig[b], tier);
    }
    fclose(f);
    if (!quiet_qc) fprintf(stderr, "[rabbitbin qc] wrote per-bin metrics: %s\n", fpath.c_str());
  }

  if (!quiet_qc)
    fprintf(stderr, "[rabbitbin qc] binning: %zu binned seqs, %d bins, "
                    "%zu marker-bearing contigs matched in %.2fs\n",
            seq2bin.size(), nbins, NE, secs(t_q));

  // ── Summary (stdout) ──────────────────────────────────────────────────────
  printf("== rabbitbin qc summary ==\n");
  printf("marker set size G       : %ld\n", G);
  printf("total bins              : %d\n", nbins);
  printf("bins with >=1 marker    : %ld\n", with_marker);
  printf("HQ (comp>%.0f, cont<%.0f)  : %ld\n", hq_comp, hq_cont, hq);
  printf("MQ (comp>=%.0f, cont<%.0f) : %ld\n", mq_comp, mq_cont, mq);
  printf("avg completeness (%%)    : %.2f\n", with_marker ? sum_comp / with_marker : 0.0);
  printf("avg contamination (%%)   : %.2f\n", with_marker ? sum_cont / with_marker : 0.0);
  printf("total time              : %.2fs (%zu threads)\n", secs(t_start), threads);
  return 0;
}

// ═══════════════════════════════════════════════════════════════════════════
// `rabbitbin refine` — fast, multithreaded DAS Tool-style consensus binning.
//
// Given several INDEPENDENT binnings of the same assembly (e.g. RabbitBin +
// MetaBAT2 + MaxBin2) and an SCG marker map, pool every input bin as a
// candidate, then iteratively pick the highest-quality bin (SCG completeness
// penalised by contamination) and remove its contigs from all other candidates,
// re-scoring as we go. The result is a non-redundant consensus binning that is
// typically better than any single input — the same idea as DAS Tool, but in
// C++ and parallel (no DIAMOND/R; markers are precomputed once).
// ═══════════════════════════════════════════════════════════════════════════

namespace rb_refine {

// Parse one binning file (CAMI bioboxes, 2-col SEQ<TAB>BIN, or rabbitbin
// members.tsv) into bins, mapping contigs through the shared interner `cn`.
// Appends to `bins` (each bin = contig ids); `src` tags the source file.
static void parse_binning(const std::string &path, rb_amber::Interner &cn,
                          std::vector<std::vector<int>> &bins,
                          std::vector<int> &bin_src, int src, size_t &nseq) {
  using namespace rb_amber;
  Mapped m;
  if (!map_file(path, m)) { cerr << "[Error!] cannot open binning: " << path << "\n"; return; }
  phmap::flat_hash_map<std::string, int> binid;   // bin name -> local bin index
  int seq_col = 0, bin_col = 1;
  bool first_line = true;
  const char *p = m.data, *end = m.data + m.len;
  const char *fb[8], *fe[8];
  while (p < end) {
    const char *nl = (const char *)memchr(p, '\n', (size_t)(end - p));
    const char *le = nl ? nl : end;
    if (le > p) {
      if (p[0] == '@') {
        int sc, bc, lc;
        if (parse_cami_header(p, le, sc, bc, lc)) { seq_col = sc; bin_col = bc; }
      } else if (p[0] == '#') {
        // comment
      } else if (first_line && (le - p >= 6) && strncmp(p, "BinNum", 6) == 0) {
        bin_col = 0; seq_col = 1; first_line = false;
      } else {
        int maxc = std::max(seq_col, bin_col) + 1;
        if (maxc > 8) maxc = 8;
        int n = split_tabs(p, le, fb, fe, maxc);
        if (n > seq_col && n > bin_col && fe[bin_col] > fb[bin_col]) {
          std::string bn(fb[bin_col], (size_t)(fe[bin_col] - fb[bin_col]));
          auto r = binid.try_emplace(bn, (int)bins.size());
          if (r.second) { bins.emplace_back(); bin_src.push_back(src); }
          int cid = cn.intern(fb[seq_col], fe[seq_col]);
          bins[r.first->second].push_back(cid);
          ++nseq;
        }
        first_line = false;
      }
    }
    if (!nl) break;
    p = nl + 1;
  }
}

} // namespace rb_refine

static int rb_cmd_refine(int ac, char *av[]) {
  using namespace rb_amber;
  using namespace rb_refine;

  std::vector<std::string> binning_paths;
  std::vector<std::string> labels_csv;
  std::string markers_path, out_path;
  size_t threads = 0;
  long set_size_opt = 0;
  double min_comp = 50.0, max_cont = 10.0;   // accept threshold (post-overlap)
  double cont_weight = 5.0;                   // score = comp - cont_weight*cont
  bool quiet_r = false;

  po::options_description desc("rabbitbin refine options", 100, 50);
  desc.add_options()
      ("help,h", "Show help")
      ("binning,i", po::value<std::vector<std::string>>(&binning_paths)->multitoken(),
       "Input binning(s): bioboxes / 2-col SEQ<TAB>BIN / members.tsv. Repeat -i or pass several. [>=2 recommended]")
      ("labels,l", po::value<std::vector<std::string>>(&labels_csv)->multitoken(), "Optional source labels (one per --binning)")
      ("markers,k", po::value<std::string>(&markers_path), "Contig->marker map (from rabbitbin_markers.sh) [required]")
      ("marker-set-size", po::value<long>(&set_size_opt)->default_value(0), "Total markers G (0=auto)")
      ("output,o", po::value<std::string>(&out_path), "Output prefix (writes <prefix>.members.tsv + <prefix>.qc.tsv) [required]")
      ("min-completeness", po::value<double>(&min_comp)->default_value(50.0), "Keep consensus bins with completeness >= this (%)")
      ("max-contamination", po::value<double>(&max_cont)->default_value(10.0), "Keep consensus bins with contamination <= this (%)")
      ("contamination-weight", po::value<double>(&cont_weight)->default_value(5.0), "Ranking score = completeness - weight*contamination")
      ("threads,t", po::value<size_t>(&threads)->default_value(0), "Threads (0=all online CPUs)")
      ("quiet,q", po::value<bool>(&quiet_r)->zero_tokens(), "Only print the summary");

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(desc).run(), vm);
  po::notify(vm);

  if (vm.count("help") || binning_paths.empty() || markers_path.empty() || out_path.empty()) {
    cerr << "\nrabbitbin refine: DAS Tool-style SCG consensus over multiple binnings\n\n"
         << "Usage: rabbitbin refine -k contigs.markers.tsv -o consensus \\\n"
         << "         -i rabbitbin.members.tsv -i metabat2.binning -i maxbin2.binning\n\n"
         << desc << "\n";
    return vm.count("help") ? 0 : 1;
  }

  long onln = sysconf(_SC_NPROCESSORS_ONLN);
  size_t hw = (onln > 0) ? (size_t)onln : 1;
  if (threads == 0) threads = hw; else threads = std::min(threads, hw);
  omp_set_num_threads((int)threads);
  auto t0 = std::chrono::steady_clock::now();
  auto secs = [&](std::chrono::steady_clock::time_point a) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  };

  // ── 1. Markers → contig interner + contig->marker lists, G ────────────────
  Interner cn;                                  // shared contig name space
  std::vector<std::vector<int>> c2mk;           // [contig id] -> marker ids
  long G = 0;
  {
    Mapped km;
    if (!map_file(markers_path, km)) { cerr << "[Error!] cannot open markers: " << markers_path << "\n"; return 1; }
    Interner mk;
    long hdr = 0; size_t hits = 0;
    const char *p = km.data, *end = km.data + km.len;
    const char *fb[64], *fe[64];
    while (p < end) {
      const char *nl = (const char *)memchr(p, '\n', (size_t)(end - p));
      const char *le = nl ? nl : end;
      if (le > p) {
        if (p[0] == '#') {
          static const char *kw = "marker_set_size"; const size_t kl = strlen(kw);
          for (const char *q = p; q + kl < le; ++q)
            if (strncmp(q, kw, kl) == 0) {
              const char *eq = (const char *)memchr(q + kl, '=', (size_t)(le - (q + kl)));
              if (eq) hdr = (long)parse_uint(eq + 1, le);
              break;
            }
        } else {
          int n = split_tabs(p, le, fb, fe, 64);
          if (n >= 2) {
            int cid = cn.intern(fb[0], fe[0]);
            if ((int)c2mk.size() <= cid) c2mk.resize(cid + 1);
            for (int c = 1; c < n; ++c)
              if (fe[c] > fb[c]) { c2mk[cid].push_back(mk.intern(fb[c], fe[c])); ++hits; }
          }
        }
      }
      if (!nl) break;
      p = nl + 1;
    }
    G = set_size_opt > 0 ? set_size_opt : (hdr > 0 ? hdr : (long)mk.names.size());
    if (G <= 0) { cerr << "[Error!] empty marker set\n"; return 1; }
    if (!quiet_r)
      fprintf(stderr, "[rabbitbin refine] markers: %zu marker-bearing contigs, "
                      "%zu hits, G=%ld in %.2fs\n", c2mk.size(), hits, G, secs(t0));
  }

  // ── 2. Parse all input binnings into candidate bins ───────────────────────
  std::vector<std::vector<int>> cand;           // candidate bin -> contig ids
  std::vector<int> cand_src;
  for (size_t s = 0; s < binning_paths.size(); ++s) {
    size_t nseq = 0, before = cand.size();
    parse_binning(binning_paths[s], cn, cand, cand_src, (int)s, nseq);
    if (!quiet_r)
      fprintf(stderr, "[rabbitbin refine] input %zu (%s): %zu bins, %zu assignments\n",
              s, binning_paths[s].c_str(), cand.size() - before, nseq);
  }
  const size_t NC = cand.size();
  const size_t NV = cn.names.size();
  if (NC == 0) { cerr << "[Error!] no candidate bins parsed\n"; return 1; }
  c2mk.resize(NV);                              // contigs only in binnings: empty markers

  // ── 3. Iterative greedy consensus (DAS Tool core) ─────────────────────────
  const double invG = 100.0 / (double)G;
  std::vector<char> used(NV, 0);
  std::vector<char> alive(NC, 1);
  // current score cache; recomputed for candidates touched by the last pick.
  std::vector<double> sc_score(NC, -1e300);
  std::vector<double> sc_comp(NC, 0.0), sc_cont(NC, 0.0);
  std::vector<int>    sc_n(NC, 0);

  auto score_cand = [&](size_t b) {
    phmap::flat_hash_map<int, int> cnt;
    int nfree = 0;
    for (int c : cand[b]) {
      if (used[c]) continue;
      ++nfree;
      for (int m : c2mk[c]) cnt[m] += 1;
    }
    long pres = 0, extra = 0;
    for (auto &mp : cnt) { if (mp.second >= 1) ++pres; if (mp.second > 1) extra += mp.second - 1; }
    double comp = (double)pres * invG, cont = (double)extra * invG;
    sc_comp[b] = comp; sc_cont[b] = cont; sc_n[b] = nfree;
    sc_score[b] = (comp >= min_comp && cont <= max_cont) ? comp - cont_weight * cont
                                                         : -1e300;
  };

  // Initial scoring (parallel).
#pragma omp parallel for schedule(dynamic, 64)
  for (size_t b = 0; b < NC; ++b) score_cand(b);

  // Inverted index contig -> candidates, to re-score only affected bins.
  std::vector<std::vector<int>> c2cand(NV);
  for (size_t b = 0; b < NC; ++b)
    for (int c : cand[b]) c2cand[c].push_back((int)b);

  std::vector<std::vector<int>> out_bins;       // accepted consensus bins
  std::vector<int> out_src;
  while (true) {
    // argmax over alive candidates.
    size_t best = NC; double best_s = -1e299;
    for (size_t b = 0; b < NC; ++b)
      if (alive[b] && sc_score[b] > best_s) { best_s = sc_score[b]; best = b; }
    if (best == NC) break;                       // nothing meets thresholds

    // Accept the winner's currently-free contigs as a consensus bin.
    std::vector<int> bin;
    bin.reserve(sc_n[best]);
    for (int c : cand[best]) if (!used[c]) bin.push_back(c);
    if (bin.empty()) { alive[best] = 0; continue; }
    out_bins.push_back(bin);
    out_src.push_back(cand_src[best]);
    alive[best] = 0;

    // Mark contigs used; collect candidates needing a re-score.
    phmap::flat_hash_set<int> touched;
    for (int c : bin) {
      used[c] = 1;
      for (int b : c2cand[c]) if (alive[b]) touched.insert(b);
    }
    std::vector<int> tv(touched.begin(), touched.end());
#pragma omp parallel for schedule(dynamic, 16)
    for (size_t i = 0; i < tv.size(); ++i) {
      score_cand((size_t)tv[i]);
      if (sc_n[tv[i]] == 0) alive[tv[i]] = 0;
    }
  }

  // ── 4. Write consensus members.tsv + per-bin qc ───────────────────────────
  std::string mpath = out_path + ".members.tsv";
  std::string qpath = out_path + ".qc.tsv";
  FILE *mf = fopen(mpath.c_str(), "w");
  if (!mf) { cerr << "[Error!] cannot write " << mpath << "\n"; return 1; }
  FILE *qf = fopen(qpath.c_str(), "w");
  if (!qf) { cerr << "[Error!] cannot write " << qpath << "\n"; fclose(mf); return 1; }
  fputs("BinNum\tSequenceName\n", mf);
  fputs("BinNum\tCompleteness\tContamination\tn_contigs\tsource\ttier\n", qf);

  long hq = 0, mq = 0;
  size_t total_contigs = 0;
  for (size_t b = 0; b < out_bins.size(); ++b) {
    // final SCG quality of the accepted bin (contigs all free at accept time)
    phmap::flat_hash_map<int, int> cnt;
    for (int c : out_bins[b]) for (int m : c2mk[c]) cnt[m] += 1;
    long pres = 0, extra = 0;
    for (auto &mp : cnt) { if (mp.second >= 1) ++pres; if (mp.second > 1) extra += mp.second - 1; }
    double comp = (double)pres * invG, cont = (double)extra * invG;
    const char *tier = (comp > 90.0 && cont < 5.0) ? "HQ"
                     : (comp >= 50.0 && cont < 10.0) ? "MQ" : "LQ";
    if (comp > 90.0 && cont < 5.0) ++hq;
    else if (comp >= 50.0 && cont < 10.0) ++mq;
    int bn = (int)b + 1;
    for (int c : out_bins[b]) { fprintf(mf, "%d\t%s\n", bn, cn.names[c].c_str()); ++total_contigs; }
    fprintf(qf, "%d\t%.2f\t%.2f\t%zu\t%d\t%s\n", bn, comp, cont, out_bins[b].size(),
            out_src[b], tier);
  }
  fclose(mf); fclose(qf);

  if (!quiet_r)
    fprintf(stderr, "[rabbitbin refine] wrote %s (%zu contigs) and %s\n",
            mpath.c_str(), total_contigs, qpath.c_str());

  printf("== rabbitbin refine summary ==\n");
  printf("input binnings         : %zu\n", binning_paths.size());
  printf("candidate bins (pooled): %zu\n", NC);
  printf("consensus bins         : %zu\n", out_bins.size());
  printf("  HQ (comp>90,cont<5)  : %ld\n", hq);
  printf("  MQ (comp>=50,cont<10): %ld\n", mq);
  printf("marker set size G      : %ld\n", G);
  printf("total time             : %.2fs (%zu threads)\n", secs(t0), threads);
  return 0;
}
