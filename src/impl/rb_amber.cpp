// RabbitBin module: rb_amber.cpp
//
// `rabbitbin amber` — a fast, multithreaded reimplementation of AMBER's genome-
// binning evaluation that reproduces AMBER's per-bin purity/completeness and the
// HQ/MQ/LQ recovered-bin counts, without pandas/matplotlib (AMBER spends minutes
// mostly in DataFrame merges + plotting).
//
// Faithful to cami_amber.binning_classes.GenomeQuery (genome binning, simple
// 1-contig-1-bin case, which is what every binner here produces):
//   confusion[bin][genome] = sum of contig lengths in `bin` belonging to `genome`
//   most_abundant_genome(bin) = argmax_genome confusion[bin][genome]   (idxmax)
//   tp_length            = confusion[bin][most_abundant_genome]
//   total_length(bin)    = sum of all contig lengths in bin (that are in the GS)
//   length_gs(genome)    = sum of contig lengths of the genome in the GS
//   Purity (bp)       = precision_bp = tp_length / total_length
//   Completeness (bp) = recall_bp    = tp_length / length_gs(most_abundant_genome)
// Recovered-bin counts (AMBER calc_num_recovered_genomes, strict >):
//   HQ: recall_bp > 0.90 AND precision_bp > 0.95
//   MQ: recall_bp > 0.70 AND precision_bp > 0.90
//   LQ: recall_bp > 0.50 AND precision_bp > 0.90
//
// Length comes from the gold standard (the query needs only SEQUENCEID->BINID),
// exactly like AMBER. Query sequences absent from the GS are dropped.
//
// Threads are used for the heavy phases: genome resolution of query contigs
// (read-only concurrent lookups into the finished GS map), the confusion build
// (bins partitioned across threads so each writes a disjoint shard, lock-free),
// and the per-bin metric pass.

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace rb_amber {

struct Mapped {
  const char *data = nullptr;
  size_t len = 0;
  int fd = -1;
  bool ok() const { return data != nullptr; }
  ~Mapped() {
    if (data && data != MAP_FAILED) munmap((void *)data, len);
    if (fd >= 0) close(fd);
  }
};

static bool map_file(const std::string &path, Mapped &m) {
  m.fd = open(path.c_str(), O_RDONLY);
  if (m.fd < 0) return false;
  struct stat st;
  if (fstat(m.fd, &st) != 0) return false;
  m.len = (size_t)st.st_size;
  if (m.len == 0) { m.data = ""; return true; }  // empty file
  void *p = mmap(nullptr, m.len, PROT_READ, MAP_PRIVATE, m.fd, 0);
  if (p == MAP_FAILED) return false;
  m.data = (const char *)p;
  return true;
}

// Parse an unsigned integer from [p,end); stop at first non-digit. Returns value.
static inline uint64_t parse_uint(const char *p, const char *end) {
  uint64_t v = 0;
  while (p < end && *p >= '0' && *p <= '9') { v = v * 10 + (uint64_t)(*p - '0'); ++p; }
  return v;
}

// Split a line [s,e) into up to maxcols tab-delimited field [begin,end) ranges.
// Returns the number of fields found.
static inline int split_tabs(const char *s, const char *e,
                             const char **fb, const char **fe, int maxcols) {
  int n = 0;
  const char *cur = s;
  while (n < maxcols) {
    const char *t = (const char *)memchr(cur, '\t', (size_t)(e - cur));
    fb[n] = cur;
    if (!t || t >= e) { fe[n] = e; ++n; break; }
    fe[n] = t;
    ++n;
    cur = t + 1;
  }
  return n;
}

// Locate the column indices for SEQUENCEID / BINID / length from a CAMI bioboxes
// "@@SEQUENCEID\tBINID\t..." header line. Returns false if not a header line.
static bool parse_cami_header(const char *s, const char *e,
                              int &seq_col, int &bin_col, int &len_col) {
  if (e - s < 2 || s[0] != '@' || s[1] != '@') return false;
  const char *p = s + 2;  // skip "@@"
  int col = 0;
  seq_col = bin_col = len_col = -1;
  while (p < e) {
    const char *t = (const char *)memchr(p, '\t', (size_t)(e - p));
    const char *fend = (t && t < e) ? t : e;
    size_t L = (size_t)(fend - p);
    auto eq = [&](const char *kw) {
      size_t kl = strlen(kw);
      return L == kl && strncmp(p, kw, kl) == 0;
    };
    if (eq("SEQUENCEID") || eq("@@SEQUENCEID")) seq_col = col;
    else if (eq("BINID")) bin_col = col;
    else if (eq("_LENGTH") || eq("LENGTH")) len_col = col;
    ++col;
    if (!t || t >= e) break;
    p = t + 1;
  }
  if (seq_col < 0) seq_col = 0;
  if (bin_col < 0) bin_col = 1;
  return true;
}

// String interner: name -> dense id, with a parallel vector of names.
struct Interner {
  phmap::flat_hash_map<std::string, int> id;
  std::vector<std::string> names;
  int intern(const char *b, const char *e) {
    std::string k(b, (size_t)(e - b));
    auto it = id.find(k);
    if (it != id.end()) return it->second;
    int v = (int)names.size();
    id.emplace(std::move(k), v);
    names.push_back(std::string(b, (size_t)(e - b)));
    return v;
  }
};

struct SeqRec { int genome; uint32_t len; };

} // namespace rb_amber

static int rb_cmd_amber(int ac, char *av[]) {
  using namespace rb_amber;

  std::string gold_path, query_path, out_path, members_path;
  size_t threads = 0;
  size_t min_length = 0;      // GS contig length filter (AMBER --min_length)
  bool quiet_amber = false;

  po::options_description desc("rabbitbin amber options", 100, 50);
  desc.add_options()
      ("help,h", "Show help")
      ("gold,g", po::value<std::string>(&gold_path), "Gold standard binning (CAMI bioboxes; needs _LENGTH) [required]")
      ("binning,i", po::value<std::string>(&query_path), "Predicted binning (CAMI bioboxes or 2-col SEQUENCEID<TAB>BINID)")
      ("members", po::value<std::string>(&members_path), "Predicted binning as rabbitbin members.tsv (BinNum<TAB>SequenceName)")
      ("output,o", po::value<std::string>(&out_path), "Write per-bin metrics TSV here (dir or file; optional)")
      ("min-length", po::value<size_t>(&min_length)->default_value(0), "Ignore GS contigs shorter than this")
      ("threads,t", po::value<size_t>(&threads)->default_value(0), "Threads (0=all online CPUs)")
      ("quiet,q", po::value<bool>(&quiet_amber)->zero_tokens(), "Only print the summary");

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(desc).run(), vm);
  po::notify(vm);

  if (query_path.empty() && !members_path.empty()) query_path = members_path;
  const bool query_is_members = !members_path.empty();

  if (vm.count("help") || gold_path.empty() || query_path.empty()) {
    cerr << "\nrabbitbin amber: fast genome-binning evaluation (AMBER-compatible)\n\n"
         << "Usage: rabbitbin amber -g gold.binning -i prediction.binning\n"
         << "       rabbitbin amber -g gold.binning --members out.members.tsv -o metrics/\n\n"
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

  // ── 1. Parse gold standard → seq2gen + genome sizes ───────────────────────
  Mapped gm;
  if (!map_file(gold_path, gm)) { cerr << "[Error!] cannot open gold: " << gold_path << "\n"; return 1; }

  Interner genomes;
  phmap::flat_hash_map<std::string, SeqRec> seq2gen;
  seq2gen.reserve(1u << 20);
  {
    int seq_col = 0, bin_col = 1, len_col = 2;
    bool header_seen = false;
    const char *p = gm.data, *end = gm.data + gm.len;
    const char *fb[8], *fe[8];
    size_t dropped_short = 0, dup = 0;
    while (p < end) {
      const char *nl = (const char *)memchr(p, '\n', (size_t)(end - p));
      const char *le = nl ? nl : end;
      if (le > p) {
        if (p[0] == '@') {
          int sc, bc, lc;
          if (parse_cami_header(p, le, sc, bc, lc)) {
            seq_col = sc; bin_col = bc;
            if (lc >= 0) len_col = lc;
            header_seen = true;
          }
        } else {
          int maxc = std::max(std::max(seq_col, bin_col), len_col) + 1;
          if (maxc > 8) maxc = 8;
          int n = split_tabs(p, le, fb, fe, maxc);
          if (n > seq_col && n > bin_col && n > len_col) {
            uint64_t L = parse_uint(fb[len_col], fe[len_col]);
            if (min_length == 0 || L >= min_length) {
              int g = genomes.intern(fb[bin_col], fe[bin_col]);
              std::string sk(fb[seq_col], (size_t)(fe[seq_col] - fb[seq_col]));
              auto r = seq2gen.try_emplace(std::move(sk), SeqRec{g, (uint32_t)L});
              if (!r.second) dup++;  // keep first (AMBER drop_duplicates SEQUENCEID)
            } else dropped_short++;
          }
        }
      }
      if (!nl) break;
      p = nl + 1;
    }
    (void)header_seen;
    // Genome sizes from the unique GS sequences (matches AMBER groupby on gs_df).
    std::vector<double> gsize(genomes.names.size(), 0.0);
    std::vector<long> gseqcnt(genomes.names.size(), 0);
    for (auto &kv : seq2gen) {
      gsize[kv.second.genome] += kv.second.len;
      gseqcnt[kv.second.genome] += 1;
    }
    if (!quiet_amber)
      fprintf(stderr, "[rabbitbin amber] gold: %zu sequences, %zu genomes "
                      "(%zu dup-seqid kept-first, %zu < min-length) in %.2fs\n",
              seq2gen.size(), genomes.names.size(), dup, dropped_short, secs(t_start));

    // ── 2. Parse query → unique seqid->bin (dedup, first wins) ──────────────
    auto t_q = std::chrono::steady_clock::now();
    Mapped qm;
    if (!map_file(query_path, qm)) { cerr << "[Error!] cannot open binning: " << query_path << "\n"; return 1; }

    Interner bins;
    phmap::flat_hash_map<std::string, int> seq2bin;
    seq2bin.reserve(1u << 20);
    size_t conflicts = 0;
    {
      int seq_col = 0, bin_col = 1;
      if (query_is_members) { bin_col = 0; seq_col = 1; }  // BinNum, SequenceName
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
            // rabbitbin members.tsv header (auto-detected or via --members): skip
            // the header row and use BinNum (col 0) / SequenceName (col 1).
            bin_col = 0; seq_col = 1;
            first_line = false;
          } else {
            int maxc = std::max(seq_col, bin_col) + 1;
            if (maxc > 8) maxc = 8;
            int n = split_tabs(qp, le, fb2, fe2, maxc);
            if (n > seq_col && n > bin_col) {
              int b = bins.intern(fb2[bin_col], fe2[bin_col]);
              std::string sk(fb2[seq_col], (size_t)(fe2[seq_col] - fb2[seq_col]));
              auto r = seq2bin.try_emplace(std::move(sk), b);
              if (!r.second && r.first->second != b) conflicts++;  // multi-bin: keep first
            }
            first_line = false;
          }
        }
        if (!nl) break;
        qp = nl + 1;
      }
    }

    // Flatten to arrays; resolve each query seq's genome+len once (read-only,
    // concurrent) and drop sequences absent from the GS.
    std::vector<std::string> uq_keys;
    std::vector<int> uq_bin;
    uq_keys.reserve(seq2bin.size());
    uq_bin.reserve(seq2bin.size());
    for (auto &kv : seq2bin) { uq_keys.push_back(kv.first); uq_bin.push_back(kv.second); }
    const size_t M = uq_keys.size();

    std::vector<int> uq_gen(M, -1);
    std::vector<uint32_t> uq_len(M, 0);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < M; ++i) {
      auto it = seq2gen.find(uq_keys[i]);
      if (it != seq2gen.end()) { uq_gen[i] = it->second.genome; uq_len[i] = it->second.len; }
    }

    const int nbins = (int)bins.names.size();
    std::vector<double> bin_total(nbins, 0.0);
    std::vector<long>   bin_seqcnt(nbins, 0);
    std::vector<phmap::flat_hash_map<int, double>> bin_gen(nbins);

    // ── 3. Confusion build, partitioned by bin across threads (lock-free) ───
    // Each thread owns bins where (bin % T == tid); it scans all records but
    // only writes its own bins, so no two threads touch the same map.
#pragma omp parallel
    {
      const int T = omp_get_num_threads();
      const int tid = omp_get_thread_num();
      for (size_t i = 0; i < M; ++i) {
        if (uq_gen[i] < 0) continue;
        int b = uq_bin[i];
        if ((b % T) != tid) continue;
        bin_total[b] += uq_len[i];
        bin_seqcnt[b] += 1;
        bin_gen[b][uq_gen[i]] += uq_len[i];
      }
    }

    // ── 4. Per-bin metrics (dominant genome, purity, completeness) ──────────
    std::vector<int>    bin_domgen(nbins, -1);
    std::vector<double> bin_tp(nbins, 0.0);
    std::vector<double> bin_prec(nbins, 0.0);
    std::vector<double> bin_recall(nbins, 0.0);
#pragma omp parallel for schedule(dynamic, 256)
    for (int b = 0; b < nbins; ++b) {
      if (bin_gen[b].empty()) continue;
      int best_g = -1; double best_bp = -1.0;
      for (auto &gp : bin_gen[b]) {
        // Tie-break: larger bp, then smaller genome id (stable, deterministic).
        if (gp.second > best_bp || (gp.second == best_bp && gp.first < best_g)) {
          best_bp = gp.second; best_g = gp.first;
        }
      }
      bin_domgen[b] = best_g;
      bin_tp[b] = best_bp;
      bin_prec[b] = (bin_total[b] > 0.0) ? best_bp / bin_total[b] : 0.0;
      double gs = gsize[best_g];
      bin_recall[b] = (gs > 0.0) ? best_bp / gs : 0.0;
    }

    // ── 5. Recovered-bin counts (AMBER strict >) + averages ─────────────────
    long hq = 0, mq = 0, lq = 0, nonempty = 0;
    double sum_prec = 0.0, sum_tp = 0.0, sum_total = 0.0;
    for (int b = 0; b < nbins; ++b) {
      if (bin_gen[b].empty()) continue;
      ++nonempty;
      sum_prec += bin_prec[b];
      sum_tp += bin_tp[b];
      sum_total += bin_total[b];
      if (bin_recall[b] > 0.90 && bin_prec[b] > 0.95) ++hq;
      if (bin_recall[b] > 0.70 && bin_prec[b] > 0.90) ++mq;
      if (bin_recall[b] > 0.50 && bin_prec[b] > 0.90) ++lq;
    }

    // Per-genome best-bin completeness (AMBER recall_df: avg over all GS genomes).
    std::vector<double> gen_best_recall(genomes.names.size(), 0.0);
    for (int b = 0; b < nbins; ++b) {
      if (bin_gen[b].empty()) continue;
      for (auto &gp : bin_gen[b]) {
        double r = (gsize[gp.first] > 0.0) ? gp.second / gsize[gp.first] : 0.0;
        if (r > gen_best_recall[gp.first]) gen_best_recall[gp.first] = r;
      }
    }
    double sum_genrecall = 0.0;
    for (double r : gen_best_recall) sum_genrecall += r;
    const size_t ngen = genomes.names.size();

    // ── 6. Optional per-bin TSV ─────────────────────────────────────────────
    if (!out_path.empty()) {
      std::string fpath = out_path;
      struct stat st;
      bool is_dir = (stat(out_path.c_str(), &st) == 0 && S_ISDIR(st.st_mode));
      if (is_dir || (!out_path.empty() && out_path.back() == '/')) {
        if (!out_path.empty() && out_path.back() != '/') fpath = out_path + "/";
        else fpath = out_path;
        boost::filesystem::create_directories(fpath);
        fpath += "metrics_per_bin.tsv";
      }
      FILE *f = fopen(fpath.c_str(), "w");
      if (!f) { cerr << "[Error!] cannot write metrics: " << fpath << "\n"; return 1; }
      fputs("BINID\tmost_abundant_genome\tPurity (bp)\tCompleteness (bp)\t"
            "total_length\ttp_length\ttotal_seq_counts\n", f);
      for (int b = 0; b < nbins; ++b) {
        if (bin_gen[b].empty()) continue;
        fprintf(f, "%s\t%s\t%.6f\t%.6f\t%.0f\t%.0f\t%ld\n",
                bins.names[b].c_str(),
                genomes.names[bin_domgen[b]].c_str(),
                bin_prec[b], bin_recall[b], bin_total[b], bin_tp[b], bin_seqcnt[b]);
      }
      fclose(f);
      if (!quiet_amber)
        fprintf(stderr, "[rabbitbin amber] wrote per-bin metrics: %s\n", fpath.c_str());
    }

    if (!quiet_amber)
      fprintf(stderr, "[rabbitbin amber] query: %zu unique seqs (%zu not in GS, "
                      "%zu multi-bin kept-first), %d bins in %.2fs\n",
              M, M - (size_t)std::count_if(uq_gen.begin(), uq_gen.end(),
                                           [](int g){ return g >= 0; }),
              conflicts, nbins, secs(t_q));

    // ── Summary (stdout) ────────────────────────────────────────────────────
    double prec_weighted = (sum_total > 0.0) ? sum_tp / sum_total : 0.0;
    printf("== rabbitbin amber summary ==\n");
    printf("genomes (gold)         : %zu\n", ngen);
    printf("predicted bins         : %ld\n", nonempty);
    printf("HQ (comp>0.90,pur>0.95): %ld\n", hq);
    printf("MQ (comp>0.70,pur>0.90): %ld\n", mq);
    printf("LQ (comp>0.50,pur>0.90): %ld\n", lq);
    printf("avg purity (bp)        : %.4f\n", nonempty ? sum_prec / nonempty : 0.0);
    printf("weighted purity (bp)   : %.4f\n", prec_weighted);
    printf("avg completeness/genome: %.4f\n", ngen ? sum_genrecall / ngen : 0.0);
    printf("total time             : %.2fs (%zu threads)\n", secs(t_start), threads);
  }
  return 0;
}
