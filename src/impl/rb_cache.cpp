// RabbitBin module: rb_cache.cpp
//
// Persistent cache for "build once, re-bin many times" (Priority 2).
//
// The expensive stages — FASTA parse, PMH sketch, depth load + rank transform,
// and the O(N^2) similarity-graph construction — produce a compact state that
// the cheap tail (edgeScore -> incidence -> label propagation -> recruit ->
// split -> output) consumes.  rb_write_cache() snapshots exactly that state to a
// little-endian binary file; rb_load_cache() restores it so a subsequent run
// skips straight to the cheap tail under different parameters.
//
// Format choice: a hand-rolled binary dump (not Boost.serialization).  The
// payload is a handful of flat POD vectors, two string vectors and two ublas
// matrices, so a direct fwrite/fread of the contiguous storage is both simplest
// and fastest, and sidesteps ublas-serialization friction.  An ABI guard (the
// sizes of StoredDistance and size_t) plus a magic + version prevent loading a
// mismatched cache.
//
// NOT cached (recomputed at load, kept derivable to keep the file small and
// avoid staleness): g_depth_unit, the small-contig unit vectors, the anynz
// cache, g_small_means (recruit rebuilds it from the raw small_depth_matrix),
// and edgeScore / incidence (the cheap tail).  Contig sequences are likewise
// not cached, so a --load-cache run emits the membership tables (the default
// output) but not per-bin FASTA / unbinned FASTA.

namespace {

template <class T> static inline void rbc_wr_pod(FILE *f, const T &v) {
  fwrite(&v, sizeof(T), 1, f);
}
template <class T> static inline bool rbc_rd_pod(FILE *f, T &v) {
  return fread(&v, sizeof(T), 1, f) == 1;
}

template <class T>
static inline void rbc_wr_vec(FILE *f, const std::vector<T> &v) {
  uint64_t n = v.size();
  fwrite(&n, sizeof(n), 1, f);
  if (n) fwrite(v.data(), sizeof(T), n, f);
}
template <class T> static inline bool rbc_rd_vec(FILE *f, std::vector<T> &v) {
  uint64_t n = 0;
  if (fread(&n, sizeof(n), 1, f) != 1) return false;
  v.resize(n);
  if (n && fread(v.data(), sizeof(T), n, f) != n) return false;
  return true;
}

static inline void rbc_wr_str(FILE *f, const std::string &s) {
  uint32_t n = (uint32_t)s.size();
  fwrite(&n, sizeof(n), 1, f);
  if (n) fwrite(s.data(), 1, n, f);
}
static inline bool rbc_rd_str(FILE *f, std::string &s) {
  uint32_t n = 0;
  if (fread(&n, sizeof(n), 1, f) != 1) return false;
  s.resize(n);
  if (n && fread(&s[0], 1, n, f) != n) return false;
  return true;
}

static inline void rbc_wr_strvec(FILE *f, const std::vector<std::string> &v) {
  uint64_t n = v.size();
  fwrite(&n, sizeof(n), 1, f);
  for (const auto &s : v) rbc_wr_str(f, s);
}
static inline bool rbc_rd_strvec(FILE *f, std::vector<std::string> &v) {
  uint64_t n = 0;
  if (fread(&n, sizeof(n), 1, f) != 1) return false;
  v.resize(n);
  for (auto &s : v) if (!rbc_rd_str(f, s)) return false;
  return true;
}

// ublas matrix<StoredDistance> with the default (row-major) layout stores its
// elements contiguously, reachable via data()[0]; dump rows, cols, then bytes.
static inline void rbc_wr_mat(FILE *f, const Matrix &m) {
  uint64_t r = m.size1(), c = m.size2();
  fwrite(&r, sizeof(r), 1, f);
  fwrite(&c, sizeof(c), 1, f);
  if (r && c) fwrite(&m.data()[0], sizeof(StoredDistance), r * c, f);
}
static inline bool rbc_rd_mat(FILE *f, Matrix &m) {
  uint64_t r = 0, c = 0;
  if (fread(&r, sizeof(r), 1, f) != 1) return false;
  if (fread(&c, sizeof(c), 1, f) != 1) return false;
  m.resize(r, c, false);
  if (r && c &&
      fread(&m.data()[0], sizeof(StoredDistance), r * c, f) != r * c)
    return false;
  return true;
}

static const char RBC_MAGIC[8] = {'R', 'B', 'C', 'A', 'C', 'H', 'E', '1'};
// v1: original (no SNV).  v2: appends an optional per-large-contig SNV block
// (feature #5) so --strain results are reusable from cache.  v1 caches still
// load (the SNV block is simply absent).
static const uint32_t RBC_VERSION = 2;

} // namespace

bool rb_write_cache(const std::string &path, const Graph &g, bool has_depth) {
  FILE *f = fopen(path.c_str(), "wb");
  if (!f) {
    cerr << "[Error!] cannot open cache for writing: " << path << "\n";
    return false;
  }
  fwrite(RBC_MAGIC, 1, 8, f);
  rbc_wr_pod(f, RBC_VERSION);
  // ABI guard: element sizes must match on load.
  uint32_t sz_sd = (uint32_t)sizeof(StoredDistance);
  uint32_t sz_st = (uint32_t)sizeof(size_t);
  rbc_wr_pod(f, sz_sd);
  rbc_wr_pod(f, sz_st);

  // ── Scalars ──────────────────────────────────────────────────────────────
  uint64_t v_nobs = nobs, v_nobs1 = nobs1, v_nds = num_depth_samples;
  rbc_wr_pod(f, v_nobs);
  rbc_wr_pod(f, v_nobs1);
  rbc_wr_pod(f, v_nds);
  unsigned long long v_seed = seed, v_ts = totalSize, v_ts1 = totalSize1;
  rbc_wr_pod(f, v_seed);
  rbc_wr_pod(f, v_ts);
  rbc_wr_pod(f, v_ts1);
  double v_sim = (double)simCutoff, v_b0 = g_pmh_baseline;
  rbc_wr_pod(f, v_sim);
  rbc_wr_pod(f, v_b0);
  uint64_t v_mc = minContig, v_msc = min_small_contig;
  rbc_wr_pod(f, v_mc);
  rbc_wr_pod(f, v_msc);
  uint8_t v_hd = has_depth ? 1 : 0;
  rbc_wr_pod(f, v_hd);

  // ── Per-contig arrays ────────────────────────────────────────────────────
  rbc_wr_strvec(f, contig_names);
  rbc_wr_strvec(f, small_contig_names);
  rbc_wr_vec(f, seq_lens);
  rbc_wr_vec(f, small_seq_lens);

  // ── Depth state ──────────────────────────────────────────────────────────
  // depth_matrix is RANKED (as left by the sketch loop); small_depth_matrix is
  // RAW (the recruit stage ranks it in place on both the normal and cache path,
  // preserving the original accumulation order).
  rbc_wr_mat(f, depth_matrix);
  rbc_wr_mat(f, small_depth_matrix);
  rbc_wr_vec(f, g_large_means);   // raw large means (split)
  rbc_wr_vec(f, g_depth_raw);     // raw depths for weighted-Jaccard metric
  rbc_wr_vec(f, g_depth_colnorm); // per-sample normalisation (may be empty)

  // ── Similarity graph topology ────────────────────────────────────────────
  rbc_wr_vec(f, g.from);
  rbc_wr_vec(f, g.to);
  rbc_wr_vec(f, g.sComp);

  // ── SNV / strain block (v2; feature #5) ──────────────────────────────────
  // Row-aligned with contig_names (nobs × num_depth_samples).  Stored so a
  // re-bin from cache reproduces <prefix>.snv.tsv / <prefix>.strain.tsv without
  // re-scanning the BAMs.  Flag = 0 means "no SNV data in this cache".
  uint8_t v_hs = (g_strain_scan && !g_snv_pi.empty()) ? 1 : 0;
  rbc_wr_pod(f, v_hs);
  if (v_hs) {
    rbc_wr_vec(f, g_snv_pi);
    rbc_wr_vec(f, g_snv_nsnv);
    rbc_wr_vec(f, g_snv_cov);
  }

  bool ok = (ferror(f) == 0);
  fclose(f);
  if (!ok) {
    cerr << "[Error!] failed while writing cache: " << path << "\n";
    return false;
  }
  return true;
}

bool rb_load_cache(const std::string &path) {
  FILE *f = fopen(path.c_str(), "rb");
  if (!f) {
    cerr << "[Error!] cannot open cache for reading: " << path << "\n";
    return false;
  }
  char magic[8];
  if (fread(magic, 1, 8, f) != 8 || memcmp(magic, RBC_MAGIC, 8) != 0) {
    cerr << "[Error!] not a RabbitBin cache file: " << path << "\n";
    fclose(f);
    return false;
  }
  uint32_t ver = 0;
  if (!rbc_rd_pod(f, ver) || ver < 1 || ver > RBC_VERSION) {
    cerr << "[Error!] unsupported cache version (" << ver << ", expected 1.."
         << RBC_VERSION << ")\n";
    fclose(f);
    return false;
  }
  uint32_t sz_sd = 0, sz_st = 0;
  if (!rbc_rd_pod(f, sz_sd) || !rbc_rd_pod(f, sz_st) ||
      sz_sd != (uint32_t)sizeof(StoredDistance) ||
      sz_st != (uint32_t)sizeof(size_t)) {
    cerr << "[Error!] cache ABI mismatch (built with a different RabbitBin)\n";
    fclose(f);
    return false;
  }

  bool ok = true;
  uint64_t v_nobs = 0, v_nobs1 = 0, v_nds = 0;
  ok &= rbc_rd_pod(f, v_nobs);
  ok &= rbc_rd_pod(f, v_nobs1);
  ok &= rbc_rd_pod(f, v_nds);
  unsigned long long v_seed = 0, v_ts = 0, v_ts1 = 0;
  ok &= rbc_rd_pod(f, v_seed);
  ok &= rbc_rd_pod(f, v_ts);
  ok &= rbc_rd_pod(f, v_ts1);
  double v_sim = 0.0, v_b0 = 0.0;
  ok &= rbc_rd_pod(f, v_sim);
  ok &= rbc_rd_pod(f, v_b0);
  uint64_t v_mc = 0, v_msc = 0;
  ok &= rbc_rd_pod(f, v_mc);
  ok &= rbc_rd_pod(f, v_msc);
  uint8_t v_hd = 0;
  ok &= rbc_rd_pod(f, v_hd);
  if (!ok) {
    cerr << "[Error!] truncated cache header: " << path << "\n";
    fclose(f);
    return false;
  }

  nobs = (size_t)v_nobs;
  nobs1 = (size_t)v_nobs1;
  num_depth_samples = (size_t)v_nds;
  g_cache_seed = v_seed;
  totalSize = v_ts;
  totalSize1 = v_ts1;
  simCutoff = (Similarity)v_sim;
  g_pmh_baseline = v_b0;
  g_cache_has_depth = (v_hd != 0);
  (void)v_mc;  // build-time minContig (informational; downstream uses CLI value)
  (void)v_msc;

  ok &= rbc_rd_strvec(f, contig_names);
  ok &= rbc_rd_strvec(f, small_contig_names);
  ok &= rbc_rd_vec(f, seq_lens);
  ok &= rbc_rd_vec(f, small_seq_lens);
  ok &= rbc_rd_mat(f, depth_matrix);
  ok &= rbc_rd_mat(f, small_depth_matrix);
  ok &= rbc_rd_vec(f, g_large_means);
  ok &= rbc_rd_vec(f, g_depth_raw);
  ok &= rbc_rd_vec(f, g_depth_colnorm);
  ok &= rbc_rd_vec(f, g_cache_from);
  ok &= rbc_rd_vec(f, g_cache_to);
  ok &= rbc_rd_vec(f, g_cache_scomp);

  // ── SNV / strain block (v2; feature #5) ──────────────────────────────────
  g_cache_has_snv = false;
  if (ok && ver >= 2) {
    uint8_t v_hs = 0;
    ok &= rbc_rd_pod(f, v_hs);
    if (ok && v_hs) {
      ok &= rbc_rd_vec(f, g_snv_pi);
      ok &= rbc_rd_vec(f, g_snv_nsnv);
      ok &= rbc_rd_vec(f, g_snv_cov);
      if (ok) {
        g_cache_has_snv = true;
        g_strain_scan = true;  // re-emit snv/strain tables from cached vectors
      }
    }
  }

  fclose(f);
  if (!ok) {
    cerr << "[Error!] truncated or corrupt cache body: " << path << "\n";
    return false;
  }
  return true;
}
