// rb_map.cpp - `rabbitbin map`: reads -> sorted + indexed BAM in one pass.
//
// Fuses align -> sort -> index in-process: the aligner hands each bam1_t to a
// sink that collects them, then we coordinate-sort (samtools-equivalent) and
// write the BAM with an on-the-fly .bai. No intermediate SAM text and no temp
// BAM on disk -- the main waste in a `bwa mem | samtools sort` pipeline.

#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <unistd.h>

#include <omp.h>

#include <htslib/sam.h>

#include "rb_map_cli.h"
#include "align/rb_align_api.h"
#include "bamsort/rb_bam_sort.h"
#ifdef RABBITBIN_ENABLE_STROBE
#include "align/rb_strobe.h"
#endif
#ifdef RABBITBIN_FUSE
#include "rabbit_depth_fuse.h"  // fused map->depth (--fused-depth)
#endif

namespace {

struct CollectCtx {
  std::vector<bam1_t *> *recs;
};

void collect_sink(void *ctx, bam1_t *b) {
  auto *c = static_cast<CollectCtx *>(ctx);
  c->recs->push_back(b);  // takes ownership; freed after writing
}

// One sample spec parsed from a --samples manifest (or built from -1/-2/-p/-o).
struct Sample {
  std::string out;          // output sorted BAM path
  std::string r1;           // reads (mate1, interleaved, or single-end)
  std::string r2;           // mate2 (empty if interleaved/single)
  bool interleaved = false; // r1 holds interleaved pairs
};

void map_usage() {
  fprintf(stderr,
      "\nrabbitbin map: reads -> sorted+indexed BAM (align|sort|index fused)\n\n"
      "Usage:\n"
      "  rabbitbin map -r ref.fa -1 r1.fq -2 r2.fq -o out.sorted.bam\n"
      "  rabbitbin map -r ref.fa -p interleaved.fq -o out.sorted.bam\n"
      "  rabbitbin map -r ref.fa --samples samples.tsv   # many samples, ONE index build\n\n"
      "Options:\n"
      "  -r, --ref FILE     Reference FASTA (required)\n"
      "  -1 FILE            Read 1 FASTQ(.gz)\n"
      "  -2 FILE            Read 2 FASTQ(.gz)\n"
      "  -p, --interleaved FILE   Interleaved paired FASTQ(.gz)\n"
      "  -o, --out FILE     Output sorted BAM (single-sample mode)\n"
      "      --samples FILE Manifest of samples to map against ONE shared index.\n"
      "                     Each non-empty/non-# line is TAB-separated:\n"
      "                       <out.bam> <reads1> [<reads2>]\n"
      "                     2 cols => reads1 is interleaved (-p); 3 cols => -1/-2.\n"
      "  -t, --threads N    Threads (default: all online CPUs)\n"
      "  -l, --level N      Output compression level 0-9 (default 6)\n"
      "  -k INT             Minimizer k (default 19)\n"
      "  -w INT             Minimizer window (default 19)\n"
      "      --band INT     Banded-DP half width (default 31)\n"
      "      --no-index     Do not write the .bai\n"
      "  -h, --help         Show this help\n");
}

// Align one sample's reads with the already-built (shared) index and coordinate
// -sort the records in memory; returns a heap-owned vector (caller writes then
// frees it) or nullptr on error. Writing is kept separate so a sample's BAM
// write can overlap the NEXT sample's alignment (see rb_cmd_map).
std::vector<bam1_t *> *align_sort_sample(const SaIndex &idx, const SaOpt &opt,
                                         sam_hdr_t *hdr, const Sample &s,
                                         int threads) {
  auto *recs = new std::vector<bam1_t *>();
  recs->reserve(1u << 22);
  CollectCtx ctx{recs};
  fprintf(stderr, "[rabbitbin map] %s: aligning -> in-memory records\n",
          s.out.c_str());
  int rc = rb_align_run(idx, opt, s.r1, s.r2, s.interleaved, threads,
                        /*out=*/nullptr, hdr, collect_sink, &ctx);
  if (rc != 0) {
    for (auto *b : *recs) bam_destroy1(b);
    delete recs;
    return nullptr;
  }
  fprintf(stderr, "[rabbitbin map] %s: sorting %zu records\n", s.out.c_str(),
          recs->size());
  rb_stable_coord_sort(*recs);
  return recs;
}

#ifdef RABBITBIN_ENABLE_STROBE
// Same as align_sort_sample but using the strobealign engine (shared mapper).
std::vector<bam1_t *> *align_sort_sample_strobe(StrobeMapper *sm,
                                                sam_hdr_t *hdr, const Sample &s,
                                                int threads) {
  auto *recs = new std::vector<bam1_t *>();
  recs->reserve(1u << 22);
  CollectCtx ctx{recs};
  fprintf(stderr, "[rabbitbin map] %s: aligning (strobe) -> in-memory records\n",
          s.out.c_str());
  int rc = rb_strobe_run(sm, s.r1, s.r2, s.interleaved, threads,
                         /*out=*/nullptr, hdr, collect_sink, &ctx);
  if (rc != 0) {
    for (auto *b : *recs) bam_destroy1(b);
    delete recs;
    return nullptr;
  }
  fprintf(stderr, "[rabbitbin map] %s: sorting %zu records\n", s.out.c_str(),
          recs->size());
  rb_stable_coord_sort(*recs);
  return recs;
}
#endif

// Parse a --samples manifest into a list of Sample specs. Returns false on a
// malformed line (and prints the offending line number).
bool parse_samples_manifest(const std::string &path,
                            std::vector<Sample> &out) {
  std::ifstream in(path);
  if (!in) {
    fprintf(stderr, "[Error!] cannot open --samples manifest: %s\n",
            path.c_str());
    return false;
  }
  std::string line;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    // strip trailing CR (tolerate CRLF manifests)
    if (!line.empty() && line.back() == '\r') line.pop_back();
    // skip blank lines and comments
    size_t first = line.find_first_not_of(" \t");
    if (first == std::string::npos || line[first] == '#') continue;

    std::vector<std::string> cols;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, '\t')) cols.push_back(field);
    if (cols.size() < 2 || cols.size() > 3 || cols[0].empty() ||
        cols[1].empty()) {
      fprintf(stderr,
              "[Error!] --samples line %d malformed (need "
              "<out.bam>\\t<reads1>[\\t<reads2>]): %s\n",
              lineno, line.c_str());
      return false;
    }
    Sample s;
    s.out = cols[0];
    s.r1 = cols[1];
    if (cols.size() == 3 && !cols[2].empty()) {
      s.r2 = cols[2];
      s.interleaved = false;
    } else {
      s.interleaved = true;  // single reads column => interleaved pairs
    }
    out.push_back(std::move(s));
  }
  if (out.empty()) {
    fprintf(stderr, "[Error!] --samples manifest %s has no samples\n",
            path.c_str());
    return false;
  }
  return true;
}

}  // namespace

int rb_cmd_map(int ac, char *av[]) {
  std::string ref, r1, r2, interleaved_path, out_path, samples_path;
  int threads = 0, level = 6, k = 19, w = 19, band = 31;
  bool no_index = false;
  bool use_rs = false;  // minimizer seeding by default (faster: smaller index)
  int read_len = 150;   // strobe engine: read-length parameter profile
#ifdef RABBITBIN_ENABLE_STROBE
  bool use_strobe = true;  // strobealign engine: faster + better binning
#else
  bool use_strobe = false;
#endif
  // --fused-depth: compute the MetaBAT depth TSV directly from each sample's
  // in-memory records (no BAM written/reread). Output is byte-identical to
  // `depth` on the written BAMs. Depth params mirror the `depth` defaults.
  std::string fused_depth_out;
  int fd_pctid = 97, fd_min_len = 1000;
  double fd_min_depth = 1.0;
  SaRsParams rs;       // defaults: s=16 wmin=5 wmax=11 maxdist=80 q=255

  enum {
    OPT_BAND = 1000,
    OPT_NOINDEX = 1001,
    OPT_SAMPLES = 1002,
    OPT_SEED = 1003,
    OPT_RS_S = 1004,
    OPT_RS_WMIN = 1005,
    OPT_RS_WMAX = 1006,
    OPT_RS_MAXDIST = 1007,
    OPT_ENGINE = 1008,
    OPT_READLEN = 1009,
    OPT_FUSED_DEPTH = 1010,
    OPT_FD_PCTID = 1011,
    OPT_FD_MINLEN = 1012,
    OPT_FD_MINDEPTH = 1013
  };
  static const struct option longopts[] = {
      {"ref", required_argument, 0, 'r'},
      {"interleaved", required_argument, 0, 'p'},
      {"out", required_argument, 0, 'o'},
      {"samples", required_argument, 0, OPT_SAMPLES},
      {"threads", required_argument, 0, 't'},
      {"level", required_argument, 0, 'l'},
      {"band", required_argument, 0, OPT_BAND},
      {"seed", required_argument, 0, OPT_SEED},
      {"engine", required_argument, 0, OPT_ENGINE},
      {"read-len", required_argument, 0, OPT_READLEN},
      {"fused-depth", required_argument, 0, OPT_FUSED_DEPTH},
      {"percent-identity", required_argument, 0, OPT_FD_PCTID},
      {"min-contig-length", required_argument, 0, OPT_FD_MINLEN},
      {"min-contig-depth", required_argument, 0, OPT_FD_MINDEPTH},
      {"rs-s", required_argument, 0, OPT_RS_S},
      {"rs-wmin", required_argument, 0, OPT_RS_WMIN},
      {"rs-wmax", required_argument, 0, OPT_RS_WMAX},
      {"rs-maxdist", required_argument, 0, OPT_RS_MAXDIST},
      {"no-index", no_argument, 0, OPT_NOINDEX},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};
  optind = 1;
  int c;
  while ((c = getopt_long(ac, av, "r:1:2:p:o:t:l:k:w:h", longopts, nullptr)) !=
         -1) {
    switch (c) {
      case 'r': ref = optarg; break;
      case '1': r1 = optarg; break;
      case '2': r2 = optarg; break;
      case 'p': interleaved_path = optarg; break;
      case 'o': out_path = optarg; break;
      case OPT_SAMPLES: samples_path = optarg; break;
      case 't': threads = atoi(optarg); break;
      case 'l': level = atoi(optarg); break;
      case 'k': k = atoi(optarg); break;
      case 'w': w = atoi(optarg); break;
      case OPT_BAND: band = atoi(optarg); break;
      case OPT_ENGINE:
        use_strobe = (std::string(optarg) == "strobe");
        break;
      case OPT_READLEN: read_len = atoi(optarg); break;
      case OPT_FUSED_DEPTH: fused_depth_out = optarg; break;
      case OPT_FD_PCTID: fd_pctid = atoi(optarg); break;
      case OPT_FD_MINLEN: fd_min_len = atoi(optarg); break;
      case OPT_FD_MINDEPTH: fd_min_depth = atof(optarg); break;
      case OPT_SEED:
        use_rs = (std::string(optarg) != "minimizer");
        break;
      case OPT_RS_S: rs.s = atoi(optarg); break;
      case OPT_RS_WMIN: rs.w_min = atoi(optarg); break;
      case OPT_RS_WMAX: rs.w_max = atoi(optarg); break;
      case OPT_RS_MAXDIST: rs.max_dist = atoi(optarg); break;
      case OPT_NOINDEX: no_index = true; break;
      case 'h': map_usage(); return 0;
      default: map_usage(); return 1;
    }
  }

  // Assemble the sample list: either a --samples manifest (many samples, one
  // shared index build) or a single sample from -1/-2/-p/-o (legacy path).
  std::vector<Sample> samples;
  if (!samples_path.empty()) {
    if (!r1.empty() || !r2.empty() || !interleaved_path.empty() ||
        !out_path.empty()) {
      fprintf(stderr,
              "[Error!] --samples is mutually exclusive with -1/-2/-p/-o\n");
      return 1;
    }
    if (!parse_samples_manifest(samples_path, samples)) return 1;
  } else {
    bool interleaved = !interleaved_path.empty();
    if (interleaved) r1 = interleaved_path;
    if (r1.empty() && optind < ac) r1 = av[optind];
    if (ref.empty() || r1.empty() || out_path.empty()) {
      map_usage();
      fprintf(stderr,
              "[Error!] --ref, reads, and --out (or --samples) are required\n");
      return 1;
    }
    Sample s;
    s.out = out_path;
    s.r1 = r1;
    s.r2 = r2;
    s.interleaved = interleaved;
    samples.push_back(std::move(s));
  }
  if (ref.empty()) {
    map_usage();
    fprintf(stderr, "[Error!] --ref is required\n");
    return 1;
  }

  long onln = sysconf(_SC_NPROCESSORS_ONLN);
  if (threads <= 0) threads = (onln > 0) ? (int)onln : 1;
  // Make OpenMP-backed parallel algorithms (e.g. the parallel stable sort in
  // rb_stable_coord_sort, __gnu_parallel) honour the requested thread count.
  omp_set_num_threads(threads);

  // Build the index ONCE and reuse it (and the header) for every sample.
  // Two engines: the vendored strobealign (default; faster + better binning)
  // or the legacy minimizer seed-and-extend (--engine minimizer).
  SaIndex idx;
  SaOpt opt;
  sam_hdr_t *hdr = nullptr;
#ifdef RABBITBIN_ENABLE_STROBE
  StrobeMapper *sm = nullptr;
  if (use_strobe) {
    sm = rb_strobe_build(ref, read_len, threads);
    if (!sm) return 1;
    hdr = rb_strobe_make_header(sm);
  } else
#endif
  {
    idx.use_randstrobe = use_rs;
    idx.rs = rs;
    if (!idx.build(ref, k, w, threads)) return 1;
    opt.band = band;
    hdr = rb_align_make_header(idx);
  }

  // Per sample: align + sort (uses all cores via the aligner pipeline) on the
  // main thread, then hand the sorted records to a SINGLE background writer so
  // the parallel BGZF write overlaps the NEXT sample's alignment. At most one
  // write is in flight, so the shared header is never touched concurrently
  // (alignment in `map` mode never uses it). The .bai is built on the MAIN
  // thread right after the write joins (htslib's index thread pool does not
  // engage from a worker thread, where it would fall back to a slow serial
  // re-read). Each sample's output is content-identical to a serial run.
  int rc = 0;
  std::thread writer;
  std::atomic<int> write_rc{0};
  std::string pending_bai;  // out path whose .bai is owed once `writer` joins
#ifdef RABBITBIN_FUSE
  std::vector<DepthColumn> fd_cols;   // fused: one depth column per sample
  std::vector<std::string> fd_names;  // fused: column labels (= sample out path)
  const bool fused = !fused_depth_out.empty();
#else
  const bool fused = false;
#endif
  for (size_t i = 0; i < samples.size(); ++i) {
    if (samples.size() > 1)
      fprintf(stderr, "[rabbitbin map] === sample %zu/%zu: %s ===\n", i + 1,
              samples.size(), samples[i].out.c_str());
    std::vector<bam1_t *> *recs = nullptr;
#ifdef RABBITBIN_ENABLE_STROBE
    if (use_strobe)
      recs = align_sort_sample_strobe(sm, hdr, samples[i], threads);
    else
#endif
      recs = align_sort_sample(idx, opt, hdr, samples[i], threads);
    if (!recs) {
      rc = 1;
      fprintf(stderr, "[Error!] sample %s failed; aborting\n",
              samples[i].out.c_str());
      break;
    }
#ifdef RABBITBIN_FUSE
    if (fused) {
      // No BAM written: accumulate this sample's depth column from the sorted
      // in-memory records, then free them. Byte-identical to writing the BAM
      // and running `depth` (same records, same per-read filter/overlap).
      fprintf(stderr, "[rabbitbin map] %s: accumulating depth (fused)\n",
              samples[i].out.c_str());
      DepthColumn col = depth_accumulate_sample(
          *recs, hdr, (float)fd_pctid / 100.0f, /*maxEdgeBases=*/75,
          /*includeEdgeBases=*/false, /*minMapQual=*/0);
      fd_cols.push_back(std::move(col));
      fd_names.push_back(samples[i].out);
      for (auto *b : *recs) bam_destroy1(b);
      delete recs;
      continue;
    }
#endif
    // Finish the previous sample's write, then index it on the main thread.
    if (writer.joinable()) writer.join();
    if (write_rc.load() != 0) {
      rc = write_rc.load();
      for (auto *b : *recs) bam_destroy1(b);
      delete recs;
      break;
    }
    if (!pending_bai.empty()) {
      rb_index_bam(pending_bai, threads);
      pending_bai.clear();
    }
    Sample s = samples[i];
    fprintf(stderr, "[rabbitbin map] %s: writing (background)\n", s.out.c_str());
    writer = std::thread([&, recs, s]() {
      int r = rb_write_sorted_bam(hdr, *recs, s.out, threads, level,
                                  /*write_index=*/false);
      rb_free_records(*recs, threads);
      delete recs;
      if (r != 0) write_rc.store(r);
    });
    if (!no_index) pending_bai = s.out;
  }
  if (writer.joinable()) writer.join();
  if (rc == 0 && write_rc.load() != 0) rc = write_rc.load();
  if (rc == 0 && !pending_bai.empty()) rb_index_bam(pending_bai, threads);

#ifdef RABBITBIN_FUSE
  if (rc == 0 && fused) {
    std::string tsv = depth_format_table(
        fd_cols, fd_names, hdr, (float)fd_pctid / 100.0f, fd_min_len,
        (float)fd_min_depth, /*maxEdgeBases=*/75, /*includeEdgeBases=*/false,
        /*intraDepthVariance=*/true);
    FILE *f = fopen(fused_depth_out.c_str(), "wb");
    if (f) {
      fwrite(tsv.data(), 1, tsv.size(), f);
      fclose(f);
      fprintf(stderr, "[rabbitbin map] fused depth -> %s (%zu bytes)\n",
              fused_depth_out.c_str(), tsv.size());
    } else {
      fprintf(stderr, "[Error!] cannot open --fused-depth output: %s\n",
              fused_depth_out.c_str());
      rc = 1;
    }
  }
#endif

  sam_hdr_destroy(hdr);
#ifdef RABBITBIN_ENABLE_STROBE
  if (sm) rb_strobe_free(sm);
#endif
  if (rc == 0) fprintf(stderr, "[rabbitbin map] done (%zu sample(s))\n",
                       samples.size());
  return rc;
}
