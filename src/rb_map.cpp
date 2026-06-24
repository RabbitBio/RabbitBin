// rb_map.cpp - `rabbitbin map`: reads -> sorted + indexed BAM in one pass.
//
// Fuses align -> sort -> index in-process: the aligner hands each bam1_t to a
// sink that collects them, then we coordinate-sort (samtools-equivalent) and
// write the BAM with an on-the-fly .bai. No intermediate SAM text and no temp
// BAM on disk -- the main waste in a `bwa mem | samtools sort` pipeline.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <string>
#include <vector>

#include <unistd.h>

#include <htslib/sam.h>

#include "rb_map_cli.h"
#include "align/rb_align_api.h"
#include "bamsort/rb_bam_sort.h"

namespace {

struct CollectCtx {
  std::vector<bam1_t *> *recs;
};

void collect_sink(void *ctx, bam1_t *b) {
  auto *c = static_cast<CollectCtx *>(ctx);
  c->recs->push_back(b);  // takes ownership; freed after writing
}

void map_usage() {
  fprintf(stderr,
      "\nrabbitbin map: reads -> sorted+indexed BAM (align|sort|index fused)\n\n"
      "Usage:\n"
      "  rabbitbin map -r ref.fa -1 r1.fq -2 r2.fq -o out.sorted.bam\n"
      "  rabbitbin map -r ref.fa -p interleaved.fq -o out.sorted.bam\n\n"
      "Options:\n"
      "  -r, --ref FILE     Reference FASTA (required)\n"
      "  -1 FILE            Read 1 FASTQ(.gz)\n"
      "  -2 FILE            Read 2 FASTQ(.gz)\n"
      "  -p, --interleaved FILE   Interleaved paired FASTQ(.gz)\n"
      "  -o, --out FILE     Output sorted BAM (required)\n"
      "  -t, --threads N    Threads (default: all online CPUs)\n"
      "  -l, --level N      Output compression level 0-9 (default 6)\n"
      "  -k INT             Minimizer k (default 15)\n"
      "  -w INT             Minimizer window (default 10)\n"
      "      --band INT     Banded-DP half width (default 31)\n"
      "      --no-index     Do not write the .bai\n"
      "  -h, --help         Show this help\n");
}

}  // namespace

int rb_cmd_map(int ac, char *av[]) {
  std::string ref, r1, r2, interleaved_path, out_path;
  int threads = 0, level = 6, k = 15, w = 10, band = 31;
  bool no_index = false;

  enum { OPT_BAND = 1000, OPT_NOINDEX = 1001 };
  static const struct option longopts[] = {
      {"ref", required_argument, 0, 'r'},
      {"interleaved", required_argument, 0, 'p'},
      {"out", required_argument, 0, 'o'},
      {"threads", required_argument, 0, 't'},
      {"level", required_argument, 0, 'l'},
      {"band", required_argument, 0, OPT_BAND},
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
      case 't': threads = atoi(optarg); break;
      case 'l': level = atoi(optarg); break;
      case OPT_BAND: band = atoi(optarg); break;
      case OPT_NOINDEX: no_index = true; break;
      case 'h': map_usage(); return 0;
      default: map_usage(); return 1;
    }
  }
  bool interleaved = !interleaved_path.empty();
  if (interleaved) r1 = interleaved_path;
  if (r1.empty() && optind < ac) r1 = av[optind];

  if (ref.empty() || r1.empty() || out_path.empty()) {
    map_usage();
    fprintf(stderr, "[Error!] --ref, reads, and --out are required\n");
    return 1;
  }
  long onln = sysconf(_SC_NPROCESSORS_ONLN);
  if (threads <= 0) threads = (onln > 0) ? (int)onln : 1;

  SaIndex idx;
  if (!idx.build(ref, k, w, threads)) return 1;
  SaOpt opt;
  opt.band = band;

  sam_hdr_t *hdr = rb_align_make_header(idx);

  // Align, collecting records in memory.
  std::vector<bam1_t *> recs;
  recs.reserve(1u << 22);
  CollectCtx ctx{&recs};
  fprintf(stderr, "[rabbitbin map] aligning -> in-memory records\n");
  int rc = rb_align_run(idx, opt, r1, r2, interleaved, threads,
                        /*out=*/nullptr, hdr, collect_sink, &ctx);
  if (rc != 0) {
    for (auto *b : recs) bam_destroy1(b);
    sam_hdr_destroy(hdr);
    return rc;
  }

  fprintf(stderr, "[rabbitbin map] sorting %zu records\n", recs.size());
  rb_stable_coord_sort(recs);

  fprintf(stderr, "[rabbitbin map] writing %s%s\n", out_path.c_str(),
          no_index ? "" : " (+ .bai)");
  rc = rb_write_sorted_bam(hdr, recs, out_path, threads, level, !no_index);

  for (auto *b : recs) bam_destroy1(b);
  sam_hdr_destroy(hdr);
  if (rc == 0) fprintf(stderr, "[rabbitbin map] done\n");
  return rc;
}
