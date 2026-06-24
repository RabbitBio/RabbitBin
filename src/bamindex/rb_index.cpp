// rb_index.cpp - `rabbitbin bai`: build a .bai index for a coordinate-sorted BAM.
//
// A BAI index is small and IO-light to build (one linear scan of the BGZF
// virtual offsets), so the win is not a novel algorithm but (a) a multithreaded
// BGZF reader to feed the index builder and (b) sharing the read pass with the
// sort writer in `rabbitbin map` (see rb_map.cpp --write-index). The standalone
// `bai` command uses htslib's sam_index_build3, which is already the reference
// implementation and produces a byte-equivalent .bai to `samtools index`.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <string>
#include <unistd.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/thread_pool.h>

#include "rb_map_cli.h"

static void bai_usage() {
  fprintf(stderr,
      "\nrabbitbin bai: build a .bai index for a coordinate-sorted BAM\n\n"
      "Usage: rabbitbin bai [options] <in.sorted.bam>\n\n"
      "Options:\n"
      "  -o, --out FILE   Output index path (default: <in.bam>.bai)\n"
      "  -t, --threads N  BGZF decompression threads (default: 4)\n"
      "  -h, --help       Show this help\n\n"
      "The BAM must be coordinate-sorted. Output is equivalent to "
      "`samtools index`.\n");
}

int rb_cmd_bai(int ac, char *av[]) {
  std::string in_path, out_path;
  int threads = 4;

  static const struct option longopts[] = {
      {"out", required_argument, 0, 'o'},
      {"threads", required_argument, 0, 't'},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};
  // Use getopt_long; reset optind so repeated subcommand parsing is clean.
  optind = 1;
  int c;
  while ((c = getopt_long(ac, av, "o:t:h", longopts, nullptr)) != -1) {
    switch (c) {
      case 'o': out_path = optarg; break;
      case 't': threads = atoi(optarg); break;
      case 'h': bai_usage(); return 0;
      default: bai_usage(); return 1;
    }
  }
  if (optind < ac) in_path = av[optind];

  if (in_path.empty()) {
    bai_usage();
    fprintf(stderr, "[Error!] an input sorted BAM is required\n");
    return 1;
  }
  if (out_path.empty()) out_path = in_path + ".bai";
  if (threads < 1) threads = 1;

  fprintf(stderr, "[rabbitbin bai] indexing %s -> %s (%d threads)\n",
          in_path.c_str(), out_path.c_str(), threads);

  // min_shift=0 selects the classic BAI format (samtools index default).
  int ret = sam_index_build3(in_path.c_str(), out_path.c_str(),
                             /*min_shift=*/0, threads);
  if (ret < 0) {
    fprintf(stderr,
            "[Error!] failed to build index (code %d). Is the BAM "
            "coordinate-sorted?\n",
            ret);
    return 1;
  }
  fprintf(stderr, "[rabbitbin bai] done\n");
  return 0;
}
