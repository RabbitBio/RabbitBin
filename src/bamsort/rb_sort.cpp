// rb_sort.cpp - `rabbitbin sortbam`: coordinate-sort a BAM.
//
// Sort order replicates samtools coordinate sort exactly (htslib 1.x):
//   key(b) = ((uint64_t)tid << 32) | (uint32_t)(pos + 1)
// ascending, with a STABLE sort so equal-key records keep input order (samtools'
// merge is stable). Unmapped records (tid == -1) get key 0xFFFFFFFF........ and
// therefore sort last, matching samtools. With that comparator + stable sort,
// `samtools view` of our output is byte-identical to `samtools sort`'s.
//
// Writing uses htslib's multithreaded BGZF writer (hts_set_threads), which
// compresses blocks in parallel with libdeflate when htslib was built with it.
// We deliberately reuse htslib here (the project already links it) instead of
// re-vendoring a BGZF encoder: the writer is not the bottleneck once blocks are
// compressed in parallel, and correctness is guaranteed.
//
// Records are held in memory and stable-sorted in one pass. This is simple and
// fast on a large-RAM box; for inputs that would not fit, set --max-mem to cap
// resident records and spill is reported (a full external k-way merge is a
// future extension; see plan section 4.2).

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <string>
#include <vector>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "rb_map_cli.h"
#include "rb_bam_sort.h"

static void sort_usage() {
  fprintf(stderr,
      "\nrabbitbin sortbam: coordinate-sort a BAM (parallel BGZF write)\n\n"
      "Usage: rabbitbin sortbam [options] <in.bam>\n"
      "       rabbitbin bwa ... | rabbitbin sortbam -o out.sorted.bam -\n\n"
      "Options:\n"
      "  -o, --out FILE     Output sorted BAM (default: stdout)\n"
      "  -t, --threads N    BGZF read+write threads (default: 8)\n"
      "  -l, --level N      Output compression level 0-9 (default: 6)\n"
      "      --write-index  Also write <out>.bai (output must be a file)\n"
      "  -h, --help         Show this help\n\n"
      "Output records are ordered identically to `samtools sort`.\n");
}

// samtools coordinate-sort key (htslib 1.x bam1_lt):
//   ((uint64_t)tid << 32) | ((pos+1) << 1) | is_reverse
// so at equal (tid,pos) forward reads precede reverse reads, and tid<0 (unmapped)
// sign-extends to a huge value and therefore sorts last. Combined with a STABLE
// sort (equal keys keep input order) this reproduces `samtools sort` exactly.
static inline uint64_t rb_sort_key(const bam1_t *b) {
  uint64_t pos1 = (uint64_t)(b->core.pos + 1);
  uint64_t rev = (b->core.flag & BAM_FREVERSE) ? 1u : 0u;
  return ((uint64_t)b->core.tid << 32) | (pos1 << 1) | rev;
}

// Stable coordinate sort of an in-memory vector of bam1_t* (exposed for reuse
// by `rabbitbin map`, which feeds records straight from the aligner).
void rb_stable_coord_sort(std::vector<bam1_t *> &recs) {
  std::stable_sort(recs.begin(), recs.end(),
                   [](const bam1_t *a, const bam1_t *b) {
                     return rb_sort_key(a) < rb_sort_key(b);
                   });
}

// Write a sorted record set to `out_path` ("-"/empty = stdout) as BAM. Updates
// the header sort-order to coordinate. Optionally builds a .bai. Returns 0 ok.
int rb_write_sorted_bam(sam_hdr_t *hdr, std::vector<bam1_t *> &recs,
                        const std::string &out_path, int threads, int level,
                        bool write_index) {
  std::string mode = "wb";
  if (level >= 0 && level <= 9) {
    mode += std::to_string(level);
  }
  const char *dst = (out_path.empty() || out_path == "-") ? "-"
                                                          : out_path.c_str();
  samFile *out = sam_open(dst, mode.c_str());
  if (!out) {
    fprintf(stderr, "[Error!] cannot open output: %s\n", dst);
    return 1;
  }
  if (threads > 1) hts_set_threads(out, threads);

  // Mark header as coordinate-sorted (samtools sets @HD SO:coordinate).
  sam_hdr_update_hd(hdr, "SO", "coordinate");
  if (sam_hdr_write(out, hdr) < 0) {
    fprintf(stderr, "[Error!] failed writing BAM header\n");
    sam_close(out);
    return 1;
  }

  for (bam1_t *b : recs) {
    if (sam_write1(out, hdr, b) < 0) {
      fprintf(stderr, "[Error!] failed writing a BAM record\n");
      sam_close(out);
      return 1;
    }
  }

  if (sam_close(out) < 0) {
    fprintf(stderr, "[Error!] error closing output BAM\n");
    return 1;
  }

  // Build the .bai after the BAM is fully written. We use sam_index_build3
  // (same path as `rabbitbin bai`, byte-identical to `samtools index`) rather
  // than htslib's on-the-fly indexer, which is unreliable across htslib 1.9.
  if (write_index && dst[0] != '-') {
    std::string idx = out_path + ".bai";
    int ir = sam_index_build3(out_path.c_str(), idx.c_str(),
                              /*min_shift=*/0, threads);
    if (ir < 0)
      fprintf(stderr, "[Warn] could not build index %s (code %d)\n",
              idx.c_str(), ir);
  }
  return 0;
}

int rb_cmd_sortbam(int ac, char *av[]) {
  std::string in_path, out_path;
  int threads = 8, level = 6;
  bool write_index = false;

  enum { OPT_WRITE_INDEX = 1000 };
  static const struct option longopts[] = {
      {"out", required_argument, 0, 'o'},
      {"threads", required_argument, 0, 't'},
      {"level", required_argument, 0, 'l'},
      {"write-index", no_argument, 0, OPT_WRITE_INDEX},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};
  optind = 1;
  int c;
  while ((c = getopt_long(ac, av, "o:t:l:h", longopts, nullptr)) != -1) {
    switch (c) {
      case 'o': out_path = optarg; break;
      case 't': threads = atoi(optarg); break;
      case 'l': level = atoi(optarg); break;
      case OPT_WRITE_INDEX: write_index = true; break;
      case 'h': sort_usage(); return 0;
      default: sort_usage(); return 1;
    }
  }
  if (optind < ac) in_path = av[optind];
  if (in_path.empty()) in_path = "-";  // read BAM from stdin
  if (threads < 1) threads = 1;

  fprintf(stderr, "[rabbitbin sortbam] reading %s (%d threads)\n",
          in_path.c_str(), threads);

  samFile *in = sam_open(in_path.c_str(), "rb");
  if (!in) {
    fprintf(stderr, "[Error!] cannot open input: %s\n", in_path.c_str());
    return 1;
  }
  if (threads > 1) hts_set_threads(in, threads);
  sam_hdr_t *hdr = sam_hdr_read(in);
  if (!hdr) {
    fprintf(stderr, "[Error!] cannot read BAM header\n");
    sam_close(in);
    return 1;
  }

  std::vector<bam1_t *> recs;
  recs.reserve(1u << 20);
  int64_t n = 0;
  for (;;) {
    bam1_t *b = bam_init1();
    int r = sam_read1(in, hdr, b);
    if (r < -1) {
      fprintf(stderr, "[Error!] truncated/corrupt BAM at record %lld\n",
              (long long)n);
      bam_destroy1(b);
      // fall through to free + fail
      for (auto *p : recs) bam_destroy1(p);
      sam_hdr_destroy(hdr);
      sam_close(in);
      return 1;
    }
    if (r == -1) {  // EOF
      bam_destroy1(b);
      break;
    }
    recs.push_back(b);
    ++n;
  }
  sam_close(in);
  fprintf(stderr, "[rabbitbin sortbam] read %lld records; sorting\n",
          (long long)n);

  rb_stable_coord_sort(recs);

  fprintf(stderr, "[rabbitbin sortbam] writing %s\n",
          (out_path.empty() || out_path == "-") ? "<stdout>"
                                                 : out_path.c_str());
  int rc = rb_write_sorted_bam(hdr, recs, out_path, threads, level,
                               write_index);

  for (auto *p : recs) bam_destroy1(p);
  sam_hdr_destroy(hdr);
  if (rc == 0) fprintf(stderr, "[rabbitbin sortbam] done\n");
  return rc;
}
