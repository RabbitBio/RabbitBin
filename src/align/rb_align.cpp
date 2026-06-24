// rb_align.cpp - `rabbitbin bwa`: align short reads to a reference (CPU).
//
// A lightweight, self-contained seed-and-extend aligner (minimizer seeds +
// banded extension; see sa_index/sa_map) wrapped as a rabbitbin subcommand.
// It is NOT a drop-in for bwa-mem output (different seeding => different ties),
// but produces valid coordinate alignments suitable for coverage/binning. The
// engine is intentionally simple so it can be optimized further in place.
//
// Reads are processed in batches: a batch is loaded, mapped in parallel with
// OpenMP, then written in input order (so output is deterministic). Paired-end
// (-1/-2 or interleaved -p) sets mate fields and the proper-pair flag.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <string>
#include <vector>

#include <omp.h>
#include <zlib.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/kseq.h>

#include "rb_map_cli.h"
#include "sa_index.h"
#include "sa_map.h"
#include "rb_align_api.h"

KSEQ_INIT(gzFile, gzread)

namespace {

struct ReadRec {
  std::string name, seq, qual;
};

// Reverse-complement a sequence string in place into `out`.
void revcomp_str(const std::string &in, std::string &out) {
  int n = (int)in.size();
  out.resize(n);
  for (int i = 0; i < n; ++i) {
    char c = in[n - 1 - i];
    switch (c) {
      case 'A': case 'a': out[i] = 'T'; break;
      case 'C': case 'c': out[i] = 'G'; break;
      case 'G': case 'g': out[i] = 'C'; break;
      case 'T': case 't': out[i] = 'A'; break;
      default: out[i] = 'N';
    }
  }
}

// Build a bam1_t for one mapped/unmapped read. Caller owns `b`.
void fill_bam(bam1_t *b, const ReadRec &r, const SaAln &a, uint16_t extra_flag,
              int mtid, int64_t mpos, int64_t isize, bool mate_rev,
              bool mate_unmapped) {
  uint16_t flag = extra_flag;
  std::string seq = r.seq, qual = r.qual;
  std::vector<uint32_t> cigar;
  int tid = -1;
  int64_t pos = -1;
  uint8_t mapq = 0;

  if (a.mapped) {
    tid = a.contig;
    pos = a.pos;
    mapq = (uint8_t)a.mapq;
    cigar = a.cigar;
    if (a.rev) {
      flag |= BAM_FREVERSE;
      std::string rs, rq;
      revcomp_str(r.seq, rs);
      seq = rs;
      if (!r.qual.empty()) {
        rq.assign(r.qual.rbegin(), r.qual.rend());
        qual = rq;
      }
    }
  } else {
    flag |= BAM_FUNMAP;
  }
  if (mate_rev) flag |= BAM_FMREVERSE;
  if (mate_unmapped) flag |= BAM_FMUNMAP;

  const char *qual_ptr = nullptr;
  std::string qbuf;
  if (!qual.empty() && qual[0] != '*') {
    // htslib expects phred values, not ASCII; convert.
    qbuf.resize(qual.size());
    for (size_t i = 0; i < qual.size(); ++i)
      qbuf[i] = (char)(qual[i] - 33);
    qual_ptr = qbuf.data();
  }

  bam_set1(b, r.name.size(), r.name.c_str(), flag, tid, pos, mapq,
           cigar.size(), cigar.empty() ? nullptr : cigar.data(), mtid, mpos,
           isize, seq.size(), seq.c_str(), qual_ptr, 0);

  if (a.mapped) {
    int32_t nm = a.nm;
    bam_aux_append(b, "NM", 'i', sizeof(int32_t), (uint8_t *)&nm);
  }
}

void align_usage() {
  fprintf(stderr,
      "\nrabbitbin bwa: align reads to a reference (CPU seed-and-extend)\n\n"
      "Usage:\n"
      "  rabbitbin bwa -r ref.fa -1 r1.fq [-2 r2.fq] [-o out.bam]\n"
      "  rabbitbin bwa -r ref.fa -p interleaved.fq -o out.bam\n"
      "  rabbitbin bwa -r ref.fa reads.fq                       # single-end to stdout SAM\n\n"
      "Options:\n"
      "  -r, --ref FILE     Reference FASTA (required)\n"
      "  -1 FILE            Read 1 FASTQ(.gz)\n"
      "  -2 FILE            Read 2 FASTQ(.gz)\n"
      "  -p, --interleaved FILE   Interleaved paired FASTQ(.gz)\n"
      "  -o, --out FILE     Output (default stdout). .bam => BAM, else SAM\n"
      "  -t, --threads N    Threads (default: all online CPUs)\n"
      "  -k INT             Minimizer k (default 15)\n"
      "  -w INT             Minimizer window (default 10)\n"
      "      --band INT     Banded-DP half width (default 31)\n"
      "  -h, --help         Show this help\n");
}

}  // namespace

// Build a SAM header (sam_hdr_t) with @SQ lines from the index. Exposed so
// `rabbitbin map` can reuse it.
sam_hdr_t *rb_align_make_header(const SaIndex &idx) {
  sam_hdr_t *h = sam_hdr_init();
  std::string hd = "@HD\tVN:1.6\tSO:unsorted\n";
  sam_hdr_add_lines(h, hd.c_str(), hd.size());
  std::string sq;
  sq.reserve(idx.contigs.size() * 48);
  for (const auto &c : idx.contigs) {
    sq += "@SQ\tSN:";
    sq += c.name;
    sq += "\tLN:";
    sq += std::to_string(c.len);
    sq += "\n";
  }
  sam_hdr_add_lines(h, sq.c_str(), sq.size());
  std::string pg = "@PG\tID:rabbitbin\tPN:rabbitbin\tVN:bwa\n";
  sam_hdr_add_lines(h, pg.c_str(), pg.size());
  return h;
}

// Core: align reads from r1 (and optional r2, or interleaved) using `idx`,
// invoking sink(b) for each finished record in input order. If sink is null,
// records are written to `out` (already opened). Returns 0 on success.
int rb_align_run(const SaIndex &idx, const SaOpt &opt, const std::string &r1,
                 const std::string &r2, bool interleaved, int threads,
                 samFile *out, sam_hdr_t *hdr, RbAlignSink sink,
                 void *sink_ctx) {
  bool paired = interleaved || !r2.empty();
  gzFile f1 = gzopen(r1.c_str(), "rb");
  if (!f1) {
    fprintf(stderr, "[Error!] cannot open reads: %s\n", r1.c_str());
    return 1;
  }
  gzFile f2 = nullptr;
  if (!r2.empty()) {
    f2 = gzopen(r2.c_str(), "rb");
    if (!f2) {
      fprintf(stderr, "[Error!] cannot open reads: %s\n", r2.c_str());
      gzclose(f1);
      return 1;
    }
  }
  kseq_t *k1 = kseq_init(f1);
  kseq_t *k2 = f2 ? kseq_init(f2) : nullptr;

  const int BATCH = 1 << 16;  // read pairs (or singles) per batch
  std::vector<ReadRec> A, B;  // mate1, mate2 (B empty if single)
  A.reserve(BATCH);
  B.reserve(BATCH);

  auto load_one = [](kseq_t *ks, ReadRec &out) -> bool {
    int r = kseq_read(ks);
    if (r < 0) return false;
    out.name.assign(ks->name.s, ks->name.l);
    out.seq.assign(ks->seq.s, ks->seq.l);
    if (ks->qual.l)
      out.qual.assign(ks->qual.s, ks->qual.l);
    else
      out.qual = "*";
    return true;
  };

  int64_t total = 0, mapped = 0;
  int rc = 0;
  for (;;) {
    A.clear();
    B.clear();
    int got = 0;
    while (got < BATCH) {
      ReadRec a;
      if (!load_one(k1, a)) break;
      if (paired) {
        ReadRec b;
        bool ok;
        if (interleaved)
          ok = load_one(k1, b);
        else
          ok = load_one(k2, b);
        if (!ok) {
          fprintf(stderr, "[Warn] odd number of paired reads; dropping last\n");
          break;
        }
        A.push_back(std::move(a));
        B.push_back(std::move(b));
      } else {
        A.push_back(std::move(a));
      }
      ++got;
    }
    if (got == 0) break;

    std::vector<SaAln> ra(A.size()), rb(paired ? B.size() : 0);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 256)
    for (int i = 0; i < (int)A.size(); ++i) {
      ra[i] = sa_align_read(idx, opt, A[i].seq.c_str(), (int)A[i].seq.size());
      if (paired)
        rb[i] = sa_align_read(idx, opt, B[i].seq.c_str(), (int)B[i].seq.size());
    }

    // Emit records in input order.
    for (int i = 0; i < (int)A.size(); ++i) {
      ++total;
      if (ra[i].mapped) ++mapped;
      if (paired) {
        ++total;
        if (rb[i].mapped) ++mapped;
        const SaAln &a1 = ra[i];
        const SaAln &a2 = rb[i];
        // mate fields
        int m1tid = a2.mapped ? a2.contig : -1;
        int64_t m1pos = a2.mapped ? a2.pos : -1;
        int m2tid = a1.mapped ? a1.contig : -1;
        int64_t m2pos = a1.mapped ? a1.pos : -1;
        bool proper = a1.mapped && a2.mapped && a1.contig == a2.contig &&
                      a1.rev != a2.rev;
        int64_t isize1 = 0, isize2 = 0;
        if (a1.mapped && a2.mapped && a1.contig == a2.contig) {
          int64_t lo = std::min(a1.pos, a2.pos);
          int64_t hi = std::max(a1.pos, a2.pos);
          // approximate end span using read lengths
          int64_t span = hi - lo + (int64_t)std::max(A[i].seq.size(), B[i].seq.size());
          isize1 = (a1.pos <= a2.pos) ? span : -span;
          isize2 = -isize1;
        }
        uint16_t f1flag = BAM_FPAIRED | BAM_FREAD1 | (proper ? BAM_FPROPER_PAIR : 0);
        uint16_t f2flag = BAM_FPAIRED | BAM_FREAD2 | (proper ? BAM_FPROPER_PAIR : 0);

        bam1_t *b1 = bam_init1();
        bam1_t *b2 = bam_init1();
        fill_bam(b1, A[i], a1, f1flag, m1tid, m1pos, isize1, a2.rev,
                 !a2.mapped);
        fill_bam(b2, B[i], a2, f2flag, m2tid, m2pos, isize2, a1.rev,
                 !a1.mapped);
        if (sink) {
          sink(sink_ctx, b1);
          sink(sink_ctx, b2);
          // sink takes ownership; do not destroy here
        } else {
          if (sam_write1(out, hdr, b1) < 0 || sam_write1(out, hdr, b2) < 0)
            rc = 1;
          bam_destroy1(b1);
          bam_destroy1(b2);
        }
      } else {
        bam1_t *b = bam_init1();
        fill_bam(b, A[i], ra[i], 0, -1, -1, 0, false, false);
        if (sink) {
          sink(sink_ctx, b);
        } else {
          if (sam_write1(out, hdr, b) < 0) rc = 1;
          bam_destroy1(b);
        }
      }
    }
    if (rc) break;
  }

  kseq_destroy(k1);
  if (k2) kseq_destroy(k2);
  gzclose(f1);
  if (f2) gzclose(f2);
  fprintf(stderr, "[rabbitbin bwa] %lld reads, %lld mapped (%.2f%%)\n",
          (long long)total, (long long)mapped,
          total ? 100.0 * mapped / total : 0.0);
  return rc;
}

int rb_cmd_bwa(int ac, char *av[]) {
  std::string ref, r1, r2, interleaved_path, out_path;
  int threads = 0, k = 15, w = 10, band = 31;

  enum { OPT_BAND = 1000 };
  static const struct option longopts[] = {
      {"ref", required_argument, 0, 'r'},
      {"interleaved", required_argument, 0, 'p'},
      {"out", required_argument, 0, 'o'},
      {"threads", required_argument, 0, 't'},
      {"band", required_argument, 0, OPT_BAND},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};
  optind = 1;
  int c;
  while ((c = getopt_long(ac, av, "r:1:2:p:o:t:k:w:h", longopts, nullptr)) !=
         -1) {
    switch (c) {
      case 'r': ref = optarg; break;
      case '1': r1 = optarg; break;
      case '2': r2 = optarg; break;
      case 'p': interleaved_path = optarg; break;
      case 'o': out_path = optarg; break;
      case 't': threads = atoi(optarg); break;
      case 'k': k = atoi(optarg); break;
      case 'w': w = atoi(optarg); break;
      case OPT_BAND: band = atoi(optarg); break;
      case 'h': align_usage(); return 0;
      default: align_usage(); return 1;
    }
  }
  // positional single-end reads
  if (r1.empty() && interleaved_path.empty() && optind < ac) r1 = av[optind];
  bool interleaved = !interleaved_path.empty();
  if (interleaved) r1 = interleaved_path;

  if (ref.empty() || r1.empty()) {
    align_usage();
    fprintf(stderr, "[Error!] --ref and reads are required\n");
    return 1;
  }
  long onln = sysconf(_SC_NPROCESSORS_ONLN);
  if (threads <= 0) threads = (onln > 0) ? (int)onln : 1;

  SaIndex idx;
  if (!idx.build(ref, k, w, threads)) return 1;

  SaOpt opt;
  opt.band = band;

  // Output: .bam -> BAM, else SAM (text); "-"/empty -> stdout SAM.
  bool bam_out = out_path.size() >= 4 &&
                 out_path.compare(out_path.size() - 4, 4, ".bam") == 0;
  const char *dst = out_path.empty() ? "-" : out_path.c_str();
  samFile *out = sam_open(dst, bam_out ? "wb" : "w");
  if (!out) {
    fprintf(stderr, "[Error!] cannot open output: %s\n", dst);
    return 1;
  }
  if (threads > 1) hts_set_threads(out, threads);
  sam_hdr_t *hdr = rb_align_make_header(idx);
  if (sam_hdr_write(out, hdr) < 0) {
    fprintf(stderr, "[Error!] failed writing header\n");
    sam_hdr_destroy(hdr);
    sam_close(out);
    return 1;
  }

  int rc = rb_align_run(idx, opt, r1, r2, interleaved, threads, out, hdr,
                        nullptr, nullptr);

  sam_hdr_destroy(hdr);
  if (sam_close(out) < 0) rc = 1;
  return rc;
}
