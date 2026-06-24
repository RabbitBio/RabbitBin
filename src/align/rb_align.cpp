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
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <map>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

#include <omp.h>
#include <zlib.h>
#if defined(USE_IGZIP)
#include <igzip_lib.h>
#endif

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/kseq.h>

#include "rb_map_cli.h"
#include "sa_index.h"
#include "sa_map.h"
#include "rb_align_api.h"

namespace {

// FASTQ byte source for kseq: PLAIN (fread), ZLIB (gzread), or ISA-L igzip when
// built with USE_IGZIP. igzip decompresses a single gzip stream ~2-4x faster
// than zlib (a DEFLATE stream's inflate is inherently sequential, so it runs on
// the pipeline's single reader thread). The decoded bytes are identical to
// zlib, so parsing -- and therefore alignment output -- is unchanged. The .gz
// vs plain choice is made from the file's magic bytes. Mirrors the FASTA reader
// facade in rabbitbin.cpp.
struct FqReader {
  enum Mode { PLAIN, ZLIB
#if defined(USE_IGZIP)
            , IGZIP
#endif
  };
  Mode mode = PLAIN;
  FILE *fp = nullptr;
  gzFile gz = nullptr;
  bool eof_ = false;
#if defined(USE_IGZIP)
  static const size_t IN_SZ = 1u << 22;     // 4 MB compressed read granularity
  static const uint32_t HDR_REQ = 1u << 16; // bytes to (re)read a gz header
  unsigned char *inbuf = nullptr;
  struct inflate_state st;
  struct isal_gzip_header hdr;
#endif

  bool open(const char *path) {
    FILE *probe = fopen(path, "rb");
    if (!probe) return false;
    unsigned char m[2] = {0, 0};
    size_t n = fread(m, 1, 2, probe);
    fclose(probe);
    bool is_gz = (n == 2 && m[0] == 0x1f && m[1] == 0x8b);
    if (!is_gz) {
      mode = PLAIN;
      fp = fopen(path, "rb");
      return fp != nullptr;
    }
#if defined(USE_IGZIP)
    mode = IGZIP;
    fp = fopen(path, "rb");
    if (!fp) return false;
    inbuf = (unsigned char *)malloc(IN_SZ);
    if (!inbuf) { fclose(fp); fp = nullptr; return false; }
    isal_gzip_header_init(&hdr);
    isal_inflate_init(&st);
    st.crc_flag = ISAL_GZIP_NO_HDR_VER;
    st.next_in = inbuf;
    st.avail_in = (uint32_t)fread(inbuf, 1, IN_SZ, fp);
    if (isal_read_gzip_header(&st, &hdr) != ISAL_DECOMP_OK) return false;
    return true;
#else
    mode = ZLIB;
    gz = gzopen(path, "r");
    if (!gz) return false;
    gzbuffer(gz, 4u << 20);
    return true;
#endif
  }

  int read(void *outv, int size) {
    char *out = (char *)outv;
    if (mode == PLAIN) return (int)fread(out, 1, (size_t)size, fp);
#if defined(USE_IGZIP)
    if (mode == IGZIP) {
      if (eof_) return 0;
      unsigned char *o = (unsigned char *)out;
      uint64_t off = 0;
      while (off < (uint64_t)size) {
        if (st.avail_in == 0) {
          if (feof(fp)) { eof_ = true; break; }
          st.next_in = inbuf;
          st.avail_in = (uint32_t)fread(inbuf, 1, IN_SZ, fp);
          if (st.avail_in == 0) { eof_ = true; break; }
        }
        st.next_out = o + off;
        st.avail_out = (uint32_t)((uint64_t)size - off);
        int r = isal_inflate(&st);
        off = (uint64_t)(st.next_out - o);
        if (r != ISAL_DECOMP_OK) { eof_ = true; break; }
        if (st.block_state == ISAL_BLOCK_FINISH) {
          if (feof(fp) && st.avail_in == 0) { eof_ = true; break; }
          if (st.avail_in == 0) {
            isal_inflate_reset(&st);
            st.next_in = inbuf;
            st.avail_in = (uint32_t)fread(inbuf, 1, IN_SZ, fp);
            if (st.avail_in == 0) { eof_ = true; break; }
          } else if (st.avail_in >= HDR_REQ) {
            unsigned char *oni = st.next_in;
            uint32_t oai = st.avail_in;
            isal_inflate_reset(&st);
            st.next_in = oni;
            st.avail_in = oai;
          } else {
            uint32_t oai = st.avail_in;
            memmove(inbuf, st.next_in, oai);
            size_t rd = feof(fp) ? 0 : fread(inbuf + oai, 1, IN_SZ - oai, fp);
            isal_inflate_reset(&st);
            st.next_in = inbuf;
            st.avail_in = oai + (uint32_t)rd;
          }
          if (isal_read_gzip_header(&st, &hdr) != ISAL_DECOMP_OK) {
            eof_ = true;
            break;
          }
        }
      }
      return (int)off;
    }
#endif
    return gzread(gz, out, (unsigned)size);  // ZLIB
  }

  void close() {
    if (gz) { gzclose(gz); gz = nullptr; }
    if (fp) { fclose(fp); fp = nullptr; }
#if defined(USE_IGZIP)
    if (inbuf) { free(inbuf); inbuf = nullptr; }
#endif
  }
};

static int fq_read(FqReader *r, void *buf, int len) { return r->read(buf, len); }

KSEQ_INIT(FqReader *, fq_read)

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
      "  -k INT             Seed k-mer (minimizer k / strobe k; default 19)\n"
      "  -w INT             Minimizer window (default 19; minimizer mode only)\n"
      "      --seed MODE    minimizer (default) | randstrobe\n"
      "      --rs-s INT     randstrobe syncmer s-mer (default 16)\n"
      "      --rs-wmin INT  randstrobe window min (default 5)\n"
      "      --rs-wmax INT  randstrobe window max (default 11)\n"
      "      --rs-maxdist INT  randstrobe max strobe gap bp (default 80)\n"
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
  FqReader f1;
  if (!f1.open(r1.c_str())) {
    fprintf(stderr, "[Error!] cannot open reads: %s\n", r1.c_str());
    return 1;
  }
  FqReader f2;
  bool have_f2 = false;
  if (!r2.empty()) {
    if (!f2.open(r2.c_str())) {
      fprintf(stderr, "[Error!] cannot open reads: %s\n", r2.c_str());
      f1.close();
      return 1;
    }
    have_f2 = true;
  }
  kseq_t *k1 = kseq_init(&f1);
  kseq_t *k2 = have_f2 ? kseq_init(&f2) : nullptr;

  const int BATCH = 8192;  // reads/pairs per batch (fine grain => load balance)

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

  // Fill one batch from the kseq stream(s); returns the number of reads/pairs.
  auto load_batch = [&](std::vector<ReadRec> &A,
                        std::vector<ReadRec> &B) -> int {
    A.clear();
    B.clear();
    int got = 0;
    while (got < BATCH) {
      ReadRec a;
      if (!load_one(k1, a)) break;
      if (paired) {
        ReadRec b;
        bool ok = interleaved ? load_one(k1, b) : load_one(k2, b);
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
    return got;
  };

  // ── Producer / consumer pipeline ──────────────────────────────────────────
  // 1 reader thread parses FASTQ batches; `n_workers` threads align + build the
  // BAM records for a whole batch each; 1 writer thread emits completed batches
  // IN INPUT ORDER. Reading, alignment, and serialization therefore overlap and
  // the workers never wait on a per-batch barrier (better core utilisation),
  // while ordered emission keeps the output byte-identical to a serial run.
  struct InBatch {
    uint64_t seq = 0;
    std::vector<ReadRec> A, B;
    int got = 0;
  };
  struct OutBatch {
    uint64_t seq = 0;
    int got = 0;
    std::vector<bam1_t *> b1, b2;
    int64_t n_total = 0, n_mapped = 0;
  };

  const int n_workers = threads > 0 ? threads : 1;
  // Cap batches in flight (reader-created .. writer-freed) to bound memory.
  const int max_inflight = std::max(4, n_workers * 3);

  std::mutex in_mtx;
  std::condition_variable in_not_empty, in_not_full;
  std::queue<InBatch *> in_q;
  int inflight = 0;
  bool reading_done = false;   // guarded by in_mtx (workers)
  bool all_produced = false;   // guarded by out_mtx (writer)
  uint64_t n_batches = 0;      // guarded by out_mtx (writer)

  std::mutex out_mtx;
  std::condition_variable out_cv;
  std::map<uint64_t, OutBatch *> out_ready;

  std::atomic<int> rc{0};
  std::atomic<bool> aborted{false};
  int64_t total = 0, mapped = 0;  // only the writer thread touches these

  auto do_abort = [&]() {
    aborted.store(true);
    { std::lock_guard<std::mutex> lk(in_mtx); }
    in_not_empty.notify_all();
    in_not_full.notify_all();
    { std::lock_guard<std::mutex> lk(out_mtx); }
    out_cv.notify_all();
  };

  // Align one input batch into records (single-threaded; runs on a worker).
  auto process_batch = [&](InBatch *ib) -> OutBatch * {
    OutBatch *ob = new OutBatch();
    ob->seq = ib->seq;
    ob->got = ib->got;
    ob->b1.assign(ib->got, nullptr);
    if (paired) ob->b2.assign(ib->got, nullptr);
    int64_t nt = 0, nm = 0;

    // Align the WHOLE batch in one batched, latency-hidden pass before building
    // records. Reads are laid out [A0,B0,A1,B1,...] (paired) or [A0,A1,...] so
    // the aligner's sub-batch pipeline keeps many reference-window prefetches in
    // flight. Per-read results are identical to aligning one at a time.
    static thread_local std::vector<const char *> seqs;
    static thread_local std::vector<int> lens;
    static thread_local std::vector<SaAln> alns;
    const int m = paired ? 2 * ib->got : ib->got;
    if ((int)seqs.size() < m) { seqs.resize(m); lens.resize(m); alns.resize(m); }
    for (int i = 0; i < ib->got; ++i) {
      int ia = paired ? 2 * i : i;
      seqs[ia] = ib->A[i].seq.c_str();
      lens[ia] = (int)ib->A[i].seq.size();
      if (paired) {
        seqs[2 * i + 1] = ib->B[i].seq.c_str();
        lens[2 * i + 1] = (int)ib->B[i].seq.size();
      }
    }
    sa_align_reads_batch(idx, opt, seqs.data(), lens.data(), m, alns.data());

    for (int i = 0; i < ib->got; ++i) {
      SaAln a1 = std::move(alns[paired ? 2 * i : i]);
      ++nt;
      if (a1.mapped) ++nm;
      if (paired) {
        SaAln a2 = std::move(alns[2 * i + 1]);
        ++nt;
        if (a2.mapped) ++nm;
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
          int64_t span =
              hi - lo + (int64_t)std::max(ib->A[i].seq.size(), ib->B[i].seq.size());
          isize1 = (a1.pos <= a2.pos) ? span : -span;
          isize2 = -isize1;
        }
        uint16_t f1flag =
            BAM_FPAIRED | BAM_FREAD1 | (proper ? BAM_FPROPER_PAIR : 0);
        uint16_t f2flag =
            BAM_FPAIRED | BAM_FREAD2 | (proper ? BAM_FPROPER_PAIR : 0);
        bam1_t *b1 = bam_init1();
        bam1_t *b2 = bam_init1();
        fill_bam(b1, ib->A[i], a1, f1flag, m1tid, m1pos, isize1, a2.rev,
                 !a2.mapped);
        fill_bam(b2, ib->B[i], a2, f2flag, m2tid, m2pos, isize2, a1.rev,
                 !a1.mapped);
        ob->b1[i] = b1;
        ob->b2[i] = b2;
      } else {
        bam1_t *b = bam_init1();
        fill_bam(b, ib->A[i], a1, 0, -1, -1, 0, false, false);
        ob->b1[i] = b;
      }
    }
    ob->n_total = nt;
    ob->n_mapped = nm;
    return ob;
  };

  std::thread reader([&]() {
    uint64_t seq = 0;
    for (;;) {
      if (aborted.load()) break;
      InBatch *ib = new InBatch();
      ib->seq = seq;
      ib->got = load_batch(ib->A, ib->B);
      if (ib->got == 0) {
        delete ib;
        break;
      }
      {
        std::unique_lock<std::mutex> lk(in_mtx);
        in_not_full.wait(
            lk, [&] { return inflight < max_inflight || aborted.load(); });
        if (aborted.load()) {
          delete ib;
          break;
        }
        ++inflight;
        in_q.push(ib);
      }
      in_not_empty.notify_one();
      ++seq;
    }
    { std::lock_guard<std::mutex> lk(in_mtx); reading_done = true; }
    in_not_empty.notify_all();
    { std::lock_guard<std::mutex> lk(out_mtx); all_produced = true; n_batches = seq; }
    out_cv.notify_all();
  });

  auto worker_fn = [&]() {
    for (;;) {
      InBatch *ib = nullptr;
      {
        std::unique_lock<std::mutex> lk(in_mtx);
        in_not_empty.wait(lk, [&] {
          return !in_q.empty() || reading_done || aborted.load();
        });
        if (in_q.empty()) {
          if (reading_done || aborted.load()) break;
          continue;
        }
        ib = in_q.front();
        in_q.pop();
      }
      if (aborted.load()) { delete ib; break; }
      OutBatch *ob = process_batch(ib);
      delete ib;  // input reads no longer needed
      {
        std::lock_guard<std::mutex> lk(out_mtx);
        out_ready[ob->seq] = ob;
      }
      out_cv.notify_one();
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(n_workers);
  for (int t = 0; t < n_workers; ++t) workers.emplace_back(worker_fn);

  std::thread writer([&]() {
    uint64_t expected = 0;
    for (;;) {
      OutBatch *ob = nullptr;
      {
        std::unique_lock<std::mutex> lk(out_mtx);
        out_cv.wait(lk, [&] {
          return out_ready.count(expected) || aborted.load() ||
                 (all_produced && expected >= n_batches);
        });
        auto it = out_ready.find(expected);
        if (it == out_ready.end()) {
          if (aborted.load() || (all_produced && expected >= n_batches)) break;
          continue;
        }
        ob = it->second;
        out_ready.erase(it);
      }
      for (int i = 0; i < ob->got; ++i) {
        if (sink) {
          sink(sink_ctx, ob->b1[i]);
          if (paired) sink(sink_ctx, ob->b2[i]);
        } else {
          if (sam_write1(out, hdr, ob->b1[i]) < 0) rc.store(1);
          bam_destroy1(ob->b1[i]);
          if (paired) {
            if (sam_write1(out, hdr, ob->b2[i]) < 0) rc.store(1);
            bam_destroy1(ob->b2[i]);
          }
        }
      }
      total += ob->n_total;
      mapped += ob->n_mapped;
      delete ob;
      ++expected;
      { std::lock_guard<std::mutex> lk(in_mtx); --inflight; }
      in_not_full.notify_one();
      if (rc.load()) { do_abort(); break; }
    }
  });

  reader.join();
  for (auto &t : workers) t.join();
  writer.join();

  // Free any batches left over after an abort (normal completion leaves none).
  while (!in_q.empty()) { delete in_q.front(); in_q.pop(); }
  for (auto &kv : out_ready) {
    if (!sink) {
      for (auto *b : kv.second->b1) if (b) bam_destroy1(b);
      for (auto *b : kv.second->b2) if (b) bam_destroy1(b);
    }
    delete kv.second;
  }

  kseq_destroy(k1);
  if (k2) kseq_destroy(k2);
  f1.close();
  if (have_f2) f2.close();
  fprintf(stderr, "[rabbitbin bwa] %lld reads, %lld mapped (%.2f%%)\n",
          (long long)total, (long long)mapped,
          total ? 100.0 * mapped / total : 0.0);
  return rc.load();
}

int rb_cmd_bwa(int ac, char *av[]) {
  std::string ref, r1, r2, interleaved_path, out_path;
  int threads = 0, k = 19, w = 19, band = 31;
  bool use_rs = false;  // minimizer seeding by default (faster: smaller index)
  SaRsParams rs;       // defaults: s=16 wmin=5 wmax=11 maxdist=80 q=255

  enum {
    OPT_BAND = 1000,
    OPT_SEED = 1001,
    OPT_RS_S = 1002,
    OPT_RS_WMIN = 1003,
    OPT_RS_WMAX = 1004,
    OPT_RS_MAXDIST = 1005
  };
  static const struct option longopts[] = {
      {"ref", required_argument, 0, 'r'},
      {"interleaved", required_argument, 0, 'p'},
      {"out", required_argument, 0, 'o'},
      {"threads", required_argument, 0, 't'},
      {"band", required_argument, 0, OPT_BAND},
      {"seed", required_argument, 0, OPT_SEED},
      {"rs-s", required_argument, 0, OPT_RS_S},
      {"rs-wmin", required_argument, 0, OPT_RS_WMIN},
      {"rs-wmax", required_argument, 0, OPT_RS_WMAX},
      {"rs-maxdist", required_argument, 0, OPT_RS_MAXDIST},
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
      case OPT_SEED:
        use_rs = (std::string(optarg) != "minimizer");
        break;
      case OPT_RS_S: rs.s = atoi(optarg); break;
      case OPT_RS_WMIN: rs.w_min = atoi(optarg); break;
      case OPT_RS_WMAX: rs.w_max = atoi(optarg); break;
      case OPT_RS_MAXDIST: rs.max_dist = atoi(optarg); break;
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
  idx.use_randstrobe = use_rs;
  idx.rs = rs;
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
