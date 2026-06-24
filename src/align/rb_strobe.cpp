// rb_strobe.cpp - RabbitBin <-> vendored strobealign engine bridge.
//
// Builds the strobemer index once (in memory) and aligns reads with strobealign's
// align_or_map_* core, converting its SAM output to bam1_t for RabbitBin's
// parallel coordinate sort + BGZF writer. See rb_strobe.h.

#include "rb_strobe.h"

#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <map>
#include <mutex>
#include <queue>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>

// strobealign engine (vendored under src/align/strobe).
#include "refs.hpp"
#include "indexparameters.hpp"
#include "index.hpp"
#include "aligner.hpp"
#include "chain.hpp"
#include "mappingparameters.hpp"
#include "statistics.hpp"
#include "insertsizedistribution.hpp"
#include "sam.hpp"
#include "aln.hpp"
#include "kseq++/kseq++.hpp"

KSEQ_INIT(gzFile, gzread)

struct StrobeMapper {
  References references;            // must outlive index (held by const ref)
  IndexParameters index_parameters;
  StrobemerIndex index;
  Aligner aligner;
  Chainer chainer;
  MappingParameters map_param;
  int read_len = 150;

  StrobeMapper(References &&refs, const IndexParameters &ip,
               const AlignmentParameters &ap, const MappingParameters &mp)
      : references(std::move(refs)),
        index_parameters(ip),
        index(references, index_parameters),
        aligner(ap),
        chainer(mp.chaining_params, index.k()),
        map_param(mp) {}
};

StrobeMapper *rb_strobe_build(const std::string &ref_fasta, int read_len,
                              int threads) {
  try {
    References refs = References::from_fasta(ref_fasta);
    fprintf(stderr,
            "[rabbitbin strobe] reference: %zu contigs, %zu bp (r=%d)\n",
            refs.size(), refs.total_length(), read_len);
    IndexParameters ip = IndexParameters::from_read_length(read_len);

    AlignmentParameters ap;
    ap.match = 2;
    ap.mismatch = 8;
    ap.gap_open = 12;
    ap.gap_extend = 1;
    ap.end_bonus = 10;

    MappingParameters mp;
    mp.r = read_len;
    mp.chaining_params = {50, 0.1f, 0.05f, 0.7f, 10000, 0.01f};
    mp.output_format = OutputFormat::SAM;
    mp.cigar_ops = CigarOps::M;  // M ops (coverage/binning friendly, like bwa)
    mp.output_unmapped = true;
    mp.use_nams = false;  // collinear chaining (default, higher accuracy)
    mp.verify();

    auto *m = new StrobeMapper(std::move(refs), ip, ap, mp);
    m->read_len = read_len;
    fprintf(stderr, "[rabbitbin strobe] building strobemer index (%d threads)\n",
            threads);
    m->index.populate(0.0002f, threads < 1 ? 1 : threads);
    fprintf(stderr, "[rabbitbin strobe] index ready\n");
    return m;
  } catch (const std::exception &e) {
    fprintf(stderr, "[Error!] strobe index build failed: %s\n", e.what());
    return nullptr;
  }
}

void rb_strobe_free(StrobeMapper *m) { delete m; }

sam_hdr_t *rb_strobe_make_header(const StrobeMapper *m) {
  sam_hdr_t *h = sam_hdr_init();
  std::string hd = "@HD\tVN:1.6\tSO:unsorted\n";
  sam_hdr_add_lines(h, hd.c_str(), hd.size());
  std::string sq;
  sq.reserve(m->references.size() * 32);
  for (size_t i = 0; i < m->references.size(); ++i) {
    sq += "@SQ\tSN:";
    sq += m->references.names[i];
    sq += "\tLN:";
    sq += std::to_string(m->references.lengths[i]);
    sq += "\n";
  }
  sam_hdr_add_lines(h, sq.c_str(), sq.size());
  std::string pg = "@PG\tID:rabbitbin\tPN:rabbitbin\tVN:strobe\n";
  sam_hdr_add_lines(h, pg.c_str(), pg.size());
  return h;
}

namespace {

struct ReadRec {
  std::string name, seq, qual;
};

struct InBatch {
  uint64_t seq = 0;
  std::vector<ReadRec> A, B;
  int got = 0;
};

struct OutBatch {
  uint64_t seq = 0;
  std::vector<bam1_t *> recs;  // finished records (in input order within batch)
  int64_t n_total = 0, n_mapped = 0;
};

klibpp::KSeq to_kseq(const ReadRec &r) {
  klibpp::KSeq k;
  k.name = r.name;
  k.seq = r.seq;
  k.qual = r.qual;
  to_uppercase(k.seq);
  return k;
}

}  // namespace

int rb_strobe_run(StrobeMapper *m, const std::string &r1, const std::string &r2,
                  bool interleaved, int threads, samFile *out, sam_hdr_t *hdr,
                  RbAlignSink sink, void *sink_ctx) {
  const bool paired = interleaved || !r2.empty();
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

  // Warm the header's RNAME->tid hash once (single-threaded) so concurrent
  // sam_parse1 lookups in the workers are read-only.
  if (m->references.size() > 0)
    (void)sam_hdr_name2tid(hdr, m->references.names[0].c_str());

  const int BATCH = 4096;
  auto load_one = [](kseq_t *ks, ReadRec &out) -> bool {
    int r = kseq_read(ks);
    if (r < 0) return false;
    out.name.assign(ks->name.s, ks->name.l);
    out.seq.assign(ks->seq.s, ks->seq.l);
    out.qual.assign(ks->qual.l ? ks->qual.s : "*", ks->qual.l ? ks->qual.l : 1);
    return true;
  };
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
        if (!ok) break;
        A.push_back(std::move(a));
        B.push_back(std::move(b));
      } else {
        A.push_back(std::move(a));
      }
      ++got;
    }
    return got;
  };

  const int n_workers = threads > 0 ? threads : 1;
  const int max_inflight = std::max(4, n_workers * 3);

  std::mutex in_mtx;
  std::condition_variable in_not_empty, in_not_full;
  std::queue<InBatch *> in_q;
  int inflight = 0;
  bool reading_done = false;

  std::mutex out_mtx;
  std::condition_variable out_cv;
  std::map<uint64_t, OutBatch *> out_ready;
  bool all_produced = false;
  uint64_t n_batches = 0;

  std::atomic<int> rc{0};
  std::atomic<bool> aborted{false};
  int64_t total = 0, mapped = 0;

  auto do_abort = [&]() {
    aborted.store(true);
    { std::lock_guard<std::mutex> lk(in_mtx); }
    in_not_empty.notify_all();
    in_not_full.notify_all();
    { std::lock_guard<std::mutex> lk(out_mtx); }
    out_cv.notify_all();
  };

  // Align one batch -> records (runs on a worker; per-worker strobealign state).
  auto process_batch = [&](InBatch *ib, int wid) -> OutBatch * {
    auto *ob = new OutBatch();
    ob->seq = ib->seq;
    // Per-worker strobealign scratch (reused across this worker's batches).
    static thread_local Sam *sam = nullptr;
    static thread_local std::string sam_string;
    static thread_local AlignmentStatistics stats;
    static thread_local InsertSizeDistribution isize;
    static thread_local std::minstd_rand rng;
    static thread_local std::vector<double> abundances;
    static thread_local bool inited = false;
    if (!inited) {
      rng.seed(1234567u + (unsigned)wid);
      abundances.assign(m->references.size(), 0.0);
      inited = true;
    }
    if (!sam) {
      sam = new Sam(sam_string, m->references, m->map_param.cigar_ops, "",
                    m->map_param.output_unmapped, false, false);
    }

    // Parse SAM text -> bam1_t with a reusable kstring.
    static thread_local kstring_t ks = {0, 0, nullptr};
    auto emit_samtext = [&](OutBatch *o) {
      // sam_string holds 1+ SAM records separated by '\n'; parse each.
      size_t start = 0;
      const std::string &s = sam_string;
      while (start < s.size()) {
        size_t nl = s.find('\n', start);
        size_t end = (nl == std::string::npos) ? s.size() : nl;
        if (end > start) {
          ks.l = 0;
          kputsn(s.data() + start, (int)(end - start), &ks);
          bam1_t *b = bam_init1();
          if (sam_parse1(&ks, hdr, b) >= 0) {
            o->recs.push_back(b);
          } else {
            bam_destroy1(b);
          }
        }
        if (nl == std::string::npos) break;
        start = nl + 1;
      }
    };

    int64_t nt = 0, nm = 0;
    for (int i = 0; i < ib->got; ++i) {
      sam_string.clear();
      klibpp::KSeq ra = to_kseq(ib->A[i]);
      if (paired) {
        klibpp::KSeq rb = to_kseq(ib->B[i]);
        align_or_map_paired(ra, rb, *sam, sam_string, stats, isize, m->aligner,
                            m->chainer, m->map_param, m->index_parameters,
                            m->references, m->index, rng, abundances);
        nt += 2;
      } else {
        align_or_map_single(ra, *sam, sam_string, stats, m->aligner, m->chainer,
                            m->map_param, m->index_parameters, m->references,
                            m->index, rng, abundances);
        nt += 1;
      }
      emit_samtext(ob);
    }
    // count mapped from flags
    for (auto *b : ob->recs)
      if (!(b->core.flag & BAM_FUNMAP)) ++nm;
    ob->n_total = nt;
    ob->n_mapped = nm;
    return ob;
  };

  std::thread reader([&]() {
    uint64_t seq = 0;
    for (;;) {
      if (aborted.load()) break;
      auto *ib = new InBatch();
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

  auto worker_fn = [&](int wid) {
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
      OutBatch *ob = process_batch(ib, wid);
      delete ib;
      {
        std::lock_guard<std::mutex> lk(out_mtx);
        out_ready[ob->seq] = ob;
      }
      out_cv.notify_one();
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(n_workers);
  for (int t = 0; t < n_workers; ++t) workers.emplace_back(worker_fn, t);

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
      for (auto *b : ob->recs) {
        if (sink) {
          sink(sink_ctx, b);
        } else {
          if (sam_write1(out, hdr, b) < 0) rc.store(1);
          bam_destroy1(b);
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

  while (!in_q.empty()) { delete in_q.front(); in_q.pop(); }
  for (auto &kv : out_ready) {
    if (!sink)
      for (auto *b : kv.second->recs) bam_destroy1(b);
    delete kv.second;
  }

  kseq_destroy(k1);
  if (k2) kseq_destroy(k2);
  gzclose(f1);
  if (f2) gzclose(f2);
  fprintf(stderr, "[rabbitbin strobe] %lld reads, %lld mapped (%.2f%%)\n",
          (long long)total, (long long)mapped,
          total ? 100.0 * mapped / total : 0.0);
  return rc.load();
}
