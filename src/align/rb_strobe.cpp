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

#include <cstring>
#include <string_view>
#include <omp.h>
#include <zlib.h>
#include <htslib/sam.h>

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

// Parallel FASTA loader: build strobealign's References from a (.gz ok) FASTA
// using all cores. The single-threaded References::from_fasta does line-by-line
// getline over millions of contigs (~20s on the 6.1M-contig CAMI3 reference);
// here the file is read once, contig-start offsets located, and each contig
// parsed in parallel. Contig ORDER is preserved (ref_id == BAM tid == @SQ order),
// and sequences are uppercased exactly as from_fasta does.
static References build_references_parallel(const std::string &fasta,
                                            int threads) {
  // Read the whole file (gzread is transparent for uncompressed input too).
  std::string raw;
  {
    gzFile fp = gzopen(fasta.c_str(), "rb");
    if (!fp) throw std::runtime_error("cannot open reference FASTA: " + fasta);
    const size_t CH = 32u << 20;
    size_t used = 0;
    for (;;) {
      raw.resize(used + CH);
      int n = gzread(fp, &raw[used], (unsigned)CH);
      if (n <= 0) break;
      used += (size_t)n;
    }
    gzclose(fp);
    raw.resize(used);
  }
  if (raw.empty() || raw[0] != '>')
    throw std::runtime_error("FASTA must begin with '>'");

  // Locate contig starts: offset 0 and every '\n' immediately followed by '>'.
  std::vector<size_t> starts;
  starts.push_back(0);
  {
    const char *base = raw.data();
    const char *p = base;
    const char *end = base + raw.size();
    while (p < end) {
      const char *nl = (const char *)memchr(p, '\n', end - p);
      if (!nl) break;
      size_t j = (size_t)(nl - base) + 1;
      if (j < raw.size() && base[j] == '>') starts.push_back(j);
      p = nl + 1;
    }
  }
  const size_t n = starts.size();
  std::vector<std::string> names(n), seqs(n);

#pragma omp parallel for num_threads(threads < 1 ? 1 : threads) schedule(dynamic, 256)
  for (long long ci = 0; ci < (long long)n; ++ci) {
    size_t s = starts[ci];
    size_t e = (ci + 1 < (long long)n) ? starts[ci + 1] : raw.size();
    const char *p = raw.data() + s;
    const char *cend = raw.data() + e;
    // Header line (after '>'), name = up to first whitespace.
    const char *hnl = (const char *)memchr(p, '\n', cend - p);
    const char *hend = hnl ? hnl : cend;
    const char *np = p + 1;  // skip '>'
    const char *ne = np;
    while (ne < hend && *ne != ' ' && *ne != '\t' && *ne != '\r' && *ne != '\v' &&
           *ne != '\f')
      ++ne;
    names[ci].assign(np, ne - np);
    // Sequence: concatenate the remaining lines, uppercased, newlines stripped.
    std::string &seq = seqs[ci];
    seq.reserve((size_t)(cend - (hnl ? hnl + 1 : cend)));
    const char *q = hnl ? hnl + 1 : cend;
    while (q < cend) {
      const char *nl = (const char *)memchr(q, '\n', cend - q);
      const char *le = nl ? nl : cend;
      for (const char *t = q; t < le; ++t) {
        char c = *t;
        if (c == '\r') continue;
        seq.push_back((char)(c & ~32));  // uppercase (matches to_uppercase)
      }
      if (!nl) break;
      q = nl + 1;
    }
  }
  return References(std::move(seqs), std::move(names));
}

StrobeMapper *rb_strobe_build(const std::string &ref_fasta, int read_len,
                              int threads) {
  try {
    References refs = build_references_parallel(ref_fasta, threads);
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

// A chunk of raw FASTQ text ending exactly on a record(-pair) boundary. The
// reader only finds boundaries (memchr over newlines) -- cheap and serial --
// while the expensive parsing (field extraction, string building) happens in the
// workers, in parallel. This removes the single-reader bottleneck that otherwise
// starves the fast strobealign workers.
struct InBatch {
  uint64_t seq = 0;
  std::string raw;   // FASTQ text (interleaved pairs, single reads, or mate1)
  std::string raw2;  // mate2 text (2-file paired only)
  bool paired = false;
  bool interleaved = false;
};

struct OutBatch {
  uint64_t seq = 0;
  std::vector<bam1_t *> recs;  // finished records (in input order within batch)
  int64_t n_total = 0, n_mapped = 0;
};

// Parse one FASTQ record (4 lines) starting at p (< end). Fills name/seq/qual
// (views into the buffer). Returns pointer just past the record, or nullptr on
// malformed input.
inline const char *parse_record(const char *p, const char *end,
                                std::string_view &name, std::string_view &seq,
                                std::string_view &qual) {
  auto getline = [&](const char *&q, std::string_view &out) -> bool {
    const char *nl = (const char *)memchr(q, '\n', end - q);
    if (!nl) return false;
    const char *e = nl;
    if (e > q && e[-1] == '\r') --e;
    out = std::string_view(q, e - q);
    q = nl + 1;
    return true;
  };
  std::string_view l0, l2;
  if (!getline(p, l0)) return nullptr;  // @name
  if (!getline(p, seq)) return nullptr;
  if (!getline(p, l2)) return nullptr;  // +
  if (!getline(p, qual)) return nullptr;
  // strip leading '@' and trailing fields from name (up to first space/tab)
  name = l0;
  if (!name.empty() && name[0] == '@') name.remove_prefix(1);
  size_t sp = name.find_first_of(" \t");
  if (sp != std::string_view::npos) name = name.substr(0, sp);
  return p;
}

inline klibpp::KSeq to_kseq(std::string_view name, std::string_view seq,
                            std::string_view qual) {
  klibpp::KSeq k;
  k.name.assign(name.data(), name.size());
  k.seq.assign(seq.data(), seq.size());
  k.qual.assign(qual.data(), qual.size());
  to_uppercase(k.seq);
  return k;
}

// Reads raw FASTQ text from a gzFile and hands out chunks that end exactly on a
// record boundary (4 lines) whose record index is a multiple of `align` (2 for
// interleaved pairs, 1 otherwise). Only newline scanning happens here (serial,
// cheap); parsing is left to the workers.
struct GzChunker {
  gzFile gz;
  std::string carry;
  bool eof = false;

  // Append complete aligned records (~target bytes worth) into `out`. Returns
  // record count, or 0 at end of input.
  int next(size_t target, int align, std::string &out) {
    char buf[1 << 20];
    while (carry.size() < target && !eof) {
      int n = gzread(gz, buf, sizeof(buf));
      if (n <= 0) { eof = true; break; }
      carry.append(buf, (size_t)n);
    }
    if (carry.empty()) return 0;
    const char *p = carry.data();
    const char *end = p + carry.size();
    const char *q = p;
    size_t lines = 0, recs = 0, best_off = 0, best_recs = 0;
    while (true) {
      const char *nl = (const char *)memchr(q, '\n', end - q);
      if (!nl) break;
      ++lines;
      q = nl + 1;
      if ((lines & 3) == 0) {  // every 4th line completes a record
        ++recs;
        if (recs % (size_t)align == 0) { best_off = q - p; best_recs = recs; }
      }
    }
    if (best_off == 0) {
      if (eof) return 0;            // no complete aligned record and no more data
      return next(target * 2, align, out);  // need a bigger window
    }
    out.assign(carry.data(), best_off);
    carry.erase(0, best_off);
    return (int)best_recs;
  }

  // Read exactly n complete records into `out` (for the 2-file paired mate).
  bool next_n(int n, std::string &out) {
    char buf[1 << 20];
    for (;;) {
      const char *p = carry.data();
      const char *end = p + carry.size();
      const char *q = p;
      size_t lines = 0, recs = 0, off = 0;
      while (recs < (size_t)n) {
        const char *nl = (const char *)memchr(q, '\n', end - q);
        if (!nl) break;
        ++lines;
        q = nl + 1;
        if ((lines & 3) == 0) { ++recs; off = q - p; }
      }
      if ((int)recs == n) {
        out.assign(carry.data(), off);
        carry.erase(0, off);
        return true;
      }
      if (eof) return false;
      int got = gzread(gz, buf, sizeof(buf));
      if (got <= 0) { eof = true; continue; }
      carry.append(buf, (size_t)got);
    }
  }
};

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
  GzChunker ch1{f1, {}, false};
  GzChunker ch2{f2, {}, false};
  const bool two_file = paired && !interleaved;  // r1 + r2 separate files
  const int chunk_align = interleaved ? 2 : 1;   // interleaved: pair boundaries
  const size_t CHUNK_BYTES = 2u << 20;           // ~2 MB raw per chunk

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
  // strobealign's Sam writes bam1_t DIRECTLY into ob->recs via its bam_sink
  // (tid == ref_id; no SAM text, no sam_parse1, no RNAME->tid lookup).
  auto process_batch = [&](InBatch *ib, int wid) -> OutBatch * {
    auto *ob = new OutBatch();
    ob->seq = ib->seq;
    static thread_local Sam *sam = nullptr;
    static thread_local std::string sam_string;  // unused in bam mode
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
    sam->set_bam_sink(&ob->recs);  // records go straight into this batch

    int64_t nt = 0, nm = 0;
    const char *p = ib->raw.data();
    const char *end = p + ib->raw.size();
    if (ib->paired && ib->interleaved) {
      std::string_view n1, s1, q1, n2, s2, q2;
      while (p < end) {
        p = parse_record(p, end, n1, s1, q1);
        if (!p) break;
        p = parse_record(p, end, n2, s2, q2);
        if (!p) break;
        klibpp::KSeq ra = to_kseq(n1, s1, q1);
        klibpp::KSeq rb = to_kseq(n2, s2, q2);
        align_or_map_paired(ra, rb, *sam, sam_string, stats, isize, m->aligner,
                            m->chainer, m->map_param, m->index_parameters,
                            m->references, m->index, rng, abundances);
        nt += 2;
      }
    } else if (ib->paired) {  // 2-file: raw=mate1, raw2=mate2
      const char *p2 = ib->raw2.data();
      const char *end2 = p2 + ib->raw2.size();
      std::string_view n1, s1, q1, n2, s2, q2;
      while (p < end && p2 < end2) {
        p = parse_record(p, end, n1, s1, q1);
        p2 = parse_record(p2, end2, n2, s2, q2);
        if (!p || !p2) break;
        klibpp::KSeq ra = to_kseq(n1, s1, q1);
        klibpp::KSeq rb = to_kseq(n2, s2, q2);
        align_or_map_paired(ra, rb, *sam, sam_string, stats, isize, m->aligner,
                            m->chainer, m->map_param, m->index_parameters,
                            m->references, m->index, rng, abundances);
        nt += 2;
      }
    } else {
      std::string_view n1, s1, q1;
      while (p < end) {
        p = parse_record(p, end, n1, s1, q1);
        if (!p) break;
        klibpp::KSeq ra = to_kseq(n1, s1, q1);
        align_or_map_single(ra, *sam, sam_string, stats, m->aligner, m->chainer,
                            m->map_param, m->index_parameters, m->references,
                            m->index, rng, abundances);
        nt += 1;
      }
    }
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
      ib->paired = paired;
      ib->interleaved = interleaved;
      int nrec = ch1.next(CHUNK_BYTES, chunk_align, ib->raw);
      if (nrec == 0) {
        delete ib;
        break;
      }
      if (two_file) {
        if (!ch2.next_n(nrec, ib->raw2)) {  // matching mate2 records
          delete ib;
          break;
        }
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

  gzclose(f1);
  if (f2) gzclose(f2);
  fprintf(stderr, "[rabbitbin strobe] %lld reads, %lld mapped (%.2f%%)\n",
          (long long)total, (long long)mapped,
          total ? 100.0 * mapped / total : 0.0);
  return rc.load();
}
