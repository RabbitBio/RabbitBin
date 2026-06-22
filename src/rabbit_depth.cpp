// rabbit_depth.cpp — summarize per-contig coverage depth from sorted BAM files.

#include <cassert>
#include <condition_variable>
#include <mutex>

#include "rabbit_depth.h"
#include "version.h"

#include "IOThreadBuffer.h"
#include "KseqReader.h"
#include "RunningStats.h"
#include "SafeOfstream.hpp"

#ifdef USE_RABBITBAM
// RabbitBAM parallel BAM reader: decompresses BGZF blocks of a single bam with
// a pool of worker threads, then hands back bam1_t records in file order via
// getBam1_t() — a drop-in for sam_read1() that removes the serial-decompress
// bottleneck of the depth stage. Record order is identical to sam_read1, so the
// emitted depth table is unchanged.
#include "BamReader.h"
#endif

ThreadBlocker tb;

static struct option long_options[] = {{"help", 0, 0, 0},
                                       {"outputDepth", 1, 0, 0},
                                       {"percentIdentity", 1, 0, 0},
                                       {"pairedContigs", 1, 0, 0},
                                       {"unmappedFastq", 1, 0, 0},
                                       {"referenceFasta", 1, 0, 0},
                                       {"shredLength", 1, 0, 0},
                                       {"shredDepth", 1, 0, 0},
                                       {"minContigLength", 1, 0, 0},
                                       {"minContigDepth", 1, 0, 0},
                                       {"minMapQual", 1, 0, 0},
                                       {"weightMapQual", 1, 0, 0},
                                       {"noIntraDepthVariance", 0, 0, 0},
                                       {"showDepth", 0, 0, 0},
                                       {"includeEdgeBases", 0, 0, 0},
                                       {"maxEdgeBases", 1, 0, 0},
                                       {"outputReadStats", 1, 0, 0},
                                       {"outputGC", 1, 0, 0},
                                       {"gcWindow", 1, 0, 0},
                                       {"outputKmers", 1, 0, 0},
                                       {"checkpoint", 0, 0, 0}};

void usage() {
  fprintf(
      stderr,
      "rabbit_depth %s %s\n"
      "Usage: rabbit_depth <options> sortedBam1 [ "
      "sortedBam2 ...]\n"
      "where options include:\n"
      "\t--outputDepth       arg  The file to put the contig by bam depth "
      "matrix (default: STDOUT)\n"
      "\t--checkpoint             The write checkpoints for every bam "
      "processed (default: off)\n"
      "\t--percentIdentity   arg  The minimum end-to-end %% identity of "
      "qualifying reads (default: 97)\n"
      "\t--pairedContigs     arg  The file to output the sparse matrix of "
      "contigs which paired reads span (default: none)\n"
      "\t--unmappedFastq     arg  The prefix to output unmapped reads from "
      "each bam file suffixed by 'bamfile.bam.fastq.gz'\n"
      "\t--noIntraDepthVariance   Do not include variance from mean depth "
      "along the contig\n"
      "\t--showDepth              Output a .depth file per bam for each contig "
      "base\n"
      "\t--minMapQual        arg  The minimum mapping quality necessary to "
      "count the read as mapped (default: 0)\n"
      "\t--weightMapQual     arg  Weight per-base depth based on the MQ of the "
      "read (i.e uniqueness) (default: 0.0 (disabled))\n"
      "\t--includeEdgeBases       When calculating depth & variance, include "
      "the 1-readlength edges (off by default)\n"
      "\t--maxEdgeBases           When calculating depth & variance, and not "
      "--includeEdgeBases, the maximum length (default:75)\n"
      "\t--referenceFasta    arg  The reference file.  (It must be the same "
      "fasta that bams used)\n"
      "\nOptions that require a --referenceFasta\n"
      "\t--outputGC          arg  The file to print the gc coverage histogram\n"
      "\t--gcWindow          arg  The sliding window size for GC calculations\n"
      "\t--outputReadStats   arg  The file to print the per read statistics\n"
      "\t--outputKmers       arg  The file to print the perfect kmer counts\n"
      "\nOptions to control shredding contigs that are under represented by "
      "the reads\n"
      "\t--shredLength       arg  The maximum length of the shreds\n"
      "\t--shredDepth        arg  The depth to generate overlapping shreds\n"
      "\t--minContigLength   arg  The mimimum length of contig to include for "
      "mapping and shredding\n"
      "\t--minContigDepth    arg  The minimum depth along contig at which to "
      "break the contig\n"
      "\n",
      RabbitBin_VERSION, RabbitBin_BUILD_DATE);
}

void abortMe(string msg) {
  cerr << msg << endl;
  exit(1);
}

bool file_exists(const string fname) {
  ifstream ifs(fname, std::ios_base::binary);
  return ifs.is_open();
}

string get_basename(const string fname) {
  string basename(fname);
  auto pos = basename.find_last_of('/');
  if (pos != string::npos)
    basename = basename.substr(pos + 1);
  pos = basename.find_last_of('.');
  if (pos != string::npos)
    basename = basename.substr(0, pos);
  return basename;
}

void printDepthTableHeader(ostream &of, const vector<string> &names,
                           bool intraDepthVariance) {
  of << "contigName\tcontigLen\ttotalAvgDepth";
  for (int bamIdx = 0; bamIdx < (int)names.size(); bamIdx++) {
    string shortName = get_basename(names[bamIdx]);
    of << "\t" << shortName;
    if (intraDepthVariance) {
      of << "\t" << shortName << "-var";
    }
  }
  of << "\n";
}

class ContigDepth {
public:
  int numSamples, contigLen;
  bool includeEdgeBases;
  int maxEdgeBases;
  double mean, variance;
  double weightedSum = 0.0;
  uint64_t totalCorrectedLength = 0;
  vector<double> means;
  vector<double> variances;
  ContigDepth() = default;
  ContigDepth(int numSamples, int contigLen, bool includeEdgeBases,
              int maxEdgeBases, bool intraDepthVariance)
      : numSamples(numSamples), contigLen(contigLen),
        includeEdgeBases(includeEdgeBases), maxEdgeBases(maxEdgeBases), mean{},
        variance{}, weightedSum{}, totalCorrectedLength{}, means{},
        variances{} {
    means.reserve(numSamples);
    if (intraDepthVariance)
      variances.resize(numSamples);
  }
  void addDepth(int averageReadSize, CountType &contigDepth) {
    if (!variances.empty())
      throw;
    auto correctedLen = getCorrectedLen(averageReadSize);
    totalCorrectedLength += correctedLen;
    double mean = contigDepth / correctedLen;
    if (mean < 0.0 || correctedLen < 0) {
      std::stringstream ss;
      ss << "WARNING: calculated a negative mean=" << mean
         << ". correctedLen=" << correctedLen << " contigDepth=" << contigDepth
         << "\nPlease report this BAM file to the RabbitBin maintainers\n";
      std::cerr << ss.str() << std::flush;
      contigDepth = 0;
      mean = 0;
    }
    if (mean > 1000000) {
      std::stringstream ss;
      ss << "WARNING: calculated a huge mean=" << mean
         << ". correctedLen=" << correctedLen << " contigDepth=" << contigDepth
         << "\nPlease report this BAM file to the RabbitBin maintainers\n";
      std::cerr << ss.str() << std::flush;
      contigDepth = 0;
      mean = 0;
    }

    weightedSum += mean * correctedLen;
    means.push_back(mean);
  }
  void addDepth(int averageReadSize, VarianceType &variance) {
    if (variances.empty())
      throw;
    auto correctedLen = getCorrectedLen(averageReadSize);
    if (variance.mean < 0.0 || variance.variance < 0.0 || correctedLen < 0) {
      std::stringstream ss;
      ss << "WARNING: received a negative mean=" << variance.mean
         << " or variance=" << variance.variance
         << ". correctedLen=" << correctedLen
         << "\nPlease report this BAM file to the RabbitBin maintainers\n";
      std::cerr << ss.str() << std::flush;
      variance.mean = 0;
      variance.variance = 0;
    }
    if (variance.mean > 1000000) {
      std::stringstream ss;
      ss << "WARNING: received a huge mean=" << variance.mean
         << " or variance=" << variance.variance
         << ". correctedLen=" << correctedLen
         << "\nPlease report this BAM file to the RabbitBin maintainers\n";
      std::cerr << ss.str() << std::flush;
      variance.mean = 0;
      variance.variance = 0;
    }
    totalCorrectedLength += correctedLen;
    variances[means.size()] = variance.variance;
    means.push_back(variance.mean);
    weightedSum += variance.mean * correctedLen;
  }
  int getCorrectedLen(int averageReadSize) {
    int x = 2 * std::min(maxEdgeBases, averageReadSize / 3);
    int correctedLen =
        (includeEdgeBases | (x >= contigLen)) ? contigLen : contigLen - x;
    assert(correctedLen > 0);
    return correctedLen;
  }
  double getTotalCorrectedDepth() const {
    double totalCorrectedDepth = 0.0;
    if (totalCorrectedLength > 0.0) {
      totalCorrectedDepth = weightedSum / totalCorrectedLength;
    }
    totalCorrectedDepth *= numSamples; // scale back up because we did not
                                       // divide by avgCorrectedLength
    return totalCorrectedDepth;
  }
  bool printDepth(ostream &of, string name, float minContigDepth) {
    auto totalAverageDepth = getTotalCorrectedDepth();
    std::string nameonly = name;
    int pos = nameonly.find_first_of(" \t");
    nameonly = nameonly.substr(0, pos);
    std::stringstream ss;
    ss << nameonly << "\t" << contigLen << "\t" << (float)totalAverageDepth;
    for (int i = 0; i < numSamples; i++) {
      if (variances.size()) {
        ss << "\t" << (float)means[i] << "\t" << (float)variances[i];
      } else {
        ss << "\t" << (float)means[i];
      }
    }
    ss << "\n";
    of << ss.str();
    return totalAverageDepth >= minContigDepth;
  }
};

void printDepthTable(ostream &of, const vector<string> &names,
                     bool intraDepthVariance, CountTypeMatrix &bamContigDepths,
                     const vector<int> &averageReadSize, float percentIdentity,
                     VarianceTypeMatrix &bamContigVariances,
                     bam_header_t *header, BoolVector &contigLengthPass,
                     BoolVector &contigDepthPass, float minContigDepth,
                     bool includeEdgeBases = false, int maxEdgeBases = 0) {

  printDepthTableHeader(of, names, intraDepthVariance);

  // Formatting ~240k rows × 40+ float fields is otherwise a serial tail. Split
  // the contig range into one CONTIGUOUS block per thread, format each block into
  // its own buffer, then emit buffers in order (preserves contig order exactly).
  const int32_t n = header->n_targets;
  const int nthreads = std::max(1, omp_get_max_threads());
  std::vector<std::string> chunks(nthreads);
#pragma omp parallel num_threads(nthreads)
  {
    const int tn = omp_get_thread_num();
    const int32_t lo = (int32_t)((int64_t)tn * n / nthreads);
    const int32_t hi = (int32_t)((int64_t)(tn + 1) * n / nthreads);
    std::ostringstream ss;
    for (int32_t contigIdx = lo; contigIdx < hi; contigIdx++) {
      if (!contigLengthPass[contigIdx])
        continue;
      int contigLen = header->target_len[contigIdx];
      ContigDepth contigDepth(names.size(), contigLen, includeEdgeBases,
                              maxEdgeBases, intraDepthVariance);
      for (int bamIdx = 0; bamIdx < (int)names.size(); bamIdx++) {
        if (intraDepthVariance)
          contigDepth.addDepth(averageReadSize[bamIdx],
                               bamContigVariances[bamIdx][contigIdx]);
        else
          contigDepth.addDepth(averageReadSize[bamIdx],
                               bamContigDepths[bamIdx].get()[contigIdx]);
      }
      contigDepthPass[contigIdx] =
          contigDepth.printDepth(ss, header->target_name[contigIdx],
                                 minContigDepth);
    }
    chunks[tn] = ss.str();
  }
  for (int t = 0; t < nthreads; t++)
    of << chunks[t];
}

void printPairedContigs(ostream &of, const int numThreads,
                        std::vector<PairedCountTypeMatrix> &bamPairedContigs,
                        bam_header_t *header) {
  of << "contigIdx\tcontigIdxMate\tAvgCoverage\n";
  for (int32_t contigIdx = 0; contigIdx < header->n_targets; contigIdx++) {
    float len = header->target_len[contigIdx];
    PairedCountType sums;
    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
      PairedCountTypeMatrix &pairedContigs = bamPairedContigs[threadNum];
      PairedCountType &pairedCounts = pairedContigs[contigIdx];
      if (pairedCounts.empty())
        continue;

      for (PairedCountType::const_iterator it = pairedCounts.begin();
           it != pairedCounts.end(); it++) {
        sums[it->first] += it->second;
      }
    }
    for (PairedCountType::const_iterator it = sums.begin(); it != sums.end();
         it++) {
      of << contigIdx << "\t" << it->first << "\t" << it->second / len << "\n";
    }
  }
}

int getGC(const char *seq, int len, bool isbam = false) {
  int count = 0, gc = 0;
  for (int i = 0; i < len; i++) {
    switch (isbam ? bam_nt16_rev_table[bam1_seqi((const uint8_t *)seq, i)]
                  : seq[i]) {
    case 'G':
    case 'C':
    case 'g':
    case 'c':
      gc++;
    case 'A':
    case 'T':
    case 'a':
    case 't':
      count++;
      break;
    }
  }
  if (count > 0) {
    return (int)((100.0 * (float)gc / (float)count) + 0.5);
  } else {
    return 0;
  }
}

int getGC(const bam1_t *b) {
  return getGC((const char *)bam1_seq(b), b->core.l_qseq, true);
}

ostream &writeReadStatsHeader(ostream &os) {
  ReadStatistics::writeHeader(os);
  os << "\tReadGC\tMappedGC\n";
  os.flush();
  return os;
}

ostream &writeReadStats(ostream &os, const bam1_t *b, ReadStatistics readstats,
                        const char *refseq, int refLen) {
  if ((b->core.flag & BAM_FUNMAP) == BAM_FUNMAP || !readstats.isValid() ||
      refseq == NULL)
    return os;
  if (b->core.pos < 0 ||
      (int64_t)b->core.pos + (int64_t)readstats.alignlen > refLen) {
    warnMalformedBamRead(
        b, "skipping read-stats output because the aligned span is outside "
           "the reference sequence");
    return os;
  }

  stringstream &ss = IOThreadBuffer::getMyBuffer(os);
  readstats.write(ss) << "\t" << getGC(b) << "\t"
                      << getGC(refseq + b->core.pos, readstats.alignlen, false)
                      << "\n";

  return os;
}

void writeUnmapedFastqFile(BamFile &myBam, string unmappedFastqFile,
                           BamHeaderT header, BoolVector &contigLengthPass,
                           BoolVector &contigDepthPass,
                           gzipFileBufPtr &bamUnmappedFastqof,
                           NameBamMap &bamReadIds) {
  cerr << "Sequestering reads to unmappedFastq files for " << myBam.getBamFile()
       << endl;
  // write any reads from low depth contigs

  cerr << "Sequestering reads on low abundance contigs for "
       << myBam.getBamName() << endl;
  for (int32_t contigIdx = 0; contigIdx < header->n_targets; contigIdx++) {
    if (contigLengthPass[contigIdx] && (!contigDepthPass[contigIdx])) {
      // for each bam find and write all reads on these newly failed contigs

      ostream os(bamUnmappedFastqof.get());
      BamUtils::writeFastqByContig(os, bamReadIds, myBam, contigIdx);
    }
  }

  cerr << "Sequestering read mates not yet output for " << myBam.getBamName()
       << endl;
  // for each bam write all orphans...

  string name = unmappedFastqFile + "-" + myBam.getBamName() + "-single.fastq";
  {

    SafeOfstream singles(name.c_str());
    ostream pairs(bamUnmappedFastqof.get());
    BamUtils::writeFastqOrphans(pairs, singles, bamReadIds, myBam);
    cerr << "Closing pair and single unmaps for " << myBam.getBamName() << endl;
    // close the two files now
    bamUnmappedFastqof.reset();
    cerr << "Closed pair and single unmaps for " << myBam.getBamName() << endl;
  }
  struct stat filestatus;
  stat(name.c_str(), &filestatus);

  if (filestatus.st_size <= 20) {
    unlink(name.c_str());
  } else {
    cerr << "Additional orphaned reads are output as singles to: " << name
         << endl;
  }

  cerr << "Freeing up memory for " << myBam.getBamName() << endl;
  bamReadIds.clear();
}

// ── Contig-sharded parallel depth (plain depth case; means only) ────────────
// The per-bam consumer is serial, so with only num_bams files the box tops out
// at ~num_bams busy cores. Coordinate-sorted BAMs let us shard each file by
// contig (tid) range and process ranges concurrently: a contig's reads live
// entirely in one shard, so each shard owns a DISJOINT set of contigs and
// writes its own per-contig depth with no locks and no cross-shard merge.
//
// Profiling showed the per-read time was dominated by (a) the RabbitBAM reader's
// mutex contention and (b) calculateVarianceContig + the per-base depthCounts
// array — and RabbitBin never reads the variance columns. So this fast path:
//   • reads raw records via header-free hts_open + bam_read1 (no reader locks,
//     no re-parsing the 6.1M-contig @SQ header per shard),
//   • computes ONLY the scalar per-contig depth (one caldepth() pass per read;
//     caldepth's return value is independent of the per-base buffer, so the mean
//     is byte-identical to the serial path), and
//   • leaves the variance columns at 0 (RabbitBin ignores them; see
//     rb_env_depth_parse_var()).
// avgRead is fixed per bam (reads are fixed length here) so the edge trim — and
// thus the depth — matches the serial path exactly.
struct DepthShard {
  int      bamIdx;
  int32_t  lo, hi;     // contig (tid) range [lo, hi)
  uint64_t startVoff;  // BGZF virtual offset of first record in the range
  bool     hasData;
};

// Parse a .bai directly to extract each non-empty contig's first record virtual
// offset (the first non-zero entry of its linear index). This SKIPS htslib's
// bam_index_load, which builds bin-hash structures for all 6.1M contigs and costs
// ~8s per file; we only read the raw .bai sequentially and grab the linear-index
// offsets we need for record-aligned shard starts. Returns false on any format
// surprise so the caller can fall back to the index path.
static bool parse_bai_first_voffs(const std::string &bamPath,
                                  std::vector<int32_t> &outTid,
                                  std::vector<uint64_t> &outVoff) {
  std::string bai = bamPath + ".bai";
  FILE *f = fopen(bai.c_str(), "rb");
  if (!f)
    return false;
  fseek(f, 0, SEEK_END);
  long sz = ftell(f);
  fseek(f, 0, SEEK_SET);
  if (sz < 8) {
    fclose(f);
    return false;
  }
  std::vector<char> buf(sz);
  if (fread(buf.data(), 1, sz, f) != (size_t)sz) {
    fclose(f);
    return false;
  }
  fclose(f);
  const char *p = buf.data();
  const char *end = p + sz;
  if (memcmp(p, "BAI\1", 4) != 0)
    return false;
  p += 4;
  int32_t n_ref;
  memcpy(&n_ref, p, 4);
  p += 4;
  for (int32_t r = 0; r < n_ref; r++) {
    if (p + 4 > end)
      return false;
    int32_t n_bin;
    memcpy(&n_bin, p, 4);
    p += 4;
    for (int32_t b = 0; b < n_bin; b++) {
      if (p + 8 > end)
        return false;
      p += 4; // bin id
      int32_t n_chunk;
      memcpy(&n_chunk, p, 4);
      p += 4;
      p += (size_t)n_chunk * 16; // skip chunk beg/end pairs
      if (p > end)
        return false;
    }
    if (p + 4 > end)
      return false;
    int32_t n_intv;
    memcpy(&n_intv, p, 4);
    p += 4;
    if ((size_t)(end - p) < (size_t)n_intv * 8)
      return false;
    uint64_t voff = 0;
    for (int32_t k = 0; k < n_intv; k++) {
      uint64_t o;
      memcpy(&o, p + (size_t)k * 8, 8);
      if (o != 0) {
        voff = o;
        break;
      }
    }
    p += (size_t)n_intv * 8;
    if (voff != 0) {
      outTid.push_back(r);
      outVoff.push_back(voff);
    }
  }
  return true;
}

static uint64_t shard_start_voff(hts_idx_t *idx, int32_t lo, int32_t hi,
                                 bool &hasData) {
  hasData = false;
  for (int32_t t = lo; t < hi; ++t) {
    hts_itr_t *it = sam_itr_queryi(idx, t, 0, HTS_POS_MAX);
    if (it) {
      if (it->n_off > 0) {
        uint64_t v = it->off[0].u;
        hts_itr_destroy(it);
        hasData = true;
        return v;
      }
      hts_itr_destroy(it);
    }
  }
  return 0;
}

static void process_depth_shard(const DepthShard &sh,
                                const std::string &bamPath,
                                const BamHeaderT &header,
                                CountType *contigDepths, float percentIdentity,
                                int maxEdgeBases, bool includeEdgeBases,
                                int avgRead, int minMapQual) {
  if (!sh.hasData)
    return;
  htsFile *fp = hts_open(bamPath.c_str(), "rb");
  if (!fp)
    return;
  BGZF *bgzf = fp->fp.bgzf;
  if (bgzf_seek(bgzf, (int64_t)sh.startVoff, SEEK_SET) < 0) {
    hts_close(fp);
    return;
  }
  bam1_t *b = bam_init1();
  DepthCounts none{};
  const int32_t n_targets = header->n_targets;
  const int edge = includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead / 3);

  while (bam_read1(bgzf, b) >= 0) {
    int32_t tid = b->core.tid;
    int32_t pos = b->core.pos;
    if (tid >= sh.hi)
      break;              // reached the next shard's contigs
    if (tid < sh.lo)
      continue;           // straggler from block-aligned seek, or unmapped(-1)
    if ((b->core.flag & BAM_FUNMAP) == 0 && (tid >= n_targets || pos < 0))
      continue;
    if ((b->core.flag & (BAM_FPAIRED | BAM_FMUNMAP)) == BAM_FPAIRED &&
        (b->core.mtid < 0 || b->core.mtid >= n_targets))
      continue;
    if (BamNameTrackingChooser::unsupportedRead(b))
      continue;
    if ((b->core.flag & BAM_FUNMAP) == BAM_FUNMAP || b->core.qual < minMapQual)
      continue;
    if (!CheckRead::checkEnd(b, header.get())) {
      b = CheckRead::fixEndClip(b, header.get());
      if (!CheckRead::checkEnd(b, header.get()))
        continue;
    }
    ReadStatistics rs;
    // One pass: caldepth's return (edge-trimmed overlap) is identical with or
    // without a per-base buffer, so this matches the serial depth exactly.
    CountType ov =
        caldepth(b, none, header->target_len[tid], NULL, edge, &rs);
    if (rs.isValid() && rs.getPctId() >= percentIdentity)
      contigDepths[tid] += ov;
  }
  bam_destroy1(b);
  hts_close(fp);
}

int main(int argc, char *argv[]) {

  // set and parse options
  float percentIdentity = (float)97 / (float)100.0;
  string outputTableFile, pairedContigsFile, unmappedFastqFile,
      referenceFastaFile;
  string outputReadStatsFile, outputGCFile, outputKmersFile;
  int gcWindow = 100;
  int shredLength = 16000, shredDepth = 5, minContigLength = 1, minMapQual = 0;
  float weightMapQual = 0.0;
  bool normalizeWeightMapQual = false;
  float minContigDepth = 0;
  bool intraDepthVariance = true;
  bool showDepth = false;
  bool includeEdgeBases = false;
  int maxEdgeBases = 75;
  bool checkpoint = false;

  while (1) {
    int option_index = 0;
    int c = getopt_long(argc, argv, "h", long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
    case 0:
      if (strcmp(long_options[option_index].name, "help") == 0) {
        usage();
        exit(0);
      } else if (strcmp(long_options[option_index].name, "percentIdentity") ==
                 0) {
        percentIdentity = atoi(optarg) / 100.0;
        cerr << "Minimum percent identity for a mapped read: "
             << percentIdentity << endl;
      } else if (strcmp(long_options[option_index].name, "outputDepth") == 0) {
        outputTableFile = optarg;
        cerr << "Output depth matrix to " << outputTableFile << endl;
      } else if (strcmp(long_options[option_index].name, "checkpoint") == 0) {
        checkpoint = true;
      } else if (strcmp(long_options[option_index].name, "outputReadStats") ==
                 0) {
        outputReadStatsFile = optarg;
        cerr << "Output Read Stats file to " << outputReadStatsFile << endl;
      } else if (strcmp(long_options[option_index].name, "outputGC") == 0) {
        outputGCFile = optarg;
        cerr << "Output GC stats file to " << outputGCFile << endl;
      } else if (strcmp(long_options[option_index].name, "gcWindow") == 0) {
        gcWindow = atoi(optarg);
        cerr << "GC sliding window: " << gcWindow << endl;
      } else if (strcmp(long_options[option_index].name, "outputKmers") == 0) {
        outputKmersFile = optarg;
        cerr << "Output Perfect Kmers file to " << outputKmersFile << endl;
      } else if (strcmp(long_options[option_index].name, "pairedContigs") ==
                 0) {
        pairedContigsFile = optarg;
        cerr << "Output pairedContigs lower triangle to " << pairedContigsFile
             << endl;
      } else if (strcmp(long_options[option_index].name, "unmappedFastq") ==
                 0) {
        unmappedFastqFile = optarg;
        cerr << "Output Unmapped Fastq to " << unmappedFastqFile << endl;
      } else if (strcmp(long_options[option_index].name, "referenceFasta") ==
                 0) {
        referenceFastaFile = optarg;
        cerr << "Reference fasta file " << referenceFastaFile << endl;
      } else if (strcmp(long_options[option_index].name, "shredLength") == 0) {
        shredLength = atoi(optarg);
        cerr << "shredLength: " << shredLength << endl;
      } else if (strcmp(long_options[option_index].name, "shredDepth") == 0) {
        shredDepth = atoi(optarg);
        cerr << "shredDepth: " << shredDepth << endl;
      } else if (strcmp(long_options[option_index].name, "minContigLength") ==
                 0) {
        minContigLength = atoi(optarg);
        cerr << "minContigLength: " << minContigLength << endl;
      } else if (strcmp(long_options[option_index].name, "minMapQual") == 0) {
        minMapQual = atoi(optarg);
        cerr << "minMapQual: " << minMapQual << endl;
      } else if (strcmp(long_options[option_index].name, "weightMapQual") ==
                 0) {
        weightMapQual = atof(optarg);
        cerr << "weightMapQual: " << weightMapQual << endl;
      } else if (strcmp(long_options[option_index].name, "minContigDepth") ==
                 0) {
        minContigDepth = atof(optarg);
        cerr << "minContigDepth: " << minContigDepth << endl;
      } else if (strcmp(long_options[option_index].name,
                        "noIntraDepthVariance") == 0) {
        intraDepthVariance = false;
        cerr << "Calculating intra contig depth variance\n";
      } else if (strcmp(long_options[option_index].name, "showDepth") == 0) {
        showDepth = true;
        cerr << "Outputing a .depth file for each bam\n";
      } else if (strcmp(long_options[option_index].name, "includeEdgeBases") ==
                 0) {
        includeEdgeBases = true;
        cerr << "Edge bases will be included in all calculations\n";
      } else if (strcmp(long_options[option_index].name, "maxEdgeBases") == 0) {
        maxEdgeBases = atoi(optarg);
        cerr << "Edge bases will be included up to " << maxEdgeBases
             << " bases\n";
      } else {
        usage();
        cerr << "Unrecognized option: " << long_options[option_index].name
             << endl;
        exit(1);
      }
      break;
    default:
      usage();
      exit(0);
    };
  }
  if (argc - optind < 1) {
    usage();
    if (argc - optind < 1)
      cerr << "You must specify one or more bam files\n\n";
    return 0;
  }
  CheckRead::_edgeBases() = maxEdgeBases;

  cerr << "rabbit_depth " << RabbitBin_VERSION << " "
       << RabbitBin_BUILD_DATE << endl;
  cerr << "Running with " << omp_get_max_threads()
       << " threads to save memory you can reduce the number of threads with "
          "the OMP_NUM_THREADS variable"
       << endl;
  cerr << "Output matrix to "
       << (outputTableFile.empty() ? "STDOUT" : outputTableFile.c_str())
       << endl;
  if (checkpoint)
    cerr << "Writing checkpoints as each bam has been processed" << endl;
  else
    cerr << "Not writing checkpoints, all depths will be in memory -- You can "
            "save memory by using --checkpoint"
         << endl;

  // assign names and allocate samfile handles
  BamFileVector bams;
  StringVector bamFilePaths;
  for (int i = optind; i < argc; i++) {
    bamFilePaths.push_back(argv[i]);
  }
  int num_bams = bamFilePaths.size();
  int largest_contig = 0;

  bool store_sequences = !outputGCFile.empty() | !outputKmersFile.empty() |
                         !outputReadStatsFile.empty();
  std::vector<string> referenceSequences;
  std::vector<int> referenceSequenceLens;
  std::vector<string> referenceSequenceNames;
  // The reference FASTA is read only to VALIDATE the contig set/lengths against
  // the BAM header (the depth math uses header->target_len). For the plain
  // sharded fast path we can skip the multi-second serial 3.8GB parse entirely.
  bool skip_fasta = pairedContigsFile.empty() && unmappedFastqFile.empty() &&
                    !checkpoint && !store_sequences && !showDepth &&
                    getenv("RABBIT_DEPTH_NO_SHARD") == NULL;
  if (!referenceFastaFile.empty() && !skip_fasta) {
    cout << "Reading reference fasta file: " << referenceFastaFile << endl;
    KseqReader reference(referenceFastaFile);
    while (reference.hasNext()) {
      auto seq = reference.getSeq();
      auto sz = seq.length();
      referenceSequenceLens.push_back(sz);
      if (sz > largest_contig)
        largest_contig = sz;
      if (store_sequences)
        referenceSequences.push_back(seq);
      referenceSequenceNames.push_back(reference.getName());
    }
    cout << "... " << referenceSequenceLens.size() << " sequences" << endl;
    if (referenceSequenceLens.empty()) {
      cerr << "ERROR: the reference was empty!: " << referenceFastaFile << endl;
      exit(1);
    }
  }

  SafeOfstream *readStats = NULL;
  if (!outputReadStatsFile.empty()) {
    readStats = new SafeOfstream(outputReadStatsFile.c_str());
    writeReadStatsHeader(*readStats);
  }

  std::vector<float> refGC;
  ReadGCStats readGCStats;

  std::vector<std::vector<uint8_t>> refGCWindows;
  if (!outputGCFile.empty() && !referenceSequences.empty()) {
    if (gcWindow <= 0) {
      cerr << "ERROR: gcWindow must be greater than 0" << endl;
      exit(1);
    }
    refGC.resize(101, 0);
    readGCStats.resize(101, RunningStats());
    refGCWindows.resize(referenceSequences.size());
    int skippedShortContigs = 0;
    for (long i = 0; i < (long)referenceSequences.size(); i++) {
      std::vector<uint8_t> &refGCs = refGCWindows[i];
      const size_t seqLen = referenceSequences[i].length();
      if (seqLen < (size_t)gcWindow) {
        skippedShortContigs++;
        continue;
      }
      const size_t numWindows = seqLen - gcWindow + 1;
      refGCs.reserve(numWindows);
      for (size_t j = 0; j < numWindows; j++) {
        int gc = getGC(referenceSequences[i].c_str() + j, gcWindow, false);
        refGCs.push_back(gc);
      }
    }
    if (skippedShortContigs > 0) {
      cerr << "Skipping GC window calculations for " << skippedShortContigs
           << " contigs shorter than gcWindow=" << gcWindow << endl;
    }
  }

  // open first bam
  bams.resize(num_bams);
  BamFile &firstBam = bams[0];
  BamHeaderT header;
  bool preprocess_complete = false;
  if (checkpoint) {
    cerr << "Opening the first bam file for a baseline header: "
         << bamFilePaths[0] << endl;
    firstBam = BamUtils::openBam(bamFilePaths[0], false);
    header = firstBam.header;
  } else if (skip_fasta) {
    // Fast (sharded) path: only the consolidated header is needed (each shard
    // opens its own handle), so parse ONE header instead of all num_bams.
    cerr << "Opening first bam for header (sharded fast path)" << endl;
    firstBam = BamUtils::openBam(bamFilePaths[0], false);
    header = firstBam.header;
    preprocess_complete = true;
  } else {
    cerr << "Opening all bam files and validating headers" << endl;
    header = BamUtils::openBamsAndConsolidateHeaders(bamFilePaths, bams, false);
    preprocess_complete = true;
  }

  assert(header.get() != NULL);
  BoolVector contigLengthPass, contigDepthPass;
  contigLengthPass.resize(header->n_targets);
  contigDepthPass.resize(header->n_targets, true);
  if (!referenceSequenceLens.empty() &&
      header->n_targets != (long)referenceSequenceLens.size()) {
    cerr << "Error: referenceFile: " << referenceFastaFile
         << " is not the same as in the bam headers! (targets: "
         << header->n_targets << " from the bam vs "
         << referenceSequenceLens.size() << " from the ref)" << endl;
    if (referenceSequenceLens.empty())
      cerr << "no reference was loaded for " << referenceFastaFile << endl;
    exit(1);
  }
  for (int32_t i = 0; i < header->n_targets; i++) {
    if ((int)header->target_len[i] >= minContigLength) {
      contigLengthPass[i] = true;
    } else {
      contigLengthPass[i] = false;
    }
    if (!referenceSequenceLens.empty() &&
        header->target_len[i] != referenceSequenceLens[i]) {
      cerr << "Error: referenceFile: " << referenceFastaFile << " contig " << i
           << " is not the same as in the bam headers (bam reports "
           << header->target_name[i] << " with " << header->target_len[i]
           << " len, reference loaded " << referenceSequenceNames[i] << " with "
           << referenceSequenceLens[i] << " len)! " << endl;
      exit(1);
    }
  }

  // make vector of unmappedFastq files, to reuse
  std::vector<gzipFileBufPtr> bamUnmappedFastqof;
  bamUnmappedFastqof.resize(num_bams);

  // allocate memory for depth and optionally variance matrixes
  CountTypeMatrix bamContigDepths;
  bamContigDepths.resize(num_bams);
  vector<int> averageReadSize;
  averageReadSize.resize(num_bams, 0);

  VarianceTypeMatrix bamContigVariances;
  if (intraDepthVariance) {
    bamContigVariances.resize(num_bams);
    if (!checkpoint) {
      VarianceType dummy;
      for (int i = 0; i < (int)num_bams; i++) {
        bamContigVariances[i].resize(header->n_targets, dummy);
      }
    }
  }

  std::vector<NameBamMap> bamReadIds;
  bamReadIds.resize(num_bams);

  int numThreads = 1;
#pragma omp parallel
  {
    if (omp_get_thread_num() == 0)
      numThreads = omp_get_num_threads();
#pragma omp for reduction(max : largest_contig)
    for (int64_t i = 0; i < (int64_t)header->n_targets; i++)
      if (header->target_len[i] > largest_contig)
        largest_contig = header->target_len[i];
  }
  // Full online-core count, captured before the per-bam cap below; the
  // contig-sharded fast path uses ALL cores, not just num_bams.
  const int g_full_threads = numThreads;
#ifdef USE_RABBITBAM
  const int rbam_total_threads = numThreads;
  int rbam_read_threads =
      std::max(1, rbam_total_threads / std::max(1, (int)num_bams));
#endif
  if (numThreads > (int)num_bams) {
    numThreads = num_bams;
    omp_set_num_threads(numThreads);
  }
  std::vector<PairedCountTypeMatrix> bamPairedContigs;
  bamPairedContigs.resize(numThreads);

  if (!pairedContigsFile.empty()) {
    cerr << "Allocating pairedContigs matrix: "
         << (numThreads * header->n_targets * sizeof(PairedCountType) / 1024 /
             1024)
         << " MB over " << numThreads << " threads" << endl;

#pragma omp parallel for schedule(static, 1)
    for (int threadNum = 0; threadNum < numThreads; threadNum++) {

      PairedCountTypeMatrix &pairedContigs = bamPairedContigs[threadNum];
      PairedCountType empty;
      pairedContigs.resize(header->n_targets, empty);
      assert((int)pairedContigs.size() == header->n_targets);
      for (int64_t i = 0; i < (int64_t)header->n_targets; i++) {
        if (!pairedContigs[i].empty())
          throw;
      }
    }
  }

  cerr << "Processing bam files with largest_contig=" << largest_contig << endl;
  bool hasAnyPairedContigs = false;

  // preallocate all forking for unmapped fastq, so we do not fork when memory
  // is tight
  if (!unmappedFastqFile.empty()) {
#pragma omp parallel for
    for (int bamIdx = 0; bamIdx < (int)num_bams; bamIdx++) {
      string name =
          unmappedFastqFile + "-" + bamFilePaths[bamIdx] + ".fastq.gz";
      bamUnmappedFastqof[bamIdx] = gzipOutputFile(name);
#pragma omp critical(OUTPUT)
      cerr << "Outputting any unmapped reads to " << name << endl;
    }
  }
  MappedKmersStats *ourMappedKmersStats = NULL;
  if (!outputKmersFile.empty() && !referenceSequences.empty()) {
    ourMappedKmersStats = new MappedKmersStats();
  }
  // Shared by both code paths (used after the consumer loop and at output).
  bool isSorted = true;
  int minAvgRead = maxEdgeBases * 3;
  vector<string> partialFiles(num_bams);

  // ── Fast path: contig-sharded means-only depth (see process_depth_shard) ──
  bool plain_mode = pairedContigsFile.empty() && unmappedFastqFile.empty() &&
                    !checkpoint && readGCStats.empty() &&
                    outputKmersFile.empty() && outputReadStatsFile.empty() &&
                    !showDepth && referenceSequences.empty() && weightMapQual == 0.0;
  bool use_sharded = plain_mode && getenv("RABBIT_DEPTH_NO_SHARD") == NULL;
  std::vector<DepthShard> shards;
  std::vector<int> sampledAvgRead(num_bams, 0);
  if (use_sharded) {
    // Many shards, balanced by COMPRESSED BYTES (≈ read count ≈ decode work),
    // not contig length: high-coverage contigs hold far more reads per base, so
    // length-balanced shards leave decode badly skewed (→ idle cores). The .bai
    // gives each contig's start virtual offset; coffset = voff>>16 is its
    // position in the compressed file, which we sample to cut equal-byte ranges.
    double shard_mult = 24.0;
    if (const char *e = getenv("RABBIT_DEPTH_SHARD_MULT"))
      shard_mult = atof(e) > 0 ? atof(e) : shard_mult;
    const int shards_per_bam =
        std::max(1, (int)std::lround((double)g_full_threads * shard_mult / num_bams));
    std::atomic<bool> all_indexed{true};
    std::vector<std::vector<DepthShard>> perBam(num_bams);
    omp_set_num_threads(std::min(g_full_threads, (int)num_bams));
#pragma omp parallel for schedule(dynamic, 1)
    for (int bi = 0; bi < (int)num_bams; bi++) {
      // Read the .bai directly for every non-empty contig's first record voff
      // (skips htslib's ~8s/file bin-hash construction).
      std::vector<int32_t> ntid;
      std::vector<uint64_t> nvoff;
      bool ok = parse_bai_first_voffs(bamFilePaths[bi], ntid, nvoff);
      if (!ok || ntid.empty()) {
        all_indexed = false;
        continue;
      }
      uint64_t firstVoff = nvoff.front();
      // Sample read length header-free: seek to the first record and scan.
      {
        htsFile *fp = hts_open(bamFilePaths[bi].c_str(), "rb");
        if (fp) {
          bgzf_seek(fp->fp.bgzf, (int64_t)firstVoff, SEEK_SET);
          bam1_t *sb = bam_init1();
          int64_t srs = 0, src = 0;
          while (src < 50000 && bam_read1(fp->fp.bgzf, sb) >= 0) {
            if ((sb->core.flag & BAM_FUNMAP) == 0) {
              srs += sb->core.l_qseq;
              src++;
            }
          }
          bam_destroy1(sb);
          hts_close(fp);
          sampledAvgRead[bi] = src ? (int)(srs / src) : 0;
        }
      }
      // Select up to shards_per_bam boundaries at equal compressed-byte spacing
      // among the (many) non-empty contigs.
      std::vector<DepthShard> sv;
      const uint64_t maxc = nvoff.back() >> 16;
      std::vector<size_t> bidx{0};
      size_t si = 0;
      for (int k = 1; k < shards_per_bam; k++) {
        uint64_t target = (uint64_t)((double)k / shards_per_bam * maxc);
        while (si < nvoff.size() && (nvoff[si] >> 16) < target)
          si++;
        if (si < ntid.size() && si > bidx.back())
          bidx.push_back(si);
      }
      for (size_t bb = 0; bb < bidx.size(); bb++) {
        size_t s0 = bidx[bb];
        int32_t lo = (bb == 0) ? 0 : ntid[s0];
        int32_t hi = (bb + 1 < bidx.size()) ? ntid[bidx[bb + 1]]
                                            : header->n_targets;
        if (hi <= lo)
          continue;
        uint64_t sv0 = (bb == 0) ? firstVoff : nvoff[s0];
        sv.push_back(DepthShard{bi, lo, hi, sv0, true});
      }
      perBam[bi] = std::move(sv);
    }
    if (!all_indexed) {
      use_sharded = false;
      cerr << "Some bams lack a .bai index; using the per-bam depth path"
           << endl;
    } else {
      for (auto &sv : perBam)
        for (auto &sh : sv)
          shards.push_back(sh);
    }
  }

  if (use_sharded) {
    cerr << "Contig-sharded depth (means-only): " << shards.size()
         << " shards over " << g_full_threads << " threads" << endl;
    for (int bi = 0; bi < (int)num_bams; bi++) {
      bamContigDepths[bi].reset(new CountType[header->n_targets]());
      if (intraDepthVariance)
        bamContigVariances[bi].assign(header->n_targets, VarianceType());
      averageReadSize[bi] = sampledAvgRead[bi];
      if (minAvgRead > averageReadSize[bi])
        minAvgRead = averageReadSize[bi];
    }
    omp_set_num_threads(g_full_threads);
#pragma omp parallel for schedule(dynamic, 1)
    for (int s = 0; s < (int)shards.size(); s++) {
      process_depth_shard(shards[s], bamFilePaths[shards[s].bamIdx], header,
                          bamContigDepths[shards[s].bamIdx].get(),
                          percentIdentity, maxEdgeBases, includeEdgeBases,
                          averageReadSize[shards[s].bamIdx], minMapQual);
    }
    // Output reads the per-bam mean from VarianceType.mean. Fill it from the
    // scalar depth sum using calculateVarianceContig's exact formula
    // (mean = sum(edge-trimmed coverage) / adjustedContigLength); variance stays
    // 0 (RabbitBin ignores it). This makes the mean depth byte-identical to the
    // serial path without ever building the per-base array.
    if (intraDepthVariance) {
#pragma omp parallel for schedule(static)
      for (int bi = 0; bi < (int)num_bams; bi++) {
        const int edge =
            includeEdgeBases ? 0 : std::min(maxEdgeBases, averageReadSize[bi] / 3);
        const CountType *cd = bamContigDepths[bi].get();
        auto &vars = bamContigVariances[bi];
        for (int32_t t = 0; t < header->n_targets; t++) {
          int32_t cl = header->target_len[t];
          int32_t adj = (cl > 2 * edge + 1) ? (cl - 2 * edge) : cl;
          vars[t].mean = (adj > 2) ? (float)((double)cd[t] / adj) : 0.0f;
          vars[t].variance = 0.0f;
        }
      }
    }
  } else {
  DepthCounts noDepthCounts{};
  vector<DepthCounts> workingDepthCounts(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1)
  for (int threadIdx = 0; threadIdx < workingDepthCounts.size(); threadIdx++)
    workingDepthCounts[threadIdx].resetBaseCounts(largest_contig + 1,
                                                  weightMapQual > 0.0);

  std::mutex m;
  std::condition_variable cv;
  CountTypeMatrix threadWorkingContigDepths(omp_get_max_threads());
  VarianceTypeMatrix threadWorkingVariance(omp_get_max_threads());

#pragma omp parallel for schedule(static, 1)
  for (int bamIdx = 0; bamIdx < (int)num_bams; bamIdx++) {

    int threadNum = omp_get_thread_num();
    string baseName = get_basename(bamFilePaths[bamIdx]);

    if (checkpoint && bamIdx == 0) {
      // block processing all bam file alignments until preprocessing of contig
      // name 2 id map
      std::lock_guard<std::mutex> lk(m);
#pragma omp critical(OUTPUT)
      cerr << "Thread " << threadNum << " preprocessing contig2id map" << endl;
      BamUtils::preprocessContigNames(header);
      preprocess_complete = true;
      cv.notify_all();
    }

    string checkpoint_depth = partialFiles[bamIdx] =
        outputTableFile + ".partial-" + baseName + ".data";
    if (checkpoint && file_exists(checkpoint_depth)) {
#pragma omp critical(OUTPUT)
      cerr << "partial output for bam " << bamIdx << ": " << baseName
           << " already checkpointed." << endl;
      continue;
    }

    BamFile &myBam = bams[bamIdx];
#pragma omp critical(OUTPUT)
    cerr << "Thread " << threadNum
         << " opening and reading the header for file: " << bamFilePaths[bamIdx]
         << endl;
    if (bamIdx > 0)
      myBam = BamUtils::openBam(bamFilePaths[bamIdx], false, header);
#pragma omp critical(OUTPUT)
    cerr << "Thread " << threadNum
         << " opened the file: " << bamFilePaths[bamIdx] << endl;
    if (checkpoint && !preprocess_complete) {
#pragma omp critical(OUTPUT)
      cerr << "Thread " << threadNum
           << " waiting for preprocessing to complete before reading "
              "alignments in : "
           << bamFilePaths[bamIdx] << endl;
      std::unique_lock<std::mutex> lk(m);
      if (!preprocess_complete) {
        cv.wait(lk, [&preprocess_complete] { return preprocess_complete; });
      }
    }
    ostream *unmappedFastq = NULL;
    if (!unmappedFastqFile.empty()) {
      unmappedFastq = new std::ostream(bamUnmappedFastqof[bamIdx].get());
    }

#pragma omp critical(OUTPUT)
    cerr << "Thread " << threadNum << " processing bam " << bamIdx << ": "
         << myBam.getBamName() << endl;

    // initialize and allocate memory structures
    if (!checkpoint)
      bamContigDepths[bamIdx].swap(threadWorkingContigDepths[threadNum]);
    if (!bamContigDepths[bamIdx])
      bamContigDepths[bamIdx].reset(new CountType[header->n_targets]);
    else
      memset(bamContigDepths[bamIdx].get(), 0,
             header->n_targets * sizeof(CountType));
    if (!bamContigDepths[bamIdx]) {
      cerr << "Could not allocate enough memory to track depth per contig"
           << endl;
      exit(1);
    }

    if (intraDepthVariance) {
      assert(bamContigVariances.size() > bamIdx);
      if (!checkpoint)
        bamContigVariances[bamIdx].swap(threadWorkingVariance[threadNum]);
      if (bamContigVariances[bamIdx].empty()) {
        VarianceType dummy{};
        bamContigVariances[bamIdx].resize(header->n_targets, dummy);
      } else {
        for (auto &x : bamContigVariances[bamIdx])
          x.reset();
      }
    }
    CountType *contigDepths = bamContigDepths[bamIdx].get();
    PairedCountTypeMatrix &pairedContigs = bamPairedContigs[threadNum];
    NameBamMap &readIds = bamReadIds[bamIdx];
    NameBamMap *tempMates = NULL; // new BamNameTrackingChooser();

    int lastMinPos = -1;

    DepthCounts &depthCounts = workingDepthCounts[omp_get_thread_num()];
    std::shared_ptr<std::ofstream> depthFile;
    if (intraDepthVariance || !readGCStats.empty()) {
      depthCounts.resetBaseCounts(largest_contig + 1, weightMapQual > 0.0);
      if (showDepth) {
        depthFile.reset(
            new SafeOfstream(string(myBam.getFilePath() + ".depth").c_str()));
      }
    }

    MappedKmersStats *mappedKmersStats = NULL;
    if (!outputKmersFile.empty() && !referenceSequences.empty()) {
      mappedKmersStats = new MappedKmersStats();
    }
    bam1_t *b = myBam.getBamCache().getBam(),
           *lastBam = myBam.getBamCache().getBam();
    int bytesRead = 0, lastTid = -1, lastPos = 0;
    int64_t readSizes = 0, readCounts = 0, readsWellMapped = 0;
    int &avgRead = averageReadSize[bamIdx];

    CheckRead *check = new CheckRead(avgRead, contigLengthPass, header);
    readIds.setTrackNamer(check); // BamNameMap now manages this memory

    // read the bam file
#ifdef USE_RABBITBAM
    // Per-file parallel reader. Records arrive in file order (== sam_read1
    // order), and b->core.tid is decoded against this file's own @SQ table, the
    // same table sam_read1 would use, so all downstream tid indexing into the
    // consolidated header is identical.
    // 3rd arg (single_parser/is_tgs) MUST be true: it tells BamReader to prime
    // un_comp for the serial getBam1_t() consumer. With false, un_comp is left
    // uninitialised and getBam1_t() dereferences garbage.
    BamReader *rbamReader =
        new BamReader(bamFilePaths[bamIdx], rbam_read_threads, true);
#endif
    while (true) {
      std::swap(b, lastBam);
#ifdef USE_RABBITBAM
      if (!rbamReader->getBam1_t(b))
        break;
      bytesRead = 1;  // success sentinel (matches sam_read1 >= 0)
#elif defined(LEGACY_SAMTOOLS)
      bytesRead = samread(myBam, b);
      if (bytesRead <= 0)
        break;
#else
      bytesRead = sam_read1(myBam, header.get(), b);
      if (bytesRead < 0)
        break;
#endif

      int32_t tid = b->core.tid;
      int32_t pos = b->core.pos;
      readSizes += b->core.l_qseq;
      readCounts++;
      if (tid >= 0) {
        if (lastTid > tid || (lastTid == tid && lastPos > pos)) {
#pragma omp critical(BAM_WARN_UNSORTED)
          cerr << "ERROR: the bam file '" << myBam.getBamName()
               << "' is not sorted!" << endl;
          isSorted = false;
          break;
        }
      }

      if (!isSorted) // to show all unsorted bams at once
        break;

      assert(tid == -1 || lastTid <= tid); // ensure this is a sorted bam!
      if ((b->core.flag & BAM_FUNMAP) == 0 &&
          (tid < 0 || tid >= header->n_targets || pos < 0)) {
        warnMalformedBamRead(
            b, "mapped read has an invalid target id or negative position");
        continue;
      }
      if ((b->core.flag & (BAM_FPAIRED | BAM_FMUNMAP)) == BAM_FPAIRED &&
          (b->core.mtid < 0 || b->core.mtid >= header->n_targets)) {
        warnMalformedBamRead(
            b, "paired mapped read has a mate target id outside the header");
        continue;
      }
      if (check->unsupportedRead(b)) {
        continue;
      }

      bool printed = false;
      if (unmappedFastq != NULL) {
        if (check->failedMapping(b)) {
          // write, no name tracking
          BamUtils::writeFastqOrStorePair(*unmappedFastq, b, readIds);
          printed = true;
        }
      }
      if ((b->core.flag & BAM_FUNMAP) == BAM_FUNMAP ||
          b->core.qual < minMapQual) {
        // exclude this read from counts
        continue;
      }

      if (mappedKmersStats != NULL && !referenceSequences.empty()) {
        mappedKmersStats->addKmerStats(b, referenceSequences);
      }

      avgRead =
          readSizes / readCounts; // running avg.  recalculate after every read
      if (depthCounts && lastTid != tid) {
        // calculate statistics for the previous contig

        if (lastTid >= 0) {
          if (intraDepthVariance)
            bamContigVariances[bamIdx][lastTid] = calculateVarianceContig(
                header.get(), lastTid, depthCounts,
                includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead / 3),
                weightMapQual, normalizeWeightMapQual);
          if (!readGCStats.empty() && depthCounts &&
              contigDepths[lastTid] / header->target_len[lastTid] >
                  minContigDepth) {
            addGCCounts(readGCStats, refGCWindows[lastTid], gcWindow,
                        depthCounts);
          }
          if (showDepth)
            *depthFile << getContigDepthByBase(
                header.get(), lastTid, depthCounts,
                includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead / 3),
                weightMapQual);
          while ((tid >= 0 && lastTid < tid) ||
                 (tid < 0 && lastTid < header->n_targets)) {
            if (++lastTid == tid || lastTid >= header->n_targets) {
              lastTid--;
              lastPos = 0;
              break;
            }
            if (intraDepthVariance)
              bamContigVariances[bamIdx][lastTid] = VarianceType();
            if (showDepth) {
              depthCounts.resetBaseCounts(header->target_len[lastTid],
                                          weightMapQual > 0.0);
              *depthFile << getContigDepthByBase(
                  header.get(), lastTid, depthCounts,
                  includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead / 3),
                  weightMapQual);
            }
          }
        }
      } // finish stats for previous contigs

      if (depthCounts && lastTid != tid && tid >= 0) {
        depthCounts.resetBaseCounts(header->target_len[tid],
                                    weightMapQual > 0.0);
      }
      if (lastTid != tid || lastMinPos < pos - 2500) {
        if (lastTid != tid) {
          lastMinPos = -1;
        } else {
          lastMinPos = pos - 2000;
        }
        if (tempMates != NULL)
          tempMates->erase(tid, lastMinPos);
      }
      lastTid = tid;
      lastPos = pos;

      if ((b->core.flag & BAM_FUNMAP) == 0) {
        // check the end soft clip
        if (!CheckRead::checkEnd(b, header.get())) {
          b = CheckRead::fixEndClip(b, header.get());
          if (!CheckRead::checkEnd(b, header.get())) {
            warnMalformedBamRead(
                b, "skipping mapped read because it still extends beyond the "
                   "reference after end clipping");
            continue;
          }
        }
        // This read is mapped, apply its depth
        ReadStatistics rs;
        const char *refseq = referenceSequences.empty() || b->core.tid < 0 ||
                                     ((b->core.flag & BAM_FUNMAP) == BAM_FUNMAP)
                                 ? NULL
                                 : referenceSequences[b->core.tid].data();
        const int refLen =
            (refseq != NULL && b->core.tid >= 0) ? header->target_len[tid] : -1;
        CountType overlapAdjusted,
            overlapRaw = caldepth(b, noDepthCounts, refLen, refseq, 0,
                                  &rs); // just quickly calculate the overlap
        if (readStats != NULL && !referenceSequences.empty() &&
            b->core.tid >= 0 && rs.isValid()) {
          writeReadStats(*readStats, b, rs,
                         referenceSequences[b->core.tid].data(),
                         header->target_len[b->core.tid]);
        }

        bool failedPctId = !rs.isValid() || rs.getPctId() < percentIdentity;
        // validate this read is a good match
        if (failedPctId) {
          if (unmappedFastq != NULL && !printed) {
            // write & track the name
            if ((b->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
              string baseName = BamUtils::getBaseName(b);
              if (tempMates != NULL) {
                NameBamMap::iterator it = tempMates->find(baseName);
                if (it != tempMates->end()) {
                  BamUtils::writeFastq(*unmappedFastq, b, it->second);
                  tempMates->erase(it);
                  readIds.erase(baseName);
                  printed = true;
                }
              }
            }
            if (!printed) {
              BamUtils::writeFastqOrStorePair(*unmappedFastq, b, readIds);
              printed = true;
            }
          }
          continue;
        }

        if (depthCounts || !includeEdgeBases) {
          // now apply the overlap to baseCounts and adjust for edges
          overlapAdjusted = caldepth(
              b, depthCounts, header->target_len[tid], refseq,
              includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead / 3));
        } else {
          overlapAdjusted = overlapRaw;
        }
        contigDepths[tid] += overlapAdjusted;
        readsWellMapped++;

        // check for edge read
        if (unmappedFastq != NULL && !printed && check->edgeRead(b)) {
          // print out any edge reads, no name tracking
          BamUtils::writeFastqOrStorePair(*unmappedFastq, b, readIds);
          printed = true;
        }

        if (unmappedFastq != NULL &&
            (b->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {

          string baseName = BamUtils::getBaseName(b);
          // check for mate edge effects
          if (!printed && check->edgeMateRead(b)) {
            // mate is edge and will be printed, not name tracking

            if (tempMates != NULL) {
              NameBamMap::iterator it = tempMates->find(baseName);
              if (it != tempMates->end()) {
                BamUtils::writeFastq(*unmappedFastq, b, it->second);
                tempMates->erase(it);
                readIds.erase(baseName);
                printed = true;
              }
            }
            if (!printed) {
              BamUtils::writeFastqOrStorePair(*unmappedFastq, b, readIds);
              printed = true;
            }
          }
          if (!printed) {
            // print if the pair is already stored to print
            bool stored = BamUtils::writeFastqOrStorePair(*unmappedFastq, b,
                                                          readIds, true);
            printed = true;
            if (!stored && tempMates != NULL) {
              tempMates->insert(baseName, b);
            }
          }
        } // if unmappedFastq && unmapped

        if (!pairedContigsFile.empty() &&
            (b->core.flag & (BAM_FPAIRED | BAM_FMUNMAP)) == BAM_FPAIRED &&
            tid >= 0 && b->core.mtid >= 0) {
          assert(tid < (int)pairedContigs.size());
          PairedCountType &pairedCounts = pairedContigs[tid];
          PairedCountType::iterator it = pairedCounts.find(b->core.mtid);
          if (it != pairedCounts.end()) {
            it->second += overlapRaw;
          } else {
            pairedCounts.insert(
                it, PairedCountType::value_type(b->core.mtid, overlapRaw));
          }

          hasAnyPairedContigs = true;
        }
      } // if mapped
      if (tempMates != NULL)
        delete tempMates;

    } // while read the bam file

#ifdef USE_RABBITBAM
    delete rbamReader;  // joins decompression worker threads
#endif

    if (unmappedFastq != NULL)
      delete unmappedFastq;

    if (mappedKmersStats != NULL && ourMappedKmersStats != NULL) {
#pragma omp critical(MAPPED_KMERS_STATS)
      { *ourMappedKmersStats += *mappedKmersStats; }
      delete mappedKmersStats;
    }

#pragma omp critical(OUTPUT)
    cerr << "Thread " << threadNum << " finished reading bam " << bamIdx << ": "
         << myBam.getBamName() << endl;

    if (!isSorted) // skip
      continue;

    myBam.getBamCache().putBam(b);
    myBam.getBamCache().putBam(lastBam);

    // calculate the statistics for the last contig
    if (intraDepthVariance) {
      while (lastTid >= 0 && lastTid < header->n_targets) {
        assert(bamContigVariances.size() > bamIdx &&
               bamContigVariances[bamIdx].size() > lastTid);
        bamContigVariances[bamIdx][lastTid] = calculateVarianceContig(
            header.get(), lastTid, depthCounts,
            includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead / 3),
            weightMapQual, normalizeWeightMapQual);
        if (!readGCStats.empty() && depthCounts &&
            contigDepths[lastTid] / header->target_len[lastTid] >
                minContigDepth) {
          addGCCounts(readGCStats, refGCWindows[lastTid], gcWindow,
                      depthCounts);
        }

        if (showDepth)
          *depthFile << getContigDepthByBase(
              header.get(), lastTid, depthCounts,
              includeEdgeBases ? 0 : std::min(maxEdgeBases, avgRead / 3),
              weightMapQual);
        lastTid++;
        lastPos = 0;
        if (lastTid < header->n_targets) {
          depthCounts.resetBaseCounts(header->target_len[lastTid],
                                      weightMapQual > 0.0);
        }
      }
    } // intraDepthVariance

    if (checkpoint) {
#pragma omp critical(OUTPUT)
      cerr << "Writing partial output for bam " << bamIdx << ": "
           << myBam.getBamName() << " as " << checkpoint_depth << endl;
      string tmp = checkpoint_depth + ".tmp";
      ofstream out(tmp);
      assert(contigDepths == bamContigDepths[bamIdx].get());

      for (int i = 0; i < header->n_targets; i++) {
        if (intraDepthVariance) {
          auto &bcv = bamContigVariances[bamIdx];
          auto &cv = bcv[i];
          out.write((const char *)&cv, sizeof(VarianceType));
        } else {
          out.write((const char *)(contigDepths + i), sizeof(*contigDepths));
        }
      }
      threadWorkingContigDepths[threadNum].swap(bamContigDepths[bamIdx]);
      if (intraDepthVariance) {
        auto &bcv = bamContigVariances[bamIdx];
        threadWorkingVariance[threadNum].swap(bcv);
      }
      out.close();
      if (rename(tmp.c_str(), checkpoint_depth.c_str()) != 0) {
        cerr << "Could not rename " << tmp << " to " << checkpoint_depth << "! "
             << strerror(errno);
        exit(1);
      }
    } // checkpoint

#pragma omp critical(OUTPUT)
    cerr << "Thread " << threadNum << " finished: " << myBam.getBamName()
         << " with " << readCounts << " reads and " << readsWellMapped
         << " readsWellMapped (" << 100. * readsWellMapped / readCounts << "%)"
         << endl;
    if (readStats != NULL) {
      IOThreadBuffer::flush(*readStats);
    }
    if (!unmappedFastqFile.empty()) {
      writeUnmapedFastqFile(myBam, unmappedFastqFile, header, contigLengthPass,
                            contigDepthPass, bamUnmappedFastqof[bamIdx],
                            bamReadIds[bamIdx]);
    }
    if (bamIdx > 0) {
#pragma omp critical(OUTPUT)
      cerr << "Closing bam " << bamIdx << ": " << myBam.getBamFile() << endl;
      myBam.close();
    }
    if (minAvgRead > avgRead)
      minAvgRead = avgRead;
  } // foreach bamIdx
  workingDepthCounts.clear();
  } // end per-bam (non-sharded) path
  if (!isSorted) {
    cerr << "Please execute 'samtools sort' on unsorted input bam files and "
            "try again!"
         << endl;
    exit(1);
  }

  if (readStats != NULL) {
    IOThreadBuffer::close(*readStats);
    delete readStats;
    readStats = NULL;
  }

  if (ourMappedKmersStats != NULL && !outputKmersFile.empty()) {
    SafeOfstream outkmers(outputKmersFile.c_str());
    ourMappedKmersStats->writeHeader(outkmers);
    ourMappedKmersStats->write(outkmers);
    delete ourMappedKmersStats;
    ourMappedKmersStats = NULL;
  }

  // output the matrix
  streambuf *buf;
  ofstream of;
  if (!outputTableFile.empty() || outputTableFile.compare("-") != 0) {
    of.open(outputTableFile.c_str());
    buf = of.rdbuf();
  } else {
    buf = cout.rdbuf();
  }
  ostream out(buf);
  if (checkpoint) {
    cerr << "All partial depths are written, combining partials into "
         << outputTableFile << endl;
    printDepthTableHeader(out, bamFilePaths, intraDepthVariance);
    vector<ifstream> partials;
    partials.reserve(partialFiles.size());
    for (auto &checkpoint_file : partialFiles) {
      partials.emplace_back(checkpoint_file);
    }
    // read many contigs at a time for each partial file for better IO
    // performance
    const int batch_size = 40960;
    std::vector<ContigDepth> batchContigDepths;
    std::vector<VarianceType> batchContigVariances;
    std::vector<CountType> batchContigDepthsCounts;
    for (int contigIdxStart = 0; contigIdxStart < header->n_targets;
         contigIdxStart += batch_size) {
      batchContigDepths.clear();
      batchContigDepths.resize(batch_size);
      if (intraDepthVariance) {
        batchContigVariances.clear();
        batchContigVariances.resize(batch_size * partials.size());
      } else {
        batchContigDepthsCounts.clear();
        batchContigDepthsCounts.resize(batch_size * partials.size());
      }
      cerr << "Wrote " << contigIdxStart << " of " << header->n_targets << ": "
           << 100.0 * contigIdxStart / header->n_targets << "%\r";
      for (int batchIdx = 0; batchIdx < batch_size; batchIdx++) {
        int contigIdx = contigIdxStart + batchIdx;
        if (contigIdx >= header->n_targets)
          break;
        if (!contigLengthPass[contigIdx])
          continue;
        auto contigLen = header->target_len[contigIdx];
        auto &batchContigDepth = batchContigDepths[batchIdx];

        batchContigDepth = ContigDepth(
            partials.size(), contigLen, includeEdgeBases,
            std::min(maxEdgeBases, minAvgRead / 3), intraDepthVariance);
      }

// read all partials for this batch of contigs into memory for processing
// this is done to reduce the number of IO calls and improve performance
#pragma omp parallel for schedule(static, 1)
      for (int bamIdx = 0; bamIdx < (int)partials.size(); bamIdx++) {
        auto &in = partials[bamIdx];
        if (!in.good())
          throw;
        for (int batchIdx = 0; batchIdx < batch_size; batchIdx++) {
          int contigIdx = contigIdxStart + batchIdx;
          if (contigIdx >= header->n_targets)
            break;
          if (!in.good())
            throw;
          if (intraDepthVariance) {
            in.read((char *)(batchContigVariances.data() + bamIdx * batch_size +
                             batchIdx),
                    sizeof(VarianceType));
          } else {
            CountType depth;
            in.read((char *)(batchContigDepthsCounts.data() +
                             bamIdx * batch_size + batchIdx),
                    sizeof(CountType));
          }
          if (!in.eof() && !in.good())
            throw;
        }
      }

#pragma omp parallel for schedule(static, 1)
      for (int batchIdx = 0; batchIdx < batch_size; batchIdx++) {
        int contigIdx = contigIdxStart + batchIdx;
        if (contigIdx >= header->n_targets)
          continue;
        if (!contigLengthPass[contigIdx])
          continue;
        auto &batchContigDepth = batchContigDepths[batchIdx];
        for (int bamIdx = 0; bamIdx < (int)partials.size(); bamIdx++) {
          if (intraDepthVariance) {
            batchContigDepth.addDepth(
                averageReadSize[bamIdx],
                batchContigVariances[bamIdx * batch_size + batchIdx]);
          } else {
            batchContigDepth.addDepth(
                averageReadSize[bamIdx],
                batchContigDepthsCounts[bamIdx * batch_size + batchIdx]);
          }
        }
      }

      /*
      int bamIdx = 0;
      for (auto &in : partials) {
        if (!in.good())
          throw;
        for (int batchIdx = 0; batchIdx < batch_size; batchIdx++) {
          int contigIdx = contigIdxStart + batchIdx;
          if (contigIdx >= header->n_targets)
            break;
          if (!contigLengthPass[contigIdx])
            continue;
          auto &batchContigDepth = batchContigDepths[batchIdx];
          if (intraDepthVariance) {
            VarianceType variance;
            in.read((char *)&variance, sizeof(VarianceType));
            batchContigDepth.addDepth(averageReadSize[bamIdx], variance);
          } else {
            CountType depth;
            in.read((char *)&depth, sizeof(CountType));
            batchContigDepth.addDepth(averageReadSize[bamIdx], depth);
          }
          if (!in.eof() && !in.good())
            throw;
        }
        bamIdx++;
      }
        */

      // print the depth for this batch of contigs to the output stream
      std::stringstream ss;
      for (int batchIdx = 0; batchIdx < batch_size; batchIdx++) {
        int contigIdx = contigIdxStart + batchIdx;
        if (contigIdx >= header->n_targets)
          break;
        if (!contigLengthPass[contigIdx])
          continue;
        auto &batchContigDepth = batchContigDepths[batchIdx];
        contigDepthPass[contigIdx] = batchContigDepth.printDepth(
            ss, header->target_name[contigIdx], minContigDepth);
      }
      out << ss.str();
    }
  } else { // not checkpoint
    cerr << "Creating depth matrix file: " << outputTableFile << endl;

    printDepthTable(out, bamFilePaths, intraDepthVariance, bamContigDepths,
                    averageReadSize, percentIdentity, bamContigVariances,
                    header.get(), contigLengthPass, contigDepthPass,
                    minContigDepth, includeEdgeBases,
                    std::min(maxEdgeBases, minAvgRead / 3));
  }
  out.flush();
  of.close();

  bamContigDepths.clear();
  bamContigVariances.clear();

  if (checkpoint)
    for (auto &checkpoint_file : partialFiles)
      unlink(checkpoint_file.c_str());

  if (!referenceFastaFile.empty() && !unmappedFastqFile.empty()) {
    string shredFileName = unmappedFastqFile + "-contigShreds.fasta";
    SafeOfstream shredsOf(shredFileName.c_str());
    cerr << "Outputing shredded contigs to " << shredFileName << endl;
    KseqReader reference(referenceFastaFile);
    int32_t contigIdx = 0;
    while (reference.hasNext()) {
      if (contigLengthPass[contigIdx] && contigDepthPass[contigIdx]) {
        shredFasta(shredsOf, reference.getName(), reference.getSeq(),
                   shredLength, shredDepth, shredLength * 0.10);
      }
      contigIdx++;
    }
  }

  if (!outputGCFile.empty() && !refGCWindows.empty()) {
    for (int contigIdx = 0; contigIdx < (int)refGCWindows.size(); contigIdx++) {
      if (contigDepthPass[contigIdx]) {
        for (int j = 0; j < (int)refGCWindows[contigIdx].size(); j++) {
          refGC[refGCWindows[contigIdx][j]]++;
        }
      }
    }
    SafeOfstream outGC(outputGCFile.c_str());
    outGC << "GC\tRef\tReads\tCoverage\tMean\tVariance\n";
    for (int i = 0; i <= 100; i++) {
      outGC << i << "\t" << refGC[i] << "\t"
            << readGCStats[i].mean() * readGCStats[i].size() << "\t"
            << (refGC[i] > 0
                    ? readGCStats[i].mean() * readGCStats[i].size() / refGC[i]
                    : 0)
            << "\t" << readGCStats[i].mean() << "\t"
            << readGCStats[i].variance() << "\n";
    }
  }

  // output pairedContigs lowerTriangle
  if (!pairedContigsFile.empty() && hasAnyPairedContigs) {
    cerr << "Creating pairedContigs matrix file: " << pairedContigsFile << endl;

    SafeOfstream of(pairedContigsFile.c_str());
    printPairedContigs(of, numThreads, bamPairedContigs, header.get());
  }
  if (!pairedContigsFile.empty() && !hasAnyPairedContigs) {
    cerr << "The data files have no paired contigs to report" << endl;
  }

  cerr << "Closing last bam file" << endl;
  // close the last bam file
  bams[0].close();
  cerr << "Finished" << endl;

  return 0;
}
