/*
 * rabbit_depth.h  –  summarize per-contig coverage depth from BAM alignments
 */

#ifndef RABBIT_DEPTH_H_
#define RABBIT_DEPTH_H_

#include "version.h"

#include <getopt.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <memory>
#include <unordered_map>

#include "BamUtils.h"
#include "CheckRead.hpp"
#include "OpenMP.h"
#include "RunningStats.h"

using namespace std;

typedef BamUtils::StringVector StringVector;

typedef uint32_t ContigIdxType;
typedef uint64_t CountType;
typedef vector<std::shared_ptr<CountType[]>> CountTypeMatrix;
typedef map<ContigIdxType, CountType> PairedCountType;
typedef vector<PairedCountType> PairedCountTypeMatrix;

class ThreadBlocker {
public:
  ThreadBlocker() : _blocked(0) {}
  bool isBlocked() const { return _blocked == omp_get_num_threads(); }
  bool mayBeBlocked() const { return _blocked > omp_get_num_threads() / 2 + 1; }
  operator int() { return _blocked; }
  void setBlockedThread(bool iAmBlocked) {
    if (iAmBlocked) {
#pragma omp atomic
      _blocked++;
    } else {
#pragma omp atomic
      _blocked--;
    }
  }

private:
  int _blocked;
};

extern ThreadBlocker tb;

class DepthCounts {
public:
  typedef int32_t BaseCountType;
  typedef std::shared_ptr<BaseCountType[]> Mem;

private:
  int _allocLen, _targetLen;
  Mem _alloc;

public:
  BaseCountType *baseCounts;
  BaseCountType *mapQualities;
  DepthCounts()
      : _allocLen(0), _targetLen(0), _alloc(), baseCounts(NULL), mapQualities(NULL) {}
  ~DepthCounts() { clear(); }
  operator bool() { return _alloc && baseCounts != NULL && _allocLen > 0; }
  void add_depth(int refpos, int len) {
    assert(refpos + len <= _targetLen);
    for (int i = refpos; i < refpos + len; i++) {
      baseCounts[i]++;
    }
  }
  void add_map_qual(int refpos, int len, int mq) {
    assert(mapQualities && refpos + len <= _targetLen);
    for (int i = refpos; i < refpos + len; i++)
      mapQualities[i] += mq;
  }
  void resetBaseCounts(int targetLen, bool mapQualitiesToo) {
    targetLen++; // protect against off-by-one errors in reported alignment
    int requiredLen = targetLen;
    if (mapQualitiesToo)
      requiredLen *= 2;
    if (_allocLen < requiredLen) {
      // round to the the next power of 2
      long power = 2;
      while (targetLen > power)
        power <<= 1;
      requiredLen = power;
      if (mapQualitiesToo)
        requiredLen *= 2;

      tb.setBlockedThread(true);
      int attempts = 0;
      while (attempts < 3) {
        if (tb.mayBeBlocked() && targetLen < 1024 * 1024) {
          // more than half the threads are blocked... lets free our memory
          // proactively
          clear();
        }
        _alloc.reset(new BaseCountType[requiredLen]);
        if (_alloc) {
          baseCounts = _alloc.get();
          if (mapQualitiesToo) {
            mapQualities = baseCounts + requiredLen / 2;
          } else {
            mapQualities = NULL;
          }
          _allocLen = requiredLen;
          break;
        }
        if (tb.isBlocked()) {
          attempts++;
          std::cerr << "WARNING: Having trouble attempting to allocate memory "
                       "for per-base depth and variance on contig len: "
                    << targetLen << ". attempt: " << attempts << std::endl;
        }
        std::cerr << "WARNING: waiting for some threads to free enough memory "
                     "for this thread to continue..."
                  << omp_get_thread_num() << std::endl;
        sleep(1);
      }
      tb.setBlockedThread(false);
      if (!_alloc || _allocLen < requiredLen) {
        std::cerr << "ERROR: Could not allocate enough memory for even one "
                     "contig of length: "
                  << targetLen << std::endl;
        exit(1);
      }
    }
    _targetLen = targetLen;
    memset(baseCounts, 0, _targetLen * sizeof(BaseCountType));
    if (mapQualities) 
      memset(mapQualities, 0, _targetLen * sizeof(BaseCountType));
  }
  void clear() {
    _alloc.reset();
    _allocLen = 0;
    _targetLen = 0;
    baseCounts = NULL;
    mapQualities = NULL;
  }
  long getAllocLen() const { return _allocLen; }
};

class VarianceType {
public:
  typedef float T;
  T mean, variance;
  VarianceType() : mean(0), variance(0) {}
  void reset() {
    mean = 0.0;
    variance = 0.0;
  }
  ~VarianceType() {}
};

typedef vector<vector<VarianceType>> VarianceTypeMatrix;
typedef boost::iostreams::filtering_streambuf<boost::iostreams::output>
    gzipFileBuf;
typedef std::shared_ptr<gzipFileBuf> gzipFileBufPtr;

class ReadStatistics {
public:
  uint32_t exactmatches, substitutions, insertions, deletions, softclips,
      hardclips, skipRef, seqstart, seqend, seqlen, alignstart, alignend,
      alignlen, nm;
  bam1_t *bam;
  bool valid;
  ReadStatistics()
      : exactmatches(0), substitutions(0), insertions(0), deletions(0),
        softclips(0), hardclips(0), seqstart(0), seqend(0), seqlen(0),
        alignstart(0), alignend(0), alignlen(0), nm(0), bam(NULL),
        valid(false) {}
  ReadStatistics(bam1_t *b, uint32_t _exactmatches, uint32_t _substitutions,
                 uint32_t _insertions, uint32_t _deletions, uint32_t _softclips,
                 uint32_t _hardclips, uint32_t _skipRef, uint32_t _nm)
      : exactmatches(_exactmatches), substitutions(_substitutions),
        insertions(_insertions), deletions(_deletions), softclips(_softclips),
        hardclips(_hardclips), skipRef(_skipRef), seqstart(0),
        seqend(b->core.l_qseq), seqlen(b->core.l_qseq), alignstart(0),
        alignend(0), alignlen(0), nm(_nm), bam(b), valid(true) {
    alignlen = calculateAlignment(b, alignstart, alignend, seqstart, seqend);
    // seqlen remains the length of the sequence...
    const char *failMsg = NULL;
    if (nm != substitutions + insertions + deletions) {
      failMsg = "nm != sub + ins + del";
    }
    if (seqlen != 0 &&
        seqlen != exactmatches + substitutions + insertions + softclips) {
      failMsg =
          "seqlen != exactmatches + substitutions + insertions + softclips";
    }
    if (alignlen != alignend - alignstart) {
      failMsg = "alignlen != alignend - alignstart";
    }
    if (alignlen - skipRef != exactmatches + substitutions + deletions) {
      failMsg =
          "alignlen - skipRef != exactmatches + substitutions + deletions";
    }
    if (seqlen != 0 &&
        insertions - deletions + softclips != seqlen - (alignlen - skipRef)) {
      failMsg =
          "insertions - deletions + softclips != seqlen - (alignlen - skipRef)";
    }
    if (seqlen != 0 && seqlen < seqend - seqstart) {
      failMsg = "non-zero seqlen < seqend - seqstart";
    }
    if (alignlen == 0 && (b->core.flag & BAM_FUNMAP) != BAM_FUNMAP) {
      failMsg = "alignment length is zero";
    }
    if (failMsg != NULL) {
      std::stringstream ss;
      ss << failMsg << ": cigar mismatch alignlen=" << alignlen
         << " seqlen=" << seqlen << " skipRef=" << skipRef
         << " exactmatches=" << exactmatches
         << " substitutions=" << substitutions
         << " deletions=" << deletions << " insertions=" << insertions
         << " softclips=" << softclips << " hardclips=" << hardclips
         << " astart=" << alignstart << " aend=" << alignend
         << " seqstart=" << seqstart << " seqend=" << seqend
         << " nm=" << nm;
      warnMalformedBamRead(b, ss.str());
      valid = false;
    }
  }

  static uint32_t calculateAlignment(bam1_t *b, uint32_t &alignstart,
                                     uint32_t &alignend, uint32_t &seqstart,
                                     uint32_t &seqend) {
    uint32_t len = 0;
    alignstart = alignend = seqstart = seqend = 0;
    if (b == NULL)
      return len;
    if ((b->core.flag & BAM_FUNMAP))
      return len;
    if (b->core.pos < 0) {
      warnMalformedBamRead(b, "mapped read has a negative reference position");
      return len;
    }

    int numCigarOperations = b->core.n_cigar;
    if (numCigarOperations == 0) {
      warnMalformedBamRead(b, "mapped read has no CIGAR operations");
      return len;
    }

    alignstart = b->core.pos;
    uint32_t *cigar = bam1_cigar(b);
    alignend = bam_endpos(b);

    seqend = b->core.l_qseq;
    int32_t op = bam_cigar_op(cigar[0]);
    int32_t oplen = bam_cigar_oplen(cigar[0]);
    int32_t optype = bam_cigar_type(cigar[0]);
    // if the first operation moves along the reference, not the sequence, add
    // to the alignstart
    if (((optype & 0x01) != 0x01) && ((optype & 0x02) == 0x02) &&
        ((op & BAM_CDEL) != BAM_CDEL))
      alignstart += oplen;
    // if the first operation moves along the sequence, not the reference, add
    // to the seqstart
    if (((optype & 0x01) == 0x01) && ((optype & 0x02) != 0x02) &&
        ((op & BAM_CINS) != BAM_CINS))
      seqstart += oplen;

    if (numCigarOperations > 1) {
      oplen = bam_cigar_oplen(cigar[numCigarOperations - 1]);
      optype = bam_cigar_type(cigar[numCigarOperations - 1]);
      // assume bam_calend does the proper thing...
      // if the last operation moves along the sequence, not the reference,
      // subtract from the seqend
      if ((optype & 0x01) == 0x01 && (optype & 0x02) != 0x02) {
        if ((uint32_t)oplen > seqend) {
          warnMalformedBamRead(
              b, "last CIGAR operation consumes more query bases than the "
                 "query sequence contains");
          return 0;
        }
        seqend -= oplen;
      }
    }
    len = alignend - alignstart;
    if (len == 0) {
      std::stringstream ss;
      ss << "zero-length alignment alignstart=" << alignstart
         << " alignend=" << alignend;
      warnMalformedBamRead(b, ss.str());
      return 0;
    }
    if (!((seqend == 0 && seqstart == 0) || seqend > seqstart)) {
      std::stringstream ss;
      ss << "invalid query span seqstart=" << seqstart
         << " seqend=" << seqend;
      warnMalformedBamRead(b, ss.str());
      return 0;
    }
    return len;
  }

  bool isValid() const { return valid; }

  static std::ostream &writeHeader(std::ostream &os) {
    os << "ReadName\tReadLen\tAlignedLen\tPctId\tMappedTID\tMappedPos\tExactMat"
          "ches\tSubstitutions\tInsertions\tDeletions\tSoftClips\tHardClips\tPc"
          "tId2\tPctId3\tPctId4\tPctId5\tPctId6\tPctId7\tPctId8\tNM";
    return os;
  }

  ostream &write(ostream &os) const {
    if (bam == NULL || !valid || (bam->core.flag & BAM_FUNMAP) == BAM_FUNMAP)
      return os;

    os << bam1_qname(bam) << "\t" << seqlen << "\t" << alignlen << "\t"
       << getPctId() << "\t" << (int)bam->core.tid << "\t" << (int)bam->core.pos
       << "\t" << exactmatches << "\t" << substitutions << "\t" << insertions
       << "\t" << deletions << "\t" << softclips << "\t" << hardclips << "\t"
       << getPctId2() << "\t" << getPctId3() << "\t" << getPctId4() << "\t"
       << getPctId5() << "\t" << getPctId6() << "\t" << getPctId7() << "\t"
       << getPctId8() << "\t" << nm;
    return os;
  }

  float getPctId() const {
    uint32_t denom =
        exactmatches + substitutions + insertions + deletions;
    return denom > 0 ? ((float)exactmatches) / denom : 0.0f;
  }
  float getPctId2() const {
    uint32_t denom = exactmatches + substitutions;
    return denom > 0 ? ((float)exactmatches) / denom : 0.0f;
  }
  float getPctId3() const {
    uint32_t denom = exactmatches + substitutions + insertions;
    return denom > 0 ? ((float)exactmatches) / denom : 0.0f;
  }
  float getPctId4() const {
    uint32_t denom = exactmatches + substitutions + deletions;
    return denom > 0 ? ((float)exactmatches) / denom : 0.0f;
  }
  float getPctId5() const {
    uint32_t denom = alignlen > skipRef ? alignlen - skipRef : 0;
    return denom > 0 ? ((float)exactmatches) / denom : 0.0f;
  }
  float getPctId6() const {
    return seqlen > 0 ? ((float)exactmatches) / seqlen : 1;
  }
  float getPctId7() const {
    uint32_t alignedBases = alignlen > skipRef ? alignlen - skipRef : 0;
    uint32_t denom = alignedBases > seqlen ? alignedBases : seqlen;
    return denom > 0 ? ((float)exactmatches) / denom : 0.0f;
  }
  float getPctId8() const {
    uint32_t alignedBases = alignlen > skipRef ? alignlen - skipRef : 0;
    uint32_t denom = alignedBases < seqlen ? alignedBases : seqlen;
    return denom > 0 ? ((float)exactmatches) / denom : 0.0f;
  }
};

gzipFileBufPtr gzipOutputFile(string fileName, int level = 1) {
  boost::filesystem::fstream file(fileName.c_str(),
                                  ios_base::out | ios_base::binary);
  gzipFileBufPtr filebuf(new gzipFileBuf());
  filebuf->push(boost::iostreams::gzip_compressor());
  filebuf->push(file);
  return filebuf;
}

typedef std::vector<RunningStats> ReadGCStats;
void addGCCounts(ReadGCStats &readGCStats,
                 const std::vector<uint8_t> &refGCWindows, int gcWindow,
                 DepthCounts depthCounts) {

  if (!depthCounts || gcWindow <= 0 || refGCWindows.empty())
    return;
  if (depthCounts.getAllocLen() < (int)(gcWindow + refGCWindows.size())) {
    std::cerr << "WARNING: skipping GC coverage statistics because the "
                 "available per-base depth buffer is too small for gcWindow="
              << gcWindow << " and " << refGCWindows.size() << " windows"
              << std::endl;
    return;
  }
  for (int j = 0; j < (int)refGCWindows.size(); j++) {
    uint8_t gc = refGCWindows[j];
    long sum = 0;
    for (int k = 0; k < gcWindow; k++) {
      sum += depthCounts.baseCounts[j + k];
    }
#pragma omp critical(pushReadGCStats)
    readGCStats[gc].push((double)sum / (double)gcWindow);
  }
}

#ifndef MAX_KMERS
#define MAX_KMERS 256
#endif
class MappedKmersStats {
public:
  typedef std::vector<int64_t> KmerStats;
  typedef std::vector<RunningStats> KmerReadStats;

private:
  int64_t kmersMapped[MAX_KMERS], kmersCount[MAX_KMERS],
      kmersUnmappedCount[MAX_KMERS];
  RunningStats readKmers[MAX_KMERS];
  int32_t tmpRead[MAX_KMERS], tmpRead2[MAX_KMERS];

public:
  MappedKmersStats() {
    for (int i = 0; i < MAX_KMERS; i++) {
      kmersMapped[i] = kmersCount[i] = kmersUnmappedCount[i] = 0;
      readKmers[i] = RunningStats();
    }
    resetRead();
  }

  virtual ~MappedKmersStats() { resetRead(); }

  void resetRead() {
    for (int i = 0; i < MAX_KMERS; i++) {
      tmpRead[i] = 0;
    }
  }

  void addKmerStats(bam1_t *b, std::vector<string> &refs) {
    // reset the tmpRead data
    resetRead();
    int seqLen = b->core.l_qseq;
    for (int i = 0; i < MAX_KMERS; i++) {
      int count = seqLen - i;
      tmpRead2[i] = count >= 0 ? count : 0;
    }

    if ((b->core.flag & BAM_FUNMAP) != BAM_FUNMAP) {
      if (b->core.tid < 0 || b->core.tid >= (int32_t)refs.size()) {
        warnMalformedBamRead(
            b, "mapped read references a contig outside the available "
               "reference sequence list");
      } else {
      calcAlignedKmers(b, refs);
      }
    }

    // record
    for (int i = 0; i < MAX_KMERS; i++) {
      if (tmpRead2[i] > 0) {
        if ((b->core.flag & BAM_FUNMAP) != BAM_FUNMAP) {
          kmersMapped[i] += tmpRead[i];
          kmersCount[i] += tmpRead2[i];
          readKmers[i].push((double)tmpRead[i] / (double)tmpRead2[i]);
        } else {
          kmersUnmappedCount[i] += tmpRead2[i];
        }
      }
    }
  }

private:
  void addKmerMatches(int matchedSize) {
    for (int _k = 0; _k < matchedSize; _k++) {
      if (_k >= MAX_KMERS)
        break;
      int k = _k + 1;
      tmpRead[_k] += matchedSize - k + 1;
    }
  }
  void calcAlignedKmers(bam1_t *b, std::vector<string> &refs) {
    // populate this read matches
    if (b->core.tid < 0 || b->core.tid >= (int32_t)refs.size()) {
      warnMalformedBamRead(
          b, "cannot calculate aligned k-mers because the target id is out of "
             "range");
      return;
    }
    uint32_t *cigar = bam1_cigar(b), seqpos = 0, refpos = b->core.pos;
    int numCigarOperations = b->core.n_cigar;
    const char *ref = refs[b->core.tid].c_str();
    const size_t refLen = refs[b->core.tid].size();
    uint8_t *seq = bam1_seq(b);
    if (b->core.pos < 0 || (size_t)b->core.pos >= refLen) {
      warnMalformedBamRead(
          b, "cannot calculate aligned k-mers because the mapped position is "
             "outside the reference sequence");
      return;
    }

    for (int i = 0; i < numCigarOperations; i++) {
      int32_t op = bam_cigar_op(cigar[i]);
      int32_t oplen = bam_cigar_oplen(cigar[i]);
      int32_t optype = bam_cigar_type(cigar[i]);
      int32_t oppos = 0;
      int32_t matchlen = 0;
      switch (op) {
      case BAM_CMATCH: // M is match or mismatch!
        // find exact matches
        if (refpos + oplen > refLen || seqpos + oplen > (uint32_t)b->core.l_qseq) {
          warnMalformedBamRead(
              b, "cannot calculate aligned k-mers because a CIGAR match "
                 "extends beyond the reference or query sequence");
          return;
        }
        while (oppos < oplen) {
          if (ref[refpos + oppos] ==
              bam_nt16_rev_table[bam1_seqi(seq, seqpos + oppos)]) {
            matchlen++;
          } else {
            if (matchlen > 0)
              addKmerMatches(matchlen);
            matchlen = 0;
          }
          oppos++;
        }
        if (matchlen > 0)
          addKmerMatches(matchlen);
        break;
      case BAM_CEQUAL: // = is an exact match
        if (refpos + oplen > refLen || seqpos + oplen > (uint32_t)b->core.l_qseq) {
          warnMalformedBamRead(
              b, "cannot calculate aligned k-mers because an exact-match "
                 "CIGAR op extends beyond the reference or query sequence");
          return;
        }
        addKmerMatches(oplen);
        break;
      default:
        break;
      }
      // increment seqpos and/or refpos
      if ((optype & 0x01) == 0x01)
        seqpos += oplen;
      if ((optype & 0x02) == 0x02)
        refpos += oplen;
    }
  }

public:
  friend MappedKmersStats operator+(const MappedKmersStats a,
                                    const MappedKmersStats b) {
    MappedKmersStats combined;
    combined.resetRead();
    for (int i = 0; i < MAX_KMERS; i++) {
      combined.kmersMapped[i] = a.kmersMapped[i] + b.kmersMapped[i];
      combined.kmersCount[i] = a.kmersCount[i] + b.kmersCount[i];
      combined.kmersUnmappedCount[i] =
          a.kmersUnmappedCount[i] + b.kmersUnmappedCount[i];
      combined.readKmers[i] = a.readKmers[i] + b.readKmers[i];
    }
    return combined;
  }
  MappedKmersStats &operator+=(const MappedKmersStats &rhs) {
    MappedKmersStats combined = *this + rhs;
    *this = combined;
    return *this;
  }
  MappedKmersStats &operator=(const MappedKmersStats &copy) {
    this->resetRead();
    for (int i = 0; i < MAX_KMERS; i++) {
      this->kmersMapped[i] = copy.kmersMapped[i];
      this->kmersCount[i] = copy.kmersCount[i];
      this->kmersUnmappedCount[i] = copy.kmersUnmappedCount[i];
      this->readKmers[i] = copy.readKmers[i];
    }
    return *this;
  }

  ostream &writeHeader(ostream &os) {
    return os << "Kmer\tMappedReadCount\tUnmappedReadCount\tMappedKmers\tNumMap"
                 "pedReads\tMean\tVariance\tSkewness\tKurtosis\n";
  }

  ostream &write(ostream &os) {
    for (int _k = 0; _k < MAX_KMERS; _k++) {
      os << (_k + 1) << "\t" << kmersCount[_k] << "\t" << kmersUnmappedCount[_k]
         << "\t" << kmersMapped[_k] << "\t" << readKmers[_k].size() << "\t"
         << readKmers[_k].mean() << "\t" << readKmers[_k].variance() << "\t"
         << readKmers[_k].skewness() << "\t" << readKmers[_k].kurtosis()
         << "\n";
    }
    return os;
  }
};

CountType calcMismatches(bam1_t *b, uint32_t seqpos, const char *ref,
                         uint32_t refpos, uint32_t oplen) {
  uint32_t oppos = 0, mismatches = 0;
  while (oppos < oplen) {
    if (toupper(ref[refpos + oppos]) !=
        bam_nt16_rev_table[bam1_seqi(bam1_seq(b), seqpos + oppos)]) {
      mismatches++;
      // fprintf(stderr, "mismatch within %dM at %d: %c %c\n", oplen, oppos,
      // ref[refpos + oppos], bam_nt16_rev_table[ bam1_seqi(bam1_seq(b),
      // seqpos+oppos) ]);
    }
    oppos++;
  }
  return mismatches;
}

static int64_t WARNING_FLAG = 0;
CountType caldepth(bam1_t *b, DepthCounts depthCounts = DepthCounts(),
                   int32_t refLen = -1, const char *ref = NULL,
                   int ignoreEdges = 0, ReadStatistics *readstats = NULL) {
  // calculate the covered bases
  // optionally increment the coveredBases in the baseCounts vector
  int len = 0;
  uint32_t *cigar = bam1_cigar(b), insertions = 0, deletions = 0, seqpos = 0,
           refpos = b->core.pos, mismatches = 0, exactmatches = 0,
           softclips = 0, hardclips = 0, skipRef = 0;
  if (refpos >= 0 && refpos >= (uint32_t)refLen && refLen > 0) {
#pragma omp critical(OUTPUT)
    cerr << "WARNING: bam has improper refpos (" << refpos << ") vs reflen ("
         << refLen << ") ref=" << b->core.tid << ": " << bam1_qname(b) << endl;
    return 0;
  }
  int numCigarOperations = b->core.n_cigar;
  uint32_t cmatches = 0;
  bool isCEQUAL = false;
  for (int i = 0; i < numCigarOperations; i++) {
    int32_t op = bam_cigar_op(cigar[i]);
    int32_t oplen = bam_cigar_oplen(cigar[i]);
    int32_t optype = bam_cigar_type(cigar[i]);
    int mismatchesInOp = 0;
    switch (op) {
    case BAM_CMATCH: // M is match or mismatch!
      cmatches += oplen;
      if (ref != NULL && refLen > 0 && refpos + oplen > (uint32_t)refLen) {
        warnMalformedBamRead(
            b, "match CIGAR op extends beyond the reference sequence");
        return 0;
      }
      if (ref != NULL && b->core.l_qseq > 0 && b->core.l_qseq >= seqpos + oplen) {
        mismatchesInOp = calcMismatches(b, seqpos, ref, refpos, oplen);
        mismatches += mismatchesInOp;
      }
      // continue to CDIFF
    case BAM_CDIFF: // X is a mismatch.  Still count depth for variance purposes
      if (op == BAM_CDIFF) {
        isCEQUAL = true;
        mismatchesInOp = oplen;
        mismatches += mismatchesInOp;
      }
      // continue to CEQUAL
    case BAM_CEQUAL: // = is an exact match
      if (op == BAM_CEQUAL) {
        isCEQUAL = true;
        if (mismatchesInOp != 0) {
          warnMalformedBamRead(
              b, "exact-match CIGAR op reported mismatches during depth "
                 "calculation");
          mismatchesInOp = 0;
        }
      }
      // note: exact matches may be incorrect at this point if ref == NULL and
      // cmatches > 0
      exactmatches +=
          oplen - mismatchesInOp; // cmatches == exactmatches + mismatches
      len += oplen;
      if (depthCounts) {
        if (depthCounts.baseCounts == NULL ||
            depthCounts.getAllocLen() < (int)(refpos + oplen)) {
          warnMalformedBamRead(
              b, "cannot apply depth because the per-base buffer is too small "
                 "for this alignment");
          return 0;
        }
        if (refLen > 0 && oplen + refpos > (uint32_t)refLen) {
          warnMalformedBamRead(
              b, "cannot apply depth because the alignment extends beyond the "
                 "reference sequence");
          return 0;
        }
        depthCounts.add_depth(refpos, oplen);

        if (depthCounts.mapQualities != NULL) {
          depthCounts.add_map_qual(refpos, oplen, b->core.qual);
        }
      }
      if (ignoreEdges > 0 && refLen > 2 * ignoreEdges) {
        if ((int)refpos < ignoreEdges) {
          // refpos is left of first ignoreEdges boundary
          // subtract the lesser of the distance to the boundary or the
          // operation length
          len -= std::min((int32_t)(ignoreEdges - refpos), oplen);
        }
        if (refLen > 0 &&
            ((uint32_t)refpos + oplen) > ((uint32_t)refLen - ignoreEdges)) {
          if (refpos + oplen > (uint32_t)refLen) {
            warnMalformedBamRead(
                b, "cannot trim edge bases because the alignment extends "
                   "beyond the reference sequence");
            return 0;
          }
          len -= std::min((int32_t)(refpos + oplen) - (refLen - ignoreEdges),
                          oplen);
        }
      }
      break;
    case BAM_CINS: // I to reference
      insertions += oplen;
      break;
    case BAM_CDEL: // D from reference
      deletions += oplen;
      break;
    case BAM_CREF_SKIP: // N skip reference bases
      skipRef += oplen;
      break;
    case BAM_CSOFT_CLIP:
      softclips += oplen;
      break;
    case BAM_CHARD_CLIP:
      hardclips += oplen;
      break;
    default:
      break;
    }
    // increment seqpos and/or refpos
    if ((optype & 0x01) == 0x01)
      seqpos += oplen;
    if ((optype & 0x02) == 0x02)
      refpos += oplen;
  }
  // use NM - (insertion + deletion) errors to calculate mismatches
  uint8_t *NM = (uint8_t *)bam_aux_get(b, "NM");
  if (NM != NULL) {
    int32_t nm = bam_aux2i(NM);
    if (nm < (int32_t)(insertions + deletions)) {
      if (!(WARNING_FLAG & 0x1)) {
#pragma omp critical(OUTPUT)
        fprintf(stderr,
                "WARNING: your aligner reports an incorrect NM field.  You "
                "should run samtools calmd! nm < ins + del: cmatch=%d nm=%d ( "
                "insert=%d + del=%d + mismatch=%d == %d) %s\n",
                cmatches, nm, insertions, deletions, mismatches,
                insertions + deletions + mismatches, bam1_qname(b));
      }
      WARNING_FLAG |= 0x1;
      nm = insertions + deletions + mismatches;
    }
    if (isCEQUAL) { // exact match from sam v1.4 has been used (= and X)
      if (cmatches > 0) {
        // now 'M' means an ambiguous base mismatch because '=' exact match and
        // 'X' exact mismatches have been reported Count any 'M' match as a
        // match then assume the aligner would have set 'X' if the ambiguous
        // base did not include the read base
        if (ref != NULL) {
          // There was a reference, these bases would have been correctly
          // counted as either matches or mismatches
        } else {
          // these bases are already counted as an exactmatch and will need to
          // be corrected by the MD and/or NM fields
          if (nm != (int32_t)(insertions + deletions + mismatches)) {
            if (nm ==
                (int32_t)(insertions + deletions + mismatches + cmatches)) {
              // CMATCH 'M' is actually a mismatch here
              // fprintf(stderr, "FIXED: Ambiguous base in maping (cmatch=%d)
              // may be the source of NM (%d) discrepency ( insert=%d + del=%d +
              // mismatch=%d == %d). corrected mismatch count to %d: %s\n",
              // cmatches, nm, insertions, deletions, mismatches,
              // insertions+deletions+mismatches, nm - (insertions + deletions),
              // bam1_qname(b));
              mismatches += cmatches;
              exactmatches -= cmatches;
            } else {
              if (!(WARNING_FLAG & 0x2)) {
#pragma omp critical(OUTPUT)
                fprintf(
                    stderr,
                    "Warning: Ambiguous base in maping (cmatch=%d) may be the "
                    "source of NM (%d) discrepency ( insert=%d + del=%d + "
                    "mismatch=%d == %d). corrected mismatch count to %d: %s\n",
                    cmatches, nm, insertions, deletions, mismatches,
                    insertions + deletions + mismatches,
                    nm - (insertions + deletions), bam1_qname(b));
              }
              WARNING_FLAG |= 0x2;
              int newmismatches = nm - (insertions + deletions);
              int64_t adjustedExactMatches =
                  (int64_t)exactmatches -
                  ((int64_t)newmismatches - (int64_t)mismatches);
              if (adjustedExactMatches < 0) {
                warnMalformedBamRead(
                    b, "NM adjustment would make the exact-match count "
                       "negative; clamping to zero");
                adjustedExactMatches = 0;
              }
              exactmatches = adjustedExactMatches;
              mismatches = newmismatches;
            }
          }
        }
      }
    } else {
      if (ref == NULL) {
        // set mismatches based on the NM field since no reference is provided
        // and exact match vs mismatch is ambiguous in the CIGAR
        int newmismatches = nm - (insertions + deletions);
        if (mismatches != 0 && mismatches != (uint32_t)newmismatches) {
          warnMalformedBamRead(
              b, "CIGAR and NM disagree on the mismatch count; using the NM "
                 "field");
        }
        // fprintf(stderr, "Info: assigning %d (was %d) mismatches to %s\n",
        // newmismatches, mismatches, bam1_qname(b));
        int64_t adjustedExactMatches =
            (int64_t)exactmatches -
            ((int64_t)newmismatches - (int64_t)mismatches);
        if (adjustedExactMatches < 0) {
          warnMalformedBamRead(
              b, "NM adjustment would make the exact-match count negative; "
                 "clamping to zero");
          adjustedExactMatches = 0;
        }
        exactmatches = adjustedExactMatches;
        mismatches = newmismatches;

      } else {
      }
    }
    if (nm != (int32_t)(insertions + deletions + mismatches)) {
      if (!(WARNING_FLAG & 0x04)) {
#pragma omp critical(OUTPUT)
        fprintf(stderr,
                "Warning: consider running calmd on your bamfile! nm (%d) != "
                "insertions (%d) + deletions (%d) + mismatches (%d) (== %d) "
                "for read %s\n",
                nm, insertions, deletions, mismatches,
                insertions + deletions + mismatches, bam1_qname(b));
      }
      WARNING_FLAG |= 0x04;
      nm = insertions + deletions + mismatches;
    }
    int32_t subs = nm - (insertions + deletions);
    if (subs != (int32_t)mismatches) {
      warnMalformedBamRead(
          b, "NM-derived substitution count still disagrees with the running "
             "mismatch count; using the NM-derived value");
      mismatches = subs;
    }
    if (subs > len) {
      len = 0;
    } else if (subs > 0) {
      len -= subs;
    }
    if (readstats != NULL) {
      *readstats = ReadStatistics(b, exactmatches, mismatches, insertions,
                                  deletions, softclips, hardclips, skipRef, nm);
    }
  } else {
    // TODO if NM == NULL -- then calculate it (with a reference, of course)!
    std::cerr << "Warning: SAM 1.3 M not 1.4 =/X and there is no NM aux field. "
                 "PercentID will be invalid!!! "
              << bam1_qname(b) << endl;
    if (readstats != NULL) {
      *readstats = ReadStatistics(b, exactmatches, mismatches, insertions,
                                  deletions, softclips, hardclips, skipRef,
                                  mismatches + insertions + deletions);
    }
  }
  if (b->core.l_qseq > 0 && seqpos != (uint32_t)b->core.l_qseq) {
    cerr << "WARNING: bam has incorrect seqpos (" << seqpos
         << ") and queryLen (" << b->core.l_qseq << "): " << bam1_qname(b)
         << endl;
  }
  if (refLen > 0 && refpos > (uint32_t)refLen) {
    cerr << "WARNING: bam has incorrect refpos (" << refpos << ") and refLen ("
         << refLen << "): " << bam1_qname(b) << endl;
  }
  return len < 0 ? 0 : len;
}

std::vector<float> getAvgMapQuals(bam_header_t *header, uint32_t contigIdx,
                                  DepthCounts depthCounts,
                                  int32_t ignoreEdges) {
  assert(depthCounts);
  uint32_t start = 0, end = header->target_len[contigIdx];
  if ((int)end > 2 * ignoreEdges) {
    start = ignoreEdges;
    end = end - ignoreEdges;
  }
  std::vector<float> avgMapQuals;
  avgMapQuals.reserve(end - start + 1);
  for (uint32_t pos = start; pos != end; pos++) {
    float avgMapQual = -1.0;
    if (depthCounts.baseCounts[pos] > 0) {
      avgMapQual = (float)depthCounts.mapQualities[pos] /
                   (float)depthCounts.baseCounts[pos];
    }
    avgMapQuals.push_back(avgMapQual);
  }
  return avgMapQuals;
}

std::vector<float> avgMapQualsToWeights(const std::vector<float> &avgMapQuals,
                                        float weightMapQual,
                                        bool normalize = false) {
  // geometric mean (avg Q score)
  // sum( -10*log10(Pincorrect) )/n translated to
  // (product(Pcorrect)^(1/n))^weightMapQual
  double totalWeight = 0.0;
  std::vector<float> weights(avgMapQuals.size(), 1.0);
  assert(weights.size() == avgMapQuals.size());
  for (int i = 0; i < (int)avgMapQuals.size(); i++) {
    float qual = avgMapQuals[i];
    float &weight = weights[i];
    if (qual > 0.0) {
      // some confidence in the mapping
      float p_incorrect = pow(10.0, (0.0 - qual) / 10.0);
      assert(p_incorrect >= 0.0);
      assert(p_incorrect <= 1.0);
      weight = pow((1.0 - p_incorrect), weightMapQual);
      assert(weight >= 0.0);
      assert(weight <= 1.0);
    } else if (qual < 0.0) {
      // no mapped reads, depth 0 is unambiguous
      weight = 1.0;
    } else if (qual == 0.0) {
      // fully ambiguous placement
      weight = 0.0; // 0.0 is the weight we want!
    }
    totalWeight += weight;
  }

  if (normalize && totalWeight > 0.0) {
    float avgWeight = (float)totalWeight / weights.size();
    for (std::vector<float>::iterator itr = weights.begin();
         itr != weights.end(); itr++) {
      *itr = *itr / avgWeight;
    }
  }
  return weights;
}

VarianceType calculateVarianceContig(bam_header_t *header, uint32_t contigIdx,
                                     DepthCounts depthCounts,
                                     int32_t ignoreEdges = 0,
                                     float weightMapQual = 0.0,
                                     bool normalize = false) {
  // returns the variance along the contig
  VarianceType x{};
  if (!depthCounts)
    return x;

  int32_t contigLength = header->target_len[contigIdx];
  std::vector<float> weights;
  bool hasWeights = depthCounts.mapQualities != NULL && weightMapQual > 0.0;
  if (hasWeights) {
    std::vector<float> avgMapQuals =
        getAvgMapQuals(header, contigIdx, depthCounts, ignoreEdges);
    weights = avgMapQualsToWeights(avgMapQuals, weightMapQual, normalize);
  }

  VarianceType::T &avgDepth = x.mean;
  VarianceType::T &variance = x.variance;
  int32_t start = 0, end = contigLength;
  if (contigLength > 2 * ignoreEdges + 1) {
    start = ignoreEdges;
    end = end - ignoreEdges;
  }
  assert(end > start);
  int32_t adjustedContigLength = end - start;
  if (adjustedContigLength <= 2)
    return x;

  assert(depthCounts.getAllocLen() >= end);
  double totalWeights = 0.0;
  for (int32_t i = start; i != end; i++) {
    assert(!hasWeights || i < weights.size());
    avgDepth += depthCounts.baseCounts[i] * (hasWeights ? weights[i] : 1.0);
    totalWeights += (hasWeights ? weights[i] : 1.0);
  }
  assert(totalWeights <= adjustedContigLength);
  avgDepth = avgDepth / adjustedContigLength;
  for (int32_t i = start; i != end; i++) {
    VarianceType::T diff = depthCounts.baseCounts[i] - avgDepth;
    variance += diff * diff * (hasWeights ? weights[i] : 1.0);
  }
  variance = variance / (VarianceType::T)(totalWeights > 2.0
                                              ? (totalWeights - 1.0)
                                              : adjustedContigLength - 1.0);
  return x;
}

std::string getContigDepthByBase(bam_header_t *header, uint32_t contigIdx,
                                 DepthCounts depthCounts,
                                 int32_t ignoreEdges = 0,
                                 float weightMapQual = 0.0) {
  assert(depthCounts);
  std::stringstream ss;
  std::string nameonly = header->target_name[contigIdx];
  int pos = nameonly.find_first_of(" \t");
  nameonly = nameonly.substr(0, pos);
  ss << nameonly;
  uint32_t start = 0, end = header->target_len[contigIdx];

  if ((int)end > 2 * ignoreEdges) {
    start = ignoreEdges;
    end = end - ignoreEdges;
  }
  assert(start < end);
  assert(depthCounts.getAllocLen() >= end);
  for (uint32_t pos = start; pos != end; pos++) {
    ss << "\t" << (int)depthCounts.baseCounts[pos];
  }
  ss << "\n";

  if (weightMapQual > 0.0 && depthCounts.mapQualities != NULL) {
    ss << nameonly << "-avgMapQual";
    std::vector<float> avgMapQuals =
        getAvgMapQuals(header, contigIdx, depthCounts, ignoreEdges);
    for (std::vector<float>::iterator itr = avgMapQuals.begin();
         itr != avgMapQuals.end(); itr++) {
      ss << "\t" << (int)(*itr + 0.5); // round
    }
    ss << "\n";
    std::vector<float> weights =
        avgMapQualsToWeights(avgMapQuals, weightMapQual);
    ss << nameonly << "-mapWeights";
    for (std::vector<float>::iterator itr = weights.begin();
         itr != weights.end(); itr++) {
      ss << "\t" << *itr;
    }
    ss << "\n";
  }
  return ss.str();
}

int shredFasta(ostream &shredsOf, const string name, const string sequence,
               int shredLength, int shredDepth, int step, int offset = 33) {
  int seqLen = sequence.length();
  if (seqLen < shredLength)
    shredLength = seqLen;
  if (step >= seqLen)
    step = seqLen - 1;

  int numShreds = 0;
  for (int i = 0; i < seqLen - (shredLength - step); i += step) {
    numShreds++;
    for (int d = 0; d < shredDepth; d++) {
      int len = min(seqLen - i, shredLength);
      shredsOf << ">" << name << "x" << d << "-" << i << "+" << len << "\n"
               << sequence.substr(i, len) << "\n";
    }
  }
  return numShreds;
}

#endif // RABBIT_DEPTH_H_
