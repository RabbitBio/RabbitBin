#ifndef RABBITBIN_H_
#define RABBITBIN_H_

#include <algorithm>
#include <cstdarg>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <deque>
#include <sys/time.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <limits>
#ifdef __APPLE__
#include <mach/mach.h>
#include <sys/sysctl.h>
#else
#include <sys/sysinfo.h>
#endif
#include <zlib.h>

#include "OpenMP.h"
#include "ProgressTracker.h"
#include "kseq.h"
#include "version.h"
KSEQ_INIT(gzFile, gzread)

#include "cuckoohash_map.hh"
#include "ranker.h"
#include "tile.h"

// force BOOST ublas optimizations
#define BOOST_UBLAS_INLINE inline
#define BOOST_UBLAS_CHECK_ENABLE 0
#define BOOST_UBLAS_USE_FAST_SAME
#define BOOST_UBLAS_TYPE_CHECK 0

#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/distributions.hpp>
#include <boost/program_options.hpp>
#include <boost/system/error_code.hpp>
#include <string>

#if (BOOST_VERSION / 100000 == 1) && (BOOST_VERSION / 100 % 1000 == 64)
#include <boost/serialization/array_wrapper.hpp>
#endif

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using std::cerr;
using std::cout;
using std::endl;
namespace po = boost::program_options;

// always use double when summing many values together to avoid precision issues

// for composition and depth_matrix calcs use float when SURPRISINGLY_NOT_FASTER is defined
//#define SURPRISINGLY_NOT_FASTER
#ifdef SURPRISINGLY_NOT_FASTER
typedef float CALC_TYPE;
#else
typedef double CALC_TYPE;
#endif

// use floats for most every value when NO_MEASUREABLE_DIFFERNCE_SMALLER_CALC is
// defined
//#define NO_MEASUREABLE_DIFFERNCE_SMALLER_CALC
#ifdef NO_MEASUREABLE_DIFFERNCE_SMALLER_CALC
typedef float Distance;
#else
typedef double Distance;
#endif

// use floats to store larger arrays and matrixes when SMALLER_STORED is defined
#define SMALLER_STORED
#ifdef SMALLER_STORED
typedef float StoredDistance;
#else
typedef double StoredDistance;
#endif

typedef Distance Similarity;
#define LOG log
#define LOG10 log10
#define SQRT sqrt
#define EXP exp
#define POW pow
#define FABS fabs

typedef std::pair<size_t, StoredDistance> Edge;

typedef boost::math::normal_distribution<StoredDistance> Normal;

static const std::string version = RabbitBin_VERSION;
static const std::string DATE = RabbitBin_BUILD_DATE;
static bool verbose = false;
static bool quiet = false;
static bool debug = false;
// Per-bin FASTA output is OFF by default: at scale it produces tens of
// thousands of files (a re-partition of the input) and forces all sequences to
// stay resident.  The binning result is always written as the membership tables
// (out.members.tsv / out.bins.tsv).  Pass --bin-fasta to also emit per-bin
// FASTA (which then keeps sequences in RAM).
static bool noBinOut = true;
static bool binFastaWanted = false;   // set by --bin-fasta → noBinOut = false
static bool noSampleDepths = false;
static Distance mergeSamplesCosign = 1.0;
static size_t min_bin_bp = 50000;
static size_t minContig = 2500; // minimum contig size for binning
static std::string inFile;
static std::string depth_file;
#ifdef RABBITBIN_FUSE
// Fused build only: compute depth in-process from sorted BAMs (no temp file) so
// BAM decompression overlaps the FASTA/sketch pass. These mirror the
// rabbit_depth CLI defaults so the in-memory depth is byte-identical.
static std::vector<std::string> fuse_bams;
static int fuse_pctid = 97;
static int fuse_min_contig_len = 1000;
static double fuse_min_contig_depth = 1.0;
static int fuse_max_edge = 75;
#endif
static bool cvExt;
static bool fullHeader = false;
static std::string outFile;
static bool onlyLabel = false;
static bool no_recruit = false;
// --no_gold: label-free multi-resolution mode. Builds the similarity graph once,
// runs label propagation under several α / edge_power settings, and selects the
// partition with the highest (gold-standard-free) modularity on the fixed
// composition graph. Lets RabbitBin auto-tune binning resolution without any
// ground truth. Default OFF — the single-resolution production path is unchanged.
static bool no_gold = false;
static size_t min_small_contig =
    1000;          // minimum contig size for small contig binning
size_t minCS = 10; // minimum cluster size for additional recruiting
static bool recruit_to_depth_centroid = false;
static Distance recruitSimFactor = 0.0;
static size_t numThreads = 0;
static Similarity calib_connected_pct = 95;
static Similarity min_edge_weight = 60;
static Similarity simCutoff = 0;
static Distance minCV = 1.0; /// TODO adjust to 0.1?
static Distance minCVSum = 1;
static bool saveCls = false;
static bool outUnbinned = false;
static size_t min_sample = 3;
static unsigned long long totalSize = 0, totalSize1 = 0;

static size_t maxEdges = 200;
static const char line_delim = '\n';
static const char tab_delim = '\t';
static const char fasta_delim = '>';
static const std::size_t buf_size = 1024 * 1024;

static std::vector<std::string> contig_names;
static std::vector<std::string> small_contig_names;
// Contig sequences are stored as VIEWS (std::string_view).  In the parallel
// mmap parse path a single-line sequence is referenced zero-copy directly in
// the (retained, never-unmapped) FASTA mmap — avoiding ~N per-contig heap
// allocations + a multi-GB memcpy that previously dominated the parse phase.
// Sequences that must be materialised (multi-line records, gzip / kseq
// fallbacks, the libdeflate streaming path) are copied once into g_seq_store
// (a deque, so element addresses are stable) and the view points there.
static std::deque<std::string>   g_seq_store;     // owned backing for gz/kseq paths
// Per-parse-thread contiguous arenas holding compacted MULTI-LINE sequences
// from the mmap path.  One big growing buffer per thread replaces N per-contig
// std::string allocations (which dominated the parse phase).  Kept alive for
// the whole run; seq views point into these.
static std::vector<std::string>  g_seq_arenas;
static std::vector<std::string_view> seqs;
static std::vector<std::string_view> small_seqs;
// Retained FASTA mmap (kept mapped for the whole run so zero-copy seq views
// stay valid).  Released by the OS at process exit.
static const char *g_fasta_mmap = nullptr;
static size_t      g_fasta_mmap_len = 0;
// Parallel length arrays – populated alongside seqs[]/small_seqs[] in every
// FASTA parse path so that size() queries don't need the full sequence in RAM.
static std::vector<size_t> seq_lens;
static std::vector<size_t> small_seq_lens;
static std::vector<StoredDistance> logSizes;

typedef std::vector<size_t> ContigVector;
typedef std::unordered_set<size_t> ContigSet;
typedef std::unordered_map<int, ContigVector> BinMap;

static size_t nobs = 0;  //# of large
static size_t nobs1 = 0; //# of small

typedef boost::numeric::ublas::matrix<StoredDistance> Matrix;
typedef boost::numeric::ublas::matrix_row<Matrix> MatrixRowType;

static Matrix depth_matrix;
static Matrix depth_var_matrix;
static Matrix depth_centroids;
static Matrix small_depth_matrix;

static size_t num_depth_samples = 0;
static unsigned long long seed = 0;

static size_t countLines(const char *f);
static size_t ncols(std::ifstream &is, int skip);
static size_t ncols(const char *f, int skip);

static double cal_depth_dist(size_t r1, size_t r2, size_t i, bool &nz);

static std::istream &safeGetline(std::istream &is, std::string &t);

static timeval t1, t2;

static void print_message(const char *format, ...) {
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stdout, format, argptr);
  cout.flush();
  va_end(argptr);
}

static void verbose_message(const char *format, ...) {
  if (verbose) {
    gettimeofday(&t2, NULL);
    int elapsed = (int)(((t2.tv_sec - t1.tv_sec) * 1000.0 +
                         (t2.tv_usec - t1.tv_usec) / 1000.0) /
                        1000.0); // seconds
    printf("[%02d:%02d:%02d] ", elapsed / 3600, (elapsed % 3600) / 60,
           elapsed % 60);
    va_list argptr;
    va_start(argptr, format);
    vfprintf(stdout, format, argptr);
    cout.flush();
    va_end(argptr);
  }
}

class Graph {
public:
  size_t n;
  std::vector<size_t> from;
  std::vector<size_t> to;
  std::vector<std::vector<size_t>>
      incs; // incidence list which has edge id instead of node id (compared to
            // adjacent list)
  std::vector<StoredDistance> sComp;
  std::vector<StoredDistance> edgeScore; // composite score (weight) of sComp and depth
  ContigSet connected_nodes;
  bool hasEdges;

  Graph(size_t num_nodes, bool hasEdges = false)
      : n(num_nodes), hasEdges(hasEdges) {
    if (hasEdges) {
      incs.resize(num_nodes);
    }
  }

  ~Graph() {}

  size_t getNodeCount() { return n; }

  size_t getEdgeCount() { return from.size(); }

  size_t getOtherNode(size_t e, size_t v) {
    assert(e < from.size() && e < to.size());
    return from[e] == v ? to[e] : from[e];
  }
};

void build_similarity_graph(Graph &g, Similarity cutoff);

static void trim_fasta_label(std::string &label) {
  size_t pos = label.find_first_of(" \t");
  if (pos != std::string::npos)
    label = label.substr(0, pos);
}

std::ostream &printFasta(std::ostream &os, string label, std::string_view seq) {
  int64_t len = seq.size();
  if (len == 0) {
    cerr << "Warning attempt to print an empty fasta!" << endl;
    return os;
  }
  os << fasta_delim << label << line_delim;
  const char *_seq = seq.data();
  const int maxWidth = 60;
  for (size_t s = 0; s < len; s += maxWidth) {
    int bytes = s + maxWidth < len ? maxWidth : len - s;
    os.write(_seq + s, bytes);
    os << line_delim;
  }
  return os;
}

// Fisher-Yates shuffle
// http://stackoverflow.com/questions/9345087/choose-m-elements-randomly-from-a-vector-containing-n-elements
template <class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random) {
  size_t left = std::distance(begin, end);
  while (num_random--) {
    bidiiter r = begin;
    std::advance(r, rand() % left);
    std::swap(*begin, *r);
    ++begin;
    --left;
  }
  return begin;
}

#ifdef __APPLE__
vm_statistics_data_t vmStats;
mach_msg_type_number_t infoCount = HOST_VM_INFO_COUNT;
#else
struct sysinfo memInfo;
#endif
double totalPhysMem = 0.;

int parseLine(char *line) {
  int i = strlen(line);
  while (*line < '0' || *line > '9')
    line++;
  line[i - 3] = '\0';
  i = atoi(line);
  return i;
}

double getTotalPhysMem() {
  if (totalPhysMem < 1) {
#ifdef __APPLE__
    kern_return_t kernReturn = host_statistics(
        mach_host_self(), HOST_VM_INFO, (host_info_t)&vmStats, &infoCount);
    if (kernReturn != KERN_SUCCESS)
      return 0;
    return (vm_page_size * (vmStats.wire_count + vmStats.active_count +
                            vmStats.inactive_count + vmStats.free_count)) /
           1024;
#else
    sysinfo(&memInfo);
    long long _totalPhysMem = memInfo.totalram;
    _totalPhysMem *= memInfo.mem_unit;
    totalPhysMem = (double)_totalPhysMem / 1024; // kb
#endif
  }
  return totalPhysMem;
}

// http://blog.csdn.net/hengshan/article/details/9201929
int getFreeMem() {
#ifdef __APPLE__
  kern_return_t kernReturn = host_statistics(mach_host_self(), HOST_VM_INFO,
                                             (host_info_t)&vmStats, &infoCount);
  if (kernReturn != KERN_SUCCESS)
    return 0;
  return (vm_page_size * vmStats.free_count) / 1024;
#else
  FILE *file = fopen("/proc/meminfo", "r");
  size_t result = 0;
  char line[128];

  while (fgets(line, 128, file) != NULL) {
    if (strncmp(line, "MemFree:", 6) == 0 ||
        strncmp(line, "Buffers:", 6) == 0 || strncmp(line, "Cached:", 6) == 0 ||
        strncmp(line, "SwapFree:", 6) == 0) {
      result += parseLine(line);
    }
  }
  fclose(file);
  return result; // Kb
#endif
}

double getUsedPhysMem() {
  return (getTotalPhysMem() - getFreeMem()) / 1024. / 1024.;
}

int cluster_by_propagation(Graph &g, std::vector<size_t> &membership,
                      std::vector<size_t> &node_order);

struct CompareEdge {
  // Total order so a top-maxEdges heap has a UNIQUE kept set / drain order
  // regardless of insertion order (required for run-to-run determinism under
  // multithreaded edge construction). Primary key: similarity (min-heap, so the
  // weakest kept edge is on top and is evicted first). Tie-break: neighbour id
  // — the heap top (evicted first) is the LARGER id, so on equal similarity the
  // smaller-id neighbour is retained deterministically.
  constexpr bool operator()(Edge const &a, Edge const &b) const noexcept {
    if (a.second != b.second) return a.second > b.second;
    return a.first < b.first;
  }
};

void promote_singleton_bins(BinMap &cls);
void output_bins(BinMap &cls);
size_t calibrate_sim_cutoff(Distance coverage = 1., bool full = false);
double cal_depth_corr(size_t r1, size_t r2, bool second_is_small = false,
                    bool first_is_centroid = false);
// bool is_small = false, bool is_centroid = false);

#endif
