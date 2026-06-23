/**
 * rabbitbin.cpp
 *
 * RabbitBin: sketch-based metagenome binning.
 * Builds a composition–abundance similarity graph from contig sketches and
 * coverage profiles, then clusters with label propagation.
 *
 * Composition similarity uses k-mer sketches and/or weighted ProbMinHash
 * (embedded in rabbit_sketch.h / probmh.h; no external sketch library).
 *
 * Pipeline: FASTA + depth TSV → sketch → graph → label propagation →
 *           recruit small contigs → optional abundance split → output bins.
 *
 * Tune via RABBIT_* environment variables (see README).
 */

// ── K-mer sketch (self-contained, no external deps) ──────────────────────
#include "rabbit_sketch.h"
#include "rabbit_invidx.h"
#ifdef RABBITBIN_FUSE
// Fused build: in-process BAM→depth (defined in rabbit_depth.cpp, linked in).
#include "rabbit_depth_fuse.h"
#endif
// ── ProbMinHash4 weighted path (RabbitSketch), selectable via RABBIT_PMH=1 ──
#include "probmh.h"
#include <queue>
#include <cstring>
#include <atomic>
#include <limits>
#include <chrono>

// ── libdeflate (fast gzip decompression) ─────────────────────────────────
#include <libdeflate.h>
#include <future>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <omp.h>
// ── mmap parallel FASTA reader ───────────────────────────────────────────
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdio>   // sysconf(_SC_NPROCESSORS_ONLN) for true online-CPU count
#ifdef HAVE_LIBNUMA
#include <numa.h>
#include <numaif.h>   // mbind / MPOL_INTERLEAVE
#endif

#include "rabbitbin.h"
#include <iomanip>

// Migrate a large read-mostly buffer to MPOL_INTERLEAVE across all NUMA nodes
// so a multi-socket all-pairs scan spreads DRAM traffic over every memory
// controller instead of saturating the first-touch node.  Page-aligned interior
// only (edge pages stay put); no-op without libnuma, on single-node systems, or
// when RABBIT_NO_NUMA_INTERLEAVE is set.  Pure placement change — bit-identical
// results.
static inline void numa_interleave_buffer(void *ptr, size_t bytes) {
#ifdef HAVE_LIBNUMA
  if (getenv("RABBIT_NO_NUMA_INTERLEAVE")) return;
  if (numa_available() < 0) return;
  const int nnodes = numa_num_configured_nodes();
  if (nnodes < 2) return;
  const size_t pg = (size_t)sysconf(_SC_PAGESIZE);
  uintptr_t base = (uintptr_t)ptr;
  uintptr_t a = (base + pg - 1) & ~(uintptr_t)(pg - 1);   // round up to page
  uintptr_t b = (base + bytes) & ~(uintptr_t)(pg - 1);    // round down to page
  if (b <= a) return;
  unsigned long maxnode = (unsigned long)nnodes + 1;
  std::vector<unsigned long> mask((maxnode + 8 * sizeof(unsigned long) - 1) /
                                      (8 * sizeof(unsigned long)),
                                  0UL);
  for (int n = 0; n < nnodes; ++n)
    mask[n / (8 * sizeof(unsigned long))] |=
        (1UL << (n % (8 * sizeof(unsigned long))));
  // MPOL_MF_MOVE physically migrates already-faulted pages (the winner array is
  // built during the parse, so it is faulted by the time we get here).
  (void)mbind((void *)a, (size_t)(b - a), MPOL_INTERLEAVE, mask.data(),
              maxnode, MPOL_MF_MOVE);
#else
  (void)ptr; (void)bytes;
#endif
}

// ── Sketch global parameters ────────────────────────────────────────────────
static int      sketch_kmer_size  = 8;     // k-mer size for sketching (8 gives best coverage)
static uint32_t sketch_size = 500;  // PMH registers (--sketch-m; CAMI-high default)
static uint32_t sketch_bits        = 2;    // bits per OPH bucket (bit-planes)

// ── Marker-guided bin splitting (Phase 2; --marker-seed <file>) ─────────────
// Seed file (MetaDecoder format): one line per single-copy marker:
//   <marker>\t<contig>\t<contig>...
// A marker present on >=2 contigs in one output bin => that bin merged >=2
// genomes (contamination); such bins are re-split with abundance KMeans
// (k = min(marker multiplicity, splitMaxK)). Empty file => disabled.
static std::string marker_seed_file;
static int         splitMaxK        = 6;  // cap on sub-clusters per split bin
static int         splitMinContigs  = 6;  // min contigs in a bin to consider
// ── Marker-FREE Phase-2 split (--split-bins) ──────────────────────────────
// Splits internally multi-modal bins using only the multi-sample abundance
// already computed (no gene prediction / HMM / seed file).  A bin is split
// into k sub-clusters (k chosen by silhouette over k=2..splitMaxK) when the
// best silhouette >= g_split_sil.  Recovers the marker-guided split accuracy
// at zero extra cost (fits the ~4s binning budget).
static bool        g_split_abundance      = true;  // ON by default; --noAbdSplit disables
static bool        g_no_split_abundance   = false; // --noAbdSplit
static double      g_split_sil      = 0.70;  // silhouette threshold (RABBIT_SPLIT_SIL)
static size_t      g_sil_sample_cap = 600;   // sample cap for O(n^2) silhouette
// Content-conserving retention of sub-min_bin_bp bins (RABBIT_SPLIT_KEEP_COHERENT,
// default ON).  The split phase used to DISCARD every bin below min_bin_bp, which
// silently deletes recoverable genome content.  That is harmless on datasets whose
// genomes are large (CAMI1/2: such bins are rare) but devastating on communities
// with many fragmented / low-abundance genomes (CAMI3 toy: 5k+ bins, real genomes
// contributing <200kb of >=2500bp contig).  Instead of an absolute-bp discard, keep
// a small bin iff it is internally coherent — its mean intra-bin depth correlation
// is at least the MEDIAN coherence of the confidently-sized (>=min_bin_bp) bins.
// The bar is self-calibrated from the same run (no per-dataset constant), so the
// rule is identical across CAMI1/2/3.  Big bins are untouched (rec is monotone in
// added bins), so this cannot lower CAMI1/2 recovery.
static bool        g_split_keep_coherent = true;
// Percentile of confident-bin coherences used as the retention bar (0..100).
// A small bin is kept iff its coherence >= this percentile of the >=min_bin_bp
// bins' coherences.  Lower = keep more partial-genome bins.  Default 10: "at
// least as internally coherent as the bottom decile of confidently-sized bins".
static double      g_split_coh_pct  = 10.0;
// Conservative split guard (RABBIT_SPLIT_GUARD, default ON).  A large, internally
// coherent bin is usually ONE near-complete genome.  The plain silhouette split
// will happily cut it into two equally-coherent halves (silhouette only measures
// separability, not whether the pieces are different GENOMES), dropping both
// halves below the HQ completeness bar — a net loss of a near-complete MAG.  The
// guard vetoes a proposed split of a >= g_split_guard_bp bin UNLESS the split
// genuinely IMPROVES coherence: min(child_coh) >= parent_coh + g_split_guard_margin.
// Separating two co-binned genomes makes the children tighter than the mixed
// parent (passes); cutting one genome in half leaves children ~= parent (vetoed).
// Data-agnostic: the bar is the bin's OWN coherence, not a per-dataset constant;
// g_split_guard_bp is a genome-scale prior (bacterial genomes are Mb-scale).
//
// DEFAULT OFF.  Empirically (CAMI3 toy, 20 samples) the default abundance split
// is already strongly net-beneficial (no-split 221 HQ → split 233 HQ), and this
// guard HURTS (228 HQ): depth coherence cannot tell a beneficial STRAIN split
// from a harmful genome halving — co-binned strains share near-identical depth
// profiles, so the parent stays coherent and the (correct) strain split is
// wrongly vetoed.  The "shattered large genomes" that motivated this were a
// min-bin-size=42 artifact; in the default 50k config those genomes are already
// recovered.  Kept as an opt-in knob (RABBIT_SPLIT_GUARD=1) for experimentation;
// a useful version would need a COMPOSITION (k-mer) criterion, not depth.
static bool        g_split_guard       = false;
static size_t      g_split_guard_bp    = 1000000;  // only guard bins >= 1 Mb
static double      g_split_guard_margin = 0.05;     // required coherence gain
// Raw per-sample mean depths of small contigs, snapshotted BEFORE small_depth_matrix is
// rank-transformed in place during recruitment, so marker_guided_split sees the
// same coverage values as the large-contig depth_matrix (means, not ranks).
static std::vector<float> g_small_means;  // nobs1 × num_depth_samples (means)
static std::vector<float> g_large_means;  // nobs  × num_depth_samples (means, pre-rank)

// ── Abundance edge-similarity metric (RABBIT_DEPTH_SIM) ────────────────────
// 0 = corr  : Spearman (Pearson on rank-transformed depth)        [default]
// 1 = wjac  : weighted Jaccard  sum_s min(a_s,b_s)/sum_s max(a_s,b_s)
//             on RAW per-sample mean depths (magnitude-aware).  Robust at LOW
//             sample counts where the rank vector is too coarse for Spearman.
// g_depth_raw holds the raw per-sample means snapshotted BEFORE depth_matrix is
// rank-transformed in place; required by the wjac path.  g_depth_colnorm[i] =
// per-sample scaling (1/mean_i) so unequal library sizes don't dominate the
// weighted Jaccard (RABBIT_DEPTH_WJAC_NORM=1, default ON).
// Default = fuse: min(corr, wjac) conjunctive depth weight.  Universally >= corr
// on real benchmarks (CAMI-high +10~16, human_gut +22~28, marine/plant/
// strain_madness tie) at no measurable runtime cost.  Set RABBIT_DEPTH_SIM=corr
// to restore the legacy Spearman-only metric.
static int                g_depth_sim       = 2;   // 0=corr, 1=wjac, 2=fuse
// Fusion combiner (RABBIT_DEPTH_FUSE): how the pattern term (max(corr,0)) and the
// magnitude term (wjac) are combined into one depth weight. All are conjunctive
// ("edge supported only if BOTH shape and magnitude agree"), a dataset-agnostic
// principle: 1=min(c,w) [default, most universal], 0=geomean sqrt(c·w), 2=product.
// min is the most robust empirically: >= corr on every real benchmark tested.
static int                g_depth_fuse      = 1;
static bool               g_depth_wjac_norm = true;
// Dissimilarity hard-cut for the wjac path (RABBIT_DEPTH_WJAC_CUT, default 0=off):
// drop edges whose weighted Jaccard < cut.  Plays the role the negative-corr gate
// plays for the correlation path — needed to separate strains whose composition
// is near-identical (low-sample, low-composition regime).
static double             g_depth_wjac_cut  = 0.0;
// Hybrid (RABBIT_DEPTH_WJAC_GATE=1): use weighted Jaccard for the positive edge
// weight BUT keep correlation's signed negative gate (corr < g_neg_depth_thr →
// hard-cut).  Correlation's negative tail is what separates strains whose
// abundance is anti-correlated across samples; wjac's [0,1] range cannot express
// "definitely different", so on its own it over-merges strains.
static bool               g_depth_wjac_gate = false;
static std::vector<float> g_depth_raw;             // nobs × num_depth_samples (raw, pre-rank)
static std::vector<double> g_depth_colnorm;        // num_depth_samples (empty = no norm)
// Centered + L2-normalised rank vectors (nobs × num_depth_samples).  Pearson on
// ranks == dot product of these unit rows, so the edge-weight corr term becomes a
// single length-S FMA dot product instead of the per-pair Welford (3 divisions /
// sample + 2 sqrt / pair).  Built once after depth_matrix is rank-transformed.
static std::vector<float>  g_depth_unit;           // nobs × num_depth_samples (unit ranks)

// ── Cache-mode staging (Priority 2) ────────────────────────────────────────
// rb_load_cache() populates these graph arrays (the Graph object is not yet
// constructed at load time); the do-block in main() moves them into Graph g.
// g_cache_seed / g_cache_has_depth carry the build-time seed + depth flag so a
// --load-cache run reproduces the cached run unless the user overrides --seed.
static std::vector<size_t>         g_cache_from;
static std::vector<size_t>         g_cache_to;
static std::vector<StoredDistance> g_cache_scomp;
static unsigned long long          g_cache_seed = 0;
static bool                        g_cache_has_depth = false;

// Minimum similarity floor (×1000) used during auto-calibration.
// Sketch Jaccard on short metagenomic contigs typically spans 0–0.3.
static const size_t RB_SIM_FLOOR = 50;

// Per-contig KMV sketches (indexed 0..nobs-1)
static std::vector<rabbit_sketch::KmerSketch *> g_sketches;

// Per-cluster centroid sketches (indexed by cluster slot, built on demand)
static std::vector<rabbit_sketch::KmerSketch *> g_centroids;

// Contiguous signature flat array: g_sig_flat[i * g_sig_nw * g_sig_np ... ]
// Laid out as: sketch i, plane p, word w → g_sig_flat[i*(nw*np) + p*nw + w]
// Populated after all getSignature() calls. Enables cache-friendly sequential
// access in the all-pairs inner loop without pointer indirection.
static std::vector<uint64_t> g_sig_flat;
static uint32_t g_sig_nw = 0;  // nwords_ per sketch
static uint32_t g_sig_np = 0;  // bit-planes per sketch (b_)
static uint32_t g_sig_m  = 0;  // m_ (bucket count)

// ── Weighted ProbMinHash4 path (RABBIT_PMH=1) ─────────────────────────────
// When enabled, the large-contig similarity graph + sim-cutoff auto-calibration
// use a frequency-weighted ProbMinHash4 sketch (small k, k-mer counts as weights)
// instead of the OPH b-bit signature.  This estimates the weighted Jaccard
// Σmin(wA,wB)/Σmax(wA,wB) over the small-k composition spectrum.
//
// Similarity = (# registers with the same winning element) / m, read from a
// per-contig contiguous winners array (m × uint64).  The downstream centroid /
// leftover / small-contig steps (off by default) keep using KmerSketch.
static bool     g_pmh_mode = true;        // RABBIT_PMH (default on; RABBIT_PMH=0 disables)
static int      g_pmh_k    = 4;           // PMH k-mer size (env RABBIT_PMHK, default 4)
static uint32_t g_pmh_m    = 0;           // = sketch_size (registers)
static uint64_t g_pmh_seed = 42;
// Set to true when PMH winners were built on-the-fly during the kseq streaming
// pass (streaming producer-consumer mode).  Signals sketch build loop to skip
// the PMH phase and the sketch setup to skip g_win_flat re-allocation.
static bool     g_pmh_built_streaming = false;
static std::vector<uint32_t> g_win_flat;  // nobs × g_pmh_m winner identities (folded 64→32-bit, for SIMD sim kernel)
static std::vector<uint64_t> g_win64_flat; // nobs × g_pmh_m full 64-bit winners (for inverted-index key lookup)
// 16-bit packed winners (default kernel input): each 32-bit folded winner is
// XOR-folded to 16 bits.  Halves the (post-seq-drop dominant) winner array and
// lets the SIMD match kernel compare 32 lanes/AVX-512 op.  The extra 1/65536
// chance collision is absorbed by the baseline correction (re-estimated on the
// packed data).  g_win_bits selects which array the similarity kernel reads.
static std::vector<uint16_t> g_win16;
static int g_win_bits = 32;               // 32 = g_win_flat (uint32); 16 = g_win16
static inline uint16_t fold_win16(uint32_t w) {
  if (w == 0) return 0;                    // preserve the unset-winner sentinel
  uint16_t h = (uint16_t)(w ^ (w >> 16));
  return h ? h : (uint16_t)1;              // never alias the 0 sentinel
}

// Baseline-corrected winner-match estimator.  At small k all contigs share
// nearly the same k-mer universe, so unrelated contigs already collide on a
// large fraction b0 of registers (the "globally common" winners).  This packs
// every real composition difference into the narrow band [b0, 1], crushing the
// dynamic range that label propagation needs.  We estimate b0 empirically (the
// median winner-match over random — i.e. mostly unrelated — pairs) and rescale
//   s' = (s - b0) / (1 - b0)
// so the unrelated mode maps to ~0 and identical composition to ~1.
static bool   g_pmh_base_on  = true;      // env RABBIT_PMH_BASE (default 1)
static double g_pmh_baseline = 0.0;       // estimated chance-collision baseline b0
// Precomputed reciprocals to replace per-pair divisions in the O(N²) hot loop
// with multiplications.  g_inv_pmh_m = 1/m (m constant); g_inv_one_minus_b0 =
// 1/(1-b0) (b0 constant after baseline estimation).  Refreshed when either set.
static double g_inv_pmh_m        = 0.0;   // 1.0 / g_pmh_m
static double g_inv_one_minus_b0 = 1.0;   // 1.0 / (1.0 - g_pmh_baseline)

// Edge-weight blend of composition (sComp) vs abundance correlation:
//   w = alpha*sComp + (1-alpha)*corr     (alpha = RABBIT_W_COMP, default 0.5)
static double g_w_comp = 0.5;

// Graph precision (CAMI-high tuned defaults; override via env):
//   RABBIT_MUTUAL_KNN     default on (RABBIT_MUTUAL_KNN=0 disables)
//   RABBIT_NEG_DEPTH      default -0.3 (set <-0.99 to disable depth gate)
//   RABBIT_EDGE_POWER=<p> raise composite edge weight to power p before LPA
static bool   g_mutual_knn  = true;
static double g_neg_depth_thr = -0.3;
static double g_edge_power  = 1.0;     // applied to composite w; 1.0 = no-op
// GC-content normalization for PMH k-mer weights (RABBIT_GC_NORM=1):
//   Instead of using raw k-mer frequency as weight, use enrichment ratio:
//     weight(v) = observed_freq(v) / expected_freq(v | composition)
//
//   RABBIT_GC_NORM=1  per-base independence model: expected(ABCDEF) = P(A)P(B)P(C)P(D)P(E)P(F)
//                     Removes GC-content bias.
//   RABBIT_GC_NORM=2  first-order Markov (dinucleotide) model:
//                       expected(ABCDEF) = P(A)·P(B|A)·P(C|B)·P(D|C)·P(E|D)·P(F|E)
//                     P(b|a) estimated from the contig's own dinucleotide counts.
//                     Removes both GC and dinucleotide-context bias (CpG suppression etc.)
//   Enrichment ratios are capped at RABBIT_GC_NORM_CAP (default 20).
static int    g_gc_norm     = 1;    // RABBIT_GC_NORM: 0=off, 1=per-base, 2=dinucleotide
static double g_gc_norm_cap = 20.0;

// Env helpers: CAMI-high best (C_w050_mut). Unset = default; *=0 disables when noted.
static const char *rb_getenv(const char *key) {
  const char *e = getenv(key);
  return (e && e[0]) ? e : nullptr;
}

static bool rb_env_pmh_on() {
  const char *e = rb_getenv("RABBIT_PMH");
  return !e || e[0] != '0';
}
static bool rb_env_mutual_knn_on() {
  const char *e = rb_getenv("RABBIT_MUTUAL_KNN");
  return !e || e[0] != '0';
}
static int rb_env_gc_norm() {
  const char *e = rb_getenv("RABBIT_GC_NORM");
  return e ? std::atoi(e) : 1;
}
static double rb_env_neg_depth_thr() {
  const char *e = rb_getenv("RABBIT_NEG_DEPTH");
  return e ? std::atof(e) : -0.3;
}
static double rb_env_w_comp() {
  const char *e = rb_getenv("RABBIT_W_COMP");
  double v = e ? std::atof(e) : 0.5;
  return (v < 0.0) ? 0.0 : (v > 1.0 ? 1.0 : v);
}
static double rb_env_split_sil() {
  const char *e = rb_getenv("RABBIT_SPLIT_SIL");
  return e ? std::atof(e) : 0.70;
}
static bool rb_env_split_keep_coherent() {
  const char *e = rb_getenv("RABBIT_SPLIT_KEEP_COHERENT");
  return !e || e[0] != '0';
}
static double rb_env_split_coh_pct() {
  const char *e = rb_getenv("RABBIT_SPLIT_COH_PCT");
  if (!e) return 10.0;
  double v = std::atof(e);
  return v < 0.0 ? 0.0 : (v > 100.0 ? 100.0 : v);
}
static bool rb_env_split_guard() {
  const char *e = rb_getenv("RABBIT_SPLIT_GUARD");
  return e && e[0] != '0';            // default OFF; opt-in via RABBIT_SPLIT_GUARD=1
}
static size_t rb_env_split_guard_bp() {
  const char *e = rb_getenv("RABBIT_SPLIT_GUARD_BP");
  if (!e) return 1000000;
  long v = std::atol(e);
  return v < 0 ? 0 : (size_t)v;
}
static double rb_env_split_guard_margin() {
  const char *e = rb_getenv("RABBIT_SPLIT_GUARD_MARGIN");
  if (!e) return 0.05;
  return std::atof(e);
}
static bool rb_env_depth_prefilter_on() {
  const char *e = rb_getenv("RABBIT_DEPTH_PREFILTER");
  return !e || e[0] != '0';
}
static bool rb_env_depth_parse_var() {
  const char *e = rb_getenv("RABBIT_DEPTH_PARSE_VAR");
  return e && e[0] == '1';
}
// RABBIT_DEPTH_SIM = fuse (default) | corr | wjac
static int rb_env_depth_sim() {
  const char *e = rb_getenv("RABBIT_DEPTH_SIM");
  if (!e) return 2;                                          // default: fuse
  if (e[0] == 'c' || e[0] == 'C' || e[0] == '0') return 0;  // correlation only
  if (e[0] == 'w' || e[0] == 'W' || e[0] == '1') return 1;  // weighted Jaccard
  return 2;                                                  // fuse corr+wjac
}
static int rb_env_depth_fuse() {
  const char *e = rb_getenv("RABBIT_DEPTH_FUSE");
  if (!e) return 1;                                          // default: min
  if (e[0] == 'g' || e[0] == 'G' || e[0] == '0') return 0;  // geometric mean
  if (e[0] == 'p' || e[0] == 'P' || e[0] == '2') return 2;  // product
  return 1;                                                  // min
}
static bool rb_env_depth_wjac_norm_on() {
  const char *e = rb_getenv("RABBIT_DEPTH_WJAC_NORM");
  return !e || e[0] != '0';
}
static double rb_env_depth_wjac_cut() {
  const char *e = rb_getenv("RABBIT_DEPTH_WJAC_CUT");
  return e ? std::atof(e) : 0.0;
}
static bool rb_env_depth_wjac_gate() {
  const char *e = rb_getenv("RABBIT_DEPTH_WJAC_GATE");
  return e && e[0] == '1';
}

// Weighted-Jaccard abundance similarity on RAW per-sample depths (g_depth_raw),
// returns a value in [0,1] ready to blend as the depth term.  Defined here (top
// of TU) so the graph edge-weight loops below can call it; data is filled in the
// sketch loop before the in-place Spearman rank transform.
static inline double cal_depth_wjac(size_t r1, size_t r2) {
  if (g_depth_raw.empty() || num_depth_samples == 0) return 0.0;
  const float *a = g_depth_raw.data() + r1 * num_depth_samples;
  const float *b = g_depth_raw.data() + r2 * num_depth_samples;
  const bool norm = !g_depth_colnorm.empty();
  double mn = 0.0, mx = 0.0;
  for (size_t i = 0; i < num_depth_samples; ++i) {
    double x = a[i], y = b[i];
    if (norm) { const double s = g_depth_colnorm[i]; x *= s; y *= s; }
    if (x < y) { mn += x; mx += y; } else { mn += y; mx += x; }
  }
  if (mx <= 0.0) return 0.0;
  return mn / mx;
}

// Fast Pearson-on-ranks for two LARGE contigs via the precomputed unit vectors:
// corr(i,j) == Σ_k u_i[k]·u_j[k] (g_depth_unit holds centered+L2-normalised rank
// rows).  Pure FMA dot product — no per-element division, no per-pair sqrt — and
// auto-vectorises.  Falls back to the Welford cal_depth_corr when the unit cache
// is absent (e.g. centroid/small-contig variants build their own unit vectors).
static inline double depth_corr_fast(size_t i, size_t j) {
  if (g_depth_unit.empty()) return cal_depth_corr(i, j);
  const float *a = g_depth_unit.data() + i * num_depth_samples;
  const float *b = g_depth_unit.data() + j * num_depth_samples;
  const size_t S = num_depth_samples;
  size_t k = 0;
  double s = 0.0;
#if defined(__AVX512F__)
  // Two float FMA accumulators (16 lanes each) for ILP, then one horizontal
  // reduce — same pattern as k4_cosine_sim.  Unit rows are L2-normalised so the
  // dot lies in [-1,1]; float accumulation error (~S·1e-7) is far below the
  // edge-weight granularity.  S < 16 (e.g. CAMI-high, S=5) skips this entirely
  // and uses the scalar double tail below, unchanged.
  __m512 acc0 = _mm512_setzero_ps(), acc1 = _mm512_setzero_ps();
  for (; k + 32 <= S; k += 32) {
    acc0 = _mm512_fmadd_ps(_mm512_loadu_ps(a + k),      _mm512_loadu_ps(b + k),      acc0);
    acc1 = _mm512_fmadd_ps(_mm512_loadu_ps(a + k + 16), _mm512_loadu_ps(b + k + 16), acc1);
  }
  for (; k + 16 <= S; k += 16)
    acc0 = _mm512_fmadd_ps(_mm512_loadu_ps(a + k), _mm512_loadu_ps(b + k), acc0);
  if (k > 0) s = (double)_mm512_reduce_add_ps(_mm512_add_ps(acc0, acc1));
#elif defined(__AVX2__)
  // AVX2 fallback: two 8-lane accumulators (mul+add; no FMA dependency so it
  // works on plain AVX2 hosts).  S < 8 skips this and uses the scalar tail.
  __m256 acc0 = _mm256_setzero_ps(), acc1 = _mm256_setzero_ps();
  for (; k + 16 <= S; k += 16) {
    acc0 = _mm256_add_ps(acc0, _mm256_mul_ps(_mm256_loadu_ps(a + k),     _mm256_loadu_ps(b + k)));
    acc1 = _mm256_add_ps(acc1, _mm256_mul_ps(_mm256_loadu_ps(a + k + 8), _mm256_loadu_ps(b + k + 8)));
  }
  for (; k + 8 <= S; k += 8)
    acc0 = _mm256_add_ps(acc0, _mm256_mul_ps(_mm256_loadu_ps(a + k), _mm256_loadu_ps(b + k)));
  if (k > 0) {
    __m256 v8 = _mm256_add_ps(acc0, acc1);
    __m128 lo = _mm_add_ps(_mm256_castps256_ps128(v8), _mm256_extractf128_ps(v8, 1));
    lo = _mm_add_ps(lo, _mm_movehl_ps(lo, lo));
    lo = _mm_add_ss(lo, _mm_shuffle_ps(lo, lo, 0x1));
    s = (double)_mm_cvtss_f32(lo);
  }
#endif
  for (; k < S; ++k) s += (double)a[k] * (double)b[k];
  return s;
}

// Vectorised float dot of two unit (centered+L2-normalised) rows of length S.
// Used by the abundance-first prune in the O(N²/2) graph pass, evaluated for
// EVERY candidate pair, so vectorising it is the single biggest pair-pass win on
// many-sample data.  float accumulation differs from the scalar reference by
// ~S·1e-7, far inside the prune's 1e-4 conservative safety margin, so pruning
// decisions (and thus the graph) are unchanged.
static inline float unit_dot_f(const float *__restrict__ a,
                               const float *__restrict__ b, uint32_t S) {
  uint32_t k = 0;
  float r = 0.0f;
#if defined(__AVX512F__)
  __m512 acc0 = _mm512_setzero_ps(), acc1 = _mm512_setzero_ps();
  for (; k + 32 <= S; k += 32) {
    acc0 = _mm512_fmadd_ps(_mm512_loadu_ps(a + k),      _mm512_loadu_ps(b + k),      acc0);
    acc1 = _mm512_fmadd_ps(_mm512_loadu_ps(a + k + 16), _mm512_loadu_ps(b + k + 16), acc1);
  }
  for (; k + 16 <= S; k += 16)
    acc0 = _mm512_fmadd_ps(_mm512_loadu_ps(a + k), _mm512_loadu_ps(b + k), acc0);
  if (k > 0) r = _mm512_reduce_add_ps(_mm512_add_ps(acc0, acc1));
#elif defined(__AVX2__)
  __m256 acc0 = _mm256_setzero_ps(), acc1 = _mm256_setzero_ps();
  for (; k + 16 <= S; k += 16) {
    acc0 = _mm256_add_ps(acc0, _mm256_mul_ps(_mm256_loadu_ps(a + k),     _mm256_loadu_ps(b + k)));
    acc1 = _mm256_add_ps(acc1, _mm256_mul_ps(_mm256_loadu_ps(a + k + 8), _mm256_loadu_ps(b + k + 8)));
  }
  for (; k + 8 <= S; k += 8)
    acc0 = _mm256_add_ps(acc0, _mm256_mul_ps(_mm256_loadu_ps(a + k), _mm256_loadu_ps(b + k)));
  if (k > 0) {
    __m256 v8 = _mm256_add_ps(acc0, acc1);
    __m128 lo = _mm_add_ps(_mm256_castps256_ps128(v8), _mm256_extractf128_ps(v8, 1));
    lo = _mm_add_ps(lo, _mm_movehl_ps(lo, lo));
    lo = _mm_add_ss(lo, _mm_shuffle_ps(lo, lo, 0x1));
    r = _mm_cvtss_f32(lo);
  }
#endif
  for (; k < S; ++k) r += a[k] * b[k];
  return r;
}

// Depth term for edge blending, dispatching on g_depth_sim.  Sets ok=false to
// hard-cut the edge (non-finite, or strongly negative correlation).
static inline double depth_edge_term(size_t i, size_t j, bool &ok) {
  ok = true;
  if (g_depth_sim == 1) {                        // weighted Jaccard (magnitude)
    double wj = cal_depth_wjac(i, j);
    if (!std::isfinite(wj)) { ok = false; return 0.0; }
    if (g_depth_wjac_cut > 0.0 && wj < g_depth_wjac_cut) { ok = false; return 0.0; }
    // Hybrid: keep correlation's signed negative gate for strain separation.
    if (g_depth_wjac_gate && num_depth_samples > 1 && g_neg_depth_thr > -0.99) {
      double corr = depth_corr_fast(i, j);
      if (std::isfinite(corr) && corr < g_neg_depth_thr) { ok = false; return 0.0; }
    }
    return wj < 0.0 ? 0.0 : (wj > 1.0 ? 1.0 : wj);
  }
  double corr = depth_corr_fast(i, j);           // Spearman (default)
  if (!std::isfinite(corr)) { ok = false; return 0.0; }
  if (g_neg_depth_thr > -0.99 && corr < g_neg_depth_thr) { ok = false; return 0.0; }
  double cp = corr < 0.0 ? 0.0 : corr;
  if (g_depth_sim == 2 && num_depth_samples > 1) {
    // Fusion: require BOTH shape (corr) and magnitude (wjac) to agree.  Keeps
    // corr's signed gate (above) AND its specificity, while letting magnitude
    // refine same-organism edges.  Conjunctive → cross-organism pairs whose
    // magnitude coincides but pattern doesn't are suppressed.
    double wj = cal_depth_wjac(i, j);
    if (!std::isfinite(wj)) wj = 0.0;
    wj = wj < 0.0 ? 0.0 : (wj > 1.0 ? 1.0 : wj);
    switch (g_depth_fuse) {
      case 1:  return cp < wj ? cp : wj;          // min
      case 2:  return cp * wj;                    // product
      default: return std::sqrt(cp * wj);         // geometric mean
    }
  }
  return cp;
}

// IDF (inverse document frequency) normalization for PMH k-mer weights.
// Controlled by RABBIT_IDF_NORM=1 (requires g_gc_norm>0 and k==4).
//
// After building per-contig GC-normalized frequency vectors, computes a
// global IDF weight for each canonical 4-mer:
//   IDF(v) = log2((N+1) / (df(v)+0.5))
// where df(v) = number of contigs in which canonical k-mer v is present.
// Final weight = GC_norm_weight(v) × IDF(v), down-weighting ubiquitous
// k-mers and amplifying organism-specific ones — similar to TF-IDF in NLP.
// Stored exact k=4 frequencies (size nobs×256) are freed after rebuild.
static bool              g_idf_norm = false;
static std::vector<float> g_k4freq_flat; // nobs × 256, non-canonical entries=0

// Exact k=4 cosine similarity mode (RABBIT_EXACT_COS=1, requires g_gc_norm>0 and k==4).
//
// Replaces PMH Jaccard with exact dot-product of L2-normalised GC-norm k=4
// frequency vectors (size nobs×256, non-canonical entries=0).  More accurate
// than sketch-based approximation on normalised k-mer frequency profiles.
// Storage: ~30 MB for N=30k contigs.  Comparison cost: 256 FMAs per pair.
static bool              g_exact_cos_cmp = false;
static std::vector<float> g_k4cosine_flat; // nobs × 256, L2-normalised GC-norm freqs

// Fast Jaccard on raw pointers (no virtual dispatch, no null checks).
// a, b: pointers to the start of two signature rows in g_sig_flat.
// Layout: plane p, word w  →  a[p*nw + w].
#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)
#include <immintrin.h>
// AVX-512 VPOPCNTDQ path: process 8 × uint64 per cycle.
// Handles the common case np∈{1,2} explicitly (unrolled inner-p loop).
static inline double jaccard_raw(const uint64_t* __restrict__ a,
                                  const uint64_t* __restrict__ b,
                                  uint32_t nw, uint32_t np, uint32_t m) {
    const uint32_t rem  = m & 63;
    const uint64_t lastmask = (rem == 0) ? ~0ULL : ((1ULL << rem) - 1);
    const uint32_t full = (rem == 0) ? nw : (nw - 1);  // words without tail
    __m512i vsum = _mm512_setzero_si512();

    uint32_t w = 0;
    if (np == 1) {
        for (; w + 8 <= full; w += 8) {
            __m512i va = _mm512_loadu_si512((const __m512i*)(a + w));
            __m512i vb = _mm512_loadu_si512((const __m512i*)(b + w));
            __m512i vmatch = _mm512_ternarylogic_epi64(va, vb, _mm512_set1_epi64(-1LL), 0x09);
            // ternary: A XNOR B = ~(A^B): LUT value for (NOT (A XOR B)) = 0x09 → wrong
            // Correct: XNOR = ~(A XOR B). Use: 0x09? No.
            // A XOR B = 0x66, NOT = ~0x66 = 0x99. So ternarylogic with C=all-ones: result[i] = ~(a[i]^b[i])
            // ternarylogic(A,B,C,0x09) = A AND (NOT B): wrong.
            // Correct ternarylogic for ~(A^B): A XNOR B = (A AND B) OR (NOT A AND NOT B)
            //   truth table: A=0,B=0 → 1; A=0,B=1 → 0; A=1,B=0 → 0; A=1,B=1 → 1 → 1001 = 0x09? No.
            // Actually ternarylogic indexes the truth table where bit0 = A=0,B=0,C=0;
            // For XNOR(A,B): (A==B) = 1 when both same.
            // truth table for XNOR(A,B): ABC: 000→1, 001→1, 010→0, 011→0, 100→0, 101→0, 110→1, 111→1
            // = binary 11000011 = 0xC3 (but indexed from LSB: bit0 for ABC=000, etc.)
            // Hmm, ternarylogic actually only uses 2 inputs for XNOR; C is ignored here.
            // _mm512_ternarylogic_epi64(A,B,C,imm8): bit i of result = imm8[4*(A_bit_i)+2*(B_bit_i)+(C_bit_i)]
            // For XNOR(A,B) ignoring C: bit = imm8[4*a + 2*b + c]
            // a=0,b=0: both c=0,1 → 1,1 → bits 0,1 = 11
            // a=0,b=1: both c=0,1 → 0,0 → bits 2,3 = 00
            // a=1,b=0: both c=0,1 → 0,0 → bits 4,5 = 00
            // a=1,b=1: both c=0,1 → 1,1 → bits 6,7 = 11
            // = 0b11000011 = 0xC3
            vmatch = _mm512_ternarylogic_epi64(va, vb, va, 0xC3);  // ~(A^B)
            vsum = _mm512_add_epi64(vsum, _mm512_popcnt_epi64(vmatch));
        }
    } else if (np == 2) {
        const uint64_t* a1 = a + nw;
        const uint64_t* b1 = b + nw;
        for (; w + 8 <= full; w += 8) {
            __m512i va0 = _mm512_loadu_si512((const __m512i*)(a  + w));
            __m512i vb0 = _mm512_loadu_si512((const __m512i*)(b  + w));
            __m512i va1v= _mm512_loadu_si512((const __m512i*)(a1 + w));
            __m512i vb1v= _mm512_loadu_si512((const __m512i*)(b1 + w));
            // match = ~(a0^b0) & ~(a1^b1)  (bucket matches iff both planes agree)
            __m512i xn0 = _mm512_ternarylogic_epi64(va0, vb0, va0, 0xC3);
            __m512i xn1 = _mm512_ternarylogic_epi64(va1v, vb1v, va1v, 0xC3);
            __m512i vmatch = _mm512_and_si512(xn0, xn1);
            vsum = _mm512_add_epi64(vsum, _mm512_popcnt_epi64(vmatch));
        }
    } else {
        for (; w + 8 <= full; w += 8) {
            __m512i vmatch = _mm512_set1_epi64(-1LL);
            for (uint32_t p = 0; p < np; ++p) {
                __m512i va = _mm512_loadu_si512((const __m512i*)(a + p*nw + w));
                __m512i vb = _mm512_loadu_si512((const __m512i*)(b + p*nw + w));
                vmatch = _mm512_and_si512(vmatch,
                            _mm512_ternarylogic_epi64(va, vb, va, 0xC3));
            }
            vsum = _mm512_add_epi64(vsum, _mm512_popcnt_epi64(vmatch));
        }
    }
    // Horizontal sum of vsum
    uint64_t matches = (uint64_t)_mm512_reduce_add_epi64(vsum);

    // Scalar tail (remaining full words)
    for (; w < full; ++w) {
        uint64_t msk = ~0ULL;
        for (uint32_t p = 0; p < np; ++p)
            msk &= ~(a[p*nw + w] ^ b[p*nw + w]);
        matches += (uint64_t)__builtin_popcountll(msk);
    }
    // Last partial word (if any)
    if (rem != 0) {
        uint64_t msk = ~0ULL;
        for (uint32_t p = 0; p < np; ++p)
            msk &= ~(a[p*nw + full] ^ b[p*nw + full]);
        matches += (uint64_t)__builtin_popcountll(msk & lastmask);
    }
    return (double)matches / (double)m;
}
#else
// Generic scalar fallback (also compiled on non-AVX-512 hosts).
static inline double jaccard_raw(const uint64_t* __restrict__ a,
                                  const uint64_t* __restrict__ b,
                                  uint32_t nw, uint32_t np, uint32_t m) {
    const uint32_t full = (m % 64 == 0) ? nw : (nw - 1);
    const uint32_t rem  = m & 63;
    const uint64_t lastmask = (rem == 0) ? ~0ULL : ((1ULL << rem) - 1);
    uint64_t matches = 0;
    for (uint32_t w = 0; w < nw; ++w) {
        uint64_t msk = ~0ULL;
        for (uint32_t p = 0; p < np; ++p) {
            const size_t off = (size_t)p * nw + w;
            msk &= ~(a[off] ^ b[off]);
        }
        if (w >= full) msk &= lastmask;
        matches += (uint64_t)__builtin_popcountll(msk);
    }
    return (double)matches / (double)m;
}
#endif

// ── Weighted ProbMinHash4 similarity (winner-identity positional match) ───
// Two registers "collide" iff the SAME element won register k in both sketches
// (winners encode element identity); collision fraction = weighted Jaccard.
// Winners are folded to 32-bit identities (collision prob ≈ 2^-32 per register,
// negligible bias) so the kernel scans half the bytes and vectorises 16-wide.
// Raw winner-match fraction (uncorrected) in [0, 1].
// Returns the raw number of matching non-zero winner registers (numerator of
// the winner-match fraction).  Callers multiply by a precomputed 1/m to avoid
// a per-pair division in the hot loop.
// Env-gated coarse phase timer: prints wall ms since first call to stderr.
static inline void rb_phase(const char *label) {
  static const bool on = (getenv("RB_TIMING") != nullptr);
  if (!on) return;
  static std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  double ms = std::chrono::duration<double, std::milli>(
                  std::chrono::steady_clock::now() - t0).count();
  long rss_mb = 0, hwm_mb = 0;
  if (FILE *f = fopen("/proc/self/status", "r")) {
    char ln[256];
    while (fgets(ln, sizeof(ln), f)) {
      long kb;
      if (sscanf(ln, "VmRSS: %ld kB", &kb) == 1) rss_mb = kb / 1024;
      else if (sscanf(ln, "VmHWM: %ld kB", &kb) == 1) hwm_mb = kb / 1024;
    }
    fclose(f);
  }
  fprintf(stderr, "[RB_TIMING] %-28s %8.1f ms   RSS=%ldMB peak=%ldMB\n",
          label, ms, rss_mb, hwm_mb);
}

// Env-gated sequence-integrity self-check: order-independent 64-bit hash over
// every (contig name + stored sequence bytes).  Used to PROVE that a change to
// the sequence-storage path keeps the loaded bytes identical (the partition is
// a deterministic function of these bytes).  Prints to stderr when RB_SEQHASH.
static inline void rb_seq_integrity_check(const char *label) {
  if (getenv("RB_SEQHASH") == nullptr) return;
  auto fnv = [](const char *p, size_t n, uint64_t h) {
    for (size_t i = 0; i < n; ++i) { h ^= (uint8_t)p[i]; h *= 1099511628211ULL; }
    return h;
  };
  uint64_t acc = 0;
  uint64_t totbp = 0;
  for (size_t i = 0; i < contig_names.size() && i < seqs.size(); ++i) {
    uint64_t h = 1469598103934665603ULL;
    h = fnv(contig_names[i].data(), contig_names[i].size(), h);
    h = fnv(seqs[i].data(), seqs[i].size(), h);
    acc ^= h; totbp += seqs[i].size();
  }
  for (size_t i = 0; i < small_contig_names.size() && i < small_seqs.size(); ++i) {
    uint64_t h = 1469598103934665603ULL;
    h = fnv(small_contig_names[i].data(), small_contig_names[i].size(), h);
    h = fnv(small_seqs[i].data(), small_seqs[i].size(), h);
    acc ^= h; totbp += small_seqs[i].size();
  }
  fprintf(stderr, "[RB_SEQHASH] %-12s seqhash=%016llx large=%zu small=%zu totbp=%llu\n",
          label, (unsigned long long)acc, seqs.size(), small_seqs.size(),
          (unsigned long long)totbp);
}

static inline uint32_t pmh_match_count(const uint32_t *__restrict__ a,
                                       const uint32_t *__restrict__ b, uint32_t m) {
  uint32_t cnt = 0;
  uint32_t k = 0;
#if defined(__AVX512F__)
  // Count lanes where a==b AND a!=0 (a==0 is the unset-winner sentinel).
  // Four independent vector accumulators: each iteration produces a compare
  // mask and folds it into the accumulator with a single masked vpaddd, so the
  // per-step kmask→GPR(KMOVW)→POPCNT→add dependency chain of the scalar-popcount
  // variant is removed.  Counts stay in vector lanes until one horizontal
  // reduction at the end.  Result is the identical integer match count.
  const __m512i vzero = _mm512_setzero_si512();
  const __m512i vone  = _mm512_set1_epi32(1);
  __m512i acc0 = vzero, acc1 = vzero, acc2 = vzero, acc3 = vzero;
  for (; k + 64 <= m; k += 64) {
    __m512i va0 = _mm512_loadu_si512((const __m512i *)(a + k));
    __m512i vb0 = _mm512_loadu_si512((const __m512i *)(b + k));
    __m512i va1 = _mm512_loadu_si512((const __m512i *)(a + k + 16));
    __m512i vb1 = _mm512_loadu_si512((const __m512i *)(b + k + 16));
    __m512i va2 = _mm512_loadu_si512((const __m512i *)(a + k + 32));
    __m512i vb2 = _mm512_loadu_si512((const __m512i *)(b + k + 32));
    __m512i va3 = _mm512_loadu_si512((const __m512i *)(a + k + 48));
    __m512i vb3 = _mm512_loadu_si512((const __m512i *)(b + k + 48));
    __mmask16 h0 = _mm512_cmpeq_epi32_mask(va0, vb0) & _mm512_cmpneq_epi32_mask(va0, vzero);
    __mmask16 h1 = _mm512_cmpeq_epi32_mask(va1, vb1) & _mm512_cmpneq_epi32_mask(va1, vzero);
    __mmask16 h2 = _mm512_cmpeq_epi32_mask(va2, vb2) & _mm512_cmpneq_epi32_mask(va2, vzero);
    __mmask16 h3 = _mm512_cmpeq_epi32_mask(va3, vb3) & _mm512_cmpneq_epi32_mask(va3, vzero);
    acc0 = _mm512_mask_add_epi32(acc0, h0, acc0, vone);
    acc1 = _mm512_mask_add_epi32(acc1, h1, acc1, vone);
    acc2 = _mm512_mask_add_epi32(acc2, h2, acc2, vone);
    acc3 = _mm512_mask_add_epi32(acc3, h3, acc3, vone);
  }
  for (; k + 16 <= m; k += 16) {
    __m512i va = _mm512_loadu_si512((const __m512i *)(a + k));
    __m512i vb = _mm512_loadu_si512((const __m512i *)(b + k));
    __mmask16 h = _mm512_cmpeq_epi32_mask(va, vb) & _mm512_cmpneq_epi32_mask(va, vzero);
    acc0 = _mm512_mask_add_epi32(acc0, h, acc0, vone);
  }
  acc0 = _mm512_add_epi32(_mm512_add_epi32(acc0, acc1),
                          _mm512_add_epi32(acc2, acc3));
  cnt = (uint32_t)_mm512_reduce_add_epi32(acc0);
#elif defined(__AVX2__)
  // AVX2 fallback: compare 8×uint32 per step.  cmpeq yields 0xFFFFFFFF (= -1) in
  // matching lanes; AND with NOT(a==0) excludes the unset-winner sentinel, then
  // subtract (−(−1)=+1) into an int32 accumulator.  Bit-identical match count.
  const __m256i vzero = _mm256_setzero_si256();
  __m256i acc = vzero;
  for (; k + 8 <= m; k += 8) {
    __m256i va = _mm256_loadu_si256((const __m256i *)(a + k));
    __m256i vb = _mm256_loadu_si256((const __m256i *)(b + k));
    __m256i eq = _mm256_cmpeq_epi32(va, vb);            // -1 where a==b
    __m256i nz = _mm256_cmpeq_epi32(va, vzero);         // -1 where a==0
    acc = _mm256_sub_epi32(acc, _mm256_andnot_si256(nz, eq));
  }
  __m128i s128 = _mm_add_epi32(_mm256_castsi256_si128(acc),
                               _mm256_extracti128_si256(acc, 1));
  s128 = _mm_add_epi32(s128, _mm_shuffle_epi32(s128, _MM_SHUFFLE(1, 0, 3, 2)));
  s128 = _mm_add_epi32(s128, _mm_shuffle_epi32(s128, _MM_SHUFFLE(2, 3, 0, 1)));
  cnt = (uint32_t)_mm_cvtsi128_si32(s128);
#endif
  for (; k < m; ++k)
    cnt += (a[k] != 0 && a[k] == b[k]) ? 1u : 0u;
  return cnt;
}

// 16-bit winner-match count (mirrors pmh_match_count for the packed array):
// counts lanes where a==b AND a!=0 (0 = unset-winner sentinel).
static inline uint32_t pmh_match_count16(const uint16_t *__restrict__ a,
                                         const uint16_t *__restrict__ b,
                                         uint32_t m) {
  uint32_t cnt = 0, k = 0;
#if defined(__AVX512BW__)
  const __m512i vz = _mm512_setzero_si512();
  for (; k + 32 <= m; k += 32) {
    __m512i va = _mm512_loadu_si512((const void *)(a + k));
    __m512i vb = _mm512_loadu_si512((const void *)(b + k));
    __mmask32 eq = _mm512_cmpeq_epi16_mask(va, vb);
    __mmask32 nz = _mm512_cmpneq_epi16_mask(va, vz);
    cnt += (uint32_t)__builtin_popcount((unsigned)(eq & nz));
  }
#elif defined(__AVX2__)
  const __m256i vz = _mm256_setzero_si256();
  __m256i acc = vz;
  for (; k + 16 <= m; k += 16) {
    __m256i va = _mm256_loadu_si256((const __m256i *)(a + k));
    __m256i vb = _mm256_loadu_si256((const __m256i *)(b + k));
    __m256i eq = _mm256_cmpeq_epi16(va, vb);
    __m256i nz = _mm256_cmpeq_epi16(va, vz);         // -1 where a==0
    acc = _mm256_sub_epi16(acc, _mm256_andnot_si256(nz, eq));
  }
  // horizontal sum of 16 int16 lanes (counts ≤ 16 per lane block, no overflow
  // across the ≤ m/16 iterations because we widen below)
  alignas(32) int16_t tmp[16];
  _mm256_store_si256((__m256i *)tmp, acc);
  for (int t = 0; t < 16; ++t) cnt += (uint32_t)tmp[t];
#endif
  for (; k < m; ++k)
    cnt += (a[k] != 0 && a[k] == b[k]) ? 1u : 0u;
  return cnt;
}

// Apply the baseline correction + clamp to a raw winner-match count → sim.
// Shared by the 32-bit and 16-bit similarity paths so both stretch identically.
static inline double pmh_sv_from_count(uint32_t cnt) {
  double j = (double)cnt * g_inv_pmh_m;
  if (g_pmh_base_on && g_pmh_baseline > 0.0 && g_pmh_baseline < 1.0) {
    j = (j - g_pmh_baseline) * g_inv_one_minus_b0;
    if (j < 0.0) j = 0.0;
  }
  if (j > 1.0 - 1e-6) j = 1.0 - 1e-6;
  return j;
}

// Winner-match fraction in [0,1].  Kept for callers that want it directly;
// uses the precomputed reciprocal when m matches the global sketch size.
static inline double pmh_match_frac(const uint32_t *__restrict__ a,
                                    const uint32_t *__restrict__ b, uint32_t m) {
  const uint32_t cnt = pmh_match_count(a, b, m);
  return (m == g_pmh_m && g_inv_pmh_m > 0.0) ? (cnt * g_inv_pmh_m)
                                             : (double)cnt / (double)m;
}

static inline double pmh_sim_rows(const uint32_t *__restrict__ a,
                                  const uint32_t *__restrict__ b, uint32_t m) {
  // count→fraction via precomputed 1/m (no division in the hot loop).
  double j = (double)pmh_match_count(a, b, m) * g_inv_pmh_m;
  // Baseline correction: stretch [b0, 1] → [0, 1] so composition differences
  // regain dynamic range (analogous to b-bit MinHash's (p-2^-b)/(1-2^-b)).
  // Multiply by precomputed 1/(1-b0) instead of dividing.
  if (g_pmh_base_on && g_pmh_baseline > 0.0 && g_pmh_baseline < 1.0) {
    j = (j - g_pmh_baseline) * g_inv_one_minus_b0;
    if (j < 0.0) j = 0.0;
  }
  // Clamp into (0, 1): downstream label propagation evaluates LOG(1 - sim).
  if (j > 1.0 - 1e-6) j = 1.0 - 1e-6;
  return j;
}

// Estimate the chance-collision baseline b0 = median raw winner-match over
// random (mostly unrelated) contig pairs.  With ~10^3 genomes the chance two
// random contigs share a genome is <1%, so the median is firmly the unrelated
// mode.  Cheap (sampled), called once before calibration.
// Forward declaration needed by estimate_pmh_baseline (defined below graph_sim)
static inline double k4_cosine_sim(const float * __restrict__ a,
                                   const float * __restrict__ b);

static double estimate_pmh_baseline(size_t nobs) {
  if (nobs < 4) return 0.0;
  const size_t NP = 30000;
  std::vector<float> vals;
  vals.reserve(NP);
  std::mt19937_64 rng(0x9E3779B97F4A7C15ULL);
  std::uniform_int_distribution<size_t> uni(0, nobs - 1);
  for (size_t s = 0; s < NP; ++s) {
    size_t i = uni(rng), j = uni(rng);
    if (i == j) continue;
    if (g_exact_cos_cmp && !g_k4cosine_flat.empty()) {
      const float *a = g_k4cosine_flat.data() + i * 256;
      const float *b = g_k4cosine_flat.data() + j * 256;
      vals.push_back((float)k4_cosine_sim(a, b));
    } else if (g_win_bits == 16) {
      if (g_pmh_m == 0) continue;
      const uint16_t *a = g_win16.data() + i * g_pmh_m;
      const uint16_t *b = g_win16.data() + j * g_pmh_m;
      vals.push_back((float)((double)pmh_match_count16(a, b, g_pmh_m) * g_inv_pmh_m));
    } else {
      if (g_pmh_m == 0) continue;
      const uint32_t *a = g_win_flat.data() + i * g_pmh_m;
      const uint32_t *b = g_win_flat.data() + j * g_pmh_m;
      vals.push_back((float)pmh_match_frac(a, b, g_pmh_m));
    }
  }
  if (vals.empty()) return 0.0;
  std::nth_element(vals.begin(), vals.begin() + vals.size() / 2, vals.end());
  return (double)vals[vals.size() / 2];
}

// Exact cosine similarity from L2-normalised k=4 GC-norm frequency vectors.
// a,b: pointers to 256-float L2-normalised rows in g_k4cosine_flat.
// Dot product of L2-normalised vectors = cosine similarity ∈ [0,1].
static inline double k4_cosine_sim(const float * __restrict__ a,
                                   const float * __restrict__ b) {
  double s = 0.0;
#if defined(__AVX512F__)
  __m512 acc = _mm512_setzero_ps();
  for (int v = 0; v < 256; v += 16)
    acc = _mm512_fmadd_ps(_mm512_loadu_ps(a + v), _mm512_loadu_ps(b + v), acc);
  s = (double)_mm512_reduce_add_ps(acc);
#elif defined(__AVX2__)
  // AVX2 fallback: two 8-lane accumulators over the fixed 256-float rows.
  __m256 acc0 = _mm256_setzero_ps(), acc1 = _mm256_setzero_ps();
  for (int v = 0; v < 256; v += 16) {
    acc0 = _mm256_add_ps(acc0, _mm256_mul_ps(_mm256_loadu_ps(a + v),     _mm256_loadu_ps(b + v)));
    acc1 = _mm256_add_ps(acc1, _mm256_mul_ps(_mm256_loadu_ps(a + v + 8), _mm256_loadu_ps(b + v + 8)));
  }
  __m256 v8 = _mm256_add_ps(acc0, acc1);
  __m128 lo = _mm_add_ps(_mm256_castps256_ps128(v8), _mm256_extractf128_ps(v8, 1));
  lo = _mm_add_ps(lo, _mm_movehl_ps(lo, lo));
  lo = _mm_add_ss(lo, _mm_shuffle_ps(lo, lo, 0x1));
  s = (double)_mm_cvtss_f32(lo);
#else
  for (int v = 0; v < 256; ++v) s += (double)a[v] * (double)b[v];
#endif
  return (s < 0.0) ? 0.0 : (s > 1.0) ? 1.0 : s;
}

// Unified pairwise similarity dispatcher used by calibration + the all-pairs
// graph build.  Selects the weighted-ProbMinHash winners path (RABBIT_PMH=1),
// exact k=4 cosine similarity (RABBIT_EXACT_COS=1), or the default OPH path.
static inline double graph_sim(size_t i, size_t j) {
  if (g_exact_cos_cmp && !g_k4cosine_flat.empty()) {
    const float *a = g_k4cosine_flat.data() + i * 256;
    const float *b = g_k4cosine_flat.data() + j * 256;
    double s = k4_cosine_sim(a, b);
    // Apply same baseline correction as PMH mode
    if (g_pmh_base_on && g_pmh_baseline > 0.0 && g_pmh_baseline < 1.0)
      s = (s - g_pmh_baseline) / (1.0 - g_pmh_baseline);
    return (s < 0.0) ? 0.0 : (s > 1.0) ? 1.0 : s;
  }
  if (g_pmh_mode) {
    if (g_win_bits == 16) {
      const uint16_t *a = g_win16.data() + (size_t)i * g_pmh_m;
      const uint16_t *b = g_win16.data() + (size_t)j * g_pmh_m;
      return pmh_sv_from_count(pmh_match_count16(a, b, g_pmh_m));
    }
    const uint32_t *a = g_win_flat.data() + (size_t)i * g_pmh_m;
    const uint32_t *b = g_win_flat.data() + (size_t)j * g_pmh_m;
    return pmh_sim_rows(a, b, g_pmh_m);
  }
  const size_t st = (size_t)g_sig_nw * g_sig_np;
  return jaccard_raw(g_sig_flat.data() + i * st, g_sig_flat.data() + j * st,
                     g_sig_nw, g_sig_np, g_sig_m);
}

// Build a frequency-weighted ProbMinHash4 sketch for one contig and copy its
// winner-identity array into out_winners[0..m).  k-mer counts (canonical) are
// used as element weights, so the resulting jaccard_weighted ≈ weighted Jaccard
// over the small-k composition spectrum.  `scratch` is a caller
// owned reusable buffer (avoids per-contig reallocation in the parallel loop).
// out_winners32: folded 32-bit winners (for SIMD pmh_match_frac kernel)
// out_winners64: full 64-bit winners (for inverted-index key lookup), may be null
// out_keys:      getWinnerIndexKeys() output (64-bit, for index building), may be null
static void build_pmh_winners(const char *seq, size_t len, int k, uint32_t m,
                              uint64_t seed, uint32_t *out_winners32,
                              std::vector<uint64_t> &scratch,
                              uint64_t *out_winners64 = nullptr,
                              std::vector<uint64_t> *out_keys = nullptr,
                              float *out_k4freq = nullptr) {
  Sketch::ProbMinHash4 pmh(m, k, seed);
  if (len >= (size_t)k && k >= 1 && k <= 32) {
    const uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    const int hishift = 2 * (k - 1);
    uint64_t fwd = 0, rev = 0;
    int valid = 0;
    // Small-k fast path: the canonical-code universe is 4^k, so for k<=8
    // (<=65536 buckets) a direct frequency histogram is far cheaper than
    // sorting every k-mer. Exact-equivalent to the sort/run-length path.
    if (k <= 8) {
      const size_t NB = (size_t)1 << (2 * k);
      scratch.assign(NB, 0);
      size_t total = 0;
      // Base + dinucleotide frequency counts for normalisation (A=0,C=1,G=2,T=3)
      uint64_t base_cnt[4]    = {};
      uint64_t dinuc_cnt[4][4] = {};
      uint8_t  prev_base = 255;
      for (size_t i = 0; i < len; ++i) {
        const uint8_t b = rabbit_sketch::RB_BASE_ENC_LUT[(uint8_t)seq[i]];
        if (b == 255) { valid = 0; fwd = 0; rev = 0; prev_base = 255; continue; }
        base_cnt[b]++;
        if (prev_base != 255) dinuc_cnt[prev_base][b]++;
        prev_base = b;
        fwd = ((fwd << 2) | b) & mask;
        rev = (rev >> 2) | ((uint64_t)(3u ^ b) << hishift);
        if (++valid >= k) { ++scratch[(fwd <= rev) ? fwd : rev]; ++total; }
      }
      if (total) {
        if (!g_gc_norm) {
          // Original path: raw normalised frequency as weight.
          const double inv = 1.0 / (double)total;
          for (uint64_t v = 0; v < NB; ++v)
            if (scratch[v]) pmh.addHash(v, (double)scratch[v] * inv);
        } else if (g_gc_norm == 1) {
          // Mode 1: per-base independence model.
          // E(v) = P(fwd) + P(rev_comp)   where P = product of per-base frequencies.
          const size_t vb = base_cnt[0]+base_cnt[1]+base_cnt[2]+base_cnt[3];
          double p[4];
          for (int b = 0; b < 4; ++b) p[b] = vb > 0 ? (double)base_cnt[b]/(double)vb : 0.25;
          const double inv = 1.0 / (double)total;
          const double cap = g_gc_norm_cap;
          if (out_k4freq && k == 4) std::fill(out_k4freq, out_k4freq + 256, 0.0f);
          for (uint64_t v = 0; v < NB; ++v) {
            if (!scratch[v]) continue;
            double p_fwd = 1.0, p_rev = 1.0;
            uint64_t code = v;
            for (int pos = 0; pos < k; ++pos) {
              const uint8_t b = (uint8_t)(code & 3u);
              p_fwd *= p[b]; p_rev *= p[3u - b];
              code >>= 2;
            }
            const double obs_f = (double)scratch[v] * inv;
            const double exp_f = p_fwd + p_rev;
            const double w = (exp_f > 1e-14) ? std::min(obs_f / exp_f, cap) : 0.0;
            if (w > 0.0) {
              pmh.addHash(v, w);
              if (out_k4freq && k == 4) out_k4freq[v] = (float)w;
            }
          }
        } else {
          // Mode 2: first-order Markov (dinucleotide) model.
          // P(ABCDEF) = P(A) × P(B|A) × P(C|B) × P(D|C) × P(E|D) × P(F|E)
          // P(b|a) estimated from contig's own dinucleotide counts (+ pseudocount=1).
          // Complement mapping A(0)↔T(3), C(1)↔G(2), i.e. comp = 3-b.
          const size_t vb = base_cnt[0]+base_cnt[1]+base_cnt[2]+base_cnt[3];
          double p_mono[4];
          for (int b = 0; b < 4; ++b)
            p_mono[b] = vb > 0 ? (base_cnt[b] + 1.0) / (vb + 4.0) : 0.25;
          // Conditional P(b|a) with Laplace pseudocount=1
          double p_cond[4][4];
          for (int a = 0; a < 4; ++a) {
            double rs = 4.0;  // pseudocount sum
            for (int b = 0; b < 4; ++b) rs += dinuc_cnt[a][b];
            for (int b = 0; b < 4; ++b)
              p_cond[a][b] = (dinuc_cnt[a][b] + 1.0) / rs;
          }
          const double inv = 1.0 / (double)total;
          const double cap = g_gc_norm_cap;
          for (uint64_t v = 0; v < NB; ++v) {
            if (!scratch[v]) continue;
            // Decode k bases: bases[0]=lsb=first in sequence
            uint8_t bases[8];   // max k=8
            { uint64_t code = v; for (int i = 0; i < k; ++i) { bases[i]=(uint8_t)(code&3u); code>>=2; } }
            // P(fwd)
            double p_fwd = p_mono[bases[0]];
            for (int i = 1; i < k; ++i) p_fwd *= p_cond[bases[i-1]][bases[i]];
            // P(rev_comp): reverse the complement sequence
            // rev_comp[i] = 3 - bases[k-1-i]
            double p_rev = p_mono[3u - bases[k-1]];
            for (int i = 1; i < k; ++i)
              p_rev *= p_cond[3u - bases[k-i]][3u - bases[k-1-i]];
            const double exp_f = p_fwd + p_rev;
            const double obs_f = (double)scratch[v] * inv;
            const double w = (exp_f > 1e-16) ? std::min(obs_f / exp_f, cap) : 0.0;
            if (w > 0.0) pmh.addHash(v, w);
          }
        }
      }
    } else {
      scratch.clear();
      for (size_t i = 0; i < len; ++i) {
        const uint8_t b = rabbit_sketch::RB_BASE_ENC_LUT[(uint8_t)seq[i]];
        if (b == 255) { valid = 0; fwd = 0; rev = 0; continue; }
        fwd = ((fwd << 2) | b) & mask;
        rev = (rev >> 2) | ((uint64_t)(3u ^ b) << hishift);
        if (++valid >= k)
          scratch.push_back((fwd <= rev) ? fwd : rev);
      }
      if (!scratch.empty()) {
        std::sort(scratch.begin(), scratch.end());
        const double inv = 1.0 / (double)scratch.size();
        size_t i = 0;
        while (i < scratch.size()) {
          size_t j = i + 1;
          while (j < scratch.size() && scratch[j] == scratch[i]) ++j;
          pmh.addHash(scratch[i], (double)(j - i) * inv);
          i = j;
        }
      }
    }
  }
  // Extract full-64-bit index keys (before folding).
  if (out_keys) pmh.getWinnerIndexKeys(*out_keys);

  const uint64_t *w64 = pmh.getWinners();
  for (uint32_t r = 0; r < m; ++r) {
    uint64_t h = w64[r];
    if (out_winners64) out_winners64[r] = h;
    // Fold 64→32 for the fast SIMD pairwise kernel.
    out_winners32[r] = (uint32_t)(h ^ (h >> 32));
  }
}

// OPH inverted index built inline during sketch construction (RABBIT_GRAPH_INDEX=1).
// Nullptr when not in use (default: all-pairs path).
static std::unique_ptr<rabbit_invidx::InvertedIndex> g_inv_idx;

// PMH-winner inverted index (RABBIT_PMH=1, default ON).
// Key = ((uint64_t)register_pos << 32) | folded_winner_32bit.
// Collision count between two contigs = PMH match numerator → raw_sim = count/m.
// Built inline during sketch loop, used for both calibration and graph build.
// Eliminates both O(N²) calibration sweeps and the all-pairs graph scan.
static std::unique_ptr<rabbit_invidx::InvertedIndex> g_pmh_idx;

// ── OMP reduction operators (must be visible before parallel loops) ───────
#pragma omp declare reduction (merge_size_t : std::vector<size_t> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()) )
#pragma omp declare reduction (merge_dist : std::vector<Distance> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()) )
#pragma omp declare reduction (merge_storeddist : std::vector<StoredDistance> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()) )

// Precomputed per-contig "has any abundance sample > minCV" flag.
// is_nz(i,j) == (g_anynz[i] || g_anynz[j]); collapses the per-pair O(num_depth_samples)
// matrix scan (run N²/2 times in the hot graph loop) into two byte loads.
// Empty when not yet built (is_nz then falls back to the scalar scan).
static std::vector<uint8_t> g_anynz;

// ── Forward declarations for helpers used by build_similarity_graph ────────────────
bool is_nz(size_t r1, size_t r2);
static void build_anynz_cache();
std::istream &safeGetline(std::istream &is, std::string &t);
size_t   calibrate_sim_cutoff(Distance coverage, bool full);
static Distance calibrate_sim_cutoff_fused(Distance coverage);
static Distance gen_fused_calib_graph(Graph &g, Distance coverage);

// ═══════════════════════════════════════════════════════════════════════════
// filesize  (unchanged)
// ═══════════════════════════════════════════════════════════════════════════
std::ifstream::pos_type filesize(const char *filename) {
  std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
  return in.tellg();
}

// ═══════════════════════════════════════════════════════════════════════════
// FASTA decompression strategy
//
// Two modes selected automatically based on compressed file size:
//
//  Small file  (compressed < RABBIT_LIBDEFLATE_MAXGB, default 2 GB):
//    libdeflate one-shot — 2× faster than zlib, but loads the entire
//    decompressed buffer into RAM (up to ~4× file size).  Safe when the
//    machine has enough free memory.
//
//  Large file  (compressed ≥ threshold):
//    Streaming zlib (gzread) — zero extra memory overhead, processes the
//    file chunk-by-chunk.  Mandatory for multi-GB assemblies to avoid OOM.
//
// Override threshold: RABBIT_LIBDEFLATE_MAXGB=<N>  (0 = always stream)
// ═══════════════════════════════════════════════════════════════════════════
struct FastaRecord {
  std::string name;  // header up to first whitespace
  const char* seq;   // pointer into the decompressed buffer (NOT NUL-terminated)
  size_t       seq_len;
};

// Parse FASTA in a flat memory buffer.  Calls cb(name, seq_ptr, seq_len) for
// each record.  Uses memchr for fast line scanning and bulk append for
// parse_fasta_buf  –  parse a FASTA in-memory buffer and call cb per record.
// Uses memchr for fast delimiter finding and bulk append per line segment for
// multi-line sequences, avoiding character-by-character iteration.
// cb signature: cb(const std::string& name, const char* seq_ptr, size_t seq_len)
template<typename CB>
static void parse_fasta_buf(const char* buf, size_t len, CB&& cb) {
  const char* p   = buf;
  const char* end = buf + len;

  while (p < end) {
    // Find next record header '>'
    const char* gt = (const char*)memchr(p, '>', (size_t)(end - p));
    if (!gt) break;
    p = gt + 1;

    // Header: name until whitespace, skip rest of line
    const char* name_start = p;
    while (p < end && (uint8_t)*p > ' ') ++p;
    std::string name(name_start, p);
    const char* nl = (const char*)memchr(p, '\n', (size_t)(end - p));
    if (!nl) break;
    p = nl + 1;

    // Sequence region: up to the next '>' or end of buffer
    const char* next_gt = (const char*)memchr(p, '>', (size_t)(end - p));
    const char* seq_end = next_gt ? next_gt : end;

    // Fast path: single-line sequence (no embedded '\n')
    const char* inner_nl = (const char*)memchr(p, '\n', (size_t)(seq_end - p));
    if (!inner_nl) {
      cb(name, p, (size_t)(seq_end - p));
      p = seq_end;
      continue;
    }

    // Multi-line: compact newlines using bulk memcpy per line segment
    std::string seq;
    seq.reserve((size_t)(seq_end - p));
    const char* seg = p;
    while (seg < seq_end) {
      const char* seg_nl = (const char*)memchr(seg, '\n', (size_t)(seq_end - seg));
      const char* seg_end2 = seg_nl ? seg_nl : seq_end;
      size_t seg_len = (size_t)(seg_end2 - seg);
      if (seg_len > 0 && seg[seg_len - 1] == '\r') --seg_len;  // strip \r
      if (seg_len > 0) seq.append(seg, seg_len);
      if (!seg_nl) break;
      seg = seg_nl + 1;
    }
    cb(name, seq.c_str(), seq.size());
    p = seq_end;
  }
}

// ─── Parallel mmap FASTA reader ──────────────────────────────────────────────
// Borrows RabbitFX's producer-consumer spirit but uses mmap instead of
// fread chunks: the OS pages are shared across all threads so there is
// zero per-thread I/O setup cost.
//
// Algorithm (mirrors RabbitTClust's sketchByFile × N threads, but for one file):
//   1. mmap the whole file (MAP_PRIVATE | MAP_POPULATE on small files).
//   2. Find N split positions — scan for '\n>' to align splits to record starts.
//   3. N OMP threads each parse [split[t], split[t+1]) independently into a
//      thread-local ContigRecord vector. For PMH mode winners are also built
//      inline (no alloc per contig — just reuse the thread-local scratch).
//   4. Serial merge: concatenate per-thread vectors in order → stable global
//      index, then populate g_win_flat from the per-thread winner arrays.
//
// Returns true on success.  Falls back transparently (returns false) on mmap
// failure so the caller can use the kseq fallback.

struct ContigRec {
  std::string  name;
  // Sequence location.  Exactly one of:
  //   seq_view_ptr != nullptr        → single-line: zero-copy slice of the mmap
  //   arena_tid    >= 0              → multi-line:  [arena_off,+len) in g_seq_arenas[tid]
  const char  *seq_view_ptr = nullptr;
  int          arena_tid = -1;
  size_t       arena_off = 0;
  size_t       len = 0;
  bool         is_small = false;
  std::vector<uint32_t> winners; // PMH winners (size pmh_m), empty if !stream_pmh
  std::vector<float>    k4freq;  // exact GC-norm k=4 weights [256], for IDF rebuild
};

static bool parse_fasta_mmap_parallel(
    const std::string   &path,
    int                  nthreads,
    size_t               minContig_arg,
    size_t               min_small_contig_arg,
    bool                 fullHeader_arg,
    bool                 stream_pmh_arg,
    bool                 no_store_seqs_arg,
    bool                 collect_tiny_arg,
    uint32_t             pmh_m_arg,
    int                  pmh_k_arg,
    // outputs  ──────────────────────────────────────────────────────────
    std::vector<ContigRec> &large_out,   // large contigs, file order
    std::vector<ContigRec> &small_out,   // small contigs, file order
    std::vector<std::pair<std::string,std::string>> &tiny_out, // <min_small_contig (name,seq)
    size_t               &num_seqs_out
) {
  // ── 1. mmap ─────────────────────────────────────────────────────────────
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0) return false;

  struct stat st{};
  if (fstat(fd, &st) < 0) { ::close(fd); return false; }
  const size_t fsz = (size_t)st.st_size;
  if (fsz == 0) { ::close(fd); return false; }

  // Use MAP_POPULATE only for files < 4 GB to avoid long fault storms on huge files.
  int mflags = MAP_PRIVATE;
  // NOTE: MAP_POPULATE is intentionally NOT used.  The parse reads each page
  // exactly once, sequentially, and (for multi-line records) MADV_DONTNEEDs it
  // immediately afterwards, so prefaulting the whole file would only inflate
  // peak RSS by the full file size before the incremental release can run.
  // MADV_SEQUENTIAL (below) already drives aggressive read-ahead.
  if (const char *ev = getenv("RABBIT_MMAP_POPULATE"))
    if (ev[0] == '1' && fsz < 4ULL * 1024 * 1024 * 1024) mflags |= MAP_POPULATE;
  void *mptr = mmap(nullptr, fsz, PROT_READ, mflags, fd, 0);
  ::close(fd);
  if (mptr == MAP_FAILED) return false;

  // Advise the kernel we'll scan sequentially.
  madvise(mptr, fsz, MADV_SEQUENTIAL);
  rb_phase("  mmap done");

  const char *buf = static_cast<const char *>(mptr);
  // Retain the mapping for the whole run: single-line sequences are stored as
  // zero-copy string_views into this buffer (see ContigRec::seq_view_ptr), so
  // it must NOT be unmapped after parsing.  The OS reclaims it at exit; since
  // the pages are read-only and file-backed they stay shared with the page
  // cache and cost little resident memory.
  g_fasta_mmap     = buf;
  g_fasta_mmap_len = fsz;

  // ── 2. Find N chunk boundaries ──────────────────────────────────────────
  // Each boundary is the offset of the '>' that starts a record.
  // Thread t owns bytes [splits[t], splits[t+1]).
  std::vector<size_t> splits(nthreads + 1);
  splits[0] = 0;
  splits[nthreads] = fsz;

  for (int t = 1; t < nthreads; ++t) {
    size_t approx = (fsz / nthreads) * t;
    // Scan forward to the next '\n>' (record boundary)
    const char *p = buf + approx;
    const char *end = buf + fsz;
    bool found = false;
    while (p < end - 1) {
      if (*p == '\n' && *(p + 1) == '>') {
        splits[t] = (size_t)(p + 1 - buf); // points at '>'
        found = true;
        break;
      }
      ++p;
    }
    if (!found) {
      // No boundary found — collapse remaining threads into last chunk
      for (int tt = t; tt < nthreads; ++tt) splits[tt] = fsz;
      break;
    }
  }

  // Compute actual number of non-empty threads
  int active_threads = 1;
  for (int t = 1; t <= nthreads; ++t)
    if (splits[t] > splits[t - 1]) active_threads = t;

  // ── 3. Parallel parse per chunk ─────────────────────────────────────────
  std::vector<std::vector<ContigRec>> tl_large(nthreads);
  std::vector<std::vector<ContigRec>> tl_small(nthreads);
  std::vector<std::vector<std::pair<std::string,std::string>>> tl_tiny(nthreads);
  std::vector<size_t> tl_numseqs(nthreads, 0);
  std::vector<std::vector<uint64_t>> tl_scratch(nthreads); // PMH scratch
  // Per-thread reusable multi-line compaction scratch, used ONLY when sequences
  // are not retained (no_store_seqs): each record is compacted here, sketched,
  // then overwritten by the next record, so the persistent arena never grows.
  std::vector<std::string> tl_seqbuf(nthreads);
  // Per-thread compacted-sequence arenas (multi-line records).  Each thread
  // appends only to its own arena, so there is no cross-thread contention; the
  // arenas are retained globally so seq views stay valid for the whole run.
  g_seq_arenas.assign(nthreads, std::string());

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
  for (int t = 0; t < nthreads; ++t) {
    if (splits[t] >= splits[t + 1]) continue;
    const char *p   = buf + splits[t];
    const char *end = buf + splits[t + 1];
    auto &my_large  = tl_large[t];
    auto &my_small  = tl_small[t];
    auto &my_tiny   = tl_tiny[t];
    auto &my_scratch= tl_scratch[t];
    auto &my_seqbuf = tl_seqbuf[t];
    size_t &my_nseq = tl_numseqs[t];

    // Reserve the arena once to the chunk's byte span: sequence content is the
    // vast majority of a FASTA chunk, so this avoids std::string's doubling
    // growth (which would over-allocate up to ~2× and inflate peak RSS) while
    // never reallocating during the parse.
    if (!no_store_seqs_arg)
      g_seq_arenas[t].reserve((size_t)(splits[t + 1] - splits[t]));

    // Incremental mmap page release: this thread reads its file region exactly
    // once, sequentially.  For multi-line records the bytes are compacted into
    // the arena and the mmap pages are never touched again, so we MADV_DONTNEED
    // the consumed range as we advance.  This collapses the parse-phase peak RSS
    // from "full mmap + arenas" to roughly "arenas only" on multi-line FASTAs
    // (the CAMI gold-standard assemblies).  Release is disabled the moment a
    // single-line record retains a zero-copy mmap view (its pages must survive);
    // page-aligned bounds keep the boundary pages shared with neighbour threads
    // intact.  No-op for single-line assemblies (which keep using mmap views and
    // never double-store, so their peak is already just the mmap).
    const long  _pg       = sysconf(_SC_PAGESIZE);
    const size_t REL_GRAN = 8u << 20;        // release in ≥8 MB aligned spans
    const char  *rel_base = p;               // start of not-yet-released region
    bool         can_release = true;         // false once an mmap view is kept

    // Fast FASTA parser — same logic as parse_fasta_buf
    while (p < end) {
      const char *gt = (const char *)memchr(p, '>', (size_t)(end - p));
      if (!gt) break;
      p = gt + 1;

      // Name EXTENT only (no allocation yet): the std::string is materialised
      // after we know the record is kept, so discarded tiny contigs (millions of
      // them on fragmented assemblies) never pay for a name allocation.
      const char *name_start = p;
      while (p < end && (uint8_t)*p > ' ') ++p;
      const char *name_end = p;
      const char *nl = (const char *)memchr(p, '\n', (size_t)(end - p));
      if (!nl) break;
      p = nl + 1;

      // Sequence body (up to next '>' or chunk end)
      const char *next_gt = (const char *)memchr(p, '>', (size_t)(end - p));
      const char *seq_end = next_gt ? next_gt : end;

      // Compact sequence (may be multi-line).
      // Single-line records (the overwhelmingly common case in assembly FASTA)
      // are referenced ZERO-COPY: the sequence is a contiguous slice of the
      // mmap, so we keep only a (ptr,len) view and never allocate/memcpy it.
      // Multi-line records still need newline removal → an owned compacted copy.
      std::string &my_arena = g_seq_arenas[t];
      const char *view_ptr = nullptr;  // non-null ⇒ zero-copy view into mmap
      int         arena_tid = -1;      // >=0 ⇒ compacted copy in my_arena
      size_t      arena_off = 0;
      size_t      seq_len;
      const char *seq_region = p;      // sequence bytes start (single-line view base)

      // Fast tiny-prune: seq_end - seq_region is an UPPER BOUND on seq_len
      // (it includes newlines which only reduce the count after removal).
      // If even the raw byte span is below min_small_contig, the record is
      // definitely tiny — skip the inner_nl scan, multi-line copy, and name
      // materialisation entirely.  On assemblies like CAMI2 plant (3.13M tiny
      // contigs out of 3.44M total), this saves the dominant memchr + arena
      // append + resize-undo for 91% of records.
      if (!collect_tiny_arg &&
          (size_t)(seq_end - seq_region) < (size_t)min_small_contig_arg) {
        p = seq_end;
        ++my_nseq;
        // MADV_DONTNEED: pages up to current p are consumed for this thread.
        if (can_release) {
          uintptr_t a = ((uintptr_t)rel_base + (uintptr_t)_pg - 1) & ~((uintptr_t)_pg - 1);
          uintptr_t b = (uintptr_t)p & ~((uintptr_t)_pg - 1);
          if (b > a && (b - a) >= REL_GRAN) {
            madvise((void *)a, (size_t)(b - a), MADV_DONTNEED);
            rel_base = (const char *)b;
          }
        }
        continue;
      }
      // A record is effectively single-line when the only newlines are the
      // trailing one(s) before the next '>'.  In that case the sequence is a
      // contiguous mmap slice and we reference it zero-copy.  (Note: the first
      // '\n' is the END of the sequence line, not evidence of multi-line — the
      // old `!memchr('\n')` test almost never fired since every record but the
      // last ends in a newline.)  Multi-line records are compacted once into
      // the thread arena (no per-contig allocation).
      const char *inner_nl = (const char *)memchr(p, '\n', (size_t)(seq_end - p));
      bool single_line;
      if (!inner_nl) {
        single_line = true;
        seq_len = (size_t)(seq_end - p);
      } else {
        single_line = true;
        for (const char *q = inner_nl + 1; q < seq_end; ++q)
          if (*q != '\n' && *q != '\r') { single_line = false; break; }
        seq_len = (size_t)(inner_nl - p);
      }
      if (single_line) {
        if (seq_len > 0 && p[seq_len - 1] == '\r') --seq_len;
        // view_ptr set below, only once the record is known to be kept.
      } else {
        // Multi-line: compact (strip newlines).  When sequences are not retained
        // (no_store_seqs), compact into a reusable per-thread scratch that is
        // overwritten each record, so the persistent arena never grows — the
        // dominant RAM saving on multi-line assemblies.  Otherwise append to the
        // retained arena so stored views stay valid.
        std::string &dst = no_store_seqs_arg ? my_seqbuf : my_arena;
        if (no_store_seqs_arg) dst.clear();
        arena_off = dst.size();
        const char *seg = p;
        while (seg < seq_end) {
          const char *seg_nl = (const char *)memchr(seg, '\n', (size_t)(seq_end - seg));
          const char *seg_e  = seg_nl ? seg_nl : seq_end;
          size_t slen = (size_t)(seg_e - seg);
          if (slen > 0 && seg[slen - 1] == '\r') --slen;
          if (slen > 0) dst.append(seg, slen);
          if (!seg_nl) break;
          seg = seg_nl + 1;
        }
        seq_len   = dst.size() - arena_off;
        if (!no_store_seqs_arg) arena_tid = t;
      }

      p = seq_end;
      ++my_nseq;

      // Early-skip discarded tiny contigs.  When unbinned output is off the tiny
      // (name,seq) pair is never consumed downstream, so skip the name allocation
      // and the full-sequence copy entirely — the dominant per-record cost on
      // assemblies with millions of sub-min_small_contig fragments.
      const bool is_large = seq_len >= minContig_arg;
      const bool is_small = !is_large && seq_len >= min_small_contig_arg;
      if (!is_large && !is_small && !collect_tiny_arg) {
        if (arena_tid >= 0) my_arena.resize(arena_off);  // undo multi-line tiny copy
        continue;
      }

      // Record is kept: materialise the name and (single-line) the zero-copy view.
      std::string name(name_start, name_end);
      if (!fullHeader_arg) {
        size_t sp = name.find_first_of(" \t");
        if (sp != std::string::npos) name.resize(sp);
      }
      // Single-line records keep a zero-copy mmap view; that pins the mmap pages
      // so page release must stop.  When the view is NOT retained (no_store_seqs)
      // we only read it now for sketching, so pages can still be released.
      if (single_line) { view_ptr = seq_region; if (!no_store_seqs_arg) can_release = false; }
      // Pointer to raw sequence bytes (mmap view, arena slice, or reusable
      // no-store scratch) for the k-mer scan.  Valid here because the chosen
      // buffer is not appended again before this use.
      const char *seq_bytes = view_ptr ? view_ptr
                              : ((no_store_seqs_arg ? my_seqbuf : my_arena).data() + arena_off);

      if (is_large) {
        ContigRec rec;
        rec.name = std::move(name);
        rec.len  = seq_len;
        rec.is_small = false;
        if (stream_pmh_arg) {
          rec.winners.resize(pmh_m_arg, 0u);
          float *k4out = nullptr;
          if ((g_idf_norm || g_exact_cos_cmp) && pmh_k_arg == 4 && g_gc_norm >= 1) {
            rec.k4freq.resize(256, 0.0f);
            k4out = rec.k4freq.data();
          }
          build_pmh_winners(seq_bytes, seq_len,
                            pmh_k_arg, pmh_m_arg, /*seed=*/42u,
                            rec.winners.data(), my_scratch,
                            /*out64=*/nullptr, /*out_keys=*/nullptr, k4out);
        }
        if (!no_store_seqs_arg) {
          rec.seq_view_ptr = view_ptr;   // zero-copy (single-line) or nullptr
          rec.arena_tid    = arena_tid;  // arena slice (multi-line) or -1
          rec.arena_off    = arena_off;
        }
        my_large.push_back(std::move(rec));
      } else if (is_small) {
        ContigRec rec;
        rec.name = std::move(name);
        rec.len  = seq_len;
        rec.is_small = true;
        if (!no_store_seqs_arg) {
          rec.seq_view_ptr = view_ptr;
          rec.arena_tid    = arena_tid;
          rec.arena_off    = arena_off;
        }
        my_small.push_back(std::move(rec));
      } else {
        // Tiny contig, kept only because collect_tiny is set (unbinned output):
        // materialise (name, seq) text for emission.
        my_tiny.emplace_back(std::move(name), std::string(seq_bytes, seq_len));
      }

      // Drop the mmap pages consumed so far (page-aligned, conservative bounds
      // so the boundary pages shared with neighbour threads are untouched).
      if (can_release) {
        uintptr_t a = ((uintptr_t)rel_base + (uintptr_t)_pg - 1) & ~((uintptr_t)_pg - 1);
        uintptr_t b = (uintptr_t)p & ~((uintptr_t)_pg - 1);
        if (b > a && (b - a) >= REL_GRAN) {
          madvise((void *)a, (size_t)(b - a), MADV_DONTNEED);
          rel_base = (const char *)b;
        }
      }
    }

    // Tail release: free this thread's remaining (fully consumed) region.
    if (can_release) {
      uintptr_t a = ((uintptr_t)rel_base + (uintptr_t)_pg - 1) & ~((uintptr_t)_pg - 1);
      uintptr_t b = (uintptr_t)end & ~((uintptr_t)_pg - 1);
      if (b > a) madvise((void *)a, (size_t)(b - a), MADV_DONTNEED);
    }
  }

  rb_phase("  parallel parse loop done");
  // NOTE: do NOT munmap(mptr, fsz) — zero-copy seq views point into it.

  // ── 4. Parallel merge in file order ─────────────────────────────────────
  // Prefix-sum the per-thread record counts to compute each thread's output
  // slice, resize the output vectors once, then move each thread's records into
  // its disjoint slice in parallel.  This yields byte-for-byte the same ordering
  // as a serial concatenation (thread blocks in ascending order, records within
  // a block in file order) while removing the O(total) serial copy.
  num_seqs_out = 0;
  std::vector<size_t> off_large(nthreads + 1, 0);
  std::vector<size_t> off_small(nthreads + 1, 0);
  std::vector<size_t> off_tiny(nthreads + 1, 0);
  for (int t = 0; t < nthreads; ++t) {
    num_seqs_out += tl_numseqs[t];
    off_large[t + 1] = off_large[t] + tl_large[t].size();
    off_small[t + 1] = off_small[t] + tl_small[t].size();
    off_tiny[t + 1]  = off_tiny[t]  + tl_tiny[t].size();
  }
  large_out.resize(off_large[nthreads]);
  small_out.resize(off_small[nthreads]);
  tiny_out.resize(off_tiny[nthreads]);
#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
  for (int t = 0; t < nthreads; ++t) {
    std::move(tl_large[t].begin(), tl_large[t].end(), large_out.begin() + off_large[t]);
    std::move(tl_small[t].begin(), tl_small[t].end(), small_out.begin() + off_small[t]);
    std::move(tl_tiny[t].begin(),  tl_tiny[t].end(),  tiny_out.begin()  + off_tiny[t]);
  }
  rb_phase("  parse serial merge done");
  return true;
}

// ─── Decompression backend for the streaming reader ─────────────────────────
// A tiny façade over three byte sources, chosen at open() time:
//   • PLAIN  : fopen + fread            (uncompressed input)
//   • ZLIB   : gzopen + gzread          (.gz, default)
//   • IGZIP  : ISA-L streaming inflate  (.gz, when built with USE_IGZIP)
// ISA-L's igzip decompresses a single gzip stream ~2–4× faster than zlib, which
// is the bottleneck on hundred-GB .gz assemblies (a single DEFLATE stream's
// inflate is inherently sequential, so it lives on one producer thread).  When
// ISA-L is not available the class compiles to the zlib path unchanged, so the
// streaming reader's behaviour and results are identical with or without ISA-L.
// read() returns the number of bytes produced (0 = EOF, <0 = error); on a
// corrupt/truncated stream the ISA-L path returns the bytes decoded so far and
// then reports EOF, matching zlib's lenient gzread behaviour.
#if defined(USE_IGZIP)
#include <igzip_lib.h>
#endif
struct RbGzReader {
  enum Mode { PLAIN, ZLIB
#if defined(USE_IGZIP)
            , IGZIP
#endif
  };
  Mode    mode  = PLAIN;
  FILE   *fp    = nullptr;
  gzFile  gz    = nullptr;
  bool    eof_  = false;
#if defined(USE_IGZIP)
  static const size_t   IN_SZ   = 1u << 22;   // 4 MB compressed read granularity
  static const uint32_t HDR_REQ = 1u << 16;   // bytes needed to (re)read a gz header
  unsigned char        *inbuf   = nullptr;
  struct inflate_state  st;
  struct isal_gzip_header hdr;
#endif

  bool open(const char *path, bool is_gz) {
    if (!is_gz) { mode = PLAIN; fp = fopen(path, "rb"); return fp != nullptr; }
#if defined(USE_IGZIP)
    mode = IGZIP;
    fp = fopen(path, "rb");
    if (!fp) return false;
    inbuf = (unsigned char *)malloc(IN_SZ);
    if (!inbuf) { fclose(fp); fp = nullptr; return false; }
    isal_gzip_header_init(&hdr);
    isal_inflate_init(&st);
    st.crc_flag = ISAL_GZIP_NO_HDR_VER;
    st.next_in  = inbuf;
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

  int64_t read(char *out, size_t size) {
    if (mode == PLAIN) return (int64_t)fread(out, 1, size, fp);
#if defined(USE_IGZIP)
    if (mode == IGZIP) {
      if (eof_) return 0;
      unsigned char *o = (unsigned char *)out;
      uint64_t off = 0;
      while (off < size) {
        if (st.avail_in == 0) {
          if (feof(fp)) { eof_ = true; break; }
          st.next_in  = inbuf;
          st.avail_in = (uint32_t)fread(inbuf, 1, IN_SZ, fp);
          if (st.avail_in == 0) { eof_ = true; break; }
        }
        st.next_out  = o + off;
        st.avail_out = (uint32_t)(size - off);
        int r = isal_inflate(&st);
        off = (uint64_t)(st.next_out - o);
        if (r != ISAL_DECOMP_OK) { eof_ = true; break; }  // tolerate truncation
        if (st.block_state == ISAL_BLOCK_FINISH) {
          if (feof(fp) && st.avail_in == 0) { eof_ = true; break; }
          // A concatenated gzip member may follow (pigz / cat a.gz b.gz …).
          if (st.avail_in == 0) {
            isal_inflate_reset(&st);
            st.next_in  = inbuf;
            st.avail_in = (uint32_t)fread(inbuf, 1, IN_SZ, fp);
            if (st.avail_in == 0) { eof_ = true; break; }
          } else if (st.avail_in >= HDR_REQ) {
            unsigned char *oni = st.next_in; uint32_t oai = st.avail_in;
            isal_inflate_reset(&st);
            st.next_in = oni; st.avail_in = oai;
          } else {
            uint32_t oai = st.avail_in;
            memmove(inbuf, st.next_in, oai);
            size_t rd = feof(fp) ? 0 : fread(inbuf + oai, 1, IN_SZ - oai, fp);
            isal_inflate_reset(&st);
            st.next_in = inbuf; st.avail_in = oai + (uint32_t)rd;
          }
          if (isal_read_gzip_header(&st, &hdr) != ISAL_DECOMP_OK) { eof_ = true; break; }
        }
      }
      return (int64_t)off;
    }
#endif
    return (int64_t)gzread(gz, out, (unsigned)size);  // ZLIB
  }

  void close() {
    if (gz) { gzclose(gz); gz = nullptr; }
    if (fp) { fclose(fp);  fp = nullptr; }
#if defined(USE_IGZIP)
    if (inbuf) { free(inbuf); inbuf = nullptr; }
#endif
  }
};

// ─── Streaming parallel gzip FASTA reader (large .gz) ────────────────────────
// For .gz assemblies too large for the one-shot libdeflate path (default
// > RABBIT_LIBDEFLATE_MAXGB compressed), the historical fallback was a single
// threaded kseq loop: one thread inflated AND parsed AND compacted AND copied
// AND sketched every record, so throughput collapsed to a single core on
// hundred-GB assemblies.
//
// A single gzip stream's inflate is inherently sequential, so it stays on ONE
// dedicated producer thread; everything else (FASTA parsing, newline removal,
// sequence storage, PMH sketching) moves onto N worker threads.  The producer
// cuts the decompressed byte stream into record-aligned blocks — every block
// both starts and ends on a '>' record boundary, so no record ever straddles
// two blocks — and feeds them through a bounded queue (backpressure caps the
// in-flight decompressed memory).  Workers drain the queue in parallel, parse
// each block into thread-local ContigRec lists (sequences compacted into the
// per-worker arenas, identical storage model to the mmap path) and build the
// PMH sketches inline.  Per-block results are stored by block id and merged in
// file order, so the resulting contig indexing is byte-for-byte identical to
// the serial kseq path — and therefore so is every downstream result.  This
// touches NEITHER the <threshold libdeflate path (used by CAMI-scale .gz) NOR
// the uncompressed mmap path, so their behaviour and results are unchanged.
//
// Output shape mirrors parse_fasta_mmap_parallel() so the caller can share the
// same merge / IDF-rebuild / exact-cosine post-processing.  Returns false (→
// caller uses the kseq fallback) only if the file cannot be opened.
static bool parse_fasta_stream_parallel(
    const std::string   &path,
    int                  nthreads,
    size_t               minContig_arg,
    size_t               min_small_contig_arg,
    bool                 fullHeader_arg,
    bool                 stream_pmh_arg,
    bool                 no_store_seqs_arg,
    bool                 collect_tiny_arg,
    uint32_t             pmh_m_arg,
    int                  pmh_k_arg,
    bool                 input_is_gz,    // .gz → zlib/igzip; else → plain fread
    // outputs  ──────────────────────────────────────────────────────────
    std::vector<ContigRec> &large_out,
    std::vector<ContigRec> &small_out,
    std::vector<std::pair<std::string,std::string>> &tiny_out,
    size_t               &num_seqs_out)
{
  RbGzReader fgz;
  if (!fgz.open(path.c_str(), input_is_gz)) { fgz.close(); return false; }

  // One producer thread (inflate) + the rest parse.  At least one parser.
  const int n_consumers = std::max(1, nthreads - 1);

  // Per-consumer persistent sequence arenas (retained for the whole run; the
  // merged seq views point into them).  Indexed by consumer slot.
  g_seq_arenas.assign(n_consumers, std::string());

  struct Block { std::vector<char> data; size_t id; };
  struct BlockOut {
    std::vector<ContigRec> large, small;
    std::vector<std::pair<std::string,std::string>> tiny;
    size_t nseq = 0;
  };

  // ── Bounded block queue (producer → consumers) ──────────────────────────
  std::mutex                 qmtx;
  std::condition_variable    q_notempty, q_notfull;
  std::deque<Block>          queue;
  bool                       producer_done = false;
  const size_t MAX_INFLIGHT = (size_t)std::max(4, n_consumers * 3);

  // ── Per-block parsed results (results[id]; file order == block id order) ─
  std::mutex                 rmtx;
  std::deque<BlockOut>       results;

  // ── Producer: inflate, then split on the last "\n>" so every emitted block
  //    is a whole number of records (the trailing partial record is carried
  //    over to the next block). ───────────────────────────────────────────
  std::thread producer([&]() {
    const size_t TARGET = 64u << 20;   // ~64 MB decompressed per block
    const size_t CHUNK  = 4u << 20;    // gzread granularity
    std::vector<char> carry;           // partial trailing record (starts at '>')
    size_t block_id = 0;
    bool   eof = false;

    while (true) {
      std::vector<char> buf;
      buf.swap(carry);                 // begins with '>' (or empty on first pass)
      if (buf.capacity() < TARGET + CHUNK) buf.reserve(TARGET + CHUNK);

      // Accumulate until we have >= TARGET bytes AND a record boundary, or EOF.
      // A single record larger than TARGET simply grows the block until the
      // next boundary appears (or EOF), so no record is ever split.
      for (;;) {
        bool boundary = buf.size() >= TARGET &&
            std::string_view(buf.data(), buf.size())
                .rfind(std::string_view("\n>", 2)) != std::string_view::npos;
        if (boundary || eof) break;
        size_t old = buf.size();
        buf.resize(old + CHUNK);
        int64_t got = fgz.read(buf.data() + old, CHUNK);
        if (got <= 0) { buf.resize(old); eof = true; }  // 0 = EOF, <0 = error
        else          buf.resize(old + (size_t)got);
      }

      if (buf.empty()) break;          // nothing left to emit

      size_t emit_len;
      if (eof) {
        emit_len = buf.size();         // last block: all records are complete
      } else {
        size_t cut = std::string_view(buf.data(), buf.size())
                         .rfind(std::string_view("\n>", 2));
        // Non-EOF exit implies a boundary was found above.
        emit_len = cut + 1;            // keep the '\n'; the next '>' starts carry
        carry.assign(buf.begin() + emit_len, buf.end());
        buf.resize(emit_len);
      }

      Block blk;
      blk.id   = block_id++;
      blk.data = std::move(buf);
      {
        std::unique_lock<std::mutex> lk(qmtx);
        q_notfull.wait(lk, [&] { return queue.size() < MAX_INFLIGHT; });
        queue.push_back(std::move(blk));
      }
      q_notempty.notify_one();

      if (eof) break;
    }

    {
      std::lock_guard<std::mutex> lk(qmtx);
      producer_done = true;
    }
    q_notempty.notify_all();
  });

  // ── Consumers: parse whole-record blocks into thread-local ContigRecs ────
  auto consumer = [&](int cidx) {
    std::string          &arena = g_seq_arenas[cidx];
    std::vector<uint64_t> scratch;   // build_pmh_winners scratch (reused)

    for (;;) {
      Block blk;
      {
        std::unique_lock<std::mutex> lk(qmtx);
        q_notempty.wait(lk, [&] { return !queue.empty() || producer_done; });
        if (queue.empty()) break;    // producer_done && drained
        blk = std::move(queue.front());
        queue.pop_front();
      }
      q_notfull.notify_one();

      BlockOut out;
      const char *p   = blk.data.data();
      const char *end = p + blk.data.size();

      while (p < end) {
        const char *gt = (const char *)memchr(p, '>', (size_t)(end - p));
        if (!gt) break;
        p = gt + 1;

        const char *name_start = p;
        while (p < end && (uint8_t)*p > ' ') ++p;
        const char *name_end = p;
        const char *nl = (const char *)memchr(p, '\n', (size_t)(end - p));
        if (!nl) break;              // header not terminated (block boundary)
        p = nl + 1;

        const char *next_gt = (const char *)memchr(p, '>', (size_t)(end - p));
        const char *seq_end = next_gt ? next_gt : end;

        ++out.nseq;

        // Fast tiny-prune: the raw byte span (incl. newlines) is an upper bound
        // on seq_len, so a sub-threshold span is definitely tiny — skip the
        // compaction + name materialisation entirely (unless tiny output is on).
        if (!collect_tiny_arg &&
            (size_t)(seq_end - p) < (size_t)min_small_contig_arg) {
          p = seq_end;
          continue;
        }

        // Compact the (possibly multi-line) sequence into this worker's arena,
        // stripping '\n' and trailing '\r'.  Recorded as an (arena_tid,off,len)
        // slice; rolled back below if the record turns out to be discarded.
        const size_t off = arena.size();
        const char  *seg = p;
        while (seg < seq_end) {
          const char *seg_nl = (const char *)memchr(seg, '\n', (size_t)(seq_end - seg));
          const char *seg_e  = seg_nl ? seg_nl : seq_end;
          size_t slen = (size_t)(seg_e - seg);
          if (slen > 0 && seg[slen - 1] == '\r') --slen;
          if (slen > 0) arena.append(seg, slen);
          if (!seg_nl) break;
          seg = seg_nl + 1;
        }
        const size_t seq_len = arena.size() - off;
        p = seq_end;

        const bool is_large = seq_len >= minContig_arg;
        const bool is_small = !is_large && seq_len >= min_small_contig_arg;

        if (!is_large && !is_small) {
          if (collect_tiny_arg) {
            std::string name(name_start, name_end);
            if (!fullHeader_arg) {
              size_t sp = name.find_first_of(" \t");
              if (sp != std::string::npos) name.resize(sp);
            }
            out.tiny.emplace_back(std::move(name),
                                  std::string(arena.data() + off, seq_len));
          }
          arena.resize(off);          // discard tiny bytes
          continue;
        }

        std::string name(name_start, name_end);
        if (!fullHeader_arg) {
          size_t sp = name.find_first_of(" \t");
          if (sp != std::string::npos) name.resize(sp);
        }
        const char *seq_bytes = arena.data() + off;  // stable until next append

        if (is_large) {
          ContigRec rec;
          rec.name     = std::move(name);
          rec.len      = seq_len;
          rec.is_small = false;
          if (stream_pmh_arg) {
            rec.winners.resize(pmh_m_arg, 0u);
            float *k4out = nullptr;
            if ((g_idf_norm || g_exact_cos_cmp) && pmh_k_arg == 4 && g_gc_norm >= 1) {
              rec.k4freq.resize(256, 0.0f);
              k4out = rec.k4freq.data();
            }
            build_pmh_winners(seq_bytes, seq_len,
                              pmh_k_arg, pmh_m_arg, /*seed=*/42u,
                              rec.winners.data(), scratch,
                              /*out64=*/nullptr, /*out_keys=*/nullptr, k4out);
          }
          if (!no_store_seqs_arg) {
            rec.arena_tid = cidx;     // arena slice retained
            rec.arena_off = off;
          } else {
            arena.resize(off);        // sketch built; drop the stored bytes
            rec.arena_tid = -1;
          }
          out.large.push_back(std::move(rec));
        } else {                       // is_small (always stored, matches mmap)
          ContigRec rec;
          rec.name      = std::move(name);
          rec.len       = seq_len;
          rec.is_small  = true;
          rec.arena_tid = cidx;
          rec.arena_off = off;
          out.small.push_back(std::move(rec));
        }
      }

      {
        std::lock_guard<std::mutex> lk(rmtx);
        if (results.size() <= blk.id) results.resize(blk.id + 1);
        results[blk.id] = std::move(out);
      }
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(n_consumers);
  for (int c = 0; c < n_consumers; ++c) workers.emplace_back(consumer, c);
  producer.join();
  for (auto &w : workers) w.join();
  fgz.close();
  rb_phase("  streaming parallel parse loop done");

  // ── Merge per-block results in file (block id) order ────────────────────
  num_seqs_out = 0;
  size_t n_large = 0, n_small = 0, n_tiny = 0;
  for (auto &b : results) { n_large += b.large.size(); n_small += b.small.size(); n_tiny += b.tiny.size(); }
  large_out.reserve(n_large);
  small_out.reserve(n_small);
  tiny_out.reserve(n_tiny);
  for (auto &b : results) {
    num_seqs_out += b.nseq;
    for (auto &r : b.large) large_out.push_back(std::move(r));
    for (auto &r : b.small) small_out.push_back(std::move(r));
    for (auto &r : b.tiny)  tiny_out.push_back(std::move(r));
  }
  rb_phase("  streaming parse merge done");
  return true;
}

// Decompress a gzip file with libdeflate into out_buf.
// Returns true on success; out_buf is NUL-terminated.
// ── Parallel multi-stream gzip decompression ──────────────────────────────
// A gzip file can contain multiple independent concatenated streams (e.g.,
// files created by `cat a.gz b.gz > ab.gz` or pigz -p1 output).  Each stream
// is independently decompressible, so we can dispatch N streams to N threads
// for an up-to-N× speedup.
//
// Algorithm:
//   1. Scan the compressed buffer for valid gzip magic+header bytes to locate
//      stream boundaries.  A valid header has magic {0x1f,0x8b}, CM=8, and
//      FLG bits 5–7 = 0.  False-positive candidates are accepted here but
//      rejected in step 3.
//   2. Read ISIZE from each stream's footer (last 4 bytes before the next
//      stream starts).  Pre-allocate the output buffer and output offsets.
//   3. Decompress each stream in parallel with OpenMP.  libdeflate validates
//      the CRC32, so any false-positive boundary causes a decompressor error.
//   4. If all streams succeed, compact any gaps caused by ISIZE over-estimates
//      and return.  If any stream fails, fall back to a single-threaded
//      sequential decompress of the entire file.
// ─────────────────────────────────────────────────────────────────────────────
static bool libdeflate_read_gz(const std::string& path,
                                std::vector<char>& out_buf) {
  FILE* f = fopen(path.c_str(), "rb");
  if (!f) return false;
  fseek(f, 0, SEEK_END);
  const size_t csz = (size_t)ftell(f);
  fseek(f, 0, SEEK_SET);
  std::vector<uint8_t> cbuf(csz);
  if (fread(cbuf.data(), 1, csz, f) != csz) { fclose(f); return false; }
  fclose(f);
  if (csz < 18) return false;

  // ── Step 1: Find candidate stream starts ─────────────────────────────
  // Three-byte filter: magic {0x1f,0x8b}, CM=8, valid FLG (bits 5-7 clear).
  // This gives ~1 false positive per 134M bytes in random data.  False
  // positives are safely rejected by the CRC check in step 3.
  std::vector<size_t> starts;
  for (size_t i = 0; i + 10 <= csz; i++) {
    if (cbuf[i]   == 0x1f && cbuf[i+1] == 0x8b &&
        cbuf[i+2] == 0x08 && (cbuf[i+3] & 0xE0) == 0) {
      starts.push_back(i);
    }
  }
  if (starts.empty()) return false;
  const int nstreams = (int)starts.size();

  // ── Step 2: Validate ISIZE for every candidate stream ────────────────
  // ISIZE (uncompressed size mod 2^32) is the last 4 bytes of each stream.
  // For a true stream, ISIZE / in_sz falls in a narrow range: typical FASTA
  // gzip compresses 3–5×, so ISIZE ≈ 3–5 × in_sz.  We allow [0.5×, 12×].
  //
  // If ANY candidate stream fails this check the scan has hit false positives
  // (magic bytes occurring by chance inside compressed data).  In that case
  // we cannot safely determine stream boundaries, so we skip the parallel
  // attempt entirely and fall through to single-threaded decompression.
  // This is the common case for standard single-stream gzip; the parallel
  // path only activates for files genuinely composed of N independent streams
  // (e.g., produced by "cat a.gz b.gz …" or multi-stream pigz output).
  std::vector<size_t> out_sizes(nstreams);
  size_t total_out = 0;
  {
    bool any_bad = (nstreams < 2);      // single-stream: no parallel benefit
    for (int i = 0; i < nstreams && !any_bad; i++) {
      size_t end   = (i + 1 < nstreams) ? starts[i+1] : csz;
      size_t in_sz = end - starts[i];
      uint32_t isize = 0;
      if (end >= 4) std::memcpy(&isize, cbuf.data() + end - 4, 4);
      // Plausible decompression ratio: 0.5× – 12× of compressed size
      const size_t lo = in_sz / 2, hi = in_sz * 12;
      if ((size_t)isize < lo || (size_t)isize > hi) {
        any_bad = true;
      } else {
        out_sizes[i] = (size_t)isize;
        total_out   += out_sizes[i];
      }
    }
    if (any_bad) goto single_thread_fallback;
  }

  {
  std::vector<size_t> out_off(nstreams);
  out_off[0] = 0;
  for (int i = 1; i < nstreams; i++)
    out_off[i] = out_off[i-1] + out_sizes[i-1];

  out_buf.resize(total_out + 1);

  // ── Step 3: Parallel decompress ───────────────────────────────────────
  std::vector<size_t> actual_outs(nstreams, 0);
  std::vector<int>    stream_ok(nstreams, 1);  // int for OMP atomic

  #pragma omp parallel for schedule(static,1) \
      num_threads(std::min(nstreams, omp_get_max_threads()))
  for (int i = 0; i < nstreams; i++) {
    size_t in_start = starts[i];
    size_t in_end   = (i + 1 < nstreams) ? starts[i+1] : csz;
    struct libdeflate_decompressor* dc = libdeflate_alloc_decompressor();
    size_t actual = 0;
    libdeflate_result r = libdeflate_gzip_decompress(
        dc, cbuf.data() + in_start, in_end - in_start,
        out_buf.data() + out_off[i], out_sizes[i], &actual);
    libdeflate_free_decompressor(dc);
    if (r == LIBDEFLATE_SUCCESS) {
      actual_outs[i] = actual;
    } else {
      stream_ok[i] = 0;
    }
  }

  // Check all streams succeeded
  bool all_ok = true;
  for (int i = 0; i < nstreams; i++) if (!stream_ok[i]) { all_ok = false; break; }

  if (all_ok) {
    // ── Step 4: Compact – close gaps from ISIZE over-estimates ──────────
    // Each stream's data is at out_off[i]; actual length is actual_outs[i].
    // We memmove each stream's block backward to pack them contiguously.
    // All moves are backward (src >= dst), valid for memmove.
    size_t write_pos = actual_outs[0];  // stream 0 already at offset 0
    for (int i = 1; i < nstreams; i++) {
      if (out_off[i] != write_pos)
        std::memmove(out_buf.data() + write_pos,
                     out_buf.data() + out_off[i], actual_outs[i]);
      write_pos += actual_outs[i];
    }
    out_buf.resize(write_pos + 1);
    out_buf[write_pos] = '\0';
    return true;
  }
  } // end parallel-path scope

  // ── Fallback: single-threaded sequential decompress ───────────────────
  // Handles: single-stream files, false-positive stream boundaries detected
  // by the ISIZE ratio check or failed CRC validation in parallel step.
  out_buf.clear();
  single_thread_fallback:
  {
    uint32_t isize_field;
    std::memcpy(&isize_field, cbuf.data() + csz - 4, 4);
    size_t out_sz = (size_t)isize_field;
    if (out_sz < csz) out_sz = csz * 5;
    out_buf.resize(out_sz + 1);
    struct libdeflate_decompressor* dc = libdeflate_alloc_decompressor();
    size_t actual = 0;
    libdeflate_result r = libdeflate_gzip_decompress(
        dc, cbuf.data(), csz, out_buf.data(), out_sz, &actual);
    libdeflate_free_decompressor(dc);
    if (r != LIBDEFLATE_SUCCESS) { out_buf.clear(); return false; }
    out_buf.resize(actual + 1);
    out_buf[actual] = '\0';
    return true;
  }
}

// ═══════════════════════════════════════════════════════════════════════════
// Async depth_matrix (abundance/depth) pre-parser
// ─────────────────────────────────────────────────────────────────────────
// Parses the depth/abundance file into a hash map (name → coverage values)
// independently of FASTA parsing.  This lets the two I/O phases run in
// parallel: the main thread decompresses and parses FASTA while this runs
// on a background thread.  After FASTA is done the caller merges the result.
// ─────────────────────────────────────────────────────────────────────────
struct RawDepthEntry {
  std::vector<float> means;
  std::vector<float> vars;
};

// parse_depth_async – reads depth_file and returns a name→RawDepthEntry map.
// Called via std::async before FASTA decompression starts.
//
// Parallel mmap parser: the depth matrix can be huge (e.g. 244 MB, 1.47 M rows
// × 43 float columns for CAMI human_gut) and the old line-by-line stringstream
// + stod-with-exceptions path was a ~4 s SERIAL bottleneck. We now mmap the
// file, split it into per-thread line-aligned chunks, and parse each field with
// strtof (no stringstream, no exceptions, no locale). Field separators are tabs;
// strtof stops cleanly at the tab/newline that follows each number.
// Core parser over an in-memory TSV buffer [buf, buf+fsz). Shared by the file
// path (parse_depth_async, which mmaps) and the fused in-process path (which
// passes the depth string returned by compute_depth_tsv_inmem). Identical
// parsing for both, so binning results do not depend on whether depth came from
// a file or memory.
static phmap::flat_hash_map<std::string, RawDepthEntry>
parse_depth_buf(const char* buf, size_t fsz, bool cvExt, int num_depth_samples,
                bool fullHeader, int nthreads) {
  using Map = phmap::flat_hash_map<std::string, RawDepthEntry>;
  Map result;
  if (nthreads < 1) nthreads = 1;
  if (!buf || fsz == 0) return result;

  const char* data_end = buf + fsz;
  // Skip header line: data starts after the first newline.
  const char* nl0 = (const char*)memchr(buf, '\n', fsz);
  const char* data_start = nl0 ? nl0 + 1 : data_end;

  // Split [data_start, data_end) into nthreads ranges, each starting on a line
  // boundary so no row is parsed twice or split across threads.
  std::vector<const char*> bounds(nthreads + 1);
  bounds[0] = data_start;
  bounds[nthreads] = data_end;
  for (int t = 1; t < nthreads; ++t) {
    const char* approx =
        data_start + (size_t)((double)(data_end - data_start) * t / nthreads);
    if (approx >= data_end) { bounds[t] = data_end; continue; }
    const char* nl = (const char*)memchr(approx, '\n', (size_t)(data_end - approx));
    bounds[t] = nl ? nl + 1 : data_end;
  }

  std::vector<std::vector<std::pair<std::string, RawDepthEntry>>> tl(nthreads);
  const bool depth_prefilter = rb_env_depth_prefilter_on();
  const size_t min_depth_len = min_small_contig;
  // MetaBAT depth files interleave mean/var pairs; vars are unused in the main
  // pipeline (depth_var_matrix is never filled). Skip strtof on var columns unless
  // RABBIT_DEPTH_PARSE_VAR=1 (for cal_depth_dist revival / debugging).
  const bool parse_var = cvExt ? false : rb_env_depth_parse_var();

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
  for (int t = 0; t < nthreads; ++t) {
    const char* p   = bounds[t];
    const char* end = bounds[t + 1];
    if (p >= end) continue;
    auto& out = tl[t];
    out.reserve((size_t)((end - p) / 64));  // rough row-size estimate

    while (p < end) {
      const char* eol = (const char*)memchr(p, '\n', (size_t)(end - p));
      const char* line_end = eol ? eol : end;
      const char* q = p;
      p = eol ? eol + 1 : end;
      if (line_end <= q) continue;  // empty line

      // name: up to first tab
      const char* tabp = (const char*)memchr(q, tab_delim, (size_t)(line_end - q));
      if (!tabp) continue;
      std::string name(q, tabp);
      if (!fullHeader) trim_fasta_label(name);
      const char* c = tabp + 1;

      // contigLen: parse before depth columns so short rows can be skipped early.
      tabp = (const char*)memchr(c, tab_delim, (size_t)(line_end - c));
      if (!tabp) continue;
      char* ep = nullptr;
      // contigLen may be written in float / scientific notation (jgi depth files
      // emit megabase lengths as e.g. "2.52094e+06").  strtod consumes the whole
      // field (ep == tabp) where strtoul would stop at the '.' and wrongly drop
      // the row — silently discarding the longest contigs (often near-complete
      // genomes), which is why CAMI-high regressed.  strtod fixes that while
      // keeping the "entire field consumed" malformed-row guard.
      double contig_len_d = strtod(c, &ep);
      if (ep != tabp) continue;
      unsigned long contig_len = (unsigned long)(contig_len_d < 0 ? 0 : contig_len_d);
      c = tabp + 1;
      if (depth_prefilter && contig_len < min_depth_len) continue;

      // totalAvgDepth: skip when not in cvExt mode
      if (!cvExt) {
        tabp = (const char*)memchr(c, tab_delim, (size_t)(line_end - c));
        if (!tabp) continue;
        c = tabp + 1;
      }

      RawDepthEntry entry;
      entry.means.resize(num_depth_samples, 0.f);
      if (parse_var) entry.vars.resize(num_depth_samples, 0.f);

      auto next_field = [&](const char* cur) -> const char* {
        const char* tp = (const char*)memchr(cur, tab_delim, (size_t)(line_end - cur));
        return tp ? tp + 1 : line_end;
      };
      for (int i = 0; i < num_depth_samples; ++i) {
        if (c >= line_end) break;
        char* ep = nullptr;
        entry.means[i] = strtof(c, &ep);  // stops at the trailing tab/newline
        c = next_field(c);
        if (!cvExt) {
          if (c >= line_end) break;
          if (parse_var)
            entry.vars[i] = strtof(c, &ep);
          c = next_field(c);
        }
      }
      out.emplace_back(std::move(name), std::move(entry));
    }
  }

  size_t total = 0;
  for (auto& v : tl) total += v.size();
  result.reserve(total);
  for (auto& v : tl) {
    for (auto& kv : v) result.emplace(std::move(kv.first), std::move(kv.second));
    std::vector<std::pair<std::string, RawDepthEntry>>().swap(v);
  }
  return result;
}

// File path: mmap the depth file then parse with the shared buffer core.
static phmap::flat_hash_map<std::string, RawDepthEntry>
parse_depth_async(const std::string& depth_file, bool cvExt, int num_depth_samples, bool fullHeader,
                int nthreads) {
  using Map = phmap::flat_hash_map<std::string, RawDepthEntry>;
  int fd = ::open(depth_file.c_str(), O_RDONLY);
  if (fd < 0) return Map();
  struct stat st;
  if (fstat(fd, &st) != 0 || st.st_size <= 0) { ::close(fd); return Map(); }
  const size_t fsz = (size_t)st.st_size;
  const char* buf = (const char*)mmap(nullptr, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
  ::close(fd);
  if (buf == MAP_FAILED) return Map();
  madvise((void*)buf, fsz, MADV_SEQUENTIAL);
  Map result = parse_depth_buf(buf, fsz, cvExt, num_depth_samples, fullHeader, nthreads);
  munmap((void*)buf, fsz);
  return result;
}

// ═══════════════════════════════════════════════════════════════════════════
// Marker-guided bin splitting (Phase 2) — C++ port of post_split.py
// ═══════════════════════════════════════════════════════════════════════════
// Lightweight dependency-free KMeans (k-means++ init, Lloyd iterations, a few
// restarts, best inertia). Features are small (n ≤ few hundred, dim = num_depth_samples,
// k ≤ splitMaxK) so this is microseconds per bin.
static std::vector<int> rb_kmeans(const std::vector<std::vector<float>> &X,
                                    int k, std::mt19937 &rng) {
  const size_t n = X.size();
  const size_t d = X.empty() ? 0 : X[0].size();
  std::vector<int> best_labels(n, 0);
  double best_inertia = std::numeric_limits<double>::infinity();
  if (n == 0 || d == 0 || k <= 1) return best_labels;

  for (int init = 0; init < 10; ++init) {
    // ── k-means++ seeding ──
    std::vector<std::vector<float>> cent(k);
    std::uniform_int_distribution<size_t> u0(0, n - 1);
    cent[0] = X[u0(rng)];
    std::vector<double> d2(n, std::numeric_limits<double>::infinity());
    for (int c = 1; c < k; ++c) {
      for (size_t r = 0; r < n; ++r) {
        double s = 0;
        for (size_t i = 0; i < d; ++i) { double t = X[r][i] - cent[c-1][i]; s += t*t; }
        if (s < d2[r]) d2[r] = s;
      }
      double sum = 0; for (double v : d2) sum += v;
      size_t pick = n - 1;
      if (sum > 0) {
        std::uniform_real_distribution<double> ur(0, sum);
        double thr = ur(rng), acc = 0;
        for (size_t r = 0; r < n; ++r) { acc += d2[r]; if (acc >= thr) { pick = r; break; } }
      } else {
        pick = u0(rng);
      }
      cent[c] = X[pick];
    }
    // ── Lloyd iterations ──
    std::vector<int> labels(n, 0);
    for (int iter = 0; iter < 200; ++iter) {
      bool changed = false;
      for (size_t r = 0; r < n; ++r) {
        double bd = std::numeric_limits<double>::infinity(); int bl = 0;
        for (int c = 0; c < k; ++c) {
          double s = 0;
          for (size_t i = 0; i < d; ++i) { double t = X[r][i] - cent[c][i]; s += t*t; }
          if (s < bd) { bd = s; bl = c; }
        }
        if (labels[r] != bl) { labels[r] = bl; changed = true; }
      }
      std::vector<std::vector<double>> sum(k, std::vector<double>(d, 0.0));
      std::vector<int> cnt(k, 0);
      for (size_t r = 0; r < n; ++r) {
        cnt[labels[r]]++;
        for (size_t i = 0; i < d; ++i) sum[labels[r]][i] += X[r][i];
      }
      for (int c = 0; c < k; ++c)
        if (cnt[c]) for (size_t i = 0; i < d; ++i) cent[c][i] = (float)(sum[c][i] / cnt[c]);
      if (!changed) break;
    }
    double inertia = 0;
    for (size_t r = 0; r < n; ++r) {
      double s = 0;
      for (size_t i = 0; i < d; ++i) { double t = X[r][i] - cent[labels[r]][i]; s += t*t; }
      inertia += s;
    }
    if (inertia < best_inertia) { best_inertia = inertia; best_labels = labels; }
  }
  return best_labels;
}

// Re-split bins whose single-copy markers are duplicated (chimeric), using
// per-sample log-abundance KMeans. Operates on the final BinMap in place.

// ── split helpers (before main) ──────────────────────────────────────────
#include "impl/rb_split.cpp"

// ── cache I/O (before main; used by both the save hook and the load path) ──
#include "impl/rb_cache.cpp"

// bin subcommand (also the legacy default entry point)
// ═══════════════════════════════════════════════════════════════════════════
static int rb_cmd_bin(int ac, char *av[]) {
  po::options_description desc("RabbitBin options", 100, 50);
  desc.add_options()
      ("help,h", "Show help")
      ("assembly,a", po::value<std::string>(&inFile), "Contig FASTA assembly (gzip ok) [required]")
      ("fasta", po::value<std::string>(&inFile), "Alias for --assembly")
      ("output,o", po::value<std::string>(&outFile), "Output path prefix [required]")
      ("depth,d", po::value<std::string>(&depth_file), "Coverage depth TSV")
      ("save-cache", po::value<std::string>(&cache_save_file), "Write a re-binning cache after the similarity graph is built, then continue normally")
      ("load-cache", po::value<std::string>(&cache_load_file), "Re-bin from a cache (skips parse/sketch/depth/graph); only cheap params (min-bin-size, split-silhouette, min-edge-score, ...) take effect")
      ("min-contig,m", po::value<size_t>(&minContig)->default_value(2500), "Minimum contig length (>=1500)")
      ("min-small-contig", po::value<size_t>(&min_small_contig)->default_value(1000), "Min length for small-contig recruiting")
      ("max-posterior", po::value<Similarity>(&calib_connected_pct)->default_value(95), "Well-connected contig percent for calibration")
      ("min-edge-score", po::value<Similarity>(&min_edge_weight)->default_value(60), "Minimum edge weight (1-99)")
      ("max-edges", po::value<size_t>(&maxEdges)->default_value(200), "Max neighbors per contig")
      ("sim-cutoff", po::value<Similarity>(&simCutoff)->default_value(0), "Composition similarity cutoff x100 (0=auto)")
      ("sketch-k", po::value<int>(&sketch_kmer_size)->default_value(8), "Sketch k-mer size")
      ("sketch-m", po::value<uint32_t>(&sketch_size)->default_value(500), "Sketch size (PMH registers)")
      ("sketch-b", po::value<uint32_t>(&sketch_bits)->default_value(2), "MinHash bucket bits")
      ("no-recruit", po::value<bool>(&no_recruit)->zero_tokens(), "Disable small-contig recruiting")
      ("no_gold", po::value<bool>(&no_gold)->zero_tokens(), "Label-free multi-resolution: sweep alpha/edge_power on the reused graph, auto-select max-modularity partition (no ground truth needed)")
      ("min-recruit-cluster", po::value<size_t>(&minCS)->default_value(10), "Min cluster size for recruiting")
      ("recruit-abd-centroid", po::value<bool>(&recruit_to_depth_centroid)->default_value(false)->zero_tokens(), "Recruit using abundance centroid")
      ("recruit-cutoff", po::value<Distance>(&recruitSimFactor)->default_value(0.0), "Recruit sim factor x sim-cutoff (0=off)")
      ("depth-no-variance", po::value<bool>(&cvExt)->zero_tokens(), "Depth file has no variance columns")
      ("full-header", po::value<bool>(&fullHeader)->zero_tokens(), "Keep full FASTA headers")
      ("min-coverage,x", po::value<Distance>(&minCV)->default_value(1), "Min per-sample mean coverage")
      ("min-coverage-sum", po::value<Distance>(&minCVSum)->default_value(1), "Min total mean coverage")
      ("min-bin-size,s", po::value<size_t>(&min_bin_bp)->default_value(50000), "Min output bin size (bp)")
      ("threads,t", po::value<size_t>(&numThreads)->default_value(0), "Threads (0=all online CPUs)")
      ("labels-only,l", po::value<bool>(&onlyLabel)->zero_tokens(), "Output contig names only")
      ("save-matrix", po::value<bool>(&saveCls)->zero_tokens(), "Save membership matrix")
      ("unbinned", po::value<bool>(&outUnbinned)->zero_tokens(), "Write unbinned FASTA")
      ("no-bin-fasta", po::value<bool>(&noBinOut)->zero_tokens(), "Skip per-bin FASTA output (default)")
      ("bin-fasta", po::value<bool>(&binFastaWanted)->zero_tokens(), "Also write per-bin FASTA files (keeps sequences in RAM)")
      ("no-sample-depths", po::value<bool>(&noSampleDepths)->zero_tokens(), "Omit per-sample depths in headers")
      ("seed", po::value<unsigned long long>(&seed)->default_value(0), "Random seed (0=time)")
      ("marker-seed", po::value<std::string>(&marker_seed_file)->default_value(""), "Optional marker seed file")
      ("split-max-k", po::value<int>(&splitMaxK)->default_value(6), "Max sub-clusters per split bin")
      ("split-bins", po::value<bool>(&g_split_abundance)->zero_tokens(), "Abundance bin splitting (default ON)")
      ("no-split", po::value<bool>(&g_no_split_abundance)->zero_tokens(), "Disable abundance splitting")
      ("split-silhouette", po::value<double>(&g_split_sil)->default_value(rb_env_split_sil()), "Silhouette split threshold")
      ("auto", po::value<bool>(&g_auto_select)->zero_tokens(), "Build the graph once, sweep configs, and auto-select the partition with the most near-complete bins (needs --markers)")
      ("ensemble", po::value<bool>(&g_ensemble)->zero_tokens(), "[experimental] Consensus over the swept configs: greedily keep the highest-quality non-overlapping bins (needs --markers). Note: the swept configs share one graph and are highly correlated, so on well-tuned data this rarely beats --auto; intended for diverse/independent partitions.")
      ("markers", po::value<std::string>(&g_markers_file), "Contig->marker map for --auto/--ensemble (from scripts/rabbitbin_markers.sh)")
      ("debug", po::value<bool>(&debug)->zero_tokens(), "Debug output")
      ("quiet,q", po::value<bool>(&quiet)->zero_tokens(), "Less verbose")
#ifdef RABBITBIN_FUSE
      ("bam", po::value<std::vector<std::string>>(&fuse_bams)->multitoken(),
       "Sorted BAM(s); depth is computed in-process (overlaps FASTA parse). "
       "Mutually exclusive with --depth. BAMs may also be given positionally.")
      ("bam-list", po::value<std::string>(&fuse_bam_list),
       "File with one sorted-BAM path per line (appended to --bam)")
      ("percent-identity", po::value<int>(&fuse_pctid)->default_value(97),
       "[--bam] Min mapped-read percent identity")
      ("min-contig-length", po::value<int>(&fuse_min_contig_len)->default_value(1),
       "[--bam] Min contig length emitted into depth (rabbit_depth default; "
       "rabbitbin applies its own --min-contig / small-contig filters downstream)")
      ("min-contig-depth", po::value<double>(&fuse_min_contig_depth)->default_value(0.0),
       "[--bam] Min contig depth emitted into depth (rabbit_depth default)")
      ("max-edge-bases", po::value<int>(&fuse_max_edge)->default_value(75),
       "[--bam] Max edge bases trimmed per contig end")
#endif
      ("verbose,v", po::value<bool>(&verbose)->zero_tokens(), "Verbose progress (default ON)");

  po::variables_map vm;
#ifdef RABBITBIN_FUSE
  // Allow `rabbitbin -a asm -o out bam1 bam2 ...` (BAMs positional → fuse path).
  po::positional_options_description pos;
  pos.add("bam", -1);
  po::store(po::command_line_parser(ac, av).options(desc).positional(pos).run(), vm);
#else
  po::store(po::command_line_parser(ac, av).options(desc).positional({}).run(), vm);
#endif
  po::notify(vm);

#ifdef RABBITBIN_FUSE
  // --bam-list: append each non-blank line (a BAM path) to the BAM vector so the
  // fuse path treats it identically to repeated --bam / positional args.
  if (!fuse_bam_list.empty()) {
    std::ifstream bl(fuse_bam_list);
    if (!bl) { cerr << "[Error!] cannot open --bam-list file: " << fuse_bam_list << "\n"; return 1; }
    std::string ln;
    while (std::getline(bl, ln)) {
      size_t b = ln.find_first_not_of(" \t\r\n");
      if (b == std::string::npos) continue;
      size_t e = ln.find_last_not_of(" \t\r\n");
      fuse_bams.push_back(ln.substr(b, e - b + 1));
    }
  }
#endif

  // Whether the user explicitly passed --seed (vs the default 0). In cache-load
  // mode this decides whether to reproduce the cached run's seed or honour an
  // override supplied on the rebin command line.
  const bool user_set_seed = vm.count("seed") && !vm["seed"].defaulted();

  // Per-bin FASTA is off by default; --bin-fasta re-enables it (and so keeps
  // sequences resident for the output stage).
  if (binFastaWanted) noBinOut = false;

  // In --load-cache mode the assembly is not re-read; only the output prefix is
  // required (contig names/lengths/depth/graph all come from the cache).
  const bool need_assembly = cache_load_file.empty();
  if (vm.count("help") || (need_assembly && inFile.empty()) || outFile.empty()) {
    cerr << "\nRabbitBin: sketch-based metagenome binning "
            "(version " << version << "; " << DATE << ")\n\n";
    cerr << desc << "\n\n";
    if (!vm.count("help")) {
      if (need_assembly && inFile.empty())
        cerr << "[Error!] --assembly is required\n";
      if (outFile.empty()) cerr << "[Error!] --output is required\n";
    }
    return vm.count("help") ? 0 : 1;
  }

  // ── Cache-mode mutual exclusions ──────────────────────────────────────────
  if (!cache_load_file.empty()) {
    if (!cache_save_file.empty()) {
      cerr << "[Error!] --load-cache and --save-cache are mutually exclusive\n";
      return 1;
    }
    if (!depth_file.empty()) {
      cerr << "[Error!] --load-cache and --depth are mutually exclusive\n";
      return 1;
    }
#ifdef RABBITBIN_FUSE
    if (!fuse_bams.empty()) {
      cerr << "[Error!] --load-cache and --bam are mutually exclusive\n";
      return 1;
    }
#endif
  }

  if (quiet && verbose)  verbose = false;
  else                   verbose = true;

  if (verbose) gettimeofday(&t1, NULL);

  if (seed == 0) seed = time(0);
  srand(seed);

  if (calib_connected_pct <= 0 || calib_connected_pct >= 100) {
    cerr << "[Error!] calib_connected_pct should be > 0 and < 100\n"; return 1;
  }
  if (min_edge_weight <= 1 || min_edge_weight >= 100) {
    cerr << "[Error!] min_edge_weight should be > 1 and < 100\n"; return 1;
  }
  if (simCutoff < 0 || simCutoff >= 100) {
    cerr << "[Error!] --sim-cutoff should be >= 0 and < 100\n"; return 1;
  }
  if (minContig < 1500) {
    cerr << "[Error!] Contig length < 1500 is not allowed.\n"; return 1;
  }
  if (min_small_contig < 500) {
    cerr << "[Error!] Min small contig length < 500 is not allowed.\n"; return 1;
  }
  if (minCV < 0) {
    cerr << "[Error!] minCV should be non-negative\n"; return 1;
  }
  if (sketch_kmer_size < 1 || sketch_kmer_size > 32) {
    cerr << "[Error!] --sketch-k must be in [1, 32]\n"; return 1;
  }
  if (sketch_size < 2) {
    cerr << "[Error!] --sketch-m must be >= 2\n"; return 1;
  }
  if (g_no_split_abundance) g_split_abundance = false;  // explicit opt-out
  minCVSum = std::max(minCV, minCVSum);

  boost::filesystem::path dir(outFile);
  boost::system::error_code ec;
  if (dir.parent_path().string().length() > 0) {
    if (boost::filesystem::is_regular_file(dir.parent_path())) {
      cerr << "Cannot create directory: " << dir.parent_path().string()
           << ", which exists as a regular file.\n";
      return 1;
    }
    if (!boost::filesystem::is_directory(dir.parent_path()) &&
        !boost::filesystem::create_directory(dir.parent_path(), ec)) {
      cerr << "Cannot create directory: " << dir.parent_path().string()
           << ": " << ec << "\n";
      return 1;
    }
  }

  print_message("RabbitBin (%s) using minContig %d, minCV %2.1f, "
                "minCVSum %2.1f, calib_connected_pct %2.0f%%, min_edge_weight %2.0f, maxEdges %d, "
                "min_bin_bp %d, sketch-k %d, sketch-m %d, seed=%lld\n",
                version.c_str(), minContig, minCV, minCVSum, calib_connected_pct, min_edge_weight,
                maxEdges, min_bin_bp, sketch_kmer_size, sketch_size, seed);

  calib_connected_pct /= 100.;  min_edge_weight /= 100.;

  // Thread count: honor the user's explicit `-t N`, capping only to the number
  // of CPUs physically online on the host. We deliberately IGNORE the OMP-based
  // cap (omp_get_max_threads() echoes a stale OMP_NUM_THREADS, e.g. =32, which
  // would silently shrink `-t 64` to 32 → looks like "64 is slower than 32").
  {
    long onln = sysconf(_SC_NPROCESSORS_ONLN);
    size_t hw = (onln > 0) ? (size_t)onln : (size_t)omp_get_num_procs();
    if (hw == 0) hw = numThreads ? numThreads : 1;
    if (numThreads == 0)
      numThreads = hw;                       // auto: use all online CPUs
    else
      numThreads = std::min(numThreads, hw); // explicit -t: cap to online CPUs
  }
  // Pin the OpenMP default and override any inherited OMP_NUM_THREADS so the
  // num_threads(numThreads) clauses and bare parallel regions all use it.
  omp_set_num_threads((int)numThreads);
  omp_set_dynamic(0);  // disable runtime auto-reduction of the team size
  verbose_message("Executing with %d threads\n", numThreads);

  // ── (KmerSketch needs no LUT pre-computation) ─

  nobs = 0; nobs1 = 0;

  // Dedup + name→row index maps.  phmap::flat_hash_map (vs std::unordered_map)
  // stores entries in a flat array, so build/lookup/teardown are markedly
  // cheaper for the hundreds-of-thousands of contig keys.  had_dup_names stays
  // false unless the FASTA actually contained a repeated header (the common
  // case); when it is false the depth-merge loop can use the contig's position
  // directly instead of re-looking-up the map (idx == ci by construction).
  phmap::flat_hash_map<std::string, size_t> contigs;
  phmap::flat_hash_map<std::string, size_t> small_contigs;
  bool had_dup_names = false;

  const int nNonFeat = cvExt ? 1 : 3;
  bool has_depth = depth_file.length() > 0;
#ifdef RABBITBIN_FUSE
  const bool fuse_mode = !fuse_bams.empty();
  if (fuse_mode) {
    if (has_depth) {
      cerr << "[Error!] --bam and --depth are mutually exclusive\n";
      return 1;
    }
    // In-process depth: variance columns present (mean+var per sample), so this
    // matches a standalone rabbit_depth file with default intraDepthVariance.
    cvExt = false;
    has_depth = true;
    num_depth_samples = (int)fuse_bams.size();
  }
#endif

  // ── Cache-mode dispatch ───────────────────────────────────────────────────
  // from_cache: skip the entire expensive feature-construction region below
  // (parse / depth / sketch / similarity graph) and restore that state from the
  // cache instead, then fall through to the cheap re-binning tail.
  const bool from_cache = !cache_load_file.empty();
  if (from_cache) {
    verbose_message("Loading cache %s [%.1fGb / %.1fGb]\n", cache_load_file.c_str(),
                    getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
    if (!rb_load_cache(cache_load_file)) return 1;
    has_depth = g_cache_has_depth;
    if (!user_set_seed) seed = g_cache_seed;  // reproduce cached run by default
    srand(seed);
    g_inv_one_minus_b0 = (g_pmh_baseline > 0.0 && g_pmh_baseline < 1.0)
                             ? (1.0 / (1.0 - g_pmh_baseline)) : 1.0;
    verbose_message("Cache loaded: nobs=%zu nobs1=%zu samples=%zu edges=%zu "
                    "seed=%llu [%.1fGb / %.1fGb]\n",
                    nobs, nobs1, num_depth_samples, g_cache_from.size(), seed,
                    getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
  } else
  // ── Validate FASTA / depth files ─────────────────────────────────────────
  // depth_future: async pre-parse of the depth_matrix file into a name→values map.
  // Launched immediately after the header scan so it runs in parallel with
  // the FASTA decompression (the main bottleneck).
  {
  using DepthMap = phmap::flat_hash_map<std::string, RawDepthEntry>;
  std::future<DepthMap> depth_future;
  {
#ifdef RABBITBIN_FUSE
    if (fuse_mode) {
      verbose_message("Computing depth in-process from %zu BAM(s) "
                      "(overlaps FASTA parse) [%.1fGb / %.1fGb]\n",
                      fuse_bams.size(), getUsedPhysMem(),
                      getTotalPhysMem() / 1024 / 1024);
      // Compute depth from BAMs AND parse the resulting TSV on the background
      // thread, so the (dominant) BAM decompression runs concurrently with the
      // FASTA decompression + sketch pass below. num_depth_samples already set.
      std::vector<std::string> bams = fuse_bams;
      float pctid = (float)fuse_pctid / 100.0f;
      int mcl = fuse_min_contig_len;
      double mcd = fuse_min_contig_depth;
      int medge = fuse_max_edge;
      bool fH = fullHeader;
      bool cv = cvExt;
      int nds = num_depth_samples;
      int nt = (int)numThreads;
      depth_future = std::async(std::launch::async, [=]() -> DepthMap {
        std::string tsv = compute_depth_tsv_inmem(
            bams, pctid, mcl, (float)mcd, medge,
            /*includeEdgeBases=*/false, /*intraDepthVariance=*/true, nt);
        return parse_depth_buf(tsv.data(), tsv.size(), cv, nds, fH, nt);
      });
    } else
#endif
    if (has_depth) {
      verbose_message("Parsing abundance file header [%.1fGb / %.1fGb]\n",
                      getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
      num_depth_samples = ncols(depth_file.c_str(), 1) - nNonFeat;
      if (!cvExt) {
        if (num_depth_samples % 2 != 0) {
          cerr << "[Error!] Number of columns (excluding the first column) in "
                  "abundance data file is not even.\n";
          exit(1);
        }
        num_depth_samples /= 2;
      }
      // Launch depth_matrix pre-parse on a background thread.  It runs while
      // libdeflate decompresses and parses the FASTA below.
      depth_future = std::async(std::launch::async,
                              parse_depth_async,
                              depth_file, cvExt, num_depth_samples, fullHeader, (int)numThreads);
    }

    rb_phase("parse start");
    verbose_message("Parsing assembly file [%.1fGb / %.1fGb]\n",
                    getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
    size_t num_seqs = 0;

    // ── FASTA decompression: libdeflate (small) vs streaming zlib (large) ─
    bool parsed_ok = false;
    const bool is_gz = (inFile.size() >= 3 &&
                        inFile.compare(inFile.size()-3, 3, ".gz") == 0);

    // Threshold: compressed file must be < RABBIT_LIBDEFLATE_MAXGB GB to use
    // libdeflate one-shot.  Larger files use streaming zlib to avoid OOM.
    // Default 2 GB → decompressed peak ≈ 9-12 GB, safe on most servers.
    // Set RABBIT_LIBDEFLATE_MAXGB=0 to always use streaming zlib.
    size_t libdeflate_max_bytes = 2ULL * 1024 * 1024 * 1024;
    if (const char *ev = rb_getenv("RABBIT_LIBDEFLATE_MAXGB"))
      libdeflate_max_bytes = (size_t)(atof(ev) * 1024.0 * 1024 * 1024);
    const size_t gz_file_bytes = is_gz ? (size_t)filesize(inFile.c_str()) : 0;
    const bool use_libdeflate  = is_gz && (gz_file_bytes <= libdeflate_max_bytes);

    std::ofstream *os = NULL;
    std::string filteredFile_cls;
    if (outUnbinned) {
      filteredFile_cls = outFile + ".";
      filteredFile_cls.append("tooShort");
      if (!onlyLabel) filteredFile_cls.append(".fa");
      os = new std::ofstream(filteredFile_cls.c_str());
      if (!os->is_open() || os->fail() || !*os) {
        cerr << "[Error!] can't open the output bin file: "
             << filteredFile_cls << "\n";
        return 1;
      }
    }

    if (use_libdeflate) {
      std::vector<char> fasta_buf;
      verbose_message("Decompressing %s with libdeflate (%.0fMB) [%.1fGb / %.1fGb]\n",
                      inFile.c_str(), gz_file_bytes / 1048576.0,
                      getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
      if (libdeflate_read_gz(inFile, fasta_buf)) {
        verbose_message("Decompressed to %zuMB, parsing sequences...\n",
                        fasta_buf.size() / 1024 / 1024);
        size_t nskip = 0, nskip1 = 0;
        parse_fasta_buf(fasta_buf.data(), fasta_buf.size() - 1,
          [&](const std::string& raw_name, const char* seq_ptr, size_t seq_len) {
            num_seqs++;
            std::string name = raw_name;
            if (!fullHeader) trim_fasta_label(name);

            if (seq_len >= (size_t)minContig) {
              if (has_depth && contigs.find(name) != contigs.end()) {
                if (debug) verbose_message("Skipping duplicate contig: %s\n", name.c_str());
                nskip++; had_dup_names = true;
                return;
              }
              contigs[name] = nobs + nskip;
              contig_names.push_back(name);
              g_seq_store.emplace_back(seq_ptr, seq_len);   // owned (buffer is transient)
              seqs.push_back(std::string_view(g_seq_store.back()));
              seq_lens.push_back(seq_len);
              nobs++;
            } else if (seq_len >= (size_t)min_small_contig) {
              if (has_depth && small_contigs.find(name) != small_contigs.end()) {
                nskip1++; had_dup_names = true;
                return;
              }
              small_contigs[name] = nobs1 + nskip1;
              small_contig_names.push_back(name);
              g_seq_store.emplace_back(seq_ptr, seq_len);   // owned (buffer is transient)
              small_seqs.push_back(std::string_view(g_seq_store.back()));
              small_seq_lens.push_back(seq_len);
              nobs1++;
            } else {
              if (outUnbinned && os) {
                if (onlyLabel)
                  *os << name << line_delim;
                else
                  printFasta(*os, name, std::string(seq_ptr, seq_len));
              }
            }
          });
        verbose_message("Parsed %zu sequences from FASTA. nobs=%d nobs1=%d\n",
                        num_seqs, nobs, nobs1);
        parsed_ok = true;
      } else {
        verbose_message("libdeflate decompression failed, falling back to streaming zlib\n");
      }
    } else if (is_gz && !use_libdeflate) {
      verbose_message("Large .gz file (%.0fMB > %.0fMB threshold): using streaming zlib "
                      "[set RABBIT_LIBDEFLATE_MAXGB=<N> to adjust]\n",
                      gz_file_bytes / 1048576.0, libdeflate_max_bytes / 1048576.0);
    }

    // ── Parallel ContigRec path ──────────────────────────────────────────
    // Two producers feed the SAME thread-local-ContigRec → file-order-merge
    // machinery (so the two share all downstream code and behave identically):
    //   • uncompressed FASTA  → parse_fasta_mmap_parallel   (shared mmap, N
    //     threads each parse one chunk; OS prefetches pages).
    //   • large .gz (above the libdeflate one-shot threshold) →
    //     parse_fasta_stream_parallel (1 inflate producer + N parse workers).
    //     This replaces the old single-threaded kseq fallback for big gzip
    //     assemblies while leaving the <threshold libdeflate path (CAMI-scale)
    //     and the uncompressed mmap path untouched.
    if (!parsed_ok) {
      const bool stream_pmh_mmap = rb_env_pmh_on();
      const bool no_store_seqs_mmap = [&]() -> bool {
        if (const char *e = rb_getenv("RABBIT_NOSTORE_SEQS")) return e[0] == '1';
        // Auto-drop: the stored contig bytes (≈ the whole assembly — the single
        // largest RAM consumer) are read again ONLY to write per-bin/unbinned
        // FASTA, to feed composition-sketch recruitment, or for the RB_SEQHASH
        // debug.  When PMH winners are built during this parse and none of those
        // consumers is active, the bytes are dead weight, so we skip storing them
        // and cut peak RSS by roughly the input size.  Bit-identical results.
        return stream_pmh_mmap && noBinOut && !outUnbinned &&
               recruitSimFactor <= 0.0 && getenv("RB_SEQHASH") == nullptr;
      }();
      const int pmh_k_mmap = [&]() -> int {
        const char *e = rb_getenv("RABBIT_PMHK"); return e ? std::atoi(e) : 4;
      }();
      const uint32_t pmh_m_mmap = sketch_size;
      // Read precision-tuning globals early so build_pmh_winners (called inline
      // in parse_fasta_mmap_parallel) already sees the correct settings.
      g_gc_norm      = rb_env_gc_norm();
      g_gc_norm_cap  = [] { const char *e = rb_getenv("RABBIT_GC_NORM_CAP");
                            double v = e ? std::atof(e) : 20.0; return (v < 1.0) ? 20.0 : v; }();
      g_idf_norm     = [] { const char *e = rb_getenv("RABBIT_IDF_NORM"); return e && e[0] == '1'; }();
      g_exact_cos_cmp = [] { const char *e = rb_getenv("RABBIT_EXACT_COS");
                             return e && std::atoi(e) >= 1; }();

      if (no_store_seqs_mmap) onlyLabel = true;

      verbose_message("Parallel %s FASTA parse (%d threads) [%.1fGb / %.1fGb]\n",
                      is_gz ? "streaming gzip" : "mmap",
                      (int)numThreads, getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

      std::vector<ContigRec> mmap_large, mmap_small;
      std::vector<std::pair<std::string,std::string>> mmap_tiny;

      bool producer_ok;
      bool used_stream_reader = is_gz;  // gz always uses the streaming reader
      if (is_gz) {
        producer_ok = parse_fasta_stream_parallel(
            inFile, (int)numThreads,
            minContig, min_small_contig,
            fullHeader,
            stream_pmh_mmap, no_store_seqs_mmap,
            /*collect_tiny=*/outUnbinned,
            pmh_m_mmap, pmh_k_mmap,
            /*input_is_gz=*/true,
            mmap_large, mmap_small, mmap_tiny,
            num_seqs);
      } else {
        // Uncompressed.  Prefer the zero-copy mmap reader (fastest — no copy,
        // OS read-ahead, N threads each parse one chunk).  Two escape hatches
        // both route to the PARALLEL streaming reader (never the single-threaded
        // kseq loop), so large uncompressed assemblies always parse on N threads:
        //   • RABBIT_FORCE_STREAM=1  — force streaming (for environments where
        //     mmap behaves poorly, e.g. some network/over-committed file systems).
        //   • automatic             — mmap could not be established (NFS /
        //     restricted-fs / mmap-capped env).
        // gzread reads plain files transparently — the same zlib decode basis the
        // legacy kseq fallback used — so parsed contigs, order and results match
        // the mmap path exactly.
        const bool force_stream = [] {
          const char *e = rb_getenv("RABBIT_FORCE_STREAM"); return e && e[0] == '1';
        }();
        producer_ok = false;
        if (!force_stream) {
          producer_ok = parse_fasta_mmap_parallel(
              inFile, (int)numThreads,
              minContig, min_small_contig,
              fullHeader,
              stream_pmh_mmap, no_store_seqs_mmap,
              /*collect_tiny=*/outUnbinned,
              pmh_m_mmap, pmh_k_mmap,
              mmap_large, mmap_small, mmap_tiny,
              num_seqs);
        }
        if (!producer_ok) {
          verbose_message("Using parallel streaming reader for uncompressed FASTA "
                          "(%s) [%.1fGb / %.1fGb]\n",
                          force_stream ? "RABBIT_FORCE_STREAM=1" : "mmap unavailable",
                          getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
          mmap_large.clear(); mmap_small.clear(); mmap_tiny.clear();
          used_stream_reader = true;
          producer_ok = parse_fasta_stream_parallel(
              inFile, (int)numThreads,
              minContig, min_small_contig,
              fullHeader,
              stream_pmh_mmap, no_store_seqs_mmap,
              /*collect_tiny=*/outUnbinned,
              pmh_m_mmap, pmh_k_mmap,
              /*input_is_gz=*/false,
              mmap_large, mmap_small, mmap_tiny,
              num_seqs);
        }
      }

      if (producer_ok) {

        verbose_message("Parallel %s parsed %zu seqs (%zu large, %zu small) "
                        "[%.1fGb / %.1fGb]\n",
                        used_stream_reader ? (is_gz ? "streaming gzip" : "streaming")
                                           : "mmap",
                        num_seqs, mmap_large.size(), mmap_small.size(),
                        getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

        // Pre-allocate g_win_flat now that we know the exact nobs
        if (stream_pmh_mmap) {
          g_win_flat.resize(mmap_large.size() * pmh_m_mmap, 0u);
        }
        // Pre-allocate k4freq flat array if IDF rebuild or exact-cosine is needed
        const bool do_idf_rebuild = g_idf_norm && pmh_k_mmap == 4 && g_gc_norm == 1;
        const bool do_collect_k4  = do_idf_rebuild ||
                                    (g_exact_cos_cmp && pmh_k_mmap == 4 && g_gc_norm >= 1);
        if (do_collect_k4) {
          g_k4freq_flat.assign(mmap_large.size() * 256, 0.0f);
        }

        // Merge results into global arrays (serial, but O(nobs) and fast).
        // Sequences become string_views into the retained mmap (single-line) or
        // the per-thread arenas (multi-line); no copy here either.
        bool took_mmap_view = false;
        auto rec_view = [&](const ContigRec &rec) -> std::string_view {
          if (rec.seq_view_ptr) { took_mmap_view = true;
            return std::string_view(rec.seq_view_ptr, rec.len); }
          return std::string_view(g_seq_arenas[rec.arena_tid].data() + rec.arena_off,
                                  rec.len);
        };
        // Reserve once (exact element counts are known) so the per-record
        // push_back never reallocates these large vectors mid-merge.
        contig_names.reserve(contig_names.size() + mmap_large.size());
        seq_lens.reserve(seq_lens.size() + mmap_large.size());
        if (!no_store_seqs_mmap) seqs.reserve(seqs.size() + mmap_large.size());
        small_contig_names.reserve(small_contig_names.size() + mmap_small.size());
        small_seq_lens.reserve(small_seq_lens.size() + mmap_small.size());
        small_seqs.reserve(small_seqs.size() + mmap_small.size());
        if (has_depth) {
          contigs.reserve(contigs.size() + mmap_large.size());
          small_contigs.reserve(small_contigs.size() + mmap_small.size());
        }
        size_t nskip = 0, nskip1 = 0;
        for (auto &rec : mmap_large) {
          if (has_depth && contigs.find(rec.name) != contigs.end()) { nskip++; had_dup_names = true; continue; }
          const size_t r = nobs++;
          contigs[rec.name] = r + nskip;
          // rec.name is not used past this point → move it into the permanent
          // store instead of copying (the map keeps its own key copy).
          contig_names.push_back(std::move(rec.name));
          seq_lens.push_back(rec.len);
          if (!no_store_seqs_mmap) seqs.push_back(rec_view(rec));
          if (stream_pmh_mmap && !rec.winners.empty()) {
            uint32_t *dst = g_win_flat.data() + r * pmh_m_mmap;
            std::memcpy(dst, rec.winners.data(), pmh_m_mmap * sizeof(uint32_t));
          }
          if (do_collect_k4 && !rec.k4freq.empty()) {
            std::memcpy(g_k4freq_flat.data() + r * 256,
                        rec.k4freq.data(), 256 * sizeof(float));
          }
        }
        for (auto &rec : mmap_small) {
          if (has_depth && small_contigs.find(rec.name) != small_contigs.end()) { nskip1++; had_dup_names = true; continue; }
          small_contigs[rec.name] = nobs1 + nskip1;
          small_contig_names.push_back(std::move(rec.name));
          if (!no_store_seqs_mmap) small_seqs.push_back(rec_view(rec));
          small_seq_lens.push_back(rec.len);
          nobs1++;
        }
        // If no single-line view referenced the mmap (e.g. a fully multi-line
        // assembly, where every sequence lives in an arena), release the
        // mapping now to avoid retaining ~filesize of resident pages.
        if (!took_mmap_view && g_fasta_mmap) {
          munmap((void *)g_fasta_mmap, g_fasta_mmap_len);
          g_fasta_mmap = nullptr; g_fasta_mmap_len = 0;
        }
        if (os) {
          for (auto &p : mmap_tiny) {
            if (onlyLabel) *os << p.first << line_delim;
            else printFasta(*os, p.first, p.second);
          }
        }

        if (stream_pmh_mmap) {
          // Trim to actual nobs (nskip may have reduced it)
          g_win_flat.resize((size_t)nobs * pmh_m_mmap);
          g_win_flat.shrink_to_fit();
          g_pmh_built_streaming = true;
          verbose_message("Parallel mmap PMH sketch complete (nobs=%zu) "
                          "[%.1fGb / %.1fGb]\n",
                          nobs, getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
        }

        // IDF rebuild: compute global IDF from collected k=4 exact frequencies,
        // then rebuild g_win_flat with TF×IDF weighted PMH sketches.
        if (do_idf_rebuild && !g_k4freq_flat.empty() && nobs > 0) {
          g_k4freq_flat.resize((size_t)nobs * 256);  // trim to actual nobs

          // Step 1: compute df[v] = count of contigs containing canonical 4-mer v
          std::vector<uint32_t> df(256, 0u);
          for (size_t r = 0; r < (size_t)nobs; ++r) {
            const float *row = g_k4freq_flat.data() + r * 256;
            for (int v = 0; v < 256; ++v)
              if (row[v] > 0.0f) df[v]++;
          }

          // Step 2: idf[v] = log2((N+1)/(df[v]+0.5))
          std::vector<float> idf(256, 0.0f);
          const double N1 = (double)nobs + 1.0;
          for (int v = 0; v < 256; ++v)
            idf[v] = (df[v] > 0) ? (float)std::log2(N1 / (df[v] + 0.5)) : 0.0f;

          // Step 3: rebuild PMH sketches in parallel with TF×IDF weights
          g_win_flat.assign((size_t)nobs * pmh_m_mmap, 0u);
          std::vector<uint64_t> idf_scratch; // dummy, not used
          #pragma omp parallel for schedule(static) num_threads(numThreads) \
                  firstprivate(idf_scratch)
          for (int r = 0; r < (int)nobs; ++r) {
            const float *row = g_k4freq_flat.data() + (size_t)r * 256;
            Sketch::ProbMinHash4 pmh(pmh_m_mmap, 4, 42u);
            for (int v = 0; v < 256; ++v) {
              if (row[v] > 0.0f && idf[v] > 0.0f)
                pmh.addHash((uint64_t)v, (double)(row[v] * idf[v]));
            }
            uint32_t *dst = g_win_flat.data() + (size_t)r * pmh_m_mmap;
            const uint64_t *w64 = pmh.getWinners();
            for (uint32_t j = 0; j < pmh_m_mmap; ++j)
              dst[j] = (uint32_t)(w64[j] ^ (w64[j] >> 32));
          }

          // Free exact freq storage unless we also need it for exact cosine
          if (!g_exact_cos_cmp)
            std::vector<float>().swap(g_k4freq_flat);
          verbose_message("IDF rebuild complete (nobs=%zu) "
                          "[%.1fGb / %.1fGb]\n",
                          nobs, getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
        }

        // Exact cosine mode: L2-normalise (optionally z-score) → g_k4cosine_flat
        // RABBIT_EXACT_COS=1: cosine of raw GC-norm freqs (each centred at 1.0 avg)
        // RABBIT_EXACT_COS=2: cosine of z-score = (GC_norm_freq - 1.0), emphasising
        //                     k-mers that deviate from their per-base expected frequency.
        if (g_exact_cos_cmp && !g_k4freq_flat.empty()) {
          g_k4freq_flat.resize((size_t)nobs * 256);
          g_k4cosine_flat.resize((size_t)nobs * 256);
          const bool use_zscore = ([] { const char *e = rb_getenv("RABBIT_EXACT_COS");
                                       return e ? std::atoi(e) : 0; }() >= 2);
          #pragma omp parallel for schedule(static) num_threads(numThreads)
          for (int r = 0; r < (int)nobs; ++r) {
            const float *src = g_k4freq_flat.data() + (size_t)r * 256;
            float *dst = g_k4cosine_flat.data() + (size_t)r * 256;
            double ss = 0.0;
            for (int v = 0; v < 256; ++v) {
              double x = use_zscore ? ((double)src[v] - 1.0) : (double)src[v];
              dst[v] = (float)x;
              ss += x * x;
            }
            const float inv_norm = (ss > 1e-30) ? (float)(1.0 / std::sqrt(ss)) : 0.0f;
            for (int v = 0; v < 256; ++v) dst[v] *= inv_norm;
          }
          std::vector<float>().swap(g_k4freq_flat);
          verbose_message("Exact cosine (z-score=%s) L2-norm complete (nobs=%zu) "
                          "[%.1fGb / %.1fGb]\n",
                          use_zscore ? "yes" : "no",
                          nobs, getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
        }

        parsed_ok = true;
      } else {
        verbose_message("parallel %s parse failed, falling back to streaming kseq\n",
                        is_gz ? "streaming gzip" : "mmap");
      }
    }

    // ── streaming zlib fallback (plain FASTA, large .gz, or libdeflate error) ─
    // Borrows RabbitTClust's ThreadPool concept: one kseq reader thread produces
    // records; N-1 OMP task threads build PMH sketches concurrently (overlap
    // decompression with computation).  Only active when RABBIT_PMH=1.
    if (!parsed_ok) {
      // ── Detect streaming sketch mode early (before kseq loop) ────────────
      // We read env vars here so g_win_flat can be pre-reserved before any
      // task is spawned (prevents reallocation while tasks hold data pointers).
      const bool stream_pmh = rb_env_pmh_on();
      const bool no_store_seqs = [&]() -> bool {
        // RABBIT_NOSTORE_SEQS=1: skip in-RAM sequence storage (saves 100GB+ for
        // huge assemblies); forces label-only output since sequences are gone.
        if (const char *e = rb_getenv("RABBIT_NOSTORE_SEQS")) return e[0] == '1';
        // Auto-drop when sequences have no downstream consumer (no per-bin/
        // unbinned FASTA, no composition-sketch recruit, no RB_SEQHASH) and PMH
        // winners are produced during this streaming parse — bit-identical, and
        // removes the dominant (~assembly-sized) memory cost.
        return stream_pmh && noBinOut && !outUnbinned &&
               recruitSimFactor <= 0.0 && getenv("RB_SEQHASH") == nullptr;
      }();
      const uint32_t s_pmh_m = sketch_size;
      const int s_pmh_k = [&]() -> int {
        const char *e = rb_getenv("RABBIT_PMHK"); return e ? std::atoi(e) : 4;
      }();
      // Per-thread scratch buffers for build_pmh_winners (reused across tasks
      // on the same thread; tasks on one thread are non-concurrent → safe).
      std::vector<std::vector<uint64_t>> s_pmh_scratch;
      if (stream_pmh) {
        s_pmh_scratch.resize(numThreads);
        // Upper-bound nobs estimate: one contig per (minContig/4) bytes compressed.
        // Capped at 2M to avoid excessive pre-allocation (2M×2048×4B = 16 GB).
        const size_t est_nobs = std::min(
            gz_file_bytes > 0
                ? gz_file_bytes / std::max((size_t)minContig / 4, (size_t)200)
                : (size_t)100000,
            (size_t)2000000);
        g_win_flat.reserve(est_nobs * s_pmh_m);
        verbose_message("Streaming sketch: pre-reserved %.0fMB for ~%zu contigs\n",
                        est_nobs * s_pmh_m * 4.0 / (1 << 20), est_nobs);
        if (no_store_seqs) {
          onlyLabel = true;
          verbose_message("RABBIT_NOSTORE_SEQS=1: sequences not stored; "
                          "output will be label-only (saves ~100GB+ RAM)\n");
        }
      }

      size_t inFileSize = filesize(inFile.c_str());
      if (!is_gz) {
        inFileSize *= 4;
        verbose_message("Estimating uncompressed size of %s as 4x: %ld\n",
                        inFile.c_str(), inFileSize);
      }
      ProgressTracker in_progress(inFileSize);
      gzFile fgz = gzopen(inFile.c_str(), "r");
      if (fgz == NULL) {
        cerr << "[Error!] can't open the sequence fasta file " << inFile << "\n";
        return 1;
      }
      kseq_t *seq = kseq_init(fgz);
      size_t nskip = 0, nskip1 = 0;

      // ── Producer-consumer: kseq in single region, PMH tasks in workers ───
      // When stream_pmh=false the parallel region degenerates to 1 thread
      // (num_threads(1)), so there is zero OMP overhead for the simple path.
#pragma omp parallel num_threads(stream_pmh ? (int)numThreads : 1)
      {
#pragma omp single
        {
          while (kseq_read(seq) >= 0) {
            num_seqs++;
            if (in_progress.isStepMarker())
              verbose_message("Parsing assembly %s [%.1fGb / %.1fGb]             \r",
                              in_progress.getProgress(),
                              getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
            in_progress.track(seq->seq.l + seq->name.l + 2);
            std::string name(seq->name.s);
            if (!fullHeader) trim_fasta_label(name);
            if (seq->seq.l >= minContig) {
              if (has_depth && contigs.find(name) != contigs.end()) {
                if (debug) verbose_message("Skipping duplicate contig: %s\n", name.c_str());
                nskip++; had_dup_names = true;
                continue;
              }
              const size_t r = nobs++;
              contigs[name] = r + nskip;
              contig_names.push_back(name);
              seq_lens.push_back(seq->seq.l);
              if (!no_store_seqs) {
                g_seq_store.emplace_back(seq->seq.s, seq->seq.l);  // kseq buffer is reused
                seqs.push_back(std::string_view(g_seq_store.back()));
              }

              if (stream_pmh) {
                // Grow g_win_flat by one row in the single thread (safe: no
                // concurrent task writes beyond existing rows).
                const size_t need = (r + 1) * s_pmh_m;
                if (need > g_win_flat.capacity()) {
                  // Drain all pending tasks before any reallocation so that
                  // in-flight tasks finish writing to the old buffer first.
#pragma omp taskwait
                  g_win_flat.reserve((r + 65536) * s_pmh_m);
                }
                g_win_flat.resize(need, 0u);

                // Copy sequence for OMP task (freed inside the task).
                char *sc = (char *)malloc(seq->seq.l + 1);
                memcpy(sc, seq->seq.s, seq->seq.l + 1);
                const size_t sl  = seq->seq.l;
                const uint32_t pm = s_pmh_m;
                const int      pk = s_pmh_k;
                const size_t   ri = r;
#pragma omp task firstprivate(ri, sc, sl, pm, pk)
                {
                  const int tid = omp_get_thread_num();
                  // g_win_flat.data() is read at task execution time; safe
                  // because any reallocation is guarded by taskwait above.
                  uint32_t *wrow = g_win_flat.data() + ri * pm;
                  build_pmh_winners(sc, sl, pk, pm, /*seed=*/42u,
                                    wrow, s_pmh_scratch[tid]);
                  free(sc);
                }
              }
            } else if (seq->seq.l >= min_small_contig) {
              if (has_depth && small_contigs.find(name) != small_contigs.end()) {
                nskip1++; had_dup_names = true;
                continue;
              }
              small_contigs[name] = nobs1 + nskip1;
              small_contig_names.push_back(name);
              g_seq_store.emplace_back(seq->seq.s, seq->seq.l);  // kseq buffer is reused
              small_seqs.push_back(std::string_view(g_seq_store.back()));
              small_seq_lens.push_back(seq->seq.l);
              nobs1++;
            } else {
              if (outUnbinned && os) {
                if (onlyLabel)
                  *os << name << line_delim;
                else
                  printFasta(*os, name, std::string(seq->seq.s, seq->seq.l));
              }
            }
          } // while kseq_read
#pragma omp taskwait // wait for remaining PMH sketch tasks
        } // end omp single
      } // end omp parallel

      kseq_destroy(seq);
      gzclose(fgz);

      if (stream_pmh) {
        // Trim to actual size; shrink_to_fit releases the over-reserved capacity.
        g_win_flat.resize((size_t)nobs * s_pmh_m);
        g_win_flat.shrink_to_fit();
        g_pmh_built_streaming = true;  // signal sketch build loop to skip PMH
        verbose_message("Streaming PMH sketch complete (nobs=%zu) "
                        "[%.1fGb / %.1fGb]\n",
                        nobs, getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
      }

      verbose_message("Parsed %d sequences from FASTA. nobs=%d nobs1=%d\n",
                      num_seqs, nobs, nobs1);
    }
    if (os) { os->close(); delete os; }

    // ── Merge pre-parsed depth_matrix data (async result) into depth_matrix matrices ───────
    // depth_future was launched before FASTA decompression and should be
    // complete (or nearly so) by the time we reach here.
    if (has_depth) {
      rb_phase("fasta parse+sketch done");
      verbose_message("Merging abundance data [%.1fGb / %.1fGb]\n",
                      getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

      DepthMap depth_map = depth_future.get();   // blocks only if not yet done

      size_t r = 0, r1 = 0;
      size_t num = 0, nskip = 0;
      size_t num1 = 0, nskip1 = 0;
      size_t ignored_too_small = 0;

      depth_matrix.resize(nobs, num_depth_samples, false);
      depth_var_matrix.resize(nobs, num_depth_samples, false);
      small_depth_matrix.resize(nobs1, num_depth_samples, false);

      // Iterate large contigs in order (preserves original row semantics).
      // Fast path: when the FASTA had no duplicate headers, contigs[name] == ci
      // by construction (kept order == insertion order, nskip stayed 0), so the
      // per-contig contigs.find() is redundant and is skipped entirely.
      for (size_t ci = 0; ci < contig_names.size(); ++ci) {
        const std::string& name = contig_names[ci];
        size_t idx, row;
        if (had_dup_names) {
          auto it = contigs.find(name);
          if (it == contigs.end()) { nskip++; continue; }
          idx = it->second;
          row = idx - nskip;
          if (row >= nobs) continue;
        } else {
          idx = ci; row = ci;
        }

        auto jt = depth_map.find(name);
        if (jt == depth_map.end()) { ignored_too_small++; continue; }
        const RawDepthEntry& entry = jt->second;

        for (int i = 0; i < num_depth_samples; ++i)
          depth_matrix(row, i) = entry.means[i];
        // NOTE: depth_var_matrix is intentionally NOT filled — its only reader,
        // cal_depth_dist(), has no live call sites, so the per-element variance
        // store here was dead work.  Re-add this if cal_depth_dist is revived.
        r++;  num++;  totalSize += seq_lens[idx];
      }

      // Iterate small contigs (same fast path as the large-contig loop above)
      for (size_t ci = 0; ci < small_contig_names.size(); ++ci) {
        const std::string& name = small_contig_names[ci];
        size_t idx, row;
        if (had_dup_names) {
          auto it = small_contigs.find(name);
          if (it == small_contigs.end()) { nskip1++; continue; }
          idx = it->second;
          row = idx - nskip1;
          if (row >= nobs1) continue;
        } else {
          idx = ci; row = ci;
        }

        auto jt = depth_map.find(name);
        if (jt == depth_map.end()) { ignored_too_small++; continue; }
        const RawDepthEntry& entry = jt->second;

        for (int i = 0; i < num_depth_samples; ++i) {
          small_depth_matrix(row, i) = entry.means[i];
        }
        r1++;  num1++;  totalSize1 += small_seq_lens[idx];
      }

      // depth_map (one entry per depth-file row — often millions, incl. all the
      // tiny contigs that never get binned) is no longer needed after the merge.
      // Its single-threaded teardown (freeing ~N hash slots + per-row mean/var
      // vectors) costs hundreds of ms; hand it to a detached thread so that the
      // free overlaps the sketch + graph-build work instead of stalling here.
      // The contigs / small_contigs dedup maps are also done being read here;
      // fold their (string-key) teardown into the same detached thread.
      std::thread([dm = std::move(depth_map),
                   cm = std::move(contigs),
                   scm = std::move(small_contigs)]() mutable { }).detach();

      nobs  = r;
      nobs1 = r1;

      // Remove any empty-string sentinel entries (legacy compatibility).
      // seq_lens is filtered in sync with seqs/contig_names.
      if (!seqs.empty()) {
        seqs.erase(std::remove(seqs.begin(), seqs.end(), std::string_view{}), seqs.end());
      }
      if (!small_seqs.empty()) {
        small_seqs.erase(std::remove(small_seqs.begin(), small_seqs.end(), std::string_view{}), small_seqs.end());
      }
      contig_names.erase(std::remove(contig_names.begin(), contig_names.end(), ""), contig_names.end());
      small_contig_names.erase(std::remove(small_contig_names.begin(), small_contig_names.end(), ""), small_contig_names.end());

      depth_matrix.resize(nobs, num_depth_samples, true);
      depth_var_matrix.resize(nobs, num_depth_samples, true);
      small_depth_matrix.resize(nobs1, num_depth_samples, true);

      verbose_message("Merged %d contigs and %d coverages from %s "
                      "[%.1fGb / %.1fGb]. Ignored %d too-small contigs.\n",
                      r, num_depth_samples, depth_file.c_str(),
                      getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024,
                      ignored_too_small);
    } else {
      // No depth_matrix file - set totalSize from sequence lengths
      for (auto l : seq_lens)       totalSize  += l;
      for (auto l : small_seq_lens) totalSize1 += l;
      contigs.clear();
      small_contigs.clear();
      depth_matrix.resize(nobs, 1, false);
      depth_var_matrix.resize(nobs, 1, false);
      small_depth_matrix.resize(nobs1, 1, false);
    }
  }
  }  // end else (non-cache feature construction)

  verbose_message("Number of target contigs: %d of large (>= %d) and %d of "
                  "small ones (>=%d & <%d). \n",
                  nobs, minContig, nobs1, 1000, minContig);

  if (nobs == 0) {
    cerr << "[Error!] There were no large target contigs. Cannot proceed. "
            "Rerun with the '-d' option for more details.\n";
    return 1;
  }

  // ── Build composition sketches ────────────────────────────────────────
  // When the inverted-index graph build is active (RABBIT_GRAPH_INDEX=1) we also
  // insert each sketch's (bucket, min-value) keys directly into a per-thread
  // phmap accumulator here, so the index is fully populated by the time sketch
  // building finishes.  This eliminates the separate N × m skKeys array.
  const bool build_index = !from_cache && []() -> bool {
    const char *e = rb_getenv("RABBIT_GRAPH_INDEX");
    return e && e[0] == '1';
  }();

  verbose_message("Building composition sketches%s. nobs=%zd k=%d m=%d\n",
                  (num_depth_samples > 1 ? " + abundance ranking" : ""),
                  nobs, sketch_kmer_size, sketch_size);

  // ── Weighted ProbMinHash4 path setup must come BEFORE builder creation ──
  g_pmh_mode    = rb_env_pmh_on();
  g_pmh_k       = [] { const char *e = rb_getenv("RABBIT_PMHK"); return e ? std::atoi(e) : 4; }();
  g_pmh_base_on = [] { const char *e = rb_getenv("RABBIT_PMH_BASE"); return !e || e[0] != '0'; }();
  g_w_comp      = rb_env_w_comp();
  g_depth_sim       = rb_env_depth_sim();
  g_depth_fuse      = rb_env_depth_fuse();
  g_depth_wjac_norm = rb_env_depth_wjac_norm_on();
  g_depth_wjac_cut  = rb_env_depth_wjac_cut();
  g_depth_wjac_gate = rb_env_depth_wjac_gate();
  g_split_keep_coherent = rb_env_split_keep_coherent();
  g_split_coh_pct = rb_env_split_coh_pct();
  g_split_guard = rb_env_split_guard();
  g_split_guard_bp = rb_env_split_guard_bp();
  g_split_guard_margin = rb_env_split_guard_margin();
  g_mutual_knn  = rb_env_mutual_knn_on();
  g_neg_depth_thr = rb_env_neg_depth_thr();
  g_edge_power  = [] { const char *e = rb_getenv("RABBIT_EDGE_POWER"); double v = e ? std::atof(e) : 1.0;
                       return (v < 0.1) ? 1.0 : v; }();
  g_gc_norm     = rb_env_gc_norm();
  g_gc_norm_cap = [] { const char *e = rb_getenv("RABBIT_GC_NORM_CAP"); double v = e ? std::atof(e) : 20.0;
                       return (v < 1.0) ? 20.0 : v; }();
  // PMH winner-identity inverted index: a diagnostic-only structure (the graph
  // build uses the all-pairs SIMD kernel, not the index — see the reset at the
  // end of the sketch section).  Building it costs the full 64-bit winner array
  // + ~m hash-map insertions per contig, all immediately discarded, so it is
  // OFF unless RABBIT_PMH_INDEX is set.
  const bool pmh_index_on = g_pmh_mode && !from_cache && (getenv("RABBIT_PMH_INDEX") != nullptr);
  if (g_pmh_mode && !from_cache) {
    g_pmh_m    = sketch_size;
    g_inv_pmh_m = (g_pmh_m > 0) ? (1.0 / (double)g_pmh_m) : 0.0;
    g_pmh_seed = 42;
    if (g_pmh_built_streaming) {
      // g_win_flat already populated during the kseq streaming pass;
      // skip re-allocation and re-computation.
      verbose_message("RABBIT_PMH=1 (streaming pre-built): "
                      "skipping g_win_flat allocation (already %zu entries)\n",
                      g_win_flat.size() / g_pmh_m);
    } else {
      g_win_flat.assign((size_t)nobs * g_pmh_m, 0u);
      // Full 64-bit winners are ONLY needed to feed the (diagnostic, off-by-
      // default) PMH inverted index.  The similarity kernel uses the 32-bit
      // folded winners in g_win_flat, so skip the 64-bit array (≈8 B/register ×
      // nobs × m — e.g. ~300 MB on plant) unless the index is explicitly asked
      // for.  See pmh_index_on below.
      if (pmh_index_on) g_win64_flat.assign((size_t)nobs * g_pmh_m, 0ULL);
    }
    verbose_message("RABBIT_PMH=1: weighted ProbMinHash4 graph metric "
                    "(k=%d, m=%u, count-weighted Jaccard, baseline_corr=%d)\n",
                    g_pmh_k, g_pmh_m, (int)g_pmh_base_on);
  }

  // Construct the builders before the parallel loop so threads can call insert().
  std::unique_ptr<rabbit_invidx::InvertedIndexBuilder> idx_builder;
  if (build_index)
    idx_builder = std::make_unique<rabbit_invidx::InvertedIndexBuilder>(
        (int)numThreads);

  // PMH-winner inverted index builder (RABBIT_PMH_INDEX only; see pmh_index_on).
  // Built in the same parallel loop — no extra pass over seqs[].
  std::unique_ptr<rabbit_invidx::InvertedIndexBuilder> pmh_idx_builder;
  if (pmh_index_on)
    pmh_idx_builder = std::make_unique<rabbit_invidx::InvertedIndexBuilder>(
        (int)numThreads);

  // ── Fusion B+E: sketch update + buildSig + sig_flat copy + freeReg
  //               + abundance ranking, ALL in ONE parallel loop ─────────────
  // In weighted-ProbMinHash mode the OPH (k=21) b-bit sketch is never read
  // (graph_sim uses the PMH winners), so skip building it entirely unless an
  // OPH-consuming option is active (inverted index or composition-min recruitment).
  // This removes a full per-contig MinHash pass + 30k heap allocations.
  const bool oph_needed = !g_pmh_mode || build_index || (recruitSimFactor > 0.0);

  g_sig_nw = (sketch_size + 63) / 64;
  g_sig_np = sketch_bits;
  g_sig_m  = sketch_size;
  const size_t sig_stride = (size_t)g_sig_nw * g_sig_np;
  // g_sig_flat holds the OPH b-bit signatures; only written/read when oph_needed
  // (graph_sim uses g_win_flat in PMH mode).  Skip the alloc+zero otherwise.
  if (oph_needed) g_sig_flat.resize(nobs * sig_stride);

  // Spearman scratch buffer (one per thread, reused across iterations)
  std::vector<StoredDistance> rowMat_proto(num_depth_samples);
  std::vector<std::vector<StoredDistance>> threadRowMat(numThreads, rowMat_proto);

  // Per-thread reusable k-mer code buffer for the weighted ProbMinHash path.
  std::vector<std::vector<uint64_t>> threadPmhScratch(numThreads);
  // Per-thread key buffer for the PMH winner inverted index.
  std::vector<std::vector<uint64_t>> threadPmhKeys(numThreads);

  g_sketches.resize(nobs, nullptr);

  // Snapshot raw per-sample mean depths BEFORE the loop below rank-transforms
  // depth_matrix in place (Fusion E / Spearman). marker_guided_split needs raw means.
  // The whole block below (snapshots + per-contig sketch/rank loop) is the
  // expensive feature construction; in --load-cache mode it is restored from the
  // cache (g_large_means/g_depth_raw cached, depth_matrix already ranked) so it
  // is skipped entirely.
  if (!from_cache) {
  if ((!marker_seed_file.empty() || g_split_abundance) && num_depth_samples >= 1) {
    g_large_means.assign((size_t)nobs * num_depth_samples, 0.0f);
    for (size_t r = 0; r < nobs; ++r)
      for (size_t i = 0; i < (size_t)num_depth_samples; ++i)
        g_large_means[r * num_depth_samples + i] = (float)depth_matrix(r, i);
  }

  // Weighted-Jaccard abundance metric needs RAW per-sample depths (magnitudes),
  // so snapshot them here BEFORE the in-place Spearman rank transform below.
  if (g_depth_sim >= 1 && num_depth_samples > 1) {
    g_depth_raw.assign((size_t)nobs * num_depth_samples, 0.0f);
    for (size_t r = 0; r < nobs; ++r)
      for (size_t i = 0; i < (size_t)num_depth_samples; ++i)
        g_depth_raw[r * num_depth_samples + i] = (float)depth_matrix(r, i);
    // Per-sample normalisation: scale each sample column by 1/mean so unequal
    // library sizes / sequencing depths don't dominate the weighted Jaccard.
    if (g_depth_wjac_norm) {
      g_depth_colnorm.assign(num_depth_samples, 1.0);
      for (size_t i = 0; i < (size_t)num_depth_samples; ++i) {
        double s = 0.0;
        for (size_t r = 0; r < nobs; ++r) s += g_depth_raw[r * num_depth_samples + i];
        double mean = s / std::max<size_t>(nobs, 1);
        g_depth_colnorm[i] = (mean > 1e-12) ? 1.0 / mean : 0.0;
      }
    } else {
      g_depth_colnorm.clear();
    }
    verbose_message("Abundance edge metric: %s (per-sample norm %s)\n",
                    g_depth_sim == 2 ? "fuse(corr×wjac)" : "weighted Jaccard",
                    g_depth_wjac_norm ? "ON" : "OFF");
  }

  {
    ProgressTracker progress(nobs);
#pragma omp parallel for num_threads(numThreads) schedule(dynamic, 1)
    for (size_t r = 0; r < nobs; ++r) {
      // ── Sketch (Fusion B): update → buildSig → sig_flat copy → freeReg
      if (oph_needed) {
        g_sketches[r] = new rabbit_sketch::KmerSketch(sketch_size, sketch_kmer_size, sketch_bits);
        g_sketches[r]->update(seqs[r].data(), seqs[r].size());  // fills reg_[]

        if (build_index) {
          const int tid = omp_get_thread_num();
          std::vector<uint64_t> ks;
          g_sketches[r]->getKeys(ks);
          for (uint64_t k : ks)
            idx_builder->insert(tid, k, (uint32_t)r);
        }

        const uint64_t* sig = g_sketches[r]->getSignature();         // reg_[]→bits_[]
        std::memcpy(g_sig_flat.data() + r * sig_stride, sig,
                    sig_stride * sizeof(uint64_t));                   // bits_[]→g_sig_flat
        // The inverted-index graph build re-derives keys via getKeys() (reads
        // reg_[]) at graph time, so keep reg_[] alive when indexing.
        if (!build_index) g_sketches[r]->freeRegisters();             // free reg_[]
      }

      // ── Weighted ProbMinHash4 winners (RABBIT_PMH=1): frequency-weighted,
      //     small-k composition spectrum ───────────────────────────────────
      // Skip when g_pmh_built_streaming: winners were already written to
      // g_win_flat during the kseq streaming pass (streaming producer-consumer).
      if (g_pmh_mode && !g_pmh_built_streaming) {
        const int tid = omp_get_thread_num();
        uint32_t *wrow32 = g_win_flat.data() + (size_t)r * g_pmh_m;
        if (pmh_index_on) {
          uint64_t *wrow64 = g_win64_flat.data() + (size_t)r * g_pmh_m;
          build_pmh_winners(seqs[r].data(), seqs[r].size(), g_pmh_k, g_pmh_m,
                            g_pmh_seed, wrow32, threadPmhScratch[tid],
                            wrow64, &threadPmhKeys[tid]);
          for (uint64_t key : threadPmhKeys[tid])
            pmh_idx_builder->insert(tid, key, (uint32_t)r);
        } else {
          // Default: only the 32-bit folded winners (consumed by the kernel).
          build_pmh_winners(seqs[r].data(), seqs[r].size(), g_pmh_k, g_pmh_m,
                            g_pmh_seed, wrow32, threadPmhScratch[tid],
                            /*out64=*/nullptr, /*out_keys=*/nullptr);
        }
      }

      // ── Spearman ranking (Fusion E): rank depth_matrix[r] in-place per thread
      if (num_depth_samples > 1) {
        auto &rowMat = threadRowMat[omp_get_thread_num()];
        MatrixRowType rRow(depth_matrix, r);
        std::copy(rRow.begin(), rRow.end(), rowMat.begin());
        rank(rowMat, rowMat);
        std::copy(rowMat.begin(), rowMat.end(), rRow.begin());
      }

      if (verbose) {
        progress.track();
        if (omp_get_thread_num() == 0 && progress.isStepMarker())
          verbose_message("Building sketch %s\r", progress.getProgress());
      }
    }
  }
  }  // end if (!from_cache): snapshots + sketch/rank loop

  // Finalise the index (merge thread maps, remove singletons, CSR flatten).
  if (build_index) {
    verbose_message("Merging inverted index thread maps...\n");
    g_inv_idx = std::make_unique<rabbit_invidx::InvertedIndex>(
        idx_builder->build((int)numThreads));
    idx_builder.reset();
    verbose_message("Inverted index ready: %zu postings / %zu keys "
                    "[%.1fGb / %.1fGb]\n",
                    g_inv_idx->totalPostings, g_inv_idx->postIdx.size(),
                    getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
  }

  // The PMH-winner inverted index (RABBIT_PMH_INDEX) is diagnostic only: the
  // graph build uses the all-pairs SIMD kernel, not the index (small-k → dense
  // postings → all-pairs is faster), so free the builder + 64-bit winners now.
  // When the index is off (default) neither was ever allocated.
  if (pmh_idx_builder) {
    pmh_idx_builder.reset();  // free builder thread-maps; index not needed
    g_win64_flat.clear(); g_win64_flat.shrink_to_fit(); // not queried at run time
  }

  verbose_message("Composition sketches ready%s. [%.1fGb / %.1fGb]"
                  "                          \n",
                  (num_depth_samples > 1 ? " + Spearman" : ""),
                  getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
  if (num_depth_samples > 1) verbose_message("Calculated spearman for large contigs\n");

  rb_phase("parse+sketch done");
  if (!from_cache) rb_seq_integrity_check("post-parse");

  // Precompute per-contig abundance non-zero flags so is_nz() in the O(N²)
  // graph loop reduces to two byte loads instead of an num_depth_samples matrix scan.
  build_anynz_cache();

  // Pack the 32-bit folded winners to 16 bits for the all-pairs kernel.  Must
  // run BEFORE baseline estimation so b0 is measured on the same representation
  // the graph kernel uses.  Default ON; RABBIT_WIN_BITS=32 keeps the exact
  // 32-bit path.  LSH (reads 32-bit winner bands) and exact-cosine paths keep
  // the 32-bit array.
  {
    int wb = 16;
    if (const char *e = getenv("RABBIT_WIN_BITS")) wb = atoi(e);
    const bool lsh_on = [] { const char *e = getenv("RABBIT_LSH");
                             return e && e[0] == '1'; }();
    if (wb == 16 && g_pmh_mode && !g_exact_cos_cmp && !lsh_on &&
        !g_win_flat.empty()) {
      g_win16.resize(g_win_flat.size());
#pragma omp parallel for schedule(static)
      for (size_t t = 0; t < g_win_flat.size(); ++t)
        g_win16[t] = fold_win16(g_win_flat[t]);
      g_win_bits = 16;
      std::vector<uint32_t>().swap(g_win_flat);   // free the 32-bit array
      verbose_message("Winner array packed to 16-bit (%.2f GB freed) "
                      "[%.1fGb / %.1fGb]\n",
                      g_win16.size() * 2.0 / (1ull << 30),
                      getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
    }
  }

  // Estimate the winner-match baseline b0 before calibration so that both the
  // pSim auto-calibration and the graph edge weights see the stretched scale.
  // In --load-cache mode b0 (and g_inv_one_minus_b0) were restored from the
  // cache, and the winner array is gone, so skip re-estimation.
  if (g_pmh_mode && g_pmh_base_on && !from_cache) {
    g_pmh_baseline = estimate_pmh_baseline(nobs);
    g_inv_one_minus_b0 = (g_pmh_baseline < 1.0)
                         ? (1.0 / (1.0 - g_pmh_baseline)) : 1.0;
    verbose_message("PMH baseline b0 = %.4f (median winner-match over random "
                    "pairs); similarities rescaled to (s-b0)/(1-b0)\n",
                    g_pmh_baseline);
  }

  // Precompute centered + L2-normalised rank vectors once (depth_matrix is fully
  // rank-transformed by now).  Lets the O(E·S) edge-weight corr term be a pure
  // dot product (depth_corr_fast) instead of the per-pair Welford (3 div/sample +
  // 2 sqrt/pair).  Cost: one O(N·S) pass.  Identical Pearson-on-ranks value the
  // abdfirst prune already relies on, so edge weights match within float noise.
  if (num_depth_samples > 1) {
    const size_t S = num_depth_samples;
    g_depth_unit.assign(nobs * S, 0.0f);
#pragma omp parallel for schedule(static) num_threads(numThreads)
    for (size_t r = 0; r < nobs; ++r) {
      double mean = 0.0;
      for (size_t k = 0; k < S; ++k) mean += depth_matrix(r, k);
      mean /= (double)S;
      double ss = 0.0;
      for (size_t k = 0; k < S; ++k) { double d = (double)depth_matrix(r, k) - mean; ss += d * d; }
      if (ss > 0.0) {
        const double inv = 1.0 / std::sqrt(ss);
        float *u = g_depth_unit.data() + r * S;
        for (size_t k = 0; k < S; ++k)
          u[k] = (float)(((double)depth_matrix(r, k) - mean) * inv);
      }
    }
  }

  BinMap cls;
  do {
    std::vector<size_t> mems;
    {
      Graph g(nobs);

      // ── 1+2. Fusion D: single O(N²/2) pass for calib + graph ───────────
      // For PMH mode with auto-calibration and large N, one tiled pass
      // simultaneously accumulates calibration maxsim[] and per-contig
      // neighbor heaps, then calibrates, then emits edges.
      // Replaces 771M (calib) + 475M (graph) = 1,246M pair-sims with 475M.
      if (from_cache) {
        // Topology (from/to/sComp) was loaded from the cache; calibration is
        // already baked into simCutoff/g_pmh_baseline. Skip the O(N^2) build.
        g.from  = std::move(g_cache_from);
        g.to    = std::move(g_cache_to);
        g.sComp = std::move(g_cache_scomp);
        verbose_message("Reusing cached similarity graph: %zu edges "
                        "[%.1fGb / %.1fGb]\n", g.getEdgeCount(),
                        getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
      } else if (simCutoff < 1. && g_pmh_mode && nobs > 25000) {
        simCutoff = gen_fused_calib_graph(g, calib_connected_pct);
      } else {
        // ── Original sequential path (non-PMH, manual simCutoff, or small N) ──
        if (simCutoff < 1.) {
          if (nobs <= 25000) {
            simCutoff = calibrate_sim_cutoff(calib_connected_pct, true);
          } else {
            verbose_message("Running fused 10-round calibration (Fusion C)...\n");
            simCutoff = calibrate_sim_cutoff_fused(calib_connected_pct);
          }
        } else {
          simCutoff *= 10;
        }
        verbose_message(
            "Finished Preparing Similarity Graph Building [pSim = %2.2f] "
            "[%.1fGb / %.1fGb]                                            \n",
            simCutoff / 10., getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
        build_similarity_graph(g, simCutoff / 1000.);
      }

      // ── Save hook: persist the post-graph state, then continue normally so
      // the same invocation also emits its bins. ───────────────────────────
      if (!cache_save_file.empty()) {
        verbose_message("Writing re-binning cache to %s [%.1fGb / %.1fGb]\n",
                        cache_save_file.c_str(), getUsedPhysMem(),
                        getTotalPhysMem() / 1024 / 1024);
        if (rb_write_cache(cache_save_file, g, has_depth))
          verbose_message("Cache written (%zu edges).\n", g.getEdgeCount());
        else
          cerr << "[Warn] failed to write cache: " << cache_save_file << "\n";
      }

      // ── --no_gold: label-free multi-resolution auto-selection ─────────────
      // Build the (expensive) similarity graph ONCE, then run many CHEAP
      // (edgeScore + incidence + label-propagation) passes under different
      // α (g_w_comp) / edge_power settings, and select the partition with the
      // highest weighted modularity on the fixed composition graph — a quality
      // signal that needs NO ground truth. The selected membership then flows
      // into the normal recruit/split/output pipeline (no early exit), so the
      // final bins are a complete, auto-resolution result.
      // Trigger: --no_gold (default grid) or RABBIT_REUSE_SWEEP="a:p;a:p;..."
      // for a custom grid. RABBIT_NO_GOLD_DUMP=1 additionally writes per-config
      // and selected AMBER .binning files (for offline evaluation). When neither
      // is set, the single-resolution production path below is byte-identical.
      std::vector<size_t> membership;     // shared by the two branches below
      bool used_no_gold = false;
      const char *sweep_env = getenv("RABBIT_REUSE_SWEEP");
      const bool qc_auto = (g_auto_select || g_ensemble);
      if (no_gold || sweep_env || qc_auto) {
        const bool dump = (getenv("RABBIT_NO_GOLD_DUMP") != nullptr) ||
                          (sweep_env && !no_gold);

        // --auto/--ensemble: load the SCG marker map so each swept partition can
        // be scored by genome quality (completeness/contamination) instead of
        // the label-free modularity proxy.
        bool qc_ready = false;
        if (qc_auto) {
          if (g_markers_file.empty()) {
            cerr << "[Error!] --auto/--ensemble require --markers FILE "
                    "(generate once via scripts/rabbitbin_markers.sh)\n";
            return 1;
          }
          qc_ready = rb_load_markers_for_contigs(g_markers_file);
          if (!qc_ready)
            cerr << "[Warn] no markers matched contigs; falling back to "
                    "modularity selection\n";
        }
        using clk = std::chrono::steady_clock;
        auto ms_since = [](clk::time_point t0) {
          return std::chrono::duration<double, std::milli>(clk::now() - t0).count();
        };

        // Parse "a:p;a:p" configs; null/"1"/empty → default grid.
        std::vector<std::pair<double, double>> cfgs;
        {
          std::string spec(sweep_env ? sweep_env : "");
          if (spec.empty() || spec == "1") {
            cfgs = {{0.3, 1.0}, {0.4, 1.0}, {0.5, 1.0},
                    {0.5, 1.5}, {0.6, 1.0}};
          } else {
            size_t p = 0;
            while (p < spec.size()) {
              size_t semi = spec.find(';', p);
              std::string tok = spec.substr(p, semi == std::string::npos ? std::string::npos : semi - p);
              size_t colon = tok.find(':');
              if (colon != std::string::npos) {
                double a = std::atof(tok.substr(0, colon).c_str());
                double pw = std::atof(tok.substr(colon + 1).c_str());
                if (a >= 0.0 && a <= 1.0 && pw >= 0.1) cfgs.emplace_back(a, pw);
              }
              if (semi == std::string::npos) break;
              p = semi + 1;
            }
          }
        }
        const size_t M = cfgs.size();
        const size_t E = g.getEdgeCount();
        static constexpr StoredDistance SSCR_MAX = 1.0f - 1e-6f;

        // Recompute g.edgeScore from the (reused) g.sComp + depth correlation,
        // mirroring the production block below, under the current g_w_comp /
        // g_edge_power globals.
        auto compute_edgescore = [&]() {
          g.edgeScore.assign(E, 0.0f);
          if (has_depth) {
#pragma omp parallel for schedule(dynamic, 1)
            for (size_t e = 0; e < E; ++e) {
              size_t i = g.from[e], j = g.to[e];
              if (num_depth_samples <= 1) {
                double w = (double)g.sComp[e];
                if (g_edge_power != 1.0) w = std::pow(w, g_edge_power);
                g.edgeScore[e] = (StoredDistance)w;
              } else {
                bool depth_ok;
                double dterm = depth_edge_term(i, j, depth_ok);
                if (!depth_ok) { g.edgeScore[e] = 0.0f; continue; }
                double w = g_w_comp * (double)g.sComp[e] + (1.0 - g_w_comp) * dterm;
                if (!std::isfinite(w) || w < (double)min_edge_weight) w = 0.0;
                if (g_edge_power != 1.0 && w > 0.0) w = std::pow(w, g_edge_power);
                g.edgeScore[e] = (StoredDistance)w;
              }
            }
          } else {
            g.edgeScore = g.sComp;
            for (auto &s : g.edgeScore) if (s < min_edge_weight) s = 0.0f;
          }
        };

        // Rebuild incidence lists from the current g.edgeScore (clamp + sort).
        auto rebuild_incs = [&]() {
          g.incs.assign(nobs, {});
          for (size_t e = 0; e < E; ++e) {
            if (g.edgeScore[e] > 0) {
              if (g.edgeScore[e] > SSCR_MAX) g.edgeScore[e] = SSCR_MAX;
              g.incs[g.from[e]].push_back(e);
              g.incs[g.to[e]].push_back(e);
            }
          }
#pragma omp parallel for schedule(dynamic, 256)
          for (size_t v = 0; v < nobs; ++v) {
            auto &inc = g.incs[v];
            std::sort(inc.begin(), inc.end(), [&](size_t e1, size_t e2) {
              return g.getOtherNode(e1, v) < g.getOtherNode(e2, v);
            });
          }
        };

        auto run_lp = [&](long long sd) {
          std::vector<size_t> membership;
          std::vector<size_t> node_order(nobs);
          std::iota(node_order.begin(), node_order.end(), 0);
          std::shuffle(node_order.begin(), node_order.end(),
                       std::default_random_engine((unsigned)sd));
          cluster_by_propagation(g, membership, node_order);
          return membership;
        };

        auto write_binning = [&](const std::vector<size_t> &mem,
                                 const std::string &path) {
          FILE *f = fopen(path.c_str(), "w");
          if (!f) { fprintf(stderr, "[REUSE] cannot write %s\n", path.c_str()); return; }
          fputs("@Version:0.9.0\n@SampleID:gsa\n@@SEQUENCEID\tBINID\t_LENGTH\n", f);
          for (size_t i = 0; i < nobs; ++i) {
            if (g.incs[i].empty()) continue;  // unbinned (no positive edge)
            fprintf(f, "%s\t%zu\t%zu\n", contig_names[i].c_str(), mem[i],
                    (i < seq_lens.size() ? seq_lens[i] : 0));
          }
          fclose(f);
        };

        auto count_bins = [&](const std::vector<size_t> &mem) {
          std::unordered_set<size_t> b;
          size_t binned = 0;
          for (size_t i = 0; i < nobs; ++i)
            if (!g.incs[i].empty()) { b.insert(mem[i]); ++binned; }
          return std::make_pair(b.size(), binned);
        };

        fprintf(stderr,
                "[REUSE] graph reused: nobs=%zu edges=%zu configs=%zu seed=%lld\n",
                (size_t)nobs, E, M, (long long)seed);

        std::vector<std::vector<size_t>> mems_all(M);
        std::string base = std::string(outFile);

        // ── Label-free selection criterion: weighted modularity Q ─────────────
        // Q is evaluated on the FIXED, parameter-independent composition graph
        // (weights = g.sComp, identical for every config), so differences in Q
        // reflect ONLY how well each config's partition cuts the graph — not the
        // per-config reweighting. This is the no-ground-truth selector: pick the
        // config with the highest Q. Precompute weighted degree + total weight
        // once (O(E)); per-config Q is then O(E + N).
        std::vector<double> wdeg(nobs, 0.0);
        double Wtot = 0.0;
        for (size_t e = 0; e < E; ++e) {
          double w = (double)g.sComp[e];
          Wtot += w; wdeg[g.from[e]] += w; wdeg[g.to[e]] += w;
        }
        const double invWtot  = Wtot > 0 ? 1.0 / Wtot : 0.0;
        const double inv2Wtot = Wtot > 0 ? 1.0 / (2.0 * Wtot) : 0.0;
        auto modularity = [&](const std::vector<size_t> &mem) -> double {
          double Wintra = 0.0;
          for (size_t e = 0; e < E; ++e)
            if (mem[g.from[e]] == mem[g.to[e]]) Wintra += (double)g.sComp[e];
          std::unordered_map<size_t, double> vol;
          for (size_t i = 0; i < nobs; ++i)
            if (wdeg[i] > 0.0) vol[mem[i]] += wdeg[i];
          double sumsq = 0.0;
          for (auto &kv : vol) { double f = kv.second * inv2Wtot; sumsq += f * f; }
          return Wintra * invWtot - sumsq;
        };

        // ── Per-config cheap tail: edgeScore + incs + LP + modularity (+QC) ───
        double tail_total = 0.0;
        std::vector<double> qscore(M, 0.0);
        std::vector<long>   qc_nc(M, 0);     // near-complete bin count per config
        std::vector<double> qc_sum(M, 0.0);  // summed completeness (tie-break)
        size_t best_c = 0; double best_q = -1e300;
        size_t best_qc_c = 0; long best_nc = -1; double best_ncsum = -1.0;
        for (size_t c = 0; c < M; ++c) {
          g_w_comp     = cfgs[c].first;
          g_edge_power = cfgs[c].second;
          auto t0 = clk::now();
          compute_edgescore();
          double es_ms = ms_since(t0);
          auto t1 = clk::now();
          rebuild_incs();
          double inc_ms = ms_since(t1);
          auto t2 = clk::now();
          mems_all[c] = run_lp(seed);
          double lp_ms = ms_since(t2);
          double pass_ms = es_ms + inc_ms + lp_ms;
          tail_total += pass_ms;
          auto bc = count_bins(mems_all[c]);
          qscore[c] = modularity(mems_all[c]);
          if (qscore[c] > best_q) { best_q = qscore[c]; best_c = c; }
          if (qc_ready) {
            std::vector<char> binned(nobs, 0);
            for (size_t i = 0; i < nobs; ++i) binned[i] = !g.incs[i].empty();
            rb_qc_score_membership(mems_all[c], binned, 90.0, 5.0,
                                   qc_nc[c], qc_sum[c]);
            if (qc_nc[c] > best_nc ||
                (qc_nc[c] == best_nc && qc_sum[c] > best_ncsum)) {
              best_nc = qc_nc[c]; best_ncsum = qc_sum[c]; best_qc_c = c;
            }
          }
          if (dump) {
            char path[4096];
            snprintf(path, sizeof(path), "%s.reuse.cfg%zu_a%.2f_p%.2f.binning",
                     base.c_str(), c, cfgs[c].first, cfgs[c].second);
            write_binning(mems_all[c], path);
          }
          fprintf(stderr,
                  "[REUSE] cfg%zu a=%.2f p=%.2f | edgeScore=%.1fms incs=%.1fms "
                  "LP=%.1fms pass=%.1fms | bins=%zu binned=%zu Q=%.5f"
                  " HQ~=%ld\n",
                  c, cfgs[c].first, cfgs[c].second, es_ms, inc_ms, lp_ms,
                  pass_ms, bc.first, bc.second, qscore[c],
                  qc_ready ? qc_nc[c] : -1L);
        }

        // ── Selection: QC near-complete count (--auto/--ensemble) else Q ──────
        // Re-materialise the selected config's edgeScore + incidence so the
        // downstream collect-bins / recruit / split / output stages operate on
        // exactly the chosen partition. membership is the selected LP labelling.
        if (qc_ready) best_c = best_qc_c;
        g_w_comp     = cfgs[best_c].first;
        g_edge_power = cfgs[best_c].second;
        compute_edgescore();
        rebuild_incs();
        membership = mems_all[best_c];
        used_no_gold = true;

        // ── --ensemble: greedy consensus across ALL swept partitions ──────────
        // DAS Tool-style: collect every (config, cluster) as a candidate bin,
        // score each by SCG quality, then greedily accept the best bins whose
        // contigs do not overlap already-accepted bins. Cheap because all
        // partitions came from one graph build. Produces `membership` directly.
        if (g_ensemble && qc_ready) {
          rb_qc_ensemble(mems_all, g, membership);
        }

        if (dump) {
          char path[4096];
          snprintf(path, sizeof(path), "%s.reuse.selected.binning", base.c_str());
          write_binning(membership, path);
        }
        if (qc_ready)
          fprintf(stderr,
                  "[REUSE] SELECTED (max near-complete) = cfg%zu a=%.2f p=%.2f "
                  "HQ~=%ld%s\n",
                  best_c, cfgs[best_c].first, cfgs[best_c].second, best_nc,
                  g_ensemble ? " then ensemble-consensus" : "");
        else
          fprintf(stderr,
                  "[REUSE] SELECTED (max modularity) = cfg%zu a=%.2f p=%.2f "
                  "Q=%.5f\n",
                  best_c, cfgs[best_c].first, cfgs[best_c].second, best_q);

        // ── Pairwise edge-agreement (diversity proxy; lower = more diverse) ────
        if (dump && M >= 2) {
          double agree_sum = 0.0; size_t npair = 0;
          for (size_t a = 0; a < M; ++a)
            for (size_t b = a + 1; b < M; ++b) {
              size_t same = 0;
              const auto &ma = mems_all[a], &mb = mems_all[b];
              for (size_t e = 0; e < E; ++e) {
                bool ja = (ma[g.from[e]] == ma[g.to[e]]);
                bool jb = (mb[g.from[e]] == mb[g.to[e]]);
                if (ja == jb) ++same;
              }
              agree_sum += (double)same / (double)std::max<size_t>(E, 1);
              ++npair;
            }
          fprintf(stderr, "[REUSE] mean pairwise edge-agreement=%.4f "
                  "(1.0=identical, lower=more diverse) over %zu pairs\n",
                  agree_sum / std::max<size_t>(npair, 1), npair);
        }

        fprintf(stderr,
                "[REUSE] SUMMARY: 1 graph build reused for %zu LP passes; "
                "per-pass tail avg=%.1fms (total %.1fms). Selected partition "
                "flows into recruit/split/output below.\n",
                M, tail_total / std::max<size_t>(M, 1), tail_total);
      }  // end --no_gold multi-resolution selection

      // ── 3. Compute depth_matrix graph weights and composite scores ──────────────
      if (!used_no_gold) {
      if (has_depth) {
        verbose_message("Calculating depth_matrix graph [%.1fGb / %.1fGb]               "
                        "                           \n",
                        getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
        size_t ne = g.getEdgeCount();
        g.edgeScore.resize(ne);

#pragma omp parallel for schedule(dynamic, 1)
        for (size_t e = 0; e < ne; ++e) {
          size_t i = g.from[e], j = g.to[e];
          if (num_depth_samples <= 1) {
            double w = (double)g.sComp[e];
            if (g_edge_power != 1.0) w = std::pow(w, g_edge_power);
            g.edgeScore[e] = (StoredDistance)w;
          } else {
            // Depth term: Spearman corr (default) or weighted Jaccard
            // (RABBIT_DEPTH_SIM=wjac). depth_edge_term applies the negative-corr
            // hard-cut for the correlation path; wjac is already in [0,1].
            bool depth_ok;
            double dterm = depth_edge_term(i, j, depth_ok);
            if (!depth_ok) { g.edgeScore[e] = 0.0f; continue; }
            double sComp_v = g.sComp[e];
            double min_edge_weight_v = min_edge_weight;
            double w = g_w_comp * sComp_v + (1.0 - g_w_comp) * dterm;
            if (!std::isfinite(w) || w < min_edge_weight_v) w = 0.0;
            // Edge power: raise to p>1 to de-emphasise borderline edges.
            if (g_edge_power != 1.0 && w > 0.0) w = std::pow(w, g_edge_power);
            g.edgeScore[e] = (StoredDistance)w;
          }
        }
      } else {
        g.edgeScore = g.sComp;
        for (auto &s : g.edgeScore) if (s < min_edge_weight) s = 0.0f;
      }

      rb_phase("graph+edgescore done");

      // ── 4. Build incidence list ────────────────────────────────────────
      // Clamp edgeScore to (0, 1-eps): Boost 1.66 cdf(chi_squared, x) throws when
      // x = -2*LOG(1-edgeScore) is non-finite (LOG(0) = -inf when edgeScore >= 1.0).
      static constexpr StoredDistance SSCR_MAX = 1.0f - 1e-6f;
      g.incs.resize(nobs);
      for (size_t e = 0; e < g.getEdgeCount(); ++e) {
        if (g.edgeScore[e] > 0) {
          if (g.edgeScore[e] > SSCR_MAX) g.edgeScore[e] = SSCR_MAX;
          g.incs[g.from[e]].push_back(e);
          g.incs[g.to[e]].push_back(e);
        }
      }
      // Determinism: the parallel edge build concatenates per-thread lists in a
      // non-reproducible order, so each node's incidence list (and hence the
      // floating-point accumulation order + first-max tie-break inside the
      // serial Gauss-Seidel label propagation) would vary run-to-run. Sorting
      // each incidence list by the neighbouring node id pins a canonical order,
      // making the converged labelling reproducible at a fixed seed/threads.
#pragma omp parallel for schedule(dynamic, 256)
      for (size_t v = 0; v < nobs; ++v) {
        auto &inc = g.incs[v];
        std::sort(inc.begin(), inc.end(), [&](size_t e1, size_t e2) {
          return g.getOtherNode(e1, v) < g.getOtherNode(e2, v);
        });
      }
      // A node is "connected" iff it has ≥1 positive-weight incident edge, i.e.
      // its incidence list is non-empty — so we read that directly instead of
      // maintaining a separate unordered_set (which cost ~2·|E| hash inserts).
      size_t n_connected = 0;
      for (size_t i = 0; i < nobs; ++i) if (!g.incs[i].empty()) ++n_connected;

      // ── 5. Label propagation ───────────────────────────────────────────
      verbose_message("Starting Label Propagation. connected=%zu nobs=%d "
                      "[%.1fGb / %.1fGb]                          \n",
                      n_connected, nobs,
                      getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

      std::vector<size_t> node_order(nobs);
      std::iota(node_order.begin(), node_order.end(), 0);
      std::shuffle(node_order.begin(), node_order.end(), std::default_random_engine(seed));

      cluster_by_propagation(g, membership, node_order);
      rb_phase("label propagation done");
      }  // end !used_no_gold (single-resolution production path)

      // ── 6. Collect bins ────────────────────────────────────────────────
      for (size_t i = 0; i < nobs; ++i) {
        if (!g.incs[i].empty())
          cls[membership[i]].push_back(i);
      }
      mems = membership;
    } // graph g destroyed here

    if (no_recruit) break;

    // ── 7. Recruit lost and small contigs ─────────────────────────────────
    std::vector<size_t> leftovers;
    {
      std::unordered_set<size_t> binned;
      for (auto &kv : cls) for (auto c : kv.second) binned.insert(c);
      for (size_t i = 0; i < nobs; ++i)
        if (binned.find(i) == binned.end()) leftovers.push_back(i);
    }

    verbose_message("Calculating Spearman corr for small and leftover "
                    "contigs [%.1fGb / %.1fGb]\n",
                    getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

    // This rank-transformed copy is consumed ONLY by the abundance-centroid
    // recruit path (recruit_to_depth_centroid, default off).  depth_matrix was
    // already rank-transformed in place during the sketch loop, so building it
    // unconditionally was both a redundant second ranking and a wasted nobs×S
    // matrix in the default path.  Build it only when actually used.
    Matrix spearman;
    if (recruit_to_depth_centroid) {
      // Allocated (matching the original (nobs,S) shape) whenever the centroid
      // reader below runs, so its row access stays in-bounds even for a single
      // depth sample.  The rank fill only applies for multi-sample depth.
      spearman.resize(nobs, num_depth_samples);
      if (num_depth_samples > 1) {
#pragma omp parallel for schedule(dynamic, 1)
        for (size_t r = 0; r < nobs; ++r) {
          auto &rowMat = threadRowMat[omp_get_thread_num()];
          const MatrixRowType rRow(depth_matrix, r);
          std::copy(rRow.begin(), rRow.end(), rowMat.begin());
          rank(rowMat, rowMat);
          MatrixRowType sRow(spearman, r);
          std::copy(rowMat.begin(), rowMat.end(), sRow.begin());
        }
      }
    }
    verbose_message("Calculated %d spearman corr for small and leftover "
                    "contigs\n", nobs);

    // ── Precompute centered + L2-normalised unit depth vectors ─────────────
    // Pearson/Spearman correlation between two contigs equals the dot product
    // of their centered+normalised depth vectors (same identity used by the
    // abdfirst prune in build_similarity_graph).  Replacing the per-pair scalar
    // Welford cal_depth_corr in the recruit loops below with a dot product
    // removes the per-element divisions and the two per-pair sqrt calls.  Built
    // only for multi-sample depth (cal_depth_corr requires num_depth_samples>1);
    // the lambdas fall back to cal_depth_corr when the vectors are absent.
    const uint32_t ABD_S = (uint32_t)num_depth_samples;
    // Large-contig unit vectors already live in the global g_depth_unit (built
    // once before the graph loop); reuse it here instead of rebuilding an
    // identical nobs×S array.  Only the small-contig unit vectors are local.
    std::vector<float> unit_small;
    auto build_unit = [&](std::vector<float>& out, const Matrix& m, size_t rows) {
      out.assign(rows * ABD_S, 0.0f);
#pragma omp parallel for schedule(static)
      for (size_t r = 0; r < rows; ++r) {
        double mean = 0.0;
        for (uint32_t k = 0; k < ABD_S; ++k) mean += m(r, k);
        mean /= ABD_S;
        double ss = 0.0;
        for (uint32_t k = 0; k < ABD_S; ++k) { double d = (double)m(r, k) - mean; ss += d * d; }
        if (ss > 0.0) {
          const double inv = 1.0 / std::sqrt(ss);
          float* u = out.data() + r * ABD_S;
          for (uint32_t k = 0; k < ABD_S; ++k)
            u[k] = (float)(((double)m(r, k) - mean) * inv);
        }
      }
    };
    // Large-large correlation: reuse the global g_depth_unit rows (identical to
    // the old locally-rebuilt unit_large).  Kept as the original float scalar dot
    // (NOT depth_corr_fast) so the accumulation order — and therefore every
    // borderline recruit decision — is bit-identical to before; this commit only
    // removes the redundant nobs×S rebuild, it must not perturb results.
    auto dcorr_ll = [&](size_t a, size_t b) -> double {
      if (g_depth_unit.empty()) return cal_depth_corr(a, b);
      const float* ua = g_depth_unit.data() + a * (size_t)ABD_S;
      const float* ub = g_depth_unit.data() + b * (size_t)ABD_S;
      float c = 0.0f;
      for (uint32_t k = 0; k < ABD_S; ++k) c += ua[k] * ub[k];
      return (double)c;
    };
    auto dcorr_ls = [&](size_t a, size_t s) -> double {
      if (g_depth_unit.empty() || unit_small.empty()) return cal_depth_corr(a, s, true);
      const float* ua = g_depth_unit.data() + a * (size_t)ABD_S;
      const float* us = unit_small.data() + s * (size_t)ABD_S;
      float c = 0.0f;
      for (uint32_t k = 0; k < ABD_S; ++k) c += ua[k] * us[k];
      return (double)c;
    };

    // ── Build cluster → index mapping ─────────────────────────────────────
    std::unordered_map<size_t, size_t> cls_id_to_idx;
    std::vector<size_t> clsIds;
    if (recruitSimFactor > 0.0 || recruit_to_depth_centroid) {
      auto sz = cls.size();
      cls_id_to_idx.reserve(sz);
      clsIds.reserve(sz);
      int idx = 0;
      for (auto &[clsId, contigIds] : cls) {
        cls_id_to_idx[clsId] = idx++;
        clsIds.push_back(clsId);
      }
    }

    // ── Build abundance centroids ─────────────────────────────────────────
    if (recruit_to_depth_centroid) {
      depth_centroids.resize(cls.size(), num_depth_samples, false);
      verbose_message("Calculating centroid abundances of existing %d bins "
                      "[%.1fGb / %.1fGb]\n",
                      cls.size(), getUsedPhysMem(),
                      getTotalPhysMem() / 1024 / 1024);
      for (auto &[clsId, contigIds] : cls) {
        std::vector<double> tmp_centroid(num_depth_samples, 0.0);
        for (auto contigId : contigIds)
          for (auto i = 0; i < (int)num_depth_samples; i++)
            tmp_centroid[i] += spearman(contigId, i);
        for (auto i = 0; i < (int)num_depth_samples; i++)
          depth_centroids(cls_id_to_idx[clsId], i) = tmp_centroid[i];
      }
    }

    // ── Build centroid sketches for composition-based recruitment ──────────
    if (recruitSimFactor > 0.0) {
      auto sz = cls.size();
      g_centroids.assign(sz, nullptr);
      verbose_message("Building centroid sketches for existing %d bins "
                      "[%.1fGb / %.1fGb]\n",
                      sz, getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);

      assert(sz == clsIds.size());
#pragma omp parallel for schedule(dynamic, 1)
      for (int idx = 0; idx < (int)sz; idx++) {
        auto &clsId     = clsIds[idx];
        auto &contigIds = cls[clsId];
        if (contigIds.empty()) continue;
        // Union merge of all contig sketches in the cluster
        rabbit_sketch::KmerSketch merged(*g_sketches[contigIds[0]]);
        for (size_t ci = 1; ci < contigIds.size(); ++ci)
          merged = merged.merge(*g_sketches[contigIds[ci]]);
        g_centroids[idx] = new rabbit_sketch::KmerSketch(std::move(merged));
      }
      // Pre-build centroid signatures (single producer per object) before any
      // concurrent jaccard() in the parallel recruitment loops below.
#pragma omp parallel for schedule(dynamic)
      for (int idx = 0; idx < (int)sz; idx++)
        if (g_centroids[idx]) (void)g_centroids[idx]->getSignature();
    }

    // The reg_[] arrays in g_sketches (and g_centroids) are no longer needed
    // after signatures are built: all subsequent similarity calls use bits_[].
    // Free them now to reduce peak RSS during the recruitment/output phase.
#pragma omp parallel for schedule(static)
    for (int r = 0; r < (int)nobs; ++r)
      if (g_sketches[r]) g_sketches[r]->freeRegisters();
    for (auto *c : g_centroids)
      if (c) c->freeRegisters();

    // ── Compute within-bin mean correlation ───────────────────────────────
    verbose_message("Calculating mean corr within the %d bins which have >= "
                    "minCS(%d) contigs...          \n",
                    cls.size(), minCS);
    std::unordered_map<size_t, StoredDistance> cls_corr;
    std::vector<BinMap::iterator> largeClusters;
    ProgressTracker prog_cls(cls.size());
#pragma omp parallel
#pragma omp single
    for (auto it = cls.begin(); it != cls.end(); ++it) {
      prog_cls.track();
      if (prog_cls.isStepMarker())
        verbose_message("... %s [%.1fGb / %.1fGb]\r",
                        prog_cls.getProgress(),
                        getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
      size_t kk = it->first;
      size_t cs = it->second.size();
      if (cs >= 250) {
#pragma omp critical(APPEND_LARGE_CLS)
        largeClusters.push_back(it);
      } else if (cs >= minCS) {
#pragma omp task
        {
          double corr = 0.;
          const auto &c = it->second;
          for (size_t i = 0; i < cs; ++i)
            for (size_t j = i + 1; j < cs; ++j)
              corr += dcorr_ll(c[i], c[j]);
          StoredDistance x = corr / (cs * (cs - 1) / 2);
#pragma omp critical(CALC_MEAN_CORR)
          cls_corr[kk] = x;
        }
      }
    }

    if (largeClusters.size() > 0) {
      size_t maxSize = 0, totalSize = 0, totalCalcs = 0;
      for (auto it : largeClusters) {
        size_t cs = it->second.size();
        if (cs > maxSize) maxSize = cs;
        totalSize  += cs;
        totalCalcs += cs * (cs - 1) / 2;
      }
      verbose_message("Processing %d large clusters. maxSize=%d avgSize=%f "
                      "totalCalcs=%d\n",
                      largeClusters.size(), maxSize,
                      (double)totalSize / largeClusters.size(), totalCalcs);
      ProgressTracker prog_lg(totalCalcs);
      for (auto it : largeClusters) {
        size_t kk = it->first;
        size_t cs = it->second.size();
        double corr = 0.;
        const auto &c = it->second;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:corr)
        for (size_t i = 0; i < cs; ++i) {
          for (size_t j = i + 1; j < cs; ++j)
            corr += dcorr_ll(c[i], c[j]);
          prog_lg.track(cs - i);
          if (omp_get_thread_num() == 0 && prog_lg.isStepMarker())
            verbose_message(".... %s\r", prog_lg.getProgress());
        }
        StoredDistance x = corr / (cs * (cs - 1) / 2);
        cls_corr[kk] = x;
      }
      verbose_message("Done with large cluster corr calcs                 \n");
    }

    // ── Fast recruit path: flatten eligible clusters into contiguous centroid
    //    unit vectors so each leftover/small recruit becomes a single dot
    //    product (centroid · u_x) instead of a per-member scan over an
    //    unordered_map.  This relies on the identity
    //        mean_corr(x, cluster) = (1/cs) Σ_i u_{c[i]} · u_x = centroid · u_x ,
    //    where centroid is the mean of the member unit depth vectors.  It also
    //    drops the partial-minCS early-exit gate (the centroid yields the exact
    //    full mean directly), so it is enabled only for the default depth-corr
    //    recruit (no abd-centroid, no sketch cutoff) with multi-sample depth
    //    where the unit vectors exist; otherwise the original loops run.  The
    //    per-bin decision (recruit iff full mean >= bin threshold, and unique)
    //    is unchanged, only the inner accumulation is folded into the centroid.
    const bool fast_recruit =
        !recruit_to_depth_centroid && recruitSimFactor <= 0.0 &&
        num_depth_samples > 1 && !g_depth_unit.empty() &&
        !getenv("RABBIT_NO_FAST_RECRUIT");
    std::vector<int>   recruit_clsid;     // eligible cluster ids (cs >= minCS)
    std::vector<float> recruit_centroid;  // numElig × ABD_S mean unit vectors
    std::vector<float> recruit_thresh;    // numElig per-bin corr thresholds
    if (fast_recruit) {
      recruit_clsid.reserve(cls.size());
      for (auto &kv : cls)
        if (kv.second.size() >= minCS) recruit_clsid.push_back(kv.first);
      const size_t numElig = recruit_clsid.size();
      recruit_centroid.assign(numElig * (size_t)ABD_S, 0.0f);
      recruit_thresh.assign(numElig, std::numeric_limits<float>::max());
#pragma omp parallel for schedule(dynamic, 64)
      for (size_t e = 0; e < numElig; ++e) {
        const auto &c = cls.find(recruit_clsid[e])->second;
        const size_t cs = c.size();
        float *cen = recruit_centroid.data() + e * (size_t)ABD_S;
        for (size_t i = 0; i < cs; ++i) {
          const float *ui = g_depth_unit.data() + c[i] * (size_t)ABD_S;
          for (uint32_t k = 0; k < ABD_S; ++k) cen[k] += ui[k];
        }
        const float inv = 1.0f / (float)cs;
        for (uint32_t k = 0; k < ABD_S; ++k) cen[k] *= inv;
        auto itc = cls_corr.find(recruit_clsid[e]);
        if (itc != cls_corr.end()) recruit_thresh[e] = (float)itc->second;
      }
      verbose_message("Flattened %d eligible bins for fast recruit "
                      "[%.1fGb / %.1fGb]\n",
                      numElig, getUsedPhysMem(),
                      getTotalPhysMem() / 1024 / 1024);
    }

    // ── Recruit leftover (un-binned) contigs ──────────────────────────────
    verbose_message("Binning lost contigs over %d leftovers and %d bins...      "
                    "          \n",
                    leftovers.size(), cls.size());
    ProgressTracker lost_progress(leftovers.size());
    BinMap cls_leftovers;

    const double sim_recruit_cutoff = simCutoff * recruitSimFactor / 1000.;
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t l = 0; l < leftovers.size(); ++l) {
      lost_progress.track();
      if (verbose && omp_get_thread_num() == 0 && lost_progress.isStepMarker())
        verbose_message("Finding lost contigs %s\r", lost_progress.getProgress());

      int best_cls = -1;
      if (fast_recruit) {
        const float *ul = g_depth_unit.data() + leftovers[l] * (size_t)ABD_S;
        const size_t numElig = recruit_clsid.size();
        for (size_t e = 0; e < numElig; ++e) {
          const float *cen = recruit_centroid.data() + e * (size_t)ABD_S;
          float corr = 0.0f;
          for (uint32_t k = 0; k < ABD_S; ++k) corr += cen[k] * ul[k];
          if (corr >= recruit_thresh[e]) {
            if (best_cls > -1) { best_cls = -1; break; }
            best_cls = recruit_clsid[e];
          }
        }
      } else
      for (auto it = cls.begin(); it != cls.end(); ++it) {
        size_t kk = it->first;
        const auto &c = it->second;
        size_t cs = c.size();
        if (cs >= minCS) {
          double corr = 0;
          if (recruitSimFactor > 0.0) {
            auto cls_idx = cls_id_to_idx[kk];
            // KmerSketch Jaccard between centroid and leftover contig
            auto sComp = g_centroids[cls_idx]
                          ? (StoredDistance)g_centroids[cls_idx]->jaccard(*g_sketches[leftovers[l]])
                          : (StoredDistance)0.0;
            if (sComp < (StoredDistance)sim_recruit_cutoff)
              continue;
          }
          if (recruit_to_depth_centroid) {
            auto cls_idx = cls_id_to_idx[kk];
            corr = cal_depth_corr(cls_idx, leftovers[l], false, true);
          } else {
            size_t i = 0;
            for (; i < minCS; ++i)
              corr += dcorr_ll(c[i], leftovers[l]);
            if (corr / minCS < cls_corr[kk]) continue;
            for (; i < cs; ++i)
              corr += dcorr_ll(c[i], leftovers[l]);
            corr /= cs;
          }
          if (corr >= cls_corr[kk]) {
            if (best_cls > -1) { best_cls = -1; break; }
            best_cls = kk;
          }
        }
      }
      if (best_cls > -1) {
#pragma omp critical(ADD_LEFTOVER_CONTIGS)
        cls_leftovers[best_cls].push_back(leftovers[l]);
      }
    }

    // Release per-contig sketch memory after leftover recruitment
    if (recruitSimFactor > 0.0) {
      for (auto *p : g_sketches) { delete p; p = nullptr; }
      g_sketches.clear();
      g_sketches.shrink_to_fit();
      verbose_message("Cleaned up sketch memory [%.1fGb / %.1fGb]              "
                      "                                      \n",
                      getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
    }

    // ── Recruit small contigs ─────────────────────────────────────────────
    BinMap cls_small;
    if (nobs1 > 0) {
      verbose_message("Binning %d small contigs...                              "
                      "         \n", nobs1);

      // Snapshot raw small-contig means for marker_guided_split before the
      // in-place rank transform below corrupts small_depth_matrix for that purpose.
      if ((!marker_seed_file.empty() || g_split_abundance) && num_depth_samples >= 1) {
        g_small_means.assign(nobs1 * (size_t)num_depth_samples, 0.0f);
        for (size_t r = 0; r < nobs1; ++r)
          for (size_t i = 0; i < (size_t)num_depth_samples; ++i)
            g_small_means[r * num_depth_samples + i] = (float)small_depth_matrix(r, i);
      }

      if (num_depth_samples > 1) {
#pragma omp parallel for schedule(dynamic, 1)
        for (size_t r = 0; r < small_depth_matrix.size1(); ++r) {
          auto &rowMat = threadRowMat[omp_get_thread_num()];
          MatrixRowType rRow(small_depth_matrix, r);
          std::copy(rRow.begin(), rRow.end(), rowMat.begin());
          rank(rowMat, rowMat);
          std::copy(rowMat.begin(), rowMat.end(), rRow.begin());
        }
        verbose_message("Finished %d spearman corr calcs\n", small_depth_matrix.size1());
      }
      threadRowMat.clear();

      if (num_depth_samples > 1) build_unit(unit_small, small_depth_matrix, nobs1);

      ProgressTracker small_progress(nobs1);
      if (fast_recruit) {
        // Cache-blocked recruit: stream the flat centroid array once per block
        // of SB small contigs (each centroid row, read once, is reused across
        // SB contigs) instead of re-reading the whole centroid array per small
        // contig.  The per-contig decision is unchanged — recruit iff exactly
        // one eligible bin has centroid·u_s >= its threshold — only the loop
        // order is blocked, so the result is identical to the per-contig path.
        const size_t numElig = recruit_clsid.size();
        constexpr size_t SB = 64;
#pragma omp parallel for schedule(dynamic)
        for (size_t s0 = 0; s0 < nobs1; s0 += SB) {
          const size_t s1 = std::min(s0 + SB, nobs1);
          const size_t nb = s1 - s0;
          int     best[SB];
          uint8_t ncand[SB];
          for (size_t t = 0; t < nb; ++t) { best[t] = -1; ncand[t] = 0; }
          for (size_t e = 0; e < numElig; ++e) {
            const float *cen = recruit_centroid.data() + e * (size_t)ABD_S;
            const float  th  = recruit_thresh[e];
            const int    cid = recruit_clsid[e];
            for (size_t t = 0; t < nb; ++t) {
              if (ncand[t] >= 2) continue;        // already ambiguous -> rejected
              const float *us = unit_small.data() + (s0 + t) * (size_t)ABD_S;
              float corr = 0.0f;
              for (uint32_t k = 0; k < ABD_S; ++k) corr += cen[k] * us[k];
              if (corr >= th && ++ncand[t] == 1) best[t] = cid;
            }
          }
          for (size_t t = 0; t < nb; ++t) {
            if (ncand[t] == 1) {
#pragma omp critical(ADD_SMALL_CONTIGS)
              cls_small[best[t]].push_back((s0 + t) + nobs);
            }
          }
          small_progress.track(nb);
          if (verbose && omp_get_thread_num() == 0 && small_progress.isStepMarker())
            verbose_message("Binning small contigs %s\r", small_progress.getProgress());
        }
      } else
#pragma omp parallel for schedule(dynamic)
      for (size_t s = 0; s < nobs1; ++s) {
        small_progress.track();
        if (verbose && omp_get_thread_num() == 0 && small_progress.isStepMarker())
          verbose_message("Binning small contigs %s\r", small_progress.getProgress());

        int best_cls = -1;
        for (auto it = cls.begin(); it != cls.end(); ++it) {
          size_t kk = it->first;
          const auto &c = it->second;
          size_t cs = c.size();
          if (cs >= minCS) {
            double corr = 0;
            if (recruit_to_depth_centroid) {
              auto cls_idx = cls_id_to_idx[kk];
              corr = cal_depth_corr(cls_idx, s, true, true);
            } else {
              size_t i = 0;
              for (; i < minCS; ++i)
                corr += dcorr_ls(c[i], s);
              if (corr / minCS < cls_corr[kk]) continue;
              for (; i < cs; ++i)
                corr += dcorr_ls(c[i], s);
              corr /= cs;
            }
            if (corr >= cls_corr[kk]) {
              if (best_cls > -1) { best_cls = -1; break; }
              if (sim_recruit_cutoff > 0.0) {
                auto cls_idx = cls_id_to_idx[kk];
                // Build a temporary sketch for this small contig
                rabbit_sketch::KmerSketch small_sk(sketch_size, sketch_kmer_size, sketch_bits);
                small_sk.update(small_seqs[s].data(), small_seqs[s].size());
                auto sComp = g_centroids[cls_idx]
                              ? (StoredDistance)g_centroids[cls_idx]->jaccard(small_sk)
                              : (StoredDistance)0.0;
                if (sComp < (StoredDistance)sim_recruit_cutoff) continue;
              }
              best_cls = kk;
            }
          }
        }
        if (best_cls > -1) {
#pragma omp critical(ADD_SMALL_CONTIGS)
          cls_small[best_cls].push_back(s + nobs);
        }
      }
    }

    // Apply leftover recruits
    for (auto it = cls_leftovers.begin(); it != cls_leftovers.end(); ++it) {
      size_t kk = it->first;
      cls[kk].insert(cls[kk].end(), cls_leftovers[kk].begin(), cls_leftovers[kk].end());
    }

    // Apply small contig recruits (at most 15% of total small bases)
    unsigned long long added_sum = 0;
    for (auto it = cls_small.begin(); it != cls_small.end(); ++it) {
      size_t kk = it->first;
      for (auto c : cls_small[kk]) {
        auto idx = c - nobs;
        added_sum += small_seq_lens[idx];
      }
    }
    if (added_sum > 0) {
      Distance fraction = (Distance)added_sum / totalSize1;
      if (fraction < .15) {
        for (auto it = cls_small.begin(); it != cls_small.end(); ++it) {
          size_t kk = it->first;
          cls[kk].insert(cls[kk].end(), cls_small[kk].begin(), cls_small[kk].end());
        }
      } else {
        verbose_message("[Info] Additional binning of small contigs was "
                        "ignored since it was too excessive [%2.2f%% (%lld "
                        "bases) of small (<%d) contigs is > %2.0f%%].\n",
                        fraction * 100, added_sum, minContig, .10 * 100);
      }
    }

  } while (false);
  rb_phase("recruit done");

  // Release centroid sketches
  for (auto *p : g_centroids) delete p;
  g_centroids.clear();

  // Release per-contig sketches (if not already released above)
  for (auto *p : g_sketches) delete p;
  g_sketches.clear();

  verbose_message("Rescuing singleton large contigs\n");
  promote_singleton_bins(cls);

  // ── Phase 2: split contaminated/multi-modal bins ─────────────────────────
  // An explicit --marker-seed takes priority (marker-guided); otherwise the
  // default marker-free abundance split runs (disable with --noAbdSplit).
  if (!marker_seed_file.empty()) {
    verbose_message("Marker-guided bin splitting...\n");
    marker_guided_split(cls);
  } else if (g_split_abundance) {
    verbose_message("Abundance-guided bin splitting (marker-free)...\n");
    abundance_guided_split(cls);
  }

  rb_phase("split done");
  verbose_message("Outputting bins\n");
  output_bins(cls);
  rb_phase("output done");

  verbose_message("Finished\n");
  return 0;
}

// ── implementation modules ─────────────────────────────────────────────
#include "impl/rb_cluster.cpp"
#include "impl/rb_graph.cpp"
#include "impl/rb_abundance.cpp"
#include "impl/rb_output.cpp"

// ── depth subcommand ────────────────────────────────────────────────────────
// `rabbitbin depth --bam-list bams.txt --out depth.tsv` (or repeated --bam /
// positional BAMs) computes a MetaBAT/JGI-format depth table in-process, byte-
// identical to the standalone rabbit_depth tool. --fasta is accepted for CLI
// parity but is optional: the depth table derives contig names+lengths from the
// BAM header.
static int rb_cmd_depth(int ac, char *av[]) {
#ifndef RABBITBIN_FUSE
  cerr << "[Error!] this rabbitbin build has no in-process BAM support "
          "(libdeflate missing); use the standalone rabbit_depth tool.\n";
  return 1;
#else
  std::string fasta, out_path;
  std::vector<std::string> bams;
  std::string bam_list;
  int pctid = 97, min_contig_len = 1, max_edge = 75;
  double min_contig_depth = 0.0;
  size_t threads = 0;
  bool no_variance = false;

  po::options_description desc("rabbitbin depth options", 100, 50);
  desc.add_options()
      ("help,h", "Show help")
      ("fasta", po::value<std::string>(&fasta), "Contig FASTA (optional; for CLI parity)")
      ("assembly,a", po::value<std::string>(&fasta), "Alias for --fasta")
      ("bam", po::value<std::vector<std::string>>(&bams)->multitoken(), "Sorted BAM(s)")
      ("bam-list", po::value<std::string>(&bam_list), "File with one sorted-BAM path per line")
      ("out,o", po::value<std::string>(&out_path), "Output depth TSV (default: stdout)")
      ("output", po::value<std::string>(&out_path), "Alias for --out")
      ("threads,t", po::value<size_t>(&threads)->default_value(0), "Threads (0=all online CPUs)")
      ("percent-identity", po::value<int>(&pctid)->default_value(97), "Min mapped-read percent identity")
      ("min-contig-length", po::value<int>(&min_contig_len)->default_value(1), "Min contig length emitted")
      ("min-contig-depth", po::value<double>(&min_contig_depth)->default_value(0.0), "Min contig depth emitted")
      ("max-edge-bases", po::value<int>(&max_edge)->default_value(75), "Max edge bases trimmed per contig end")
      ("no-variance", po::value<bool>(&no_variance)->zero_tokens(), "Omit per-sample variance columns");

  po::variables_map vm;
  po::positional_options_description pos;
  pos.add("bam", -1);
  po::store(po::command_line_parser(ac, av).options(desc).positional(pos).run(), vm);
  po::notify(vm);

  if (!bam_list.empty()) {
    std::ifstream bl(bam_list);
    if (!bl) { cerr << "[Error!] cannot open --bam-list file: " << bam_list << "\n"; return 1; }
    std::string ln;
    while (std::getline(bl, ln)) {
      size_t b = ln.find_first_not_of(" \t\r\n");
      if (b == std::string::npos) continue;
      size_t e = ln.find_last_not_of(" \t\r\n");
      bams.push_back(ln.substr(b, e - b + 1));
    }
  }

  if (vm.count("help") || bams.empty()) {
    cerr << "\nrabbitbin depth: BAM(s) -> MetaBAT/JGI depth TSV\n\n"
         << "Usage: rabbitbin depth --bam-list bams.txt --out depth.tsv\n"
         << "       rabbitbin depth --fasta contigs.fa --bam s1.bam s2.bam -o depth.tsv\n\n"
         << desc << "\n";
    if (!vm.count("help") && bams.empty())
      cerr << "[Error!] at least one BAM is required (--bam / --bam-list / positional)\n";
    return vm.count("help") ? 0 : 1;
  }

  long onln = sysconf(_SC_NPROCESSORS_ONLN);
  size_t hw = (onln > 0) ? (size_t)onln : 1;
  if (threads == 0) threads = hw; else threads = std::min(threads, hw);

  fprintf(stderr, "[rabbitbin depth] %zu BAM(s), %zu threads, pctid=%d, "
                  "min-contig-length=%d, min-contig-depth=%.3g, max-edge-bases=%d\n",
          bams.size(), threads, pctid, min_contig_len, min_contig_depth, max_edge);

  std::string tsv;
  try {
    tsv = compute_depth_tsv_inmem(bams, (float)pctid / 100.0f, min_contig_len,
                                  (float)min_contig_depth, max_edge,
                                  /*includeEdgeBases=*/false,
                                  /*intraDepthVariance=*/!no_variance, (int)threads);
  } catch (const std::exception &e) {
    cerr << "[Error!] depth computation failed: " << e.what() << "\n";
    return 1;
  }

  if (out_path.empty() || out_path == "-") {
    fwrite(tsv.data(), 1, tsv.size(), stdout);
  } else {
    FILE *f = fopen(out_path.c_str(), "wb");
    if (!f) { cerr << "[Error!] cannot open output: " << out_path << "\n"; return 1; }
    fwrite(tsv.data(), 1, tsv.size(), f);
    fclose(f);
    fprintf(stderr, "[rabbitbin depth] wrote %s (%zu bytes)\n",
            out_path.c_str(), tsv.size());
  }
  return 0;
#endif
}

// ── amber subcommand (fast, multithreaded genome-binning evaluation) ────────
#include "impl/rb_amber.cpp"
// ── qc subcommand (reference-free SCG completeness/contamination) ───────────
#include "impl/rb_qc.cpp"

// main: subcommand dispatch (bin / depth / amber). A bare `rabbitbin -a ... -o`
// (no recognised subcommand) is treated as the legacy bin invocation so all
// existing scripts keep working.
// ═══════════════════════════════════════════════════════════════════════════
int main(int ac, char *av[]) {
  if (ac >= 2) {
    std::string sub = av[1];
    if (sub == "bin")   return rb_cmd_bin(ac - 1, av + 1);
    if (sub == "depth") return rb_cmd_depth(ac - 1, av + 1);
    if (sub == "amber") return rb_cmd_amber(ac - 1, av + 1);
    if (sub == "qc")    return rb_cmd_qc(ac - 1, av + 1);
    if (sub == "refine") return rb_cmd_refine(ac - 1, av + 1);
    if (sub == "--help" || sub == "-h") {
      cerr << "RabbitBin " << version << " (" << DATE << ")\n\n"
           << "Usage: rabbitbin <command> [options]\n\n"
           << "Commands:\n"
           << "  bin     Metagenome binning (default; legacy `rabbitbin -a ... -o ...` also works)\n"
           << "  depth   Summarize sorted BAM(s) into a MetaBAT/JGI depth TSV\n"
           << "  amber   Fast multithreaded genome-binning evaluation (AMBER-compatible)\n"
           << "  qc      Reference-free bin quality (SCG completeness/contamination)\n"
           << "  refine  DAS Tool-style SCG consensus over multiple binnings\n\n"
           << "Run `rabbitbin <command> --help` for command-specific options.\n";
      return 0;
    }
  }
  // Legacy / default: treat all args as bin options.
  return rb_cmd_bin(ac, av);
}
