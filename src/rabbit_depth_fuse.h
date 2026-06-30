#pragma once
// Declaration of the in-memory depth summarizer defined in rabbit_depth.cpp.
// Used by rabbitbin (built with RABBITBIN_FUSE) when BAM input is detected, to
// compute MetaBAT-format depth from sorted BAMs without a temp file, so BAM
// decompression overlaps the binning FASTA/sketch pass. The returned string is
// byte-identical to the standalone rabbit_depth depth file.
#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <htslib/sam.h>

// ── Per-contig SNV / strain extraction (feature #5) ─────────────────────────
// Collected in the SAME pass as depth (no second 88GB BAM scan).  For each
// (sample, contig) we build a transient per-position allele pileup from the
// already-filtered reads, then reduce it to three contig-level scalars; the
// per-position buffer is freed as soon as the contig's reads are consumed, so
// peak memory is O(one contig length) per worker thread.  Variants are called
// reference-free (consensus = majority allele) — exactly the intra-species
// polymorphism signal strain resolution needs, and it avoids needing the
// assembly FASTA resident at depth time.
struct SnvContigStat {
  uint32_t covered = 0;   // positions with depth >= minDepth (eligible sites)
  uint32_t nsnv = 0;      // polymorphic sites (minor-allele freq >= minAf, >=2 reads)
  double pi_sum = 0.0;    // sum over covered sites of per-site nucleotide diversity
};

// One polymorphic position, captured per (sample, contig) for strain resolution.
// maj/alt are the consensus / second-allele nucleotide indices (0..3 = A,C,G,T)
// so the same site reported by different samples can be matched + oriented by
// allele identity.  alt-allele frequency = altCount/total traces the minor
// strain's per-sample relative abundance (the DESMAN/inStrain linkage signal).
struct SnvSite {
  uint32_t pos;
  uint8_t maj;
  uint8_t alt;
  uint32_t altCount;
  uint32_t total;
};

// Output channel filled by compute_depth_tsv_inmem when SNV extraction is on.
// stats is laid out [sampleIdx * n_targets + tid]; names[tid] is the BAM
// reference name so the caller can map tid -> contig row by name.
struct SnvResult {
  bool enabled = false;
  int minDepth = 5;       // min position depth to score a site
  double minAf = 0.05;    // min minor-allele frequency to call a SNV
  int minBaseQ = 0;       // min base quality (0 = ignore quals)
  bool resolve = false;   // also collect per-site allele data (strain resolution)
  int maxSites = 96;      // cap of retained polymorphic sites per (contig,sample)
  int32_t n_targets = 0;
  int num_bams = 0;
  std::vector<std::string> names;       // [n_targets]
  std::vector<SnvContigStat> stats;     // [num_bams * n_targets]
  // Per (sample, contig) capped list of polymorphic sites, indexed b*n_targets+tid.
  // Empty unless resolve == true.  Populated only for polymorphic contigs.
  std::vector<std::vector<SnvSite>> sites;
};

// ── Fused map->depth (rabbitbin map --fused-depth) ──────────────────────────
// One sample's accumulated per-contig depth column + its sampled avg read len.
struct DepthColumn {
  std::shared_ptr<uint64_t[]> depth;  // length = header->n_targets
  int avgRead = 0;
};

// Accumulate one sample's depth column directly from its in-memory, coordinate-
// sorted bam1_t records (as produced by `map`), then the caller frees them.
DepthColumn depth_accumulate_sample(const std::vector<bam1_t *> &recs,
                                    bam_hdr_t *header, float percentIdentity,
                                    int maxEdgeBases, bool includeEdgeBases,
                                    int minMapQual);

// Format the MetaBAT depth TSV from all samples' accumulated columns.
std::string depth_format_table(const std::vector<DepthColumn> &cols,
                               const std::vector<std::string> &sampleNames,
                               bam_hdr_t *header, float percentIdentity,
                               int minContigLength, float minContigDepth,
                               int maxEdgeBases, bool includeEdgeBases,
                               bool intraDepthVariance);

// Structured per-contig depth (means only) — an in-memory alternative to the TSV
// string.  Lets the fused binning path skip formatting millions of contigs to
// text and parsing them back: ~format+parse of the depth round-trip.  Emits the
// SAME contig set as the TSV (length >= minContigLength); means[i] is one float
// per sample for contig names[i] (length lens[i]).
struct DepthMatrixOut {
  std::vector<std::string> names;
  std::vector<int32_t> lens;
  std::vector<std::vector<float>> means; // means[contig][sample]
  // Paired-end cross-contig linkage (feature: PE refinement).  Filled only when
  // collectPELink is requested.  Each entry is (compact_row_a, compact_row_b,
  // count) with a<b, summed over all BAMs.  compact_row_* indexes the same
  // dense contig order as names/lens/means above, so the caller maps a row to a
  // contig by names[row].  A read pair with mates on two different kept contigs
  // is a physical adjacency signal (insert-size scale) — empirically ~98% of
  // such links are intra-genome, so it is a high-precision "same genome" gate.
  std::vector<std::tuple<int32_t, int32_t, uint32_t>> pe_links;
};

// When outCols != nullptr, compute_depth_tsv_inmem fills *outCols directly and
// returns "" (no TSV formatting).  Otherwise it returns the MetaBAT TSV string.
std::string compute_depth_tsv_inmem(const std::vector<std::string> &bamFilePaths,
                                    float percentIdentity, int minContigLength,
                                    float minContigDepth, int maxEdgeBases,
                                    bool includeEdgeBases,
                                    bool intraDepthVariance, int numThreads,
                                    SnvResult *snv = nullptr,
                                    DepthMatrixOut *outCols = nullptr,
                                    int minMapQual = 0,
                                    int dualMapQual = 0,
                                    bool collectPELink = false);

// Generic BAM/CRAM, reference-aware depth (used for CRAM input or when a
// reference is supplied). Same output format as compute_depth_tsv_inmem.
std::string compute_depth_tsv_generic(const std::vector<std::string> &bamFilePaths,
                                      float percentIdentity, int minContigLength,
                                      float minContigDepth, int maxEdgeBases,
                                      bool includeEdgeBases,
                                      bool intraDepthVariance, int numThreads,
                                      const std::string &referenceFasta);
