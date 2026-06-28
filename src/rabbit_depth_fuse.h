#pragma once
// Declaration of the in-memory depth summarizer defined in rabbit_depth.cpp.
// Used by rabbitbin (built with RABBITBIN_FUSE) when BAM input is detected, to
// compute MetaBAT-format depth from sorted BAMs without a temp file, so BAM
// decompression overlaps the binning FASTA/sketch pass. The returned string is
// byte-identical to the standalone rabbit_depth depth file.
#include <cstdint>
#include <memory>
#include <string>
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

// Output channel filled by compute_depth_tsv_inmem when SNV extraction is on.
// stats is laid out [sampleIdx * n_targets + tid]; names[tid] is the BAM
// reference name so the caller can map tid -> contig row by name.
struct SnvResult {
  bool enabled = false;
  int minDepth = 5;       // min position depth to score a site
  double minAf = 0.05;    // min minor-allele frequency to call a SNV
  int minBaseQ = 0;       // min base quality (0 = ignore quals)
  int32_t n_targets = 0;
  int num_bams = 0;
  std::vector<std::string> names;       // [n_targets]
  std::vector<SnvContigStat> stats;     // [num_bams * n_targets]
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

std::string compute_depth_tsv_inmem(const std::vector<std::string> &bamFilePaths,
                                    float percentIdentity, int minContigLength,
                                    float minContigDepth, int maxEdgeBases,
                                    bool includeEdgeBases,
                                    bool intraDepthVariance, int numThreads,
                                    SnvResult *snv = nullptr);

// Generic BAM/CRAM, reference-aware depth (used for CRAM input or when a
// reference is supplied). Same output format as compute_depth_tsv_inmem.
std::string compute_depth_tsv_generic(const std::vector<std::string> &bamFilePaths,
                                      float percentIdentity, int minContigLength,
                                      float minContigDepth, int maxEdgeBases,
                                      bool includeEdgeBases,
                                      bool intraDepthVariance, int numThreads,
                                      const std::string &referenceFasta);
