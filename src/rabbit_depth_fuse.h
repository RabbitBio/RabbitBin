#pragma once
// Declaration of the in-memory depth summarizer defined in rabbit_depth.cpp.
// Used by rabbitbin (built with RABBITBIN_FUSE) when BAM input is detected, to
// compute MetaBAT-format depth from sorted BAMs without a temp file, so BAM
// decompression overlaps the binning FASTA/sketch pass. The returned string is
// byte-identical to the standalone rabbit_depth depth file.
#include <memory>
#include <string>
#include <vector>

#include <htslib/sam.h>

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
                                    bool intraDepthVariance, int numThreads);

// Generic BAM/CRAM, reference-aware depth (used for CRAM input or when a
// reference is supplied). Same output format as compute_depth_tsv_inmem.
std::string compute_depth_tsv_generic(const std::vector<std::string> &bamFilePaths,
                                      float percentIdentity, int minContigLength,
                                      float minContigDepth, int maxEdgeBases,
                                      bool includeEdgeBases,
                                      bool intraDepthVariance, int numThreads,
                                      const std::string &referenceFasta);
