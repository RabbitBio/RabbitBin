#pragma once
// Declaration of the in-memory depth summarizer defined in rabbit_depth.cpp.
// Used by rabbitbin (built with RABBITBIN_FUSE) when BAM input is detected, to
// compute MetaBAT-format depth from sorted BAMs without a temp file, so BAM
// decompression overlaps the binning FASTA/sketch pass. The returned string is
// byte-identical to the standalone rabbit_depth depth file.
#include <string>
#include <vector>

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
