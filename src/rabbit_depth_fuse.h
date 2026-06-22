#pragma once
// Declaration of the in-memory depth summarizer defined in rabbit_depth.cpp.
// Used by the fused rabbitbin_fuse build (compiled with RABBITBIN_FUSE) to
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
