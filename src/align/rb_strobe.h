#ifndef RB_STROBE_H_
#define RB_STROBE_H_

// Bridge between RabbitBin and the vendored strobealign alignment engine
// (src/align/strobe, MIT, Kristoffer Sahlin). strobealign is both much faster
// and produces higher-quality binning than the legacy minimizer seed-and-extend
// engine, so it is the preferred `rabbitbin map` backend.
//
// The StrobeMapper owns the in-memory reference + strobemer index + SSW aligner
// (built ONCE; reused across all samples -- a rabbitbin-specific win over
// standalone strobealign which re-reads its on-disk index per invocation). Per
// read(-pair) we call strobealign's align_or_map_* to produce SAM text, parse it
// into bam1_t via htslib, and feed RabbitBin's parallel coordinate sort + BGZF
// writer. Output records are emitted in input order.

#include <string>

#include <htslib/sam.h>

#include "rb_align_api.h"  // RbAlignSink

// Opaque handle holding References + StrobemerIndex + Aligner + Chainer + params.
struct StrobeMapper;

// Build the mapper from a reference FASTA (.gz ok). read_len drives strobealign's
// IndexParameters profile (e.g. 150). Returns nullptr on error.
StrobeMapper *rb_strobe_build(const std::string &ref_fasta, int read_len,
                             int threads);
void rb_strobe_free(StrobeMapper *m);

// Build a SAM header (@SQ from the reference). Caller owns the result.
sam_hdr_t *rb_strobe_make_header(const StrobeMapper *m);

// Align reads from r1 (+ optional r2, or interleaved) with the strobealign
// engine. Each finished bam1_t is passed to sink(sink_ctx, b) in input order
// (sink takes ownership) when sink != nullptr; otherwise written to out/hdr.
// Returns 0 on success.
int rb_strobe_run(StrobeMapper *m, const std::string &r1, const std::string &r2,
                  bool interleaved, int threads, samFile *out, sam_hdr_t *hdr,
                  RbAlignSink sink, void *sink_ctx);

#endif  // RB_STROBE_H_
