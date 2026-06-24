#ifndef RB_ALIGN_API_H_
#define RB_ALIGN_API_H_

// Internal API of the CPU aligner, used by `rabbitbin map` to drive alignment
// and receive bam1_t records via a sink callback (no SAM serialization).

#include <string>

#include <htslib/sam.h>

#include "sa_index.h"
#include "sa_map.h"

// Sink receives ownership of each finished bam1_t (must bam_destroy1 it).
typedef void (*RbAlignSink)(void *ctx, bam1_t *b);

// Build a SAM header (@SQ from index contigs). Caller owns the result.
sam_hdr_t *rb_align_make_header(const SaIndex &idx);

// Align reads. If `sink` is non-null, each record is passed to sink(sink_ctx,b)
// in input order (sink owns it); otherwise records are written to `out`/`hdr`.
int rb_align_run(const SaIndex &idx, const SaOpt &opt, const std::string &r1,
                 const std::string &r2, bool interleaved, int threads,
                 samFile *out, sam_hdr_t *hdr, RbAlignSink sink,
                 void *sink_ctx);

#endif  // RB_ALIGN_API_H_
