#ifndef RB_BAM_SORT_H_
#define RB_BAM_SORT_H_

// Shared coordinate-sort + BAM writer used by `rabbitbin sortbam` and reused by
// `rabbitbin map` (which feeds records straight from the aligner, no temp BAM).

#include <string>
#include <vector>

#include <htslib/sam.h>

// Stable coordinate sort (samtools-equivalent key + stable tie-break).
void rb_stable_coord_sort(std::vector<bam1_t *> &recs);

// Write already-sorted records as BAM to out_path ("-"/empty = stdout).
// Sets @HD SO:coordinate, optionally builds <out>.bai. Returns 0 on success.
int rb_write_sorted_bam(sam_hdr_t *hdr, std::vector<bam1_t *> &recs,
                        const std::string &out_path, int threads, int level,
                        bool write_index);

#endif  // RB_BAM_SORT_H_
