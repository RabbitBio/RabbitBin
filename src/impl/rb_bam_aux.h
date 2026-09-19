#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <htslib/sam.h>

namespace rabbit_bam {

// Depth only needs NM. Traverse ordinary scalar/string tags in one inlined
// loop instead of calling bam_aux_next for every preceding tag. Bounds and
// string termination are checked before returning a pointer. Arrays, unknown
// types, missing NM and malformed data go through HTSlib so its diagnostics,
// errno and compatibility behavior remain authoritative. No tag order or
// position is assumed, and duplicate NM tags still select the first one.
inline uint8_t *aux_nm(const bam1_t *record) {
  uint8_t *p = bam_get_aux(record);
  uint8_t *end = record->data + record->l_data;
  while (end - p >= 4) {
    size_t size;
    switch (p[2]) {
    case 'A': case 'c': case 'C': size = 1; break;
    case 's': case 'S': size = 2; break;
    case 'i': case 'I': case 'f': size = 4; break;
    case 'd': size = 8; break;
    case 'Z': case 'H': {
      const uint8_t *nul = static_cast<const uint8_t *>(
          std::memchr(p + 3, 0, static_cast<size_t>(end - (p + 3))));
      if (!nul) return bam_aux_get(record, "NM");
      size = static_cast<size_t>(nul - (p + 3)) + 1;
      break;
    }
    default:
      return bam_aux_get(record, "NM");
    }
    if (size > static_cast<size_t>(end - (p + 3)))
      return bam_aux_get(record, "NM");
    if (p[0] == 'N' && p[1] == 'M') return p + 2;
    p += 3 + size;
  }
  return bam_aux_get(record, "NM");
}

} // namespace rabbit_bam
