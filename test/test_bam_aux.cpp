#include "src/impl/rb_bam_aux.h"

#include <htslib/hts_log.h>
#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

using Bytes = std::vector<uint8_t>;
static size_t checked = 0;

static void check(const Bytes &aux) {
  // An unmapped, zero-length record with a padded query name. Extra allocated
  // bytes protect older HTSlib versions when examining truncated tag headers;
  // l_data still ends exactly at the requested malformed input boundary.
  Bytes data{'r', '0', 0, 0};
  data.insert(data.end(), aux.begin(), aux.end());
  const size_t end = data.size();
  data.resize(end + 32, 0);
  bam1_t record{};
  record.core.l_qname = 4;
  record.core.flag = BAM_FUNMAP;
  record.data = data.data();
  record.l_data = static_cast<int>(end);
  errno = E2BIG;
  uint8_t *expected = bam_aux_get(&record, "NM");
  const int expected_errno = errno;
  errno = E2BIG;
  uint8_t *actual = rabbit_bam::aux_nm(&record);
  const int actual_errno = errno;
  if (actual != expected || actual_errno != expected_errno) {
    std::cerr << "NM lookup mismatch, case " << checked
              << ", aux bytes=" << aux.size()
              << ", expected errno=" << expected_errno
              << ", actual errno=" << actual_errno << '\n';
    std::exit(1);
  }
  ++checked;
}

static void append(Bytes &out, uint64_t value, size_t width) {
  for (size_t i = 0; i < width; ++i) out.push_back(value >> (8 * i));
}

static Bytes field(std::mt19937 &rng, bool nm) {
  Bytes out{uint8_t(nm ? 'N' : 'x'), uint8_t(nm ? 'M' : 'y')};
  static const uint8_t types[] = {'A','c','C','s','S','i','I','f','d','Z','H','B'};
  static const size_t widths[] = {1,1,1,2,2,4,4,4,8};
  const size_t type = rng() % 12;
  out.push_back(types[type]);
  if (type < 9) {
    append(out, uint64_t(rng()) << 32 | rng(), widths[type]);
  } else if (type < 11) {
    for (size_t i = 0, n = rng() % 100; i < n; ++i)
      out.push_back("0123456789ABCDEF"[rng() % 16]);
    out.push_back(0);
  } else {
    // Integer and float arrays force the compatibility path, including empty
    // arrays and NM encoded as an array rather than a scalar.
    const size_t subtype = 1 + rng() % 7;
    out.push_back(types[subtype]);
    const size_t count = rng() % 10;
    append(out, count, 4);
    for (size_t i = 0; i < count; ++i) append(out, rng(), widths[subtype]);
  }
  return out;
}

int main() {
  hts_set_log_level(HTS_LOG_OFF); // Expected diagnostics from malformed cases.
  std::mt19937 rng(20260919);
  check({});
  for (size_t round = 0; round < 20000; ++round) {
    Bytes aux;
    const size_t count = rng() % 24;
    for (size_t tag = 0; tag < count; ++tag) {
      const Bytes encoded = field(rng, rng() % 5 == 0);
      aux.insert(aux.end(), encoded.begin(), encoded.end());
    }
    check(aux);
  }
  // All truncation positions, unknown types, malformed arrays, unterminated
  // matched/unmatched strings, duplicates, and a large variable-size prefix.
  for (const Bytes &encoded : std::vector<Bytes>{
         {'N','M','i',1,0,0,0}, {'x','y','i',1,0,0,0,'N','M','C',7},
         {'N','M','Z','a','b',0}, {'x','y','Z','a','b',0,'N','M','C',9},
         {'x','y','B','I',2,0,0,0,1,0,0,0,2,0,0,0,'N','M','C',1},
         {'x','y','B','I',0,0,16,0}, {'x','y','B','Q',1,0,0,0,0},
         {'x','y','Q',0,'N','M','C',1}, {'N','M','Q',0},
         {'N','M','C',3,'N','M','C',8}}) {
    for (size_t len = 0; len <= encoded.size(); ++len)
      check(Bytes(encoded.begin(), encoded.begin() + len));
  }
  Bytes large{'x','y','Z'};
  large.resize(100003, 'x');
  large.insert(large.end(), {0,'N','M','I',255,255,255,255});
  check(large);
  std::cout << checked << " NM lookups match HTSlib, including errno\n";
}
