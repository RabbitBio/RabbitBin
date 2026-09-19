#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace rabbit_abundance {

// Batch independent dot products across columns, rather than horizontally
// reducing one dot product per pair. The fixed batch is a storage layout, not
// an ISA/cache-size choice: ordinary C++ works with any compiler vector width.
// This is ONLY a rejection filter. Every survivor must still be checked with
// the original dot routine before it reaches the candidate-neighbour heaps.
class DotFilter {
public:
  static constexpr size_t width = 32;

  bool build(const float *rows, size_t count, size_t dims, float threshold,
             int threads) {
    data_.clear();
    dims_ = 0;
    if (!rows || !count || !dims || !std::isfinite(threshold)) return false;
    const double eps = std::numeric_limits<float>::epsilon();
    // Both the batch dot and the original SIMD/scalar dot have at most this
    // many rounding operations on a dependency path. The extra operations
    // cover vector-accumulator merges and horizontal reductions.
    const double nerr = (2.0 * dims + 16.0) * eps;
    if (nerr >= 0.5) return false;
    const size_t blocks = count / width + (count % width != 0);
    if (dims > std::numeric_limits<size_t>::max() / width / blocks) return false;
    data_.assign(blocks * dims * width, 0.0f);
    double max_norm_sq = 0.0;
    int bad = 0;
#pragma omp parallel for num_threads(threads) schedule(static) \
    reduction(max : max_norm_sq) reduction(| : bad)
    for (size_t b = 0; b < blocks; ++b) {
      float *dst = data_.data() + b * dims * width;
      for (size_t lane = 0; lane < width && b * width + lane < count; ++lane) {
        const float *src = rows + (b * width + lane) * dims;
        double norm_sq = 0.0;
        for (size_t k = 0; k < dims; ++k) {
          const float x = src[k];
          dst[k * width + lane] = x;
          norm_sq += (double)x * x;
          if (!std::isfinite(x)) bad = 1;
        }
        max_norm_sq = std::max(max_norm_sq, norm_sq);
      }
    }
    if (bad || !std::isfinite(max_norm_sq)) {
      std::vector<float>().swap(data_);
      return false;
    }
    // |dot_new - dot_old| <= 2*gamma_n*||a||*||b|| by Cauchy-Schwarz.
    // Use 4*gamma_n for additional room for norm accumulation and rounding the
    // bound. The absolute term also covers underflow/flush-to-zero. Round the
    // cutoff DOWN, so a float cast cannot make the rejection test stricter.
    const double gamma = nerr / (1.0 - nerr);
    const double slack = 4.0 * gamma * max_norm_sq +
        (4.0 * dims + 32.0) * std::numeric_limits<float>::min();
    cutoff_ = std::nextafter((float)((double)threshold - slack),
                            -std::numeric_limits<float>::infinity());
    dims_ = dims;
    return true;
  }

  // Returned mask includes only lanes that require the ORIGINAL dot test.
  // Padding lanes are masked by the caller, which also handles is_nz and the
  // upper-triangle boundary. Neither PMH scores nor top-k ordering change.
#if defined(__GNUC__) || defined(__clang__)
  // Keep the dot accumulators out of the large graph/heap OpenMP worker's
  // register allocation; otherwise inlining spills them on every sample.
  __attribute__((noinline))
#endif
  uint32_t candidates(const float *left, size_t block) const {
    const float *right = data_.data() + block * dims_ * width;
#if defined(__GNUC__) || defined(__clang__)
    // Keep the 32 independent sums as arithmetic values across samples.
    // An array updated by an inner SIMD loop can be stored and reloaded at
    // every sample boundary. These small value groups avoid that array
    // recurrence without selecting an ISA or changing sample order.
    DotLanes s0{}, s1{}, s2{}, s3{}, s4{}, s5{}, s6{}, s7{};
    for (size_t k = 0; k < dims_; ++k) {
      const float a = left[k];
      const float *r = right + k * width;
      s0 += a * load_lanes(r);
      s1 += a * load_lanes(r + 4);
      s2 += a * load_lanes(r + 8);
      s3 += a * load_lanes(r + 12);
      s4 += a * load_lanes(r + 16);
      s5 += a * load_lanes(r + 20);
      s6 += a * load_lanes(r + 24);
      s7 += a * load_lanes(r + 28);
    }
    // Build the 32-bit survivor mask with integer bit operations, without
    // storing the sums and visiting all 32 lanes through scalar branches.
    // Complementing '<' retains NaNs, exactly as !(sum < cutoff) below does.
    static_assert(width == 32, "dot mask width mismatch");
    const DotLanes threshold = {cutoff_, cutoff_, cutoff_, cutoff_};
    const DotFlags bits = {1u, 2u, 4u, 8u};
    DotFlags mask = (DotFlags)~(s0 < threshold) & bits;
    mask |= (DotFlags)~(s1 < threshold) & (bits << 4);
    mask |= (DotFlags)~(s2 < threshold) & (bits << 8);
    mask |= (DotFlags)~(s3 < threshold) & (bits << 12);
    mask |= (DotFlags)~(s4 < threshold) & (bits << 16);
    mask |= (DotFlags)~(s5 < threshold) & (bits << 20);
    mask |= (DotFlags)~(s6 < threshold) & (bits << 24);
    mask |= (DotFlags)~(s7 < threshold) & (bits << 28);
    return mask[0] | mask[1] | mask[2] | mask[3];
#else
    float sums[width] = {};
    for (size_t k = 0; k < dims_; ++k) {
      const float a = left[k];
#pragma omp simd
      for (size_t lane = 0; lane < width; ++lane)
        sums[lane] += a * right[k * width + lane];
    }
    uint32_t mask = 0;
    for (size_t lane = 0; lane < width; ++lane)
      if (!(sums[lane] < cutoff_)) mask |= uint32_t(1) << lane;
    return mask;
#endif
  }

  size_t bytes() const { return data_.size() * sizeof(float); }

private:
#if defined(__GNUC__) || defined(__clang__)
  // Compiler vector values work on all of its targets, including scalar
  // lowering when vector instructions are unavailable. No architecture
  // intrinsics, alignment assumptions or CPU/cache-size dispatch are used.
  using DotLanes = float __attribute__((vector_size(4 * sizeof(float))));
  using DotFlags = uint32_t __attribute__((vector_size(4 * sizeof(uint32_t))));
  static inline DotLanes load_lanes(const float *p) {
    DotLanes value;
    std::memcpy(&value, p, sizeof(value));
    return value;
  }
#endif

  size_t dims_ = 0;
  float cutoff_ = 0;
  std::vector<float> data_;
};

} // namespace rabbit_abundance
