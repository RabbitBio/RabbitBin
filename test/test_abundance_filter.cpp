#include "src/impl/rb_abundance_filter.h"

#include <iostream>
#include <random>
#include <stdexcept>

static void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}

// Independent scalar and double oracles. Exercise arbitrary unit vectors,
// zero rows, duplicate/opposite rows and thresholds on either side of an
// actual dot product, including the partial final block.
int main() {
  std::mt19937 rng(219);
  std::normal_distribution<float> normal;
  size_t checked = 0, rejected = 0;
  for (size_t dims : {3, 4, 7, 8, 15, 16, 23, 32, 46, 64, 127, 128, 129, 200, 257}) {
    const size_t count = 67;
    std::vector<float> rows(count * dims);
    for (size_t i = 1; i < count; ++i) {
      double norm = 0.0;
      for (size_t k = 0; k < dims; ++k) {
        float &x = rows[i * dims + k];
        x = normal(rng);
        norm += (double)x * x;
      }
      for (size_t k = 0; k < dims; ++k) rows[i * dims + k] /= std::sqrt(norm);
    }
    for (size_t k = 0; k < dims; ++k) {
      rows[2 * dims + k] = rows[dims + k];
      rows[3 * dims + k] = -rows[dims + k];
    }
    float boundary = 0.0f;
    for (size_t k = 0; k < dims; ++k) boundary += rows[dims + k] * rows[4 * dims + k];
    for (float threshold : {0.0f, 0.7152319f, 0.99f, -0.5f, boundary,
                             std::nextafter(boundary, -1.0f),
                             std::nextafter(boundary, 1.0f)}) {
      rabbit_abundance::DotFilter filter;
      require(filter.build(rows.data(), count, dims, threshold, 1), "filter build failed");
      for (size_t i = 0; i < count; ++i) {
        for (size_t b = 0; b * filter.width < count; ++b) {
          const uint32_t mask = filter.candidates(rows.data() + i * dims, b);
          for (size_t lane = 0; lane < filter.width && b * filter.width + lane < count; ++lane) {
            const size_t j = b * filter.width + lane;
            float scalar = 0.0f;
            double exact = 0.0;
            for (size_t k = 0; k < dims; ++k) {
              const float a = rows[i * dims + k], x = rows[j * dims + k];
              scalar += a * x;
              exact += (double)a * x;
            }
            ++checked;
            if (!(mask & (uint32_t(1) << lane))) {
              ++rejected;
              require(scalar < threshold && exact < threshold, "false-negative rejection");
            }
          }
        }
      }
    }
    rabbit_abundance::DotFilter invalid;
    rows[0] = std::numeric_limits<float>::quiet_NaN();
    require(!invalid.build(rows.data(), count, dims, 0.7f, 1), "NaN must fall back");
    rows[0] = std::numeric_limits<float>::infinity();
    require(!invalid.build(rows.data(), count, dims, 0.7f, 1), "infinity must fall back");
  }
  require(rejected > checked / 4, "filter unexpectedly does no useful work");
  std::cout << "Checked " << checked << " pairs, " << rejected
            << " conservative rejections; no false negatives\n";
}
