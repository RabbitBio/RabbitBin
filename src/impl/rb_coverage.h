#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

// Equal weight for each informative coverage dimension. Joint absence is not
// evidence of similarity: skip 0/0, but retain a zero contribution when only
// one side is absent. Profiles with no informative dimensions get zero support.
template <class Left, class Right>
inline double rb_mean_coverage_ratio(std::size_t dimensions, Left left, Right right) {
  double sum = 0.0;
  std::size_t informative = 0;
  for (std::size_t s = 0; s < dimensions; ++s) {
    const double x = left(s), y = right(s);
    if (!std::isfinite(x) || !std::isfinite(y) || x < 0.0 || y < 0.0)
      return 0.0;
    const double maximum = std::max(x, y);
    if (maximum == 0.0) continue;
    sum += std::min(x, y) / maximum;
    ++informative;
  }
  return informative ? sum / static_cast<double>(informative) : 0.0;
}
