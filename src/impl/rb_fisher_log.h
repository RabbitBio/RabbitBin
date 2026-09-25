#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

// -log Q(m, x/2), for positive integer m. Scale the finite Poisson sum
// around its largest term, rather than starting at exp(-x/2), which can
// underflow even when the full tail is representable. All relative terms
// are <= 1. No probability clipping or empirical constants are required.
inline double rb_fisher_neg_log_sf(int m, double x) {
  if (m <= 0 || !(x >= 0.0))
    return std::numeric_limits<double>::quiet_NaN();
  const double u = 0.5 * x;
  if (u == 0.0 || std::isinf(u) || m == 1) return u;
  const int peak = u >= (double)(m - 1) ? m - 1 : (int)u;
  double sum = 1.0, term = 1.0;
  for (int k = peak; k > 0; --k) {
    term *= (double)k / u;
    sum += term;
  }
  term = 1.0;
  for (int k = peak + 1; k < m; ++k) {
    term *= u / (double)k;
    sum += term;
  }
  const double result = u - peak * std::log(u) +
                        std::lgamma((double)peak + 1.0) - std::log(sum);
  return std::max(0.0, result);
}

// A one-representable-double tie in log space; numerical, not a tunable
// biological/coverage threshold. The production CDF retains its old rule.
inline bool rb_fisher_log_tie(double current, double best) {
  return current >= std::nextafter(best, -std::numeric_limits<double>::infinity());
}
