#include "../../src/impl/rb_fisher_log.h"
#include <boost/math/special_functions/gamma.hpp>
#include <cassert>
#include <cstdio>

int main() {
  size_t cases = 0;
  double max_error = 0.0;
  // Covers weak scores, the mode, CDF saturation, exp(-u) underflow, and
  // the full production max-degree range. Reference uses long double.
  for (int m = 1; m <= 200; ++m) {
    for (double u : {0.0, 1e-12, 0.01, 0.5, 1.0, 10.0, 30.0, 40.0,
                     100.0, 500.0, 744.0, 746.0, 1000.0, 3000.0,
                     m * 0.5, m * 1.25, m * 13.82}) {
      const double actual = rb_fisher_neg_log_sf(m, 2.0 * u);
      const long double q = boost::math::gamma_q((long double)m, (long double)u);
      const double expected = (double)(-std::log(q));
      const double error = std::fabs(actual - expected);
      max_error = std::max(max_error, error);
      // Reference comparison tolerance is a test bound, not an algorithm knob.
      if (!std::isfinite(actual) || error > 2e-11 * (1.0 + expected)) {
        std::fprintf(stderr, "m=%d u=%.17g actual=%.17g expected=%.17g\n",
                     m, u, actual, expected);
        return 1;
      }
      ++cases;
    }
  }
  assert(std::isinf(rb_fisher_neg_log_sf(3, INFINITY)));
  assert(rb_fisher_log_tie(10.0, 10.0));
  assert(!rb_fisher_log_tie(10.0, 10.001));
  std::printf("%zu log-tail cases match long-double gamma_q; max_abs_error=%.3g\n",
              cases, max_error);
}
