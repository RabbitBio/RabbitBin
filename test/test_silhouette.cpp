#include "../src/impl/rb_silhouette.h"
#include <cstdlib>
#include <iostream>

static double score(const std::vector<double> &x, const std::vector<int> &labels,
                    int k, size_t cap = 600) {
  std::mt19937 rng(42);
  return rb_mean_silhouette(x.size(), labels, k, cap, rng,
      [&](size_t a, size_t b) { return std::abs(x[a] - x[b]); });
}

static void check(const char *name, double actual, double expected) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-12) {
    std::cerr << name << ": expected " << expected << ", got " << actual << '\n';
    std::exit(1);
  }
}

int main() {
  // The former singleton exclusion gives 0.740277..., incorrectly passing .70.
  check("singleton stays in denominator",
        score({0,1,2,3,4,10}, {0,0,0,0,0,1}, 2),
        (0.75 + 29.0/36.0 + 0.8125 + 0.75 + 7.0/12.0) / 6.0);
  check("separated pairs", score({0,1,10,11}, {0,0,1,1}, 2),
        (1.0 - 1.0/10.5 + 1.0 - 1.0/9.5) / 2.0);
  check("negative silhouettes", score({0,1,10,11}, {0,1,0,1}, 2), -0.45);
  check("coincident points", score({0,0,0,0}, {0,0,1,1}, 2), 0.0);
  check("one represented cluster", score({0,1,2,3}, {0,0,0,0}, 2), -1.0);
  check("all singleton clusters", score({0,1,2,3}, {0,1,2,3}, 4), -1.0);
  check("empty subsample", score({0,1,2,3}, {0,0,1,1}, 2, 0), -1.0);
  check("one sampled observation", score({0,1,2,3}, {0,0,1,1}, 2, 1), -1.0);
  // For two point masses, non-singleton observations score 1 and singletons 0.
  const std::vector<int> labels{0,0,0,0,1,1};
  std::vector<size_t> idx{0,1,2,3,4,5};
  std::mt19937 rng(42);
  std::shuffle(idx.begin(), idx.end(), rng);
  size_t ones = 0;
  for (size_t i = 0; i < 4; ++i) ones += labels[idx[i]] == 1;
  const double expected = ones == 0 ? -1.0 : (ones == 1 ? 0.75 : 1.0);
  check("subsample denominator", score({0,0,0,0,10,10}, labels, 2, 4), expected);
  std::cout << "Silhouette numerical regression checks passed\n";
}
