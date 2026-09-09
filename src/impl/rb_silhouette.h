#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

// Mean over every sampled observation. Singletons and zero-distance cases
// contribute zero, including to the denominator. Return -1 for a sample on
// which silhouette is undefined (fewer than two clusters, or all singletons).
template <typename Distance>
static double rb_mean_silhouette(size_t n, const std::vector<int> &labels, int k,
                                 size_t sample_cap, std::mt19937 &rng,
                                 Distance dist) {
  if (n < 3 || k < 2) return -1.0;
  std::vector<size_t> idx(n);
  for (size_t i = 0; i < n; ++i) idx[i] = i;
  if (n > sample_cap) {
    std::shuffle(idx.begin(), idx.end(), rng);
    idx.resize(sample_cap);
  }
  const size_t m = idx.size();
  std::vector<int> csz(k, 0);
  for (size_t a : idx) ++csz[labels[a]];
  const size_t represented = (size_t)std::count_if(
      csz.begin(), csz.end(), [](int count) { return count > 0; });
  if (represented < 2 || represented >= m) return -1.0;

  // Preserve the production upper-triangle distance accumulation order.
  std::vector<double> sumd(m * (size_t)k, 0.0);
  for (size_t ii = 0; ii < m; ++ii) {
    const size_t a = idx[ii];
    for (size_t jj = ii + 1; jj < m; ++jj) {
      const size_t b = idx[jj];
      const double dij = dist(a, b);
      sumd[ii * (size_t)k + (size_t)labels[b]] += dij;
      sumd[jj * (size_t)k + (size_t)labels[a]] += dij;
    }
  }

  double sil_sum = 0.0;
  for (size_t ii = 0; ii < m; ++ii) {
    const int la = labels[idx[ii]];
    if (csz[la] <= 1) continue;  // s=0; this row still counts in m
    const double *row_sumd = sumd.data() + ii * (size_t)k;
    const double ai = row_sumd[la] / (double)(csz[la] - 1);
    double bi = std::numeric_limits<double>::infinity();
    for (int c = 0; c < k; ++c) {
      if (c == la || csz[c] == 0) continue;
      const double mc = row_sumd[c] / (double)csz[c];
      if (mc < bi) bi = mc;
    }
    const double denom = std::max(ai, bi);
    if (denom > 0.0) sil_sum += (bi - ai) / denom;
  }
  return sil_sum / (double)m;
}
