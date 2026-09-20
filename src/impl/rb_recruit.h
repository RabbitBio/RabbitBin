#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

inline double rb_recruit_confidence(double candidate, double competitor) {
  if (!std::isfinite(candidate) || !std::isfinite(competitor))
    return -std::numeric_limits<double>::infinity();
  const double candidate_distance = std::max(0.0, 1.0 - candidate);
  const double competitor_distance = std::max(0.0, 1.0 - competitor);
  if (candidate_distance == 0.0)
    return competitor_distance > 0.0
               ? std::numeric_limits<double>::infinity() : 0.0;
  if (competitor_distance == 0.0)
    return -std::numeric_limits<double>::infinity();
  return std::log(competitor_distance / candidate_distance);
}

struct RbRecruitEvidence {
  // 1: source core wins; -1: a wrong core wins; 0: tied.
  int outcome = 0;
  double winner = 0.0;
  double source_counterfactual = 0.0;
  double wrong_counterfactual = 0.0;
};

struct RbRecruitThreshold {
  double value = std::numeric_limits<double>::infinity();
  double tpr = 0.0;
  double fpr = 0.0;
  double youden = 0.0;
};

// Select the largest threshold among the ROC points that maximize Youden's J
// while satisfying FPR <= max_fpr.  The ROC origin is included, so a score
// with no positive discrimination leaves the threshold at +infinity and
// disables recruitment.
inline RbRecruitThreshold rb_select_recruit_threshold(
    const std::vector<double> &positives,
    const std::vector<double> &negatives, double max_fpr) {
  RbRecruitThreshold result;
  if (positives.empty() || negatives.empty()) return result;

  struct RocPoint { double value; bool positive; };
  std::vector<RocPoint> points;
  points.reserve(positives.size() + negatives.size());
  size_t positive_count = 0, negative_count = 0;
  for (double value : positives)
    if (!std::isnan(value)) {
      points.push_back({value, true});
      ++positive_count;
    }
  for (double value : negatives)
    if (!std::isnan(value)) {
      points.push_back({value, false});
      ++negative_count;
    }
  if (positive_count == 0 || negative_count == 0) return result;
  std::sort(points.begin(), points.end(),
            [](const RocPoint &a, const RocPoint &b) {
              return a.value > b.value;
            });

  size_t true_positive = 0, false_positive = 0;
  long double best_youden_numerator = 0.0L;
  for (size_t i = 0; i < points.size();) {
    size_t next = i;
    while (next < points.size() && points[next].value == points[i].value) {
      if (points[next].positive) ++true_positive;
      else ++false_positive;
      ++next;
    }
    const double tpr = (double)true_positive / positive_count;
    const double fpr = (double)false_positive / negative_count;
    const double youden = tpr - fpr;
    // All points have the same positive_count * negative_count denominator.
    // Comparing this integer-valued numerator makes exact Youden ties
    // independent of floating-point rounding.
    const long double youden_numerator =
        (long double)true_positive * negative_count -
        (long double)false_positive * positive_count;
    // Points are visited from largest to smallest threshold.  A strict update
    // therefore implements the conservative largest-threshold tie break.
    if (fpr <= max_fpr && youden_numerator > best_youden_numerator) {
      result.value = points[i].value;
      result.tpr = tpr;
      result.fpr = fpr;
      result.youden = youden;
      best_youden_numerator = youden_numerator;
    }
    i = next;
  }
  return result;
}

// best_other is sufficient when there are exactly two cores.  If a wrong core
// wins, its runner-up is the stronger of the source core and second wrong core;
// in the two-core case the source core is the only runner-up.
inline bool rb_make_recruit_evidence(double own, double best_other,
                                     double second_other,
                                     RbRecruitEvidence &evidence) {
  if (!std::isfinite(own) || !std::isfinite(best_other)) return false;
  const double wrong_competitor = std::isfinite(second_other)
      ? std::max(own, second_other) : own;
  evidence.source_counterfactual =
      rb_recruit_confidence(own, best_other);
  evidence.wrong_counterfactual =
      rb_recruit_confidence(best_other, wrong_competitor);
  if (own > best_other) {
    evidence.outcome = 1;
    evidence.winner = evidence.source_counterfactual;
  } else if (best_other > own) {
    evidence.outcome = -1;
    evidence.winner = evidence.wrong_counterfactual;
  } else {
    evidence.outcome = 0;
    evidence.winner = 0.0;
  }
  return true;
}
