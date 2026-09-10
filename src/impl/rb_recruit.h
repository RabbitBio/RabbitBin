#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

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
