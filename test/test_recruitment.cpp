#include "../src/impl/rb_recruit.h"
#include "../src/impl/rb_coverage.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

static void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(1);
  }
}

int main() {
  auto check_ratio = [](std::vector<double> a, std::vector<double> b,
                        double expected) {
    const double score = rb_mean_coverage_ratio(a.size(),
        [&](size_t s) { return a[s]; }, [&](size_t s) { return b[s]; });
    require(std::abs(score - expected) < 1e-12,
            "unexpected per-sample coverage ratio");
  };
  check_ratio({2.0}, {1.8}, 0.9);
  check_ratio({2.0}, {0.5}, 0.25);
  check_ratio({10.0, 1.0}, {9.0, 0.5}, 0.7);
  check_ratio({1.0, 10.0}, {0.9, 5.0}, 0.7); // library rescaling cancels
  check_ratio({10.0, 9.0}, {9.0, 10.0}, 0.9); // opposite ranks still match
  check_ratio({10.0, 0.0}, {9.0, 0.0}, 0.9); // joint absence is uninformative
  check_ratio({10.0, 0.0}, {9.0, 5.0}, 0.45);
  check_ratio({0.0, 0.0}, {0.0, 0.0}, 0.0);
  check_ratio({0.0, 0.0}, {5.0, 0.0}, 0.0);
  check_ratio({std::nan("")}, {5.0}, 0.0);

  RbRecruitEvidence evidence;
  const double no_second_wrong = -std::numeric_limits<double>::infinity();

  require(rb_make_recruit_evidence(0.9, 0.4, no_second_wrong, evidence),
          "two-core evidence was rejected");
  require(evidence.outcome == 1 && evidence.winner > 0.0 &&
              evidence.wrong_counterfactual < 0.0,
          "clean two-core member did not provide fallback negative evidence");

  require(rb_make_recruit_evidence(0.4, 0.9, no_second_wrong, evidence),
          "misclassified two-core evidence was rejected");
  require(evidence.outcome == -1 && evidence.winner > 0.0 &&
              evidence.source_counterfactual < 0.0,
          "misclassified two-core member lost its winner or fallback positive");

  require(rb_make_recruit_evidence(0.8, 0.7, 0.75, evidence),
          "three-core evidence was rejected");
  require(evidence.outcome == 1 && evidence.wrong_counterfactual < 0.0,
          "strongest wrong core did not use its strongest competitor");

  require(!rb_make_recruit_evidence(
              std::nan(""), 0.7, 0.6, evidence),
          "invalid source-core score was accepted");

  std::cout << "Recruitment calibration regression checks passed\n";
}
