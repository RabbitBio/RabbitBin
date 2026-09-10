#include "../src/impl/rb_recruit.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

static void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(1);
  }
}

int main() {
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
