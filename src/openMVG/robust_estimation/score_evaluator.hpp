
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_SCORE_EVALUATOR_HPP
#define OPENMVG_ROBUST_SCORE_EVALUATOR_HPP

#include <vector>

namespace openMVG {
namespace robust{

/// Templated Functor class to evaluate a given model over a set of samples.
template<typename Kernel>
class ScorerEvaluator {
public:
  explicit ScorerEvaluator(double threshold) : threshold_(threshold) {}

  template <typename T>
  double Score(const Kernel &kernel,
    const typename Kernel::Model &model,
    const std::vector<T> &samples,
    std::vector<T> *inliers) const
  {
    double cost = 0.0;
    for (size_t j = 0; j < samples.size(); ++j) {
      const double error = kernel.Error(samples[j], model);
      if (error < threshold_) {
        cost += error;
        inliers->push_back(samples[j]);
      } else {
        cost += threshold_;
      }
    }
    return cost;
  }
private:
  double threshold_;
};

} // namespace robust
} // namespace openMVG

#endif // OPENMVG_ROBUST_SCORE_EVALUATOR_HPP
