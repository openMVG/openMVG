
// Copyright (c) 2007, 2008 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_SCORE_EVALUATOR_H_
#define OPENMVG_ROBUST_ESTIMATION_SCORE_EVALUATOR_H_

#include <vector>

namespace openMVG {
namespace robust{

/// Templated Functor class to evaluate a given model over a set of samples.
template<typename Kernel>
class ScorerEvaluator {
public:
  ScorerEvaluator(double threshold) : threshold_(threshold) {}

  template <typename T>
  double Score(const Kernel &kernel,
    const typename Kernel::Model &model,
    const std::vector<T> &samples,
    std::vector<T> *inliers) const
  {
    double cost = 0.0;
    for (size_t j = 0; j < samples.size(); ++j) {
      double error = kernel.Error(samples[j], model);
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

#endif // OPENMVG_ROBUST_ESTIMATION_SCORE_EVALUATOR_H_
