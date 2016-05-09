
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_MAX_CONSENSUS_H_
#define OPENMVG_ROBUST_ESTIMATION_MAX_CONSENSUS_H_

#include "openMVG/robust_estimation/rand_sampling.hpp"
#include <limits>
#include <vector>

namespace openMVG {
namespace robust{

/// Naive implementation of RANSAC without noise and iteration reduction options
/// Pick max_iteration times N_samples and Fit a solution.
/// Keep the solution that have the most of inlier
///  under the confidence threshold of the Scorer.
///
/// Requirements :
/// 1. The model type.
/// 2. The minimum number of samples needed to fit.
/// 3. A way to convert samples to a model.
/// 4. A way to convert samples and a model to an error.
///
/// 1. Kernel::Model
/// 2. Kernel::MINIMUM_SAMPLES
/// 3. Kernel::Fit(vector<int>, vector<Kernel::Model> *)
/// 4. Kernel::Error(Model, int) -> error
template<typename Kernel, typename Scorer>
typename Kernel::Model MaxConsensus(const Kernel &kernel,
  const Scorer &scorer,
  std::vector<size_t> *best_inliers = nullptr , size_t max_iteration = 1024) {

    const size_t min_samples = Kernel::MINIMUM_SAMPLES;
    const size_t total_samples = kernel.NumSamples();

    size_t best_num_inliers = 0;
    typename Kernel::Model best_model;

    // Test if we have sufficient points to for the kernel.
    if (total_samples < min_samples) {
      if (best_inliers) {
        best_inliers->resize(0);
      }
      return best_model;
    }

    // In this robust estimator, the scorer always works on all the data points
    // at once. So precompute the list ahead of time.
    std::vector<size_t> all_samples(total_samples);
    std::iota(all_samples.begin(), all_samples.end(), 0);

    std::vector<size_t> sample;
    for (size_t iteration = 0;  iteration < max_iteration; ++iteration) {
      UniformSample(min_samples, &all_samples, &sample);

        std::vector<typename Kernel::Model> models;
        kernel.Fit(sample, &models);

        // Compute costs for each fit.
        for (size_t i = 0; i < models.size(); ++i) {
          std::vector<size_t> inliers;
          scorer.Score(kernel, models[i], all_samples, &inliers);

          if (best_num_inliers < inliers.size()) {
            best_num_inliers = inliers.size();
            best_model = models[i];
            if (best_inliers) {
              best_inliers->swap(inliers);
            }
          }
        }
    }
    return best_model;
}

} // namespace robust
} // namespace openMVG
#endif // OPENMVG_ROBUST_ESTIMATION_MAX_CONSENSUS_H_
