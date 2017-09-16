// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_MAX_CONSENSUS_HPP
#define OPENMVG_ROBUST_ESTIMATION_MAX_CONSENSUS_HPP

#include <numeric>
#include <limits>
#include <random>
#include <vector>

#include "openMVG/robust_estimation/rand_sampling.hpp"

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
typename Kernel::Model MaxConsensus
(
  const Kernel &kernel,
  const Scorer &scorer,
  std::vector<uint32_t> *best_inliers = nullptr,
  uint32_t max_iteration = 1024
)
{
  const uint32_t min_samples = Kernel::MINIMUM_SAMPLES;
  const uint32_t total_samples = kernel.NumSamples();

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
  std::vector<uint32_t> all_samples(total_samples);
  std::iota(all_samples.begin(), all_samples.end(), 0);

  // Random number generator configuration
  std::mt19937 random_generator(std::mt19937::default_seed);

  std::vector<uint32_t> sample;
  for (uint32_t iteration = 0;  iteration < max_iteration; ++iteration) {
    UniformSample(min_samples, random_generator, &all_samples, &sample);

      std::vector<typename Kernel::Model> models;
      kernel.Fit(sample, &models);

      // Compute costs for each fit.
      for (const auto& model_it : models) {
        std::vector<uint32_t> inliers;
        scorer.Score(kernel, model_it, all_samples, &inliers);

        if (best_num_inliers < inliers.size()) {
          best_num_inliers = inliers.size();
          best_model = model_it;
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
#endif // OPENMVG_ROBUST_ESTIMATION_MAX_CONSENSUS_HPP
