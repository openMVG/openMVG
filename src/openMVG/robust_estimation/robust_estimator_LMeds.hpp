// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_LMEDS_HPP
#define OPENMVG_ROBUST_ESTIMATION_LMEDS_HPP

#include <algorithm>
#include <numeric>
#include <limits>
#include <vector>

#include "openMVG/robust_estimation/rand_sampling.hpp"
#include "openMVG/robust_estimation/robust_ransac_tools.hpp"

namespace openMVG {
namespace robust{

/// \brief Variant of RANSAC using Least Median of Squares.
/// \details Instead of using a fixed threshold to distinguish inlier/outlier,
/// find the threshold at 1-\a outlierRatio quantile of residuals and keep the
/// parameters minimizing this threshold. The final threshold
/// returned in \a outlierThreshold is a multiple of this and can be used to
/// filter out outliers.
/// LMedS : Z. Zhang. Determining The Epipolar Geometry And Its Uncertainty. A Review
/// IJCV 1998
template <typename Kernel>
  double LeastMedianOfSquares(const Kernel &kernel,
    typename Kernel::Model * model = nullptr ,
    double* outlierThreshold = nullptr ,
    double outlierRatio=0.5,
    double minProba=0.99)
{
  const size_t min_samples = Kernel::MINIMUM_SAMPLES;
  const size_t total_samples = kernel.NumSamples();

  std::vector<double> residuals(total_samples); // Array for storing residuals
  std::vector<uint32_t> vec_sample(min_samples);

  double dBestMedian = std::numeric_limits<double>::max();

  // Required number of iterations is evaluated from outliers ratio
  const uint32_t N = (min_samples<total_samples)?
    getNumSamples(minProba, outlierRatio, min_samples): 0;

  // Precompute the index [0,n] that will be used for random sampling
  std::vector<uint32_t> all_samples(total_samples);
  std::iota(all_samples.begin(), all_samples.end(), 0);

  //--
  // Random number generation
  std::mt19937 random_generator(std::mt19937::default_seed);

  for (uint32_t i=0; i < N; i++)
  {
    // Get Samples indexes
    UniformSample(min_samples, random_generator, &all_samples, &vec_sample);

    // Estimate parameters: the solutions are stored in a vector
    std::vector<typename Kernel::Model> models;
    kernel.Fit(vec_sample, &models);

    // Now test the solutions on the whole data
    for (const auto& model_it : models)
    {
      //Compute Residuals :
      for (uint32_t l = 0; l < total_samples; ++l)
      {
        residuals[l] = kernel.Error(l, model_it);
      }

      // Compute median
      const auto itMedian = residuals.begin() +
        uint32_t( total_samples*(1.-outlierRatio) );
      std::nth_element(residuals.begin(), itMedian, residuals.end());
      const double median = *itMedian;

      // Store best solution
      if (median < dBestMedian)
      {
        dBestMedian = median;
        if (model) (*model) = model_it;
      }
    }
  }

  // This array of precomputed values corresponds to the inverse
  //  cumulative function for a normal distribution. For more information
  //  consult the litterature (Robust Regression for Outlier Detection,
  //  rouseeuw-leroy). The values are computed for each 5%
  static const double ICDF[21] =
  {
    1.4e16, 15.94723940, 7.957896558, 5.287692054,
    3.947153876, 3.138344200, 2.595242369, 2.203797543,
    1.906939402, 1.672911853, 1.482602218, 1.323775627,
    1.188182950, 1.069988721, 0.9648473415, 0.8693011162,
    0.7803041458, 0.6946704675, 0.6079568319,0.5102134568,
    0.3236002672
  };

  // Evaluate the outlier threshold
  if (outlierThreshold)
  {
    const double sigma = ICDF[uint32_t((1.-outlierRatio)*20.)] *
      (1. + 5. / double(total_samples - min_samples));
    *outlierThreshold = (double)(sigma * sigma * dBestMedian * 4.);
    if (N==0) *outlierThreshold = std::numeric_limits<double>::max();
  }

  return dBestMedian;
}

} // namespace robust
} // namespace openMVG

#endif // OPENMVG_ROBUST_ESTIMATION_LMEDS_HPP
