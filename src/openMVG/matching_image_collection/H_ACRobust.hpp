
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/multiview/solver_homography_kernel.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include <limits>

namespace openMVG {
using namespace openMVG::robust;

//-- A contrario Functor to filter putative corresponding points
//--  thanks estimation of the homography matrix.
struct GeometricFilter_HMatrix_AC
{
  GeometricFilter_HMatrix_AC(
    double dPrecision = std::numeric_limits<double>::infinity(),
    size_t iteration = 4096)
    : m_dPrecision(dPrecision), m_stIteration(iteration)  {};

  /// Robust fitting of the HOMOGRAPHY matrix
  void Fit(
    const std::pair<size_t, size_t> pairIndex,
    const Mat & xA,
    const std::pair<size_t, size_t> & imgSizeA,
    const Mat & xB,
    const std::pair<size_t, size_t> & imgSizeB,
    std::vector<size_t> & vec_inliers) const
  {
    vec_inliers.clear();

    // Define the AContrario adapted Homography matrix solver
    typedef ACKernelAdaptor<
      openMVG::homography::kernel::FourPointSolver,
      openMVG::homography::kernel::AsymmetricError,
      UnnormalizerI,
      Mat3>
      KernelType;

    KernelType kernel(
      xA, imgSizeA.first, imgSizeA.second,
      xB, imgSizeB.first, imgSizeB.second,
      false); // configure as point to point error model.
    
    // Robustly estimate the Homography matrix with A Contrario ransac
    Mat3 H;
    double upper_bound_precision = m_dPrecision;
    std::pair<double,double> ACRansacOut =
      ACRANSAC(kernel, vec_inliers, m_stIteration, &H, upper_bound_precision);

    if (vec_inliers.size() < KernelType::MINIMUM_SAMPLES *2.5)  {
      vec_inliers.clear();
    }
  }

  double m_dPrecision;  //upper_bound of the precision
  size_t m_stIteration; //maximal number of used iterations
};

}; // namespace openMVG
