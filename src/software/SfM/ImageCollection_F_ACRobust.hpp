
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include <limits>


//-- A contrario Functor to filter putative corresponding points
struct GeometricFilter_FMatrix_AC
{
  GeometricFilter_FMatrix_AC(
    double dPrecision = std::numeric_limits<double>::infinity(),
    size_t iteration = 4096)
    : m_dPrecision(dPrecision), m_stIteration(iteration) {};

  /// Robust fitting of the FUNDAMENTAL matrix
  void Fit(
    const Mat & xA,
    const std::pair<size_t, size_t> & imgSizeA,
    const Mat & xB,
    const std::pair<size_t, size_t> & imgSizeB,
    std::vector<size_t> & vec_inliers) const
  {
    using namespace openMVG;
    using namespace openMVG::robust;
    vec_inliers.resize(0);
    // Define the AContrario adapted Fundamental matrix solver
    typedef ACKernelAdaptor<
      openMVG::fundamental::kernel::SevenPointSolver,
      openMVG::fundamental::kernel::SimpleError,
      UnnormalizerT,
      Mat3>
      KernelType;

    KernelType kernel(xA, imgSizeA.first, imgSizeA.second,
                      xB, imgSizeB.first, imgSizeB.second, true);

    // Robustly estimate the Fundamental matrix with A Contrario ransac
    Mat3 F;
    double upper_bound_precision = m_dPrecision;
    std::pair<double,double> ACRansacOut =
      ACRANSAC(kernel, vec_inliers, m_stIteration, &F, upper_bound_precision);
    const double & threshold = ACRansacOut.first;
    const double & NFA = ACRansacOut.second;

    if (vec_inliers.size() < KernelType::MINIMUM_SAMPLES *2.5)  {
      vec_inliers.clear();
    }
  }

  double m_dPrecision;  //upper_bound of the precision
  size_t m_stIteration; //maximal number of iteration used
};

