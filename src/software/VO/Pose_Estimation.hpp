// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POSE_ESTIMATION_VO_HPP
#define POSE_ESTIMATION_VO_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/solver_homography_kernel.hpp"

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

namespace openMVG  {
namespace VO  {

using namespace openMVG::robust;

static const size_t ACRANSAC_ITER = 64;

struct ResectionSquaredResidualError {
  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  // Return the squared error
  static double Error(const Mat34 & P, const Vec2 & pt2D, const Vec3 & pt3D){
    return (Project(P, pt3D) - pt2D).squaredNorm();
  }
};


struct Pose_Estimator
{
  bool computeH
  (
    const Mat & x1,
    const Mat & x2,
    const size_t w,
    const size_t h
  )
  {
    //-- Homography robust estimation
    using KernelType =
      ACKernelAdaptor<
        openMVG::homography::kernel::FourPointSolver,
        openMVG::homography::kernel::AsymmetricError,
        UnnormalizerI,
        Mat3>;

    KernelType kernel(
      x1, w, h,
      x2, w, h,
      false); // configure as point to point error model.

    std::vector<uint32_t> vec_inliers;
    Mat3 H;
    ACRANSAC(kernel, vec_inliers, 1024, &H, std::numeric_limits<double>::infinity(),
      true);

    // Check the homography support some point to be considered as valid
    return vec_inliers.size() > KernelType::MINIMUM_SAMPLES *2.5;
  }

  bool computeE
  (
    const Mat3 & K,
    const Mat & x1,
    const Mat & x2,
    const std::pair<size_t, size_t> & size_ima
  )
  {
    // Use the 5 point solver to estimate E
    using SolverType = openMVG::essential::kernel::FivePointKernel;
    // Define the AContrario adaptor
    using KernelType =
      ACKernelAdaptorEssential<
        SolverType,
        openMVG::fundamental::kernel::EpipolarDistanceError,
        Mat3>;

    KernelType kernel(x1, size_ima.first, size_ima.second,
      x2, size_ima.first, size_ima.second, K, K);

    std::vector<uint32_t> vec_inliers;
    Mat3 pE;
    double precision = std::numeric_limits<double>::infinity();
    // Robustly estimation of the Essential matrix and it's precision
    ACRANSAC(kernel, vec_inliers, ACRANSAC_ITER, &pE, precision, false);

    return vec_inliers.size() > 2.5 * SolverType::MINIMUM_SAMPLES;
  }

  bool computeP
  (
    const Mat3 & K,
    const Mat & pt2D,
    const Mat & pt3D
  )
  {
    using KernelType =
      ACKernelAdaptorResection_K<
        openMVG::euclidean_resection::P3PSolver,
        ResectionSquaredResidualError,
        Mat34>;

    KernelType kernel(pt2D, pt3D, K);

    std::vector<uint32_t> vec_inliers;
    Mat34 P;
    const double dPrecision = std::numeric_limits<double>::infinity();

    // Robustly estimation of the Projection matrix and it's precision
    ACRANSAC(kernel, vec_inliers, ACRANSAC_ITER, &P, dPrecision, true);

    // Test if the found model is valid
    return vec_inliers.size() > 2.5 * openMVG::euclidean_resection::P3PSolver::MINIMUM_SAMPLES;
  }
};

} // namespace VO
} // namespace openMVG


#endif // POSE_ESTIMATION_VO_HPP
