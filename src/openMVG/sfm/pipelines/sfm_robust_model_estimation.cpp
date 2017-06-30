// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"

#include <array>

#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

namespace openMVG {
namespace sfm {

bool estimate_Rt_fromE
(
  const Mat3X & x1,
  const Mat3X & x2,
  const Mat3 & E,
  const std::vector<uint32_t> & vec_inliers,
  Mat3 * R,
  Vec3 * t,
  std::vector<uint32_t> * vec_selected_points,
  std::vector<Vec3> * vec_points,
  const double positive_depth_solution_ratio
)
{
  // Accumulator to find the best solution
  std::array<uint32_t, 4> f{0, 0, 0, 0};

  // Recover plausible rotation and translation from E.
  std::vector<Mat3> Rs;  // Rotation matrix.
  std::vector<Vec3> ts;  // Translation matrix.
  MotionFromEssential(E, &Rs, &ts);

  // Find which solution is the best:
  // - count how many triangulated observations are in front of the cameras
  std::vector< std::vector<uint32_t> > vec_newInliers(Rs.size());
  std::vector< std::vector<Vec3> > vec_3D(Rs.size());

  const Mat34 P1 = HStack(Mat3::Identity(), Vec3::Zero());

  for (unsigned int i = 0; i < Rs.size(); ++i)
  {
    const Mat34 P2 = HStack(Rs[i], ts[i]);
    Vec3 X;

    for (const uint32_t & inlier_it : vec_inliers)
    {
      const Vec3
        & x1_ = x1.col(inlier_it),
        & x2_ = x2.col(inlier_it);
      TriangulateDLT(P1, x1_, P2, x2_, &X);
      // Test if X is in front of the two cameras
      const Vec3 Mc = Rs[i] * X + ts[i];
      if (x2_.dot(Mc) > 0 && x1_.dot(X)  > 0)
      {
        ++f[i];
        vec_newInliers[i].push_back(inlier_it);
        vec_3D[i].push_back(X);
      }
    }
  }

  // Check if there is a valid solution:
  const auto iter = std::max_element(f.begin(), f.end());
  if (*iter == 0)
  {
    // There is no right solution with points in front of the cameras
    return false;
  }

  // Export the best solution data
  const size_t index = std::distance(f.begin(), iter);
  if (R)
    (*R) = Rs[index];
  if (t)
    (*t) = ts[index];
  if (vec_selected_points)
    (*vec_selected_points) = vec_newInliers[index];
  if (vec_points)
    (*vec_points) = vec_3D[index];

  // Test if the best solution is good by using the ratio of the two best solution score
  std::array<uint32_t,4> f_sorted (f);
  std::sort(f_sorted.begin(), f_sorted.end());
  const double ratio = f_sorted[2] / static_cast<double>(f_sorted[3]);
  return (ratio < positive_depth_solution_ratio);
}

bool robustRelativePose
(
  const Mat3 & K1, const Mat3 & K2,
  const Mat & x1, const Mat & x2,
  RelativePose_Info & relativePose_info,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count
)
{
  // Define the AContrario adaptor
  using KernelType = robust::ACKernelAdaptorEssential<
    openMVG::essential::kernel::FivePointKernel,
    openMVG::fundamental::kernel::EpipolarDistanceError,
    Mat3>;

  KernelType kernel(x1, size_ima1.first, size_ima1.second,
                    x2, size_ima2.first, size_ima2.second,
                    K1, K2);

  // Robustly estimation of the Essential matrix and it's precision
  const std::pair<double,double> acRansacOut = robust::ACRANSAC(
    kernel, relativePose_info.vec_inliers,
    max_iteration_count, &relativePose_info.essential_matrix,
    relativePose_info.initial_residual_tolerance, false);

  relativePose_info.found_residual_precision = acRansacOut.first;

  if (relativePose_info.vec_inliers.size() <
      2.5 * openMVG::essential::kernel::FivePointKernel::MINIMUM_SAMPLES )
  {
    return false; // no sufficient coverage (the model does not support enough samples)
  }

  // estimation of the relative poses
  Mat3 R;
  Vec3 t;
  if (!estimate_Rt_fromE(
    K1.inverse() * x1.colwise().homogeneous(),
    K2.inverse() * x2.colwise().homogeneous(),
    relativePose_info.essential_matrix,
    relativePose_info.vec_inliers, &R, &t))
  {
    return false; // cannot find a valid [R|t] couple that makes the inliers in front of the camera.
  }

  // Store [R|C] for the second camera, since the first camera is [Id|0]
  relativePose_info.relativePose = geometry::Pose3(R, -R.transpose() * t);
  return true;
}

} // namespace sfm
} // namespace openMVG
