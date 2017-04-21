// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_ROBUST_MODEL_ESTIMATION_HPP
#define OPENMVG_SFM_SFM_ROBUST_MODEL_ESTIMATION_HPP

#include <limits>
#include <utility>
#include <vector>

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {
namespace sfm {

/**
* @brief Estimate the best possible Rotation/Translation from E.
*  Four are possible, keep the one with most of the point in front.
*
* @param[in] x1 bearing vectors corresponding to image observation in image 1
* @param[in] x2 bearing vectors corresponding to image observation in image 2
* @param[in] E essential matrix
* @param[in] vec_inliers selection of x1, x2 columns that are used
* @param[out] R estimated rotation
* @param[out] t estimated translation
* @param[out] vec_selected_points return the index of vec_inliers that are
*    in front of the cameras
* @param[out] vec_points return the 3D point corresponding to
*    vec_selected_points indexes
* @param[in] positive_depth_solution_ratio Pourcentage ratio threshold used
*    to discard if there is two good solution that have many points in front
*    of the cameras
*/
bool estimate_Rt_fromE
(
  const Mat3X & x1,
  const Mat3X & x2,
  const Mat3 & E,
  const std::vector<uint32_t> & vec_inliers,
  Mat3 * R = nullptr,
  Vec3 * t = nullptr,
  std::vector<uint32_t> * vec_selected_points = nullptr,
  std::vector<Vec3> * vec_points = nullptr,
  const double positive_depth_solution_ratio = 0.7
);

bool estimate_Rt_fromE
(
  const Mat3 & K1, const Mat3 & K2,
  const Mat & x1, const Mat & x2,
  const Mat3 & E, const std::vector<uint32_t> & vec_inliers,
  Mat3 * R, Vec3 * t
);

struct RelativePose_Info
{
  Mat3 essential_matrix;
  geometry::Pose3 relativePose;
  std::vector<uint32_t> vec_inliers;
  double initial_residual_tolerance;
  double found_residual_precision;

  RelativePose_Info()
    :initial_residual_tolerance(std::numeric_limits<double>::infinity()),
    found_residual_precision(std::numeric_limits<double>::infinity())
  {}
};

/**
 * @brief Estimate the Relative pose between two view from point matches and K matrices
 *  by using a robust essential matrix estimation.
 *
 * @param[in] K1 camera 1 intrinsics
 * @param[in] K2 camera 2 intrinsics
 * @param[in] x1 camera 1 image points
 * @param[in] x2 camera 2 image points
 * @param[out] relativePose_info relative pose information
 * @param[in] size_ima1 width, height of image 1
 * @param[in] size_ima2 width, height of image 2
 * @param[in] max iteration count
 */
bool robustRelativePose
(
  const Mat3 & K1, const Mat3 & K2,
  const Mat & x1, const Mat & x2,
  RelativePose_Info & relativePose_info,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count = 4096
);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_ROBUST_MODEL_ESTIMATION_HPP
