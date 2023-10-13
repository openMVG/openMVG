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

namespace openMVG { namespace cameras { struct IntrinsicBase; } }

namespace openMVG {
namespace sfm {

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
 * @param[in] intrinsics1 camera 1 intrinsics
 * @param[in] intrinsics2 camera 2 intrinsics
 * @param[in] x1 image points in image 1
 * @param[in] x2 image points in image 2
 * @param[out] relativePose_info relative pose information
 * @param[in] size_ima1 width, height of image 1
 * @param[in] size_ima2 width, height of image 2
 * @param[in] max iteration count
 */
bool robustRelativePose
(
  const cameras::IntrinsicBase * intrinsics1,
  const cameras::IntrinsicBase * intrinsics2,
  const Mat & x1,
  const Mat & x2,
  RelativePose_Info & relativePose_info,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count = 4096
);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_ROBUST_MODEL_ESTIMATION_HPP
