// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2022 Pierre MOULON, Ricardo Fabbri, Gabriel Andrade

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

struct RelativePoseTrifocal_Info
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
bool robustRelativePoseTrifocal
(
  const cameras::IntrinsicBase *cam[nviews],
  std::array<Mat, nviews> pxdatum,
  RelativePose_Info & relativePose_info,
  const size_t max_iteration_count = 1024
)
{
  /// XXX get bearing vectors datum
  using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver, 
                         Trifocal3PointPositionTangentialSolver>;
  
  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2]);

  double constexpr threshold_pix = 4;
  vector<uint32_t> inliers;
  const auto model = MaxConsensus(trifocal_kernel, 
      ScorerEvaluator<TrifocalKernel>(threshold_pix), &inliers, max_iteration_count);

  // get back the cameras
}
