// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2022 Pierre MOULON, Ricardo Fabbri, Gabriel Andrade

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -----------------------------------------------------------------------------
// This is to be in openmvg/src/openMVG/sfm/pipelines
// mimmicking sfm_robust_model_estimation.{cpp,hpp} therein
// -----------------------------------------------------------------------------

#ifndef OPENMVG_SFM_SFM_ROBUST_MODEL_ESTIMATION_TRIFOCAL_HPP
#define OPENMVG_SFM_SFM_ROBUST_MODEL_ESTIMATION_TRIFOCAL_HPP

#include <limits>
#include <array>
#include <vector>
#include "openMVG/multiview/multiview_match_constraint.hpp"
#include "openMVG/multiview/trifocal/trifocal_model.hpp"




namespace openMVG { namespace cameras { struct IntrinsicBase; } }

namespace openMVG {
namespace sfm {

  
struct RelativePoseTrifocal_Info
{
  trifocal::trifocal_model_t relativePoseTrifocal;
  std::vector<uint32_t> vec_inliers;
  double initial_residual_tolerance;
  double found_residual_precision;

  RelativePoseTrifocal_Info() :
    initial_residual_tolerance(std::numeric_limits<double>::infinity()),
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
  const cameras::IntrinsicBase *intrinsics[3],
  std::array<Mat, 3> pxdatum,
  RelativePoseTrifocal_Info & relativePoseTrifocal_info,
  double threshold_px,
  const size_t max_iteration_count, /* = 128, */ // this default is overriden in sequential_SfM*.hpp
  MultiviewMatchConstraint multiview_match_constraint
);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_ROBUST_MODEL_ESTIMATION_TRIFOCAL_HPP
