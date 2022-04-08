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
  Trifocal3PointPositionTangentialSolver::trifocal_model_t RelativePoseTrifocal;
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
  const cameras::IntrinsicBase *intrinsics[nviews],
  std::array<Mat, nviews> pxdatum,
  RelativePoseTrifocal_Info & relativePoseTrifocal_info,
  const size_t max_iteration_count = 1024
)
{
  constexpr unsigned nviews = 3, npts = 3;

  std::array<Mat, nviews> datum;

  for (unsigned v=0; v < nviews; ++v)
    for (unsigned ip=0; ip < npts; ++ip)
      datum[v].col(ip) = (*intrinsics[v])(pxdatum[v].col(ip));
        
  using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver, 
                         Trifocal3PointPositionTangentialSolver>;
  
  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2]); // perhaps pass K

  double constexpr threshold_pix = 4; // TODO: use ACRANSAC
  relativePoseTrifocal_info.RelativePoseTrifocal 
    = MaxConsensus(trifocal_kernel, 
      ScorerEvaluator<TrifocalKernel>(threshold_pix), 
      &relativePoseTrifocal_info.inliers, max_iteration_count);

  // TODO might have to re compute residual tolerance or agument the
  // MaxConsensus parameters to return that number, not just inliers.
  // Perhaps move to ACRansac since its interface already provides that.

  // chirality test is done inside the solve TODO

  // TODO important: reconstruct and reproject all inliers with orientations and
  // check that orientations match either inside ransac or as post filtering of
  // correspondences
}
