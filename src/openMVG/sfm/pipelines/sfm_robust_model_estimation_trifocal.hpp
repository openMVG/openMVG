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
#include <utility>
#include <vector>

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"
#include "openMVG/multiview/trifocal/solver_trifocal_three_point.hpp"
#include "openMVG/multiview/trifocal/three_view_kernel.hpp"
#include "openMVG/multiview/trifocal/solver_trifocal_metrics.hpp"



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
  const size_t max_iteration_count = 1024
)
{
  constexpr unsigned nviews = 3, npts = 3;

  if (!intrinsics[0] || !intrinsics[1] || !intrinsics[2])
    return false;

  std::array<Mat, nviews> datum;

  for (unsigned v=0; v < nviews; ++v)
    for (unsigned ip=0; ip < npts; ++ip)
      datum[v].col(ip) = (*intrinsics[v])(pxdatum[v].col(ip));
        
  using TrifocalKernel = trifocal::ThreeViewKernel<trifocal::Trifocal3PointPositionTangentialSolver, 
                         trifocal::NormalizedSquaredPointReprojectionOntoOneViewError>;
  
  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2]); // perhaps pass K

  // TODO: we are assuming all images have the same intrinsics
  double threshold_normalized_squared 
    = trifocal::NormalizedSquaredPointReprojectionOntoOneViewError::
    threshold_pixel_to_normalized(4.0, (double (*)[3])(double *)((dynamic_cast<const cameras::Pinhole_Intrinsic *> (intrinsics[0]))->K().data())); // TODO: use ACRANSAC
  
  threshold_normalized_squared *= threshold_normalized_squared;
  relativePoseTrifocal_info.relativePoseTrifocal 
    = MaxConsensus(trifocal_kernel, 
      robust::ScorerEvaluator<TrifocalKernel>(threshold_normalized_squared), 
      &relativePoseTrifocal_info.vec_inliers, max_iteration_count);

  if (relativePoseTrifocal_info.vec_inliers.size() <
      1.5 * TrifocalKernel::Solver::MINIMUM_SAMPLES )
  {
    return false; // no sufficient coverage (the model does not support enough samples)
  }

  // TODO might have to re compute residual tolerance or agument the
  // MaxConsensus parameters to return that number, not just inliers.
  // Perhaps move to ACRansac since its interface already provides that.

  // chirality test is done inside the solve TODO

  // TODO important: reconstruct and reproject all inliers with orientations and
  // check that orientations match either inside ransac or as post filtering of
  // correspondences
  return true;
}

#if 0
TODO ACRANSAC work in progress
/**
 * @brief Estimate the Relative pose between two view from point matches and K matrices
 *  by using a robust essential matrix estimation.
 *
 *  Uses ACRansac
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
bool robustRelativePoseTrifocal2
(
  const cameras::IntrinsicBase *intrinsics[3],
  std::array<Mat, 3> pxdatum,
  RelativePoseTrifocal_Info & relativePoseTrifocal_info,
  const size_t max_iteration_count = 1024
)
{
  constexpr unsigned nviews = 3, npts = 3;
  if (!intrinsics[0] || !intrinsics[1] || !intrinsics[2])
    return false;

  std::array<Mat, nviews> datum;

  for (unsigned v=0; v < nviews; ++v)
    for (unsigned ip=0; ip < npts; ++ip)
      datum[v].col(ip) = (*intrinsics[v])(pxdatum[v].col(ip));
        
  using TrifocalKernel = trifocal::ThreeViewKernel<trifocal::Trifocal3PointPositionTangentialSolver, 
                         trifocal::Trifocal3PointPositionTangentialSolver>;
  
  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2]); // perhaps pass K

  // TODO: we are assuming all images have the same intrinsics
  double threshold_normalized_squared 
    = trifocal::NormalizedSquaredPointReprojectionOntoOneViewError::
    threshold_pixel_to_normalized(4.0, (double (*)[3])(double *)((dynamic_cast<const cameras::Pinhole_Intrinsic *> (intrinsics[0]))->K().data()));
  
  threshold_normalized_squared *= threshold_normalized_squared;
    
  // Robustly estimation of the Model and its precision
  const auto ac_ransac_output = robust::ACRANSAC(
    trifocal_kernel, relativePoseTrifocal_info.vec_inliers,
    max_iteration_count, &relativePoseTrifocal_info.RelativePoseTrifocal,
    relativePoseTrifocal_info.initial_residual_tolerance, false);

  relativePose_info.found_residual_precision = ac_ransac_output.first;

  if (relativePose_info.vec_inliers.size() <
      1.5 * KernelType::Solver::MINIMUM_SAMPLES )
  {
    return false; // no sufficient coverage (the model does not support enough samples)
  }
    

  // TODO might have to re compute residual tolerance or agument the
  // MaxConsensus parameters to return that number, not just inliers.
  // Perhaps move to ACRansac since its interface already provides that.

  // chirality test is done inside the solve TODO

  // TODO important: reconstruct and reproject all inliers with orientations and
  // check that orientations match either inside ransac or as post filtering of
  // correspondences
  return true;
}
#endif
} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_ROBUST_MODEL_ESTIMATION_TRIFOCAL_HPP
