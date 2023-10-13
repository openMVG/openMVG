// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// -----------------------------------------------------------------------------
// This is to be in openmvg/src/openMVG/sfm/pipelines
// mimmicking sfm_robust_model_estimation.{cpp,hpp} therein
// -----------------------------------------------------------------------------


#include <utility>
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"
#include "openMVG/multiview/trifocal/solver_trifocal_three_point.hpp"
#include "openMVG/multiview/trifocal/three_view_kernel.hpp"
#include "openMVG/multiview/trifocal/solver_trifocal_metrics.hpp"
#include "openMVG/sfm/pipelines/sequential/sfm_robust_model_estimation_trifocal.hpp"


using namespace openMVG::cameras;
using namespace openMVG::geometry;

namespace openMVG {
namespace sfm {


static void
invert_intrinsics_tgt(
    const Mat3 &K,
    const double px_tgt_coords[2],
    double normalized_tgt_coords[2])
{
  const double *tp = px_tgt_coords;
  double *t = normalized_tgt_coords;
  t[1] = tp[1]/K(1,1);
  t[0] = (tp[0] - K(0,1)*t[1])/K(0,0);
  double n = hypot(t[0], t[1]);
  t[0] /= n; t[1] /= n;
}

bool robustRelativePoseTrifocal
(
  const cameras::IntrinsicBase *intrinsics[3],
  std::array<Mat, 3> pxdatum,
  RelativePoseTrifocal_Info & relativePoseTrifocal_info,
  double threshold_px,
  const size_t max_iteration_count
)
{
  OPENMVG_LOG_INFO << "npts = " << pxdatum[0].cols();
  constexpr unsigned nviews = 3;
  long int npts = pxdatum[0].cols();

  if (!intrinsics[0] || !intrinsics[1] || !intrinsics[2])
    return false;

  std::array<Mat, nviews> datum;
  for (unsigned v = 0; v < nviews; ++v) {
    datum[v].resize(4,npts);
    for (unsigned ip = 0; ip < npts; ++ip) { // workaround for inverting intrisics based on given structure
                                           // Get 3D cam coords from pxdatum ->
                                           // get eigen matrix 3x1
                                           // then convert into eigen vector and normalize it
      //OPENMVG_LOG_INFO << "datum point in pixels:" << pxdatum[v].col(ip).head(2);
      //OPENMVG_LOG_INFO << "datum tangent in pixels:" << pxdatum[v].col(ip).tail(2);
      datum[v].col(ip).head(2) = (*intrinsics[v])(pxdatum[v].col(ip).head<2>()).colwise().hnormalized();
      const cameras::Pinhole_Intrinsic *Kin = dynamic_cast<const cameras::Pinhole_Intrinsic *>(intrinsics[v]);
      assert(Kin);
      invert_intrinsics_tgt(Kin->K(), pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
      datum[v].col(ip).tail(2) = datum[v].col(ip).tail(2).normalized();
      //OPENMVG_LOG_INFO << "datum point in units:" << datum[v].col(ip).head(2);
      //OPENMVG_LOG_INFO << "datum tangent in units:" << datum[v].col(ip).tail(2);
    }
  }
  using TrifocalKernel = trifocal::ThreeViewKernel<trifocal::Trifocal3PointPositionTangentialSolver, 
                         trifocal::NormalizedSquaredPointReprojectionOntoOneViewError>;
  
  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2]); // perhaps pass K

  OPENMVG_LOG_INFO << "Initialized kernel. Calling relativePoseTrifocal";
  // TODO: we are assuming all images have the same intrinsics
  double threshold_normalized_squared 
    = trifocal::NormalizedSquaredPointReprojectionOntoOneViewError::
    threshold_pixel_to_normalized(threshold_px, (double (*)[3])(double *)((dynamic_cast<const cameras::Pinhole_Intrinsic *> (intrinsics[0]))->K().data()));
  threshold_normalized_squared *= threshold_normalized_squared;
  OPENMVG_LOG_INFO << "RANSAC threshold is " << threshold_normalized_squared;
  relativePoseTrifocal_info.found_residual_precision = threshold_px; // XXX TODO: // improve

  relativePoseTrifocal_info.relativePoseTrifocal 
    = MaxConsensus(trifocal_kernel, 
      robust::ScorerEvaluator<TrifocalKernel>(threshold_normalized_squared), 
      &relativePoseTrifocal_info.vec_inliers, max_iteration_count);

  OPENMVG_LOG_INFO << "Number of inliers " << relativePoseTrifocal_info.vec_inliers.size();

  // for Debug
  OPENMVG_LOG_INFO << "Best inlier residual: ";
  for (unsigned i=0; i < relativePoseTrifocal_info.vec_inliers.size(); ++i) {
    OPENMVG_LOG_INFO << "\tInlier " << i <<  " " << 
      trifocal::NormalizedSquaredPointReprojectionOntoOneViewError::Error( 
        relativePoseTrifocal_info.relativePoseTrifocal, 
        datum[0].col(relativePoseTrifocal_info.vec_inliers[i]), 
        datum[1].col(relativePoseTrifocal_info.vec_inliers[i]), 
        datum[2].col(relativePoseTrifocal_info.vec_inliers[i]));
  }

  if (relativePoseTrifocal_info.vec_inliers.size() <=
      TrifocalKernel::Solver::MINIMUM_SAMPLES)
  {
    OPENMVG_LOG_INFO << "No sufficient coverage";
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
