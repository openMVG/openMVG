// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"

#include <array>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/multiview/motion_from_essential.hpp"
#include "openMVG/multiview/solver_essential_eight_point.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

using namespace openMVG::cameras;
using namespace openMVG::geometry;

namespace openMVG {
namespace sfm {

bool robustRelativePose
(
  const IntrinsicBase * intrinsics1,
  const IntrinsicBase * intrinsics2,
  const Mat & x1,
  const Mat & x2,
  RelativePose_Info & relativePose_info,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  const size_t max_iteration_count
)
{
  if (!intrinsics1 || !intrinsics2)
    return false;

  // Compute the bearing vectors
  const Mat3X
    bearing1 = (*intrinsics1)(x1),
    bearing2 = (*intrinsics2)(x2);

  if (isPinhole(intrinsics1->getType())
      && isPinhole(intrinsics2->getType()))
  {
    // Define the AContrario adaptor to use the 5 point essential matrix solver.
    using KernelType = robust::ACKernelAdaptorEssential<
      openMVG::essential::kernel::FivePointSolver,
      openMVG::fundamental::kernel::EpipolarDistanceError,
      Mat3>;
    KernelType kernel(x1, bearing1, size_ima1.first, size_ima1.second,
                      x2, bearing2, size_ima2.first, size_ima2.second,
                      dynamic_cast<const cameras::Pinhole_Intrinsic*>(intrinsics1)->K(),
                      dynamic_cast<const cameras::Pinhole_Intrinsic*>(intrinsics2)->K());

    // Robustly estimation of the Model and its precision
    const auto ac_ransac_output = robust::ACRANSAC(
      kernel, relativePose_info.vec_inliers,
      max_iteration_count, &relativePose_info.essential_matrix,
      relativePose_info.initial_residual_tolerance, false);

    relativePose_info.found_residual_precision = ac_ransac_output.first;

    if (relativePose_info.vec_inliers.size() <
        2.5 * KernelType::Solver::MINIMUM_SAMPLES )
    {
      return false; // no sufficient coverage (the model does not support enough samples)
    }
  }
  else
  {
    // Define the AContrario adaptor to use the 8 point essential matrix solver.
    typedef openMVG::robust::ACKernelAdaptor_AngularRadianError<
        openMVG::EightPointRelativePoseSolver,
        // openMVG::essential::kernel::FivePointSolver,
        openMVG::AngularError,
        Mat3>
        KernelType;

    KernelType kernel(bearing1, bearing2);

    // Robustly estimate the Essential matrix with A Contrario ransac
    const double upper_bound_precision =
      (relativePose_info.initial_residual_tolerance == std::numeric_limits<double>::infinity()) ?
        std::numeric_limits<double>::infinity()
        : D2R(relativePose_info.initial_residual_tolerance);
    const auto ac_ransac_output =
      ACRANSAC(kernel, relativePose_info.vec_inliers,
        max_iteration_count, &relativePose_info.essential_matrix,
        upper_bound_precision, false);

    const double & threshold = ac_ransac_output.first;
    relativePose_info.found_residual_precision = R2D(threshold); // Degree

    if (relativePose_info.vec_inliers.size() <
        2.5 * KernelType::Solver::MINIMUM_SAMPLES )
    {
      return false; // no sufficient coverage (the model does not support enough samples)
    }
  }

  // estimation of the relative poses based on the cheirality test
  Pose3 relative_pose;
  if (!RelativePoseFromEssential(
    bearing1,
    bearing2,
    relativePose_info.essential_matrix,
    relativePose_info.vec_inliers, &relative_pose))
  {
    return false;
  }
  relativePose_info.relativePose = relative_pose;
  return true;
}

} // namespace sfm
} // namespace openMVG
