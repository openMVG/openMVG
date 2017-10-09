// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_RESECTION_KERNEL_HPP
#define OPENMVG_MULTIVIEW_RESECTION_KERNEL_HPP

#include <vector>

#include "openMVG/multiview/two_view_kernel.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {
namespace resection {
namespace kernel {

/**
 * Six-point resection
 * P Matrix estimation (Pose estimation)
 * Compute a projection matrix using linear least squares.
 * Rely on Linear Resection algorithm.
 * Work from 6 to N points.
 */
struct SixPointResectionSolver {
  enum { MINIMUM_SAMPLES = 6 };
  enum { MAX_MODELS = 1 };
  // Solve the problem of camera pose.
  // First 3d point will be translated in order to have X0 = (0,0,0,1)
  static void Solve
  (
    const Mat &pt2D,
    const Mat &pt3D,
    std::vector<Mat34> *P,
    bool bcheck = true
  );

  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  static double Error
  (
    const Mat34 & P,
    const Vec2 & pt2D,
    const Vec3 & pt3D
  );
};

//-- Usable solver for the 6pt Resection estimation
using PoseResectionKernel =
  two_view::kernel::Kernel<
    SixPointResectionSolver, // Model estimator
    SixPointResectionSolver, // Error metric
    Mat34>;

}  // namespace kernel
}  // namespace resection
}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_RESECTION_KERNEL_HPP
