// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINFINITY_COMPUTER_VISION_RESECTION_KERNEL_HPP
#define OPENMVG_LINFINITY_COMPUTER_VISION_RESECTION_KERNEL_HPP

#include <vector>

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/two_view_kernel.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {
namespace lInfinityCV {
namespace kernel {

/**
 * Six-point resection
 * P Matrix estimation (Pose estimation)
 * Rely on L1 Resection algorithm.
 * Work from 6 to N points.
 */
struct l1SixPointResectionSolver {
  enum { MINIMUM_SAMPLES = 6 };
  enum { MAX_MODELS = 1 };
  // Solve the problem of camera pose.
  // First 3d point will be translated in order to have X0 = (0,0,0,1)
  static void Solve(const Mat &pt2D, const Mat &pt3D, std::vector<Mat34> *P);

  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  static double Error(const Mat34 & P, const Vec2 & pt2D, const Vec3 & pt3D)
  {
    return ( Project(P, pt3D) - pt2D ).norm();
  }
};

//-- Usable solver for the l1 6pt Resection Estimation
using l1PoseResectionKernel =
  two_view::kernel::Kernel<
    l1SixPointResectionSolver,
    l1SixPointResectionSolver,
    Mat34>;

}  // namespace kernel
}  // namespace lInfinityCV
}  // namespace openMVG

#endif  // OPENMVG_LINFINITY_COMPUTER_VISION_RESECTION_KERNEL_HPP
