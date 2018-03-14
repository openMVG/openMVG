// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_RESECTION_P3P_KE_HPP
#define OPENMVG_MULTIVIEW_RESECTION_P3P_KE_HPP

#include "openMVG/multiview/two_view_kernel.hpp"

namespace openMVG {
namespace euclidean_resection {

struct P3PSolver_Ke {
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 4 };

  // Solve the absolute camera pose problem.
  // Use "An Efficient Algebraic Solution to the Perspective-Three-Point Problem".
  // Tong Ke, Stergios Roumeliotis. CVPR 2017
  static void Solve
  (
    const Mat & bearing_vectors,
    const Mat & X, // 3D points
    std::vector<Mat34> *models
  );

  // Compute the angular residual between the bearing vector and the 3d point projection vector
  static double Error
  (
    const Mat34 & P,
    const Vec3 & bearing_vector,
    const Vec3 & pt3D
  );
};

//-- Usable solver for robust estimation framework
using PoseResectionKernel_P3P_Ke =
  two_view::kernel::Kernel<
    P3PSolver_Ke, // Model estimator
    P3PSolver_Ke, // Error metric
    Mat34>;

}  // namespace euclidean_resection
}  // namespace openMVG

#endif // OPENMVG_MULTIVIEW_RESECTION_P3P_KE_HPP
