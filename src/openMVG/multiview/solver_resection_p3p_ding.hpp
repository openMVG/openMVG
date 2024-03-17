// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2024 Yaqing Ding
// Adapted to openMVG by Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef OPENMVG_MULTIVIEW_RESECTION_P3P_DING_HPP
#define OPENMVG_MULTIVIEW_RESECTION_P3P_DING_HPP

#include "openMVG/multiview/two_view_kernel.hpp"
#include "openMVG/multiview/solver_resection_metrics.hpp"
#include "openMVG/multiview/solver_resection_p3p_nordberg.hpp"

namespace openMVG {
namespace euclidean_resection {

struct P3PSolver_Ding {
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 4 };

  // Solve the absolute camera pose problem.
  // Use "Revisiting the P3P Problem." CVPR23
  // Y Ding, J Yang, V Larsson, C Olsson, K Åström. 
  static void Solve(const Mat &bearing_vectors,
                    const Mat &X, // 3D points
                    std::vector<Mat34> *models);
};

//-- Usable solver for robust estimation framework
using PoseResectionKernel_P3P_Ding =
    two_view::kernel::Kernel<P3PSolver_Ding, // Model estimator
                             resection::AngularReprojectionError, // Error
                                                                  // metric
                             Mat34>;

} // namespace euclidean_resection
} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_RESECTION_P3P_DING_HPP
