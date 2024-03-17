// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Michael Persson
// Adapted to openMVG by Romain Janvier and Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef OPENMVG_MULTIVIEW_RESECTION_P3P_NORDBERG_HPP
#define OPENMVG_MULTIVIEW_RESECTION_P3P_NORDBERG_HPP

#include "openMVG/multiview/two_view_kernel.hpp"
#include "openMVG/multiview/solver_resection_metrics.hpp"

namespace openMVG {
namespace euclidean_resection {

struct P3PSolver_Nordberg {
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 4 };

  static inline bool root2real(const double b, const double c, double &r1,
                               double &r2) {
    const double v = b * b - 4.0 * c;
    if (v <= 0.0) {
      r1 = r2 = -0.5 * b;
      return v >= 0.0;
    }
    const double y = std::sqrt(v);
    if (b < 0.0) {
      r1 = 0.5 * (-b + y);
      r2 = 0.5 * (-b - y);
    } else {
      r1 = 2.0 * c / (-b + y);
      r2 = 2.0 * c / (-b - y);
    }
    return true;
  };

  // Solve the absolute camera pose problem.
  // Use "Lambda Twist: An Accurate Fast Robust
  // Perspective Three Point (P3P) Solver
  // Persson, M.; Nordberg, K.
  static void Solve(const Mat &bearing_vectors,
                    const Mat &X, // 3D points
                    std::vector<Mat34> *models);
};

//-- Usable solver for robust estimation framework
using PoseResectionKernel_P3P_Nordberg =
    two_view::kernel::Kernel<P3PSolver_Nordberg, // Model estimator
                             resection::AngularReprojectionError, // Error metric
                             Mat34>;

} // namespace euclidean_resection
} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_RESECTION_P3P_NORDBERG_HPP
