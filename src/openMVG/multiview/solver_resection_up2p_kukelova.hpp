// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_RESECTION_UP2P_HPP
#define OPENMVG_MULTIVIEW_RESECTION_UP2P_HPP

#include "openMVG/multiview/two_view_kernel.hpp"
#include "openMVG/multiview/solver_resection_metrics.hpp"

namespace openMVG {
namespace euclidean_resection {

struct UP2PSolver_Kukelova {
  enum { MINIMUM_SAMPLES = 2 };
  enum { MAX_MODELS = 2 };

  // Solve the absolute camera pose problem for the p2p setting (upright camera).
  // [1] Closed-form solutions to the minimal absolute pose problems with known vertical direction
  // Zuzana Kukelova, Martin Bujnak, Tomas Pajdla
  // ACCV 2010
  static void Solve
  (
    const Mat & bearing_vectors,
    const Mat & X, // 3D points
    std::vector<Mat34> *models
  );
};

//-- Usable solver for robust estimation framework
using PoseResectionKernel_UP2P_Kukelova =
  two_view::kernel::Kernel<
    UP2PSolver_Kukelova, // Model estimator
    resection::AngularReprojectionError, // Error metric
    Mat34>;

}  // namespace euclidean_resection
}  // namespace openMVG

#endif // OPENMVG_MULTIVIEW_RESECTION_UP2P_HPP
