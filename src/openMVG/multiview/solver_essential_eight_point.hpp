// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_SPHERICAL_HPP
#define OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_SPHERICAL_HPP

#include <vector>
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {

/**
 * Eight-point algorithm for solving for the essential matrix from bearing
 * vector correspondences.
 * See page 294 in HZ Result 11.1.
 *
 */
struct EightPointRelativePoseSolver {
  enum { MINIMUM_SAMPLES = 8 };
  enum { MAX_MODELS = 1 };
  static void Solve
  (
    const Mat3X & x1,
    const Mat3X & x2,
    std::vector<Mat3> * pvec_E
  );
};

// Return the angular error between [0; PI/2]
struct AngularError
{
  static double Error
  (
    const Mat3 & model,
    const Vec3 & x1,
    const Vec3 & x2
  );
};

} // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_SPHERICAL_HPP
