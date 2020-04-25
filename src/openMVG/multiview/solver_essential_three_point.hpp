// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017, Romain JANVIER & Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_THREE_POINT_HPP
#define OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_THREE_POINT_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {

/**
 * @brief Computes the relative pose of two orthographic cameras from 3 correspondences.
 *
 * \param x1 Points in the first image. One per column.
 * \param x2 Corresponding points in the second image. One per column.
 * \param E  A list of at most 2 candidates orthographic essential matrix solutions.
 */
void ThreePointsRelativePose
(
  const Mat2X &x1,
  const Mat2X &x2,
  std::vector<Mat3> *E
);




} // namespace openMVG
#endif //OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_THREE_POINT_HPP
