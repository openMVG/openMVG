// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017, Romain JANVIER & Pierre MOULON
// Copyright (c) 2020 Pierre MOULON.

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

namespace essential {
namespace kernel {

/**
 * Three-point algorithm for solving for the essential matrix from bearing
 * vector correspondences assuming upright images.
 * Implementation of [1] section 3.3. Linear 3-point Algorithm
 * Note: this is an approximate solver, not a minimal solver
 *
 * [1] "Fast and Reliable Minimal Relative Pose Estimation under Planar Motion"
 * Sunglok Choi, Jong-Hwan Kim, 2018
 *
 * [2] Street View Goes Indoors: Automatic Pose Estimation From Uncalibrated Unordered Spherical Panoramas.
 * Mohamed Aly and Jean-Yves Bouguet.
 * IEEE Workshop on Applications of Computer Vision (WACV), Colorado, January 2012.
 *
 * Comment [2] and [1] propose both a Direct Linear Method using 3 correspondences.
 * Note that they are using gravity axis and that [1] provides more details about the fact that
 * the linear formulation is not a minimal solver since it cannot enforce the
 * "Pythagorean" identity on sin^2(t) + cos^2(t) = 1
 *
 */
struct ThreePointUprightRelativePoseSolver {
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 1 };
  static void Solve
  (
    const Mat3X & bearing_a,
    const Mat3X & bearing_b,
    std::vector<Mat3> * pvec_E
  );
};

} // namespace kernel
} // namespace essential

} // namespace openMVG
#endif //OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_THREE_POINT_HPP
