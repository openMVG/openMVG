// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017, Romain JANVIER & Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//
//
// Three point orthographic essential matrix estimation
// This solver is based on the following paper:
// Magnus Oskarsson, "Two-View Orthographic Epipolar Geometry: Minimal and Optimal Solvers",
// Journal of Mathematical Imaging and Vision, 2017
//
// Original matlab implementation can be found online
// https://github.com/hamburgerlady/ortho-gem

#ifndef OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_THREE_POINT_HPP
#define OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_THREE_POINT_HPP

#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {

/**
 * @brief Computes the relative pose of two calibrated cameras from 3 correspondences.
 *
 * \param x1 Points in the first image.  One per column.
 * \param x2 Corresponding points in the second image. One per column.
 * \param E  A list of at most 2 candidates orthographic essential matrix solutions.
 */
void ThreePointsRelativePose(const Mat2X &x1, const Mat2X &x2,
                             std::vector<Mat3> *E) {
  Mat3 E_estim;

  const Vec2 xd1 = x1.col(1) - x1.col(0);
  const Vec2 yd1 = x1.col(2) - x1.col(0);
  const Vec2 xd2 = x2.col(1) - x2.col(0);
  const Vec2 yd2 = x2.col(2) - x2.col(0);

  const double denom = xd1.x() * yd1.y() - xd1.y() * yd1.x();

  const double aac = ( xd1.y() * yd2.x() - xd2.x() * yd1.y() ) / denom;
  const double aad = ( xd1.y() * yd2.y() - xd2.y() * yd1.y() ) / denom;
  const double bbc = ( xd2.x() * yd1.x() - xd1.x() * yd2.x() ) / denom;
  const double bbd = ( xd2.y() * yd1.x() - xd1.x() * yd2.y() ) / denom;

  // cache aac * aac
  const double aac_sq = aac * aac;

  const double dd_2 = - aac_sq + aad * aad - bbc * bbc + bbd * bbd;
  const double dd_1c = 2.0 * aac * aad + 2.0 * bbc * bbd;
  const double dd_0 = aac_sq + bbc * bbc - 1.0;

  const double d4_4 = dd_1c * dd_1c + dd_2 * dd_2;
  const double d4_2 = -dd_1c * dd_1c + 2.0 * dd_0 * dd_2;
  const double d4_0 = dd_0 * dd_0;

  const double tmp = sqrt( ( d4_2 * d4_2 - 4.0 * d4_4 * d4_0 ) );
  const double dsol_1 = sqrt( - ( d4_2 + tmp ) / d4_4 / 2.0 );

  const double tmp_csol = - aac_sq + aad * aad - bbc * bbc + bbd * bbd;

  const double csol_1 = - ( tmp_csol * dsol_1 * dsol_1 + aac_sq + bbc * bbc - 1.0 ) /
                  ( 2.0 * aac * aad * dsol_1 + 2.0 * bbc * bbd * dsol_1 );
  const double asol_1 = aac * csol_1 + aad * dsol_1;
  const double bsol_1 = bbc * csol_1 + bbd * dsol_1;
  const double esol_1 = - asol_1 * x1(0,0) - bsol_1 * x1(1,0) - csol_1 * x2(0,0) - dsol_1 * x2(1,0);

  E_estim << 0.0, 0.0, asol_1,
             0.0, 0.0, bsol_1,
             csol_1, dsol_1, esol_1;
  E->push_back(E_estim);

  const double dsol_2 = sqrt( - ( d4_2 - tmp ) / d4_4 / 2.0 );
  const double csol_2 = - ( tmp_csol *  dsol_2 * dsol_2 + aac * aac + bbc * bbc - 1.0 ) /
          ( 2.0 * aac * aad * dsol_2 + 2.0 * bbc * bbd * dsol_2 );
  const double asol_2 = aac * csol_2 + aad * dsol_2;
  const double bsol_2 = bbc * csol_2 + bbd * dsol_2;
  const double esol_2 = - asol_2 * x1(0,0) - bsol_2 * x1(1,0) - csol_2 * x2(0,0) - dsol_2 * x2(1,0);

  E_estim << 0.0, 0.0, asol_2,
             0.0, 0.0, bsol_2,
             csol_2, dsol_2, esol_2;
  E->push_back(E_estim);
}

}
#endif //OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_THREE_POINT_HPP
