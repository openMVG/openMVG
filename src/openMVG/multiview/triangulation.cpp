// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {

// HZ 12.2 pag.312
void TriangulateDLT
(
  const Mat34 &P1,
  const Vec3 &x1,
  const Mat34 &P2,
  const Vec3 &x2,
  Vec4 *X_homogeneous
)
{
  // Solve:
  // [cross(x0,P0) X = 0]
  // [cross(x1,P1) X = 0]
  Mat4 design;
  design.row(0) = x1[0] * P1.row(2) - x1[2] * P1.row(0);
  design.row(1) = x1[1] * P1.row(2) - x1[2] * P1.row(1);
  design.row(2) = x2[0] * P2.row(2) - x2[2] * P2.row(0);
  design.row(3) = x2[1] * P2.row(2) - x2[2] * P2.row(1);

  Eigen::JacobiSVD<Mat4> svd( design, Eigen::ComputeFullV );
  ( *X_homogeneous ) = svd.matrixV().col( 3 );
}

void TriangulateDLT
(
  const Mat34 &P1,
  const Vec3 &x1,
  const Mat34 &P2,
  const Vec3 &x2,
  Vec3 *X_euclidean
)
{
  Vec4 X_homogeneous;
  TriangulateDLT(P1, x1, P2, x2, &X_homogeneous);
  (*X_euclidean) = X_homogeneous.hnormalized();
}

}  // namespace openMVG
