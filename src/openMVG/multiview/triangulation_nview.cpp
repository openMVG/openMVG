// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/numeric/nullspace.hpp"

#include <limits>

namespace openMVG {

void TriangulateNView
(
  const Mat3X &x,
  const std::vector<Mat34> &Ps,
  Vec4 *X
)
{
  assert(X != nullptr);
  const Mat2X::Index nviews = x.cols();
  assert(static_cast<size_t>(nviews) == Ps.size());

  Mat A = Mat::Zero(3 * nviews, 4 + nviews);
  for (Mat::Index i = 0; i < nviews; ++i)
  {
    A.block<3, 4>(3 * i, 0)     = -Ps[i];
    A.block<3,1> (3 * i, 4 + i) = x.col(i);
  }
  Vec X_and_alphas(4 + nviews);
  Nullspace(A, X_and_alphas);
  *X = X_and_alphas.head(4);
}

bool TriangulateNViewAlgebraic
(
  const Mat3X & points,
  const std::vector<Mat34>& poses,
  Vec4* X
)
{
  assert(poses.size() == points.cols());

  Mat4 AtA = Mat4::Zero();
  for (Mat3X::Index i = 0; i < points.cols(); ++i)
  {
    const Vec3 point_norm = points.col(i).normalized();
    const Mat34 cost =
        poses[i] -
        point_norm * point_norm.transpose() * poses[i];
    AtA += cost.transpose() * cost;
  }

  Eigen::SelfAdjointEigenSolver<Mat4> eigen_solver(AtA);
  *X = eigen_solver.eigenvectors().col(0);
  return eigen_solver.info() == Eigen::Success;
}

}  // namespace openMVG
