// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020, Viktor Larsson
// Copyright (c) 2020, Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_up2p_kukelova.hpp"
#include "openMVG/numeric/numeric.h"

using namespace openMVG;

namespace openMVG {
namespace euclidean_resection {

// Port from PoseLib
// https://github.com/vlarsson/PoseLib/blob/044c55022ff078f46b65195eda3f467aa85dad49/PoseLib/univariate.cc#L47
/* Solves the quadratic equation a*x^2 + b*x + c = 0 */
int solve_quadratic_real(double a, double b, double c, double roots[2]) {

  const double b2m4ac = b * b - 4 * a * c;
  if (b2m4ac < 0)
    return 0;

  const double sq = std::sqrt(b2m4ac);

  // Choose sign to avoid cancellations
  roots[0] = (b > 0) ? (2 * c) / (-b - sq) : (2 * c) / (-b + sq);
  roots[1] = c / (a * roots[0]);

  return 2;
}

// Port of up2p from PoseLib to OpenMVG
// PoseLib, a library that provides a collection of minimal solvers for
//  camera pose estimation.
// See original up2p here https://github.com/vlarsson/PoseLib/blob/master/PoseLib/up2p.cc
void UP2PSolver_Kukelova::Solve
(
  const Mat & x,
  const Mat & X,
  std::vector<Mat34> * models
)
{
  Eigen::Matrix<double, 4, 4> A;
  Eigen::Matrix<double, 4, 2> b;

  A << -x.col(0)(2), 0, x.col(0)(0), X.col(0)(0) * x.col(0)(2) - X.col(0)(2) * x.col(0)(0), 0, -x.col(0)(2), x.col(0)(1), -X.col(0)(1) * x.col(0)(2) - X.col(0)(2) * x.col(0)(1), -x.col(1)(2), 0, x.col(1)(0), X.col(1)(0) * x.col(1)(2) - X.col(1)(2) * x.col(1)(0), 0, -x.col(1)(2), x.col(1)(1), -X.col(1)(1) * x.col(1)(2) - X.col(1)(2) * x.col(1)(1);
  b << -2 * X.col(0)(0) * x.col(0)(0) - 2 * X.col(0)(2) * x.col(0)(2), X.col(0)(2) * x.col(0)(0) - X.col(0)(0) * x.col(0)(2), -2 * X.col(0)(0) * x.col(0)(1), X.col(0)(2) * x.col(0)(1) - X.col(0)(1) * x.col(0)(2), -2 * X.col(1)(0) * x.col(1)(0) - 2 * X.col(1)(2) * x.col(1)(2), X.col(1)(2) * x.col(1)(0) - X.col(1)(0) * x.col(1)(2), -2 * X.col(1)(0) * x.col(1)(1), X.col(1)(2) * x.col(1)(1) - X.col(1)(1) * x.col(1)(2);

  b = A.inverse() * b;

  const double c2 = b(3, 0);
  const double c3 = b(3, 1);

  double qq[2];
  const int sols = solve_quadratic_real(1.0, c2, c3, qq);

  models->clear();
  for (int i = 0; i < sols; ++i) {
    const double q = qq[i];
    const double q2 = q * q;
    const double inv_norm = 1.0 / (1 + q2);
    const double cq = (1 - q2) * inv_norm;
    const double sq = 2 * q * inv_norm;

    Mat3 rotation = Mat3::Identity();
    rotation(0, 0) = cq;
    rotation(0, 2) = sq;
    rotation(2, 0) = -sq;
    rotation(2, 2) = cq;

    Vec3 translation;
    translation = b.block<3, 1>(0, 0) * q + b.block<3, 1>(0, 1);
    translation *= -inv_norm;

    const Mat34 pose = HStack(rotation, translation);
    models->push_back(pose);
  }
}

} // namespace euclidean_resection
} // namespace openMVG
