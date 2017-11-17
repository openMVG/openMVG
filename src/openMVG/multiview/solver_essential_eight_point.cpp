// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_essential_eight_point.hpp"
#include "openMVG/numeric/nullspace.hpp"
#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"

namespace openMVG {

void EightPointRelativePoseSolver::Solve
(
  const Mat3X & x1,
  const Mat3X & x2,
  std::vector<Mat3> * pvec_E
)
{
  assert(8 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  MatX9 epipolar_constraint(x1.cols(), 9);
  epipolar_constraint.fill(0.0);
  fundamental::kernel::EncodeEpipolarEquation(x1, x2, &epipolar_constraint);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 9, 9>> solver
    (epipolar_constraint.transpose() * epipolar_constraint);
  const Vec9 e = solver.eigenvectors().leftCols<1>();
  Mat3 E = Map<const RMat3>(e.data());

  // Find the closest essential matrix to E in frobenius norm
  // E = UD'VT
  if (x1.cols() > 8) {
    Eigen::JacobiSVD<Mat3> USV(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Vec3 d = USV.singularValues();
    const double a = d[0];
    const double b = d[1];
    d << (a + b) / 2., (a + b) / 2., 0.0;
    E = USV.matrixU() * d.asDiagonal() * USV.matrixV().transpose();
  }
  pvec_E->emplace_back(E);
}


// Return the angular error between [0; PI/2]
double AngularError::Error
(
  const Mat3 & model,
  const Vec3 & x1,
  const Vec3 & x2
)
{
  const Vec3 Em1 = (model * x1).normalized();
  const double angleVal = (x2.transpose() * Em1);
  return std::abs(std::asin(angleVal));
}

} // namespace openMVG
