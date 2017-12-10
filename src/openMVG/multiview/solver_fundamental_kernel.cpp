
// Copyright (c) 2010 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/numeric/nullspace.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/numeric/poly.h"

namespace openMVG {
namespace fundamental {
namespace kernel {

void SevenPointSolver::Solve
(
  const Mat2X &x1, const Mat2X &x2, std::vector<Mat3> *F
)
{
  assert(7 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  Vec9 f1, f2;
  if (x1.cols() == 7) {
    // Set up the homogeneous system Af = 0 from the equations x'T*F*x = 0.
    using Mat9 = Eigen::Matrix<double, 9, 9>;
    // In the minimal solution use fixed sized matrix to let Eigen and the
    //  compiler doing the maximum of optimization.
    Mat9 epipolar_constraint = Mat::Zero(9,9);
    EncodeEpipolarEquation(x1.colwise().homogeneous(),
                           x2.colwise().homogeneous(),
                           &epipolar_constraint);
    // Find the two F matrices in the nullspace of epipolar_constraint.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 9, 9>> solver
      (epipolar_constraint.transpose() * epipolar_constraint);
    f1 = solver.eigenvectors().leftCols<2>().col(0);
    f2 = solver.eigenvectors().leftCols<2>().col(1);
  }
  else  {
    // Set up the homogeneous system Af = 0 from the equations x'T*F*x = 0.
    Mat epipolar_constraint(x1.cols(), 9);
    EncodeEpipolarEquation(x1.colwise().homogeneous(),
                           x2.colwise().homogeneous(),
                           &epipolar_constraint);
    // Find the two F matrices in the nullspace of epipolar_constraint.
    Eigen::SelfAdjointEigenSolver<Mat> solver
      (epipolar_constraint.transpose() * epipolar_constraint);
    f1 = solver.eigenvectors().leftCols<2>().col(0);
    f2 = solver.eigenvectors().leftCols<2>().col(1);
  }

  const Mat3 F1 = Map<RMat3>(f1.data());
  const Mat3 F2 = Map<RMat3>(f2.data());

  // Then, use the condition det(F) = 0 to determine F. In other words, solve
  // det(F1 + a*F2) = 0 for a.
  double  a = F1(0, 0), j = F2(0, 0),
          b = F1(0, 1), k = F2(0, 1),
          c = F1(0, 2), l = F2(0, 2),
          d = F1(1, 0), m = F2(1, 0),
          e = F1(1, 1), n = F2(1, 1),
          f = F1(1, 2), o = F2(1, 2),
          g = F1(2, 0), p = F2(2, 0),
          h = F1(2, 1), q = F2(2, 1),
          i = F1(2, 2), r = F2(2, 2);

  // Run fundamental_7point_coeffs.py to get the below coefficients.
  // The coefficients are in ascending powers of alpha, i.e. P[N]*x^N.
  const double P[4] = {
    a*e*i + b*f*g + c*d*h - a*f*h - b*d*i - c*e*g,
    a*e*r + a*i*n + b*f*p + b*g*o + c*d*q + c*h*m + d*h*l + e*i*j + f*g*k -
    a*f*q - a*h*o - b*d*r - b*i*m - c*e*p - c*g*n - d*i*k - e*g*l - f*h*j,
    a*n*r + b*o*p + c*m*q + d*l*q + e*j*r + f*k*p + g*k*o + h*l*m + i*j*n -
    a*o*q - b*m*r - c*n*p - d*k*r - e*l*p - f*j*q - g*l*n - h*j*o - i*k*m,
    j*n*r + k*o*p + l*m*q - j*o*q - k*m*r - l*n*p};

  // Solve for the roots of P[3]*x^3 + P[2]*x^2 + P[1]*x + P[0] = 0.
  double roots[3];
  const int num_roots = SolveCubicPolynomial(P, roots);

  // Build the fundamental matrix for each solution.
  for (int kk = 0; kk < num_roots; ++kk)  {
    F->emplace_back(F1 + roots[kk] * F2);
  }
}

void EightPointSolver::Solve
(
  const Mat2X &x1, const Mat2X &x2, std::vector<Mat3> *Fs
)
{
  assert(8 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  Vec9 f;
  if (x1.cols() == 8) {
    using Mat9 = Eigen::Matrix<double, 9, 9>;
    // In the minimal solution use fixed sized matrix to let Eigen and the
    //  compiler doing the maximum of optimization.
    Mat9 epipolar_constraint = Mat::Zero(9,9);
    EncodeEpipolarEquation(x1.colwise().homogeneous(),
                           x2.colwise().homogeneous(),
                           &epipolar_constraint);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 9, 9>> solver
      (epipolar_constraint.transpose() * epipolar_constraint);
    f = solver.eigenvectors().leftCols<1>();
  }
  else  {
    MatX9 epipolar_constraint(x1.cols(), 9);
    epipolar_constraint.fill(0.0);
    EncodeEpipolarEquation(x1, x2, &epipolar_constraint);
    Eigen::SelfAdjointEigenSolver<MatX9> solver
      (epipolar_constraint.transpose() * epipolar_constraint);
    f = solver.eigenvectors().leftCols<1>();
  }

  Mat3 F = Map<RMat3>(f.data());

  // Force the fundamental property if the A matrix has full rank.
  // HZ 11.1.1 pag.280
  if (x1.cols() > 8) {
    // Force fundamental matrix to have rank 2
    Eigen::JacobiSVD<Mat3> USV(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Vec3 d((Vec3() << USV.singularValues().head<2>(), 0.0).finished());
    F = USV.matrixU() * d.asDiagonal() * USV.matrixV().transpose();
  }
  Fs->push_back(F);
}

double SampsonError::Error
(
  const Mat3 &F, const Vec2 &x, const Vec2 &y
)
{
  // See page 287 equation (11.9) of HZ.
  const Vec3 F_x = F * x.homogeneous();
  const Vec3 Ft_y = F.transpose() * y.homogeneous();
  return Square(y.homogeneous().dot(F_x))
    / (F_x.head<2>().squaredNorm() + Ft_y.head<2>().squaredNorm());
}

double SymmetricEpipolarDistanceError::Error
(
  const Mat3 &F, const Vec2 &x, const Vec2 &y
)
{
  // See page 288 equation (11.10) of HZ.
  const Vec3 F_x = F * x.homogeneous();
  const Vec3 Ft_y = F.transpose() * y.homogeneous();
  return Square(y.homogeneous().dot(F_x)) *
    ( 1.0 / F_x.head<2>().squaredNorm()
      + 1.0 / Ft_y.head<2>().squaredNorm())
    / 4.0;  // The divide by 4 is to make this match the Sampson distance.
}


double EpipolarDistanceError::Error
(
  const Mat3 &F, const Vec2 &x, const Vec2 &y
)
{
  // Transfer error in image 2
  // See page 287 equation (11.9) of HZ.
  const Vec3 F_x = F * x.homogeneous();
  return Square(F_x.dot(y.homogeneous())) /  F_x.head<2>().squaredNorm();
}

}  // namespace kernel
}  // namespace fundamental
}  // namespace openMVG
