
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

#include "openMVG/multiview/solver_essential_five_point.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"

#include <iostream>

namespace openMVG {

Mat FivePointsNullspaceBasis(const Mat2X &x1, const Mat2X &x2) {
  Eigen::Matrix<double,9, 9> A;
  A.setZero();  // Make A square until Eigen supports rectangular SVD.
  fundamental::kernel::EncodeEpipolarEquation(x1, x2, &A);
  Eigen::JacobiSVD<Eigen::Matrix<double, 9, 9> > svd(A, Eigen::ComputeFullV);
  return svd.matrixV().topRightCorner<9, 4>();
}

Vec o1(const Vec &a, const Vec &b) {
  Vec res = Vec::Zero(20);

  res(coef_xx) = a(coef_x) * b(coef_x);
  res(coef_xy) = a(coef_x) * b(coef_y)
               + a(coef_y) * b(coef_x);
  res(coef_xz) = a(coef_x) * b(coef_z)
               + a(coef_z) * b(coef_x);
  res(coef_yy) = a(coef_y) * b(coef_y);
  res(coef_yz) = a(coef_y) * b(coef_z)
               + a(coef_z) * b(coef_y);
  res(coef_zz) = a(coef_z) * b(coef_z);
  res(coef_x)  = a(coef_x) * b(coef_1)
               + a(coef_1) * b(coef_x);
  res(coef_y)  = a(coef_y) * b(coef_1)
               + a(coef_1) * b(coef_y);
  res(coef_z)  = a(coef_z) * b(coef_1)
               + a(coef_1) * b(coef_z);
  res(coef_1)  = a(coef_1) * b(coef_1);

  return res;
}

Vec o2(const Vec &a, const Vec &b) {
  Vec res(20);

  res(coef_xxx) = a(coef_xx) * b(coef_x);
  res(coef_xxy) = a(coef_xx) * b(coef_y)
                + a(coef_xy) * b(coef_x);
  res(coef_xxz) = a(coef_xx) * b(coef_z)
                + a(coef_xz) * b(coef_x);
  res(coef_xyy) = a(coef_xy) * b(coef_y)
                + a(coef_yy) * b(coef_x);
  res(coef_xyz) = a(coef_xy) * b(coef_z)
                + a(coef_yz) * b(coef_x)
                + a(coef_xz) * b(coef_y);
  res(coef_xzz) = a(coef_xz) * b(coef_z)
                + a(coef_zz) * b(coef_x);
  res(coef_yyy) = a(coef_yy) * b(coef_y);
  res(coef_yyz) = a(coef_yy) * b(coef_z)
                + a(coef_yz) * b(coef_y);
  res(coef_yzz) = a(coef_yz) * b(coef_z)
                + a(coef_zz) * b(coef_y);
  res(coef_zzz) = a(coef_zz) * b(coef_z);
  res(coef_xx)  = a(coef_xx) * b(coef_1)
                + a(coef_x)  * b(coef_x);
  res(coef_xy)  = a(coef_xy) * b(coef_1)
                + a(coef_x)  * b(coef_y)
                + a(coef_y)  * b(coef_x);
  res(coef_xz)  = a(coef_xz) * b(coef_1)
                + a(coef_x)  * b(coef_z)
                + a(coef_z)  * b(coef_x);
  res(coef_yy)  = a(coef_yy) * b(coef_1)
                + a(coef_y)  * b(coef_y);
  res(coef_yz)  = a(coef_yz) * b(coef_1)
                + a(coef_y)  * b(coef_z)
                + a(coef_z)  * b(coef_y);
  res(coef_zz)  = a(coef_zz) * b(coef_1)
                + a(coef_z)  * b(coef_z);
  res(coef_x)   = a(coef_x)  * b(coef_1)
                + a(coef_1)  * b(coef_x);
  res(coef_y)   = a(coef_y)  * b(coef_1)
                + a(coef_1)  * b(coef_y);
  res(coef_z)   = a(coef_z)  * b(coef_1)
                + a(coef_1)  * b(coef_z);
  res(coef_1)   = a(coef_1)  * b(coef_1);

  return res;
}

Mat FivePointsPolynomialConstraints(const Mat &E_basis) {
  // Build the polynomial form of E (equation (8) in Stewenius et al. [1])
  Vec E[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      E[i][j] = Vec::Zero(20);
      E[i][j](coef_x) = E_basis(3 * i + j, 0);
      E[i][j](coef_y) = E_basis(3 * i + j, 1);
      E[i][j](coef_z) = E_basis(3 * i + j, 2);
      E[i][j](coef_1) = E_basis(3 * i + j, 3);
    }
  }

  // The constraint matrix.
  Mat M(10, 20);
  int mrow = 0;

  // Determinant constraint det(E) = 0; equation (19) of Nister [2].
  M.row(mrow++) = o2(o1(E[0][1], E[1][2]) - o1(E[0][2], E[1][1]), E[2][0]) +
                  o2(o1(E[0][2], E[1][0]) - o1(E[0][0], E[1][2]), E[2][1]) +
                  o2(o1(E[0][0], E[1][1]) - o1(E[0][1], E[1][0]), E[2][2]);

  // Cubic singular values constraint.
  // Equation (20).
  Vec EET[3][3];
  for (int i = 0; i < 3; ++i) {    // Since EET is symmetric, we only compute
    for (int j = 0; j < 3; ++j) {  // its upper triangular part.
      if (i <= j) {
        EET[i][j] = o1(E[i][0], E[j][0])
                  + o1(E[i][1], E[j][1])
                  + o1(E[i][2], E[j][2]);
      } else {
        EET[i][j] = EET[j][i];
      }
    }
  }

  // Equation (21).
  Vec (&L)[3][3] = EET;
  const Vec trace  = 0.5 * (EET[0][0] + EET[1][1] + EET[2][2]);
  for (const int i : {0,1,2}) {
    L[i][i] -= trace;
  }

  // Equation (23).
  for (const int i : {0,1,2}) {
    for (const int j : {0,1,2}) {
      Vec LEij = o2(L[i][0], E[0][j])
               + o2(L[i][1], E[1][j])
               + o2(L[i][2], E[2][j]);
      M.row(mrow++) = LEij;
    }
  }

  return M;
}

void FivePointsRelativePose(const Mat2X &x1,
                            const Mat2X &x2,
                            std::vector<Mat3> *Es) {
  // Step 1: Nullspace Extraction.
  const Eigen::Matrix<double, 9, 4> E_basis = FivePointsNullspaceBasis(x1, x2);

  // Step 2: Constraint Expansion.
  const Eigen::Matrix<double, 10, 20> E_constraints = FivePointsPolynomialConstraints(E_basis);

  // Step 3: Gauss-Jordan Elimination (done thanks to a LU decomposition).
  using Mat10 = Eigen::Matrix<double, 10, 10>;
  Eigen::FullPivLU<Mat10> c_lu(E_constraints.block<10, 10>(0, 0));
  const Mat10 M = c_lu.solve(E_constraints.block<10, 10>(0, 10));

  // For next steps we follow the matlab code given in Stewenius et al [1].

  // Build action matrix.

  const Mat10 & B = M.topRightCorner<10,10>();
  Mat10 At = Mat10::Zero(10,10);
  At.block<3, 10>(0, 0) = B.block<3, 10>(0, 0);
  At.row(3) = B.row(4);
  At.row(4) = B.row(5);
  At.row(5) = B.row(7);
  At(6,0) = At(7,1) = At(8,3) = At(9,6) = -1;

  Eigen::EigenSolver<Mat10> eigensolver(At);
  const auto& eigenvectors = eigensolver.eigenvectors();
  const auto& eigenvalues = eigensolver.eigenvalues();

  // Build essential matrices for the real solutions.
  Es->reserve(10);
  for (int s = 0; s < 10; ++s) {
    // Only consider real solutions.
    if (eigenvalues(s).imag() != 0) {
      continue;
    }
    Mat3 E;
    Eigen::Map<Vec9 >(E.data()) =
        E_basis * eigenvectors.col(s).tail<4>().real();
    Es->emplace_back(E.transpose());
  }
}

} // namespace openMVG
