
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
  Eigen::JacobiSVD<Eigen::Matrix<double,9, 9> > svd(A,Eigen::ComputeFullV);
  return svd.matrixV().topRightCorner<9,4>();
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
  Vec trace  = 0.5 * (EET[0][0] + EET[1][1] + EET[2][2]);
  for (int i = 0; i < 3; ++i) {
    L[i][i] -= trace;
  }

  // Equation (23).
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      Vec LEij = o2(L[i][0], E[0][j])
               + o2(L[i][1], E[1][j])
               + o2(L[i][2], E[2][j]);
      M.row(mrow++) = LEij;
    }
  }

  return M;
}

void FivePointsGaussJordan(Mat *Mp) {
  Mat &M = *Mp;

  // Gauss Elimination.
  for (int i = 0; i < 10; ++i) {
    M.row(i) /= M(i,i);
    for (int j = i + 1; j < 10; ++j) {
      M.row(j) = M.row(j) / M(j,i) - M.row(i);
    }
  }
  // Backsubstitution.
  for (int i = 9; i >= 0; --i) {
    for (int j = 0; j < i; ++j) {
      M.row(j) = M.row(j) - M(j,i) * M.row(i);
    }
  }
}

void FivePointsRelativePose(const Mat2X &x1,
                            const Mat2X &x2,
                            vector<Mat3> *Es) {
  // Step 1: Nullspace Extraction.
  Mat E_basis = FivePointsNullspaceBasis(x1, x2);

  // Step 2: Constraint Expansion.
  Mat M = FivePointsPolynomialConstraints(E_basis);

  // Step 3: Gauss-Jordan Elimination.
  FivePointsGaussJordan(&M);

  // For next steps we follow the matlab code given in Stewenius et al [1].

  // Build action matrix.
  Mat B = M.topRightCorner<10,10>();
  Mat At = Mat::Zero(10,10);
  At.row(0) = -B.row(0);
  At.row(1) = -B.row(1);
  At.row(2) = -B.row(2);
  At.row(3) = -B.row(4);
  At.row(4) = -B.row(5);
  At.row(5) = -B.row(7);
  At(6,0) = 1;
  At(7,1) = 1;
  At(8,3) = 1;
  At(9,6) = 1;

  // Compute solutions from action matrix's eigenvectors.
  Eigen::EigenSolver<Mat> es(At);
  typedef Eigen::EigenSolver<Mat>::EigenvectorsType Matc;
  Matc V = es.eigenvectors();
  Matc SOLS(4, 10);
  SOLS.row(0) = V.row(6).array() / V.row(9).array();
  SOLS.row(1) = V.row(7).array() / V.row(9).array();
  SOLS.row(2) = V.row(8).array() / V.row(9).array();
  SOLS.row(3).setOnes();

  // Get the ten candidate E matrices in vector form.
  Matc Evec = E_basis * SOLS;

  // Build essential matrices for the real solutions.
  Es->reserve(10);
  for (int s = 0; s < 10; ++s) {
    Evec.col(s) /= Evec.col(s).norm();
    bool is_real = true;
    for (int i = 0; i < 9; ++i) {
      if (Evec(i,s).imag() != 0) {
        is_real = false;
        break;
      }
    }
    if (is_real) {
      Mat3 E;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          E(i, j) = Evec(3 * i + j, s).real();
        }
      }
      Es->push_back(E);
    }
  }
}

} // namespace openMVG

