
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

#include "openMVG/multiview/solver_resection_kernel.hpp"

#include <cassert>

namespace openMVG {
namespace resection {
namespace kernel {

using namespace std;

void translate(
  const Mat3X & X,
  const Vec3 & vecTranslation,
  Mat3X * XPoints)
{
  XPoints->resize(X.rows(), X.cols());
  for (Mat3X::Index i=0; i<X.cols(); ++i)
    XPoints->col(i) = X.col(i) + vecTranslation;
}

double NullspaceRatio
(
  const Eigen::Ref<const Mat> & A,
  Eigen::Ref<Vec> nullspace
)
{
  Eigen::JacobiSVD<Mat> svd(A, Eigen::ComputeFullV);
  nullspace = svd.matrixV().col(A.cols()-1);
  return svd.singularValues()(A.cols()-2) / svd.singularValues()(0);
}

/// Setup the Direct Linear Transform.
///  Use template in order to support fixed or dynamic sized matrix.
void BuildActionMatrix
(
  Eigen::Ref<Mat> A,
  const Mat &pt2D,
  const Mat &XPoints
)
{
  const size_t n = pt2D.cols();
  for (size_t i = 0; i < n; ++i) {
    size_t row_index = i * 2;
    const Vec3 & X = XPoints.col(i);
    const Vec2 & x = pt2D.col(i);
    A(row_index,  0) =  X(0);
    A(row_index,  1) =  X(1);
    A(row_index,  2) =  X(2);
    A(row_index,  3) =  1.0;
    A(row_index,  8) = -X(0) * x(0);
    A(row_index,  9) = -X(1) * x(0);
    A(row_index, 10) = -X(2) * x(0);
    A(row_index, 11) = -1.0 * x(0);

    row_index = i * 2 + 1;
    A(row_index,  4) =  X(0);
    A(row_index,  5) =  X(1);
    A(row_index,  6) =  X(2);
    A(row_index,  7) =  1.0;
    A(row_index,  8) = -X(0) * x(1);
    A(row_index,  9) = -X(1) * x(1);
    A(row_index, 10) = -X(2) * x(1);
    A(row_index, 11) = -1.0 * x(1);
  }
  // Normalize each row
  for (size_t i = 0; i < static_cast<size_t>(A.rows()); ++i)
    A.row(i).normalize();
}

void SixPointResectionSolver::Solve(
  const Mat &pt2D,
  const Mat &pt3d,
  vector<Mat34> *Ps,
  bool bcheck)
{
  assert(2 == pt2D.rows());
  assert(3 == pt3d.rows());
  assert(6 <= pt2D.cols());
  assert(pt2D.cols() == pt3d.cols());

  //-- Translate 3D points in order to have X0 = (0,0,0,1).
  Vec3 vecTranslation = - pt3d.col(0);
  Mat4 translationMatrix = Mat4::Identity();
  translationMatrix << 1, 0, 0, vecTranslation(0),
                       0, 1, 0, vecTranslation(1),
                       0, 0, 1, vecTranslation(2),
                       0, 0, 0, 1;
  Mat3X XPoints;
  translate(pt3d, vecTranslation, &XPoints);

  const size_t n = pt2D.cols();

  using Vec12 = Eigen::Matrix<double, 12, 1>;
  Vec12 p;
  double ratio = -1.0;
  if (n==6) {
    // In the case of minimal configuration we use fixed sized matrix to let
    //  Eigen and the compiler doing the maximum of optimization.
    using Mat12 = Eigen::Matrix<double, 12, 12>;
    Mat12 A = Mat12::Zero(12, 12);
    BuildActionMatrix(A, pt2D, XPoints);
    ratio = NullspaceRatio(A, p);
  }
  else  {
    Mat A = Mat::Zero(n*2, 12);
    BuildActionMatrix(A, pt2D, XPoints);
    ratio = NullspaceRatio(A, p);
  }
  if (bcheck) {
    if (ratio > 1e-5) //Assert that at least only one solution if found by SVD
    {
      Mat34 P = Map<Mat>(p.data(),4,3).transpose();
      P = P * translationMatrix;
      P /= P(2,3);

      Mat3 K, R;
      Vec3 t;
      KRt_From_P(P,&K,&R,&t);

      //Assert point in front of the cam
      size_t cpt = 0;
      for (size_t i = 0; i < n; ++i) {
        cpt += (Depth(R, t, pt3d.col(i))>0) ? 1 : 0;
      }
      if (cpt == n) {
        Ps->push_back(P);
      }
    }
  }
  else  {
    Mat34 P = Map<Mat>(p.data(),4,3).transpose();
    P = P * translationMatrix;
    P /= P(2,3);
    Ps->push_back(P);
  }
}

}  // namespace kernel
}  // namespace resection
}  // namespace openMVG
