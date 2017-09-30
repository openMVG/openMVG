// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_kernel.hpp"

#include "openMVG/multiview/projection.hpp"

#include <cassert>

namespace openMVG {
namespace resection {
namespace kernel {

using namespace std;

void translate
(
  const Mat3X & X,
  const Vec3 & vecTranslation,
  Mat3X * XPoints
)
{
  XPoints->resize(X.rows(), X.cols());
  (*XPoints) = X.colwise() + vecTranslation;
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
  using Vec12 = Eigen::Matrix<double, 12, 1>;

  const size_t n = pt2D.cols();
  for (size_t i = 0; i < n; ++i) {
    size_t row_index = i * 2;
    const Vec3 & X = XPoints.col(i);
    const Vec2 & x = pt2D.col(i);
    A.row(row_index) =
      (Vec12() <<  X.homogeneous(),
                   Vec4::Zero(),
                   - x.x() * X.homogeneous()).finished().normalized();
    row_index = i * 2 + 1;
    A.row(row_index) =
      (Vec12() << Vec4::Zero(),
                  X.homogeneous(),
                  - x.y() * X.homogeneous()).finished().normalized();
  }
}

void SixPointResectionSolver::Solve
(
  const Mat &pt2D,
  const Mat &pt3d,
  vector<Mat34> *Ps,
  bool bcheck
)
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
    Mat A = Mat::Zero(n * 2, 12);
    BuildActionMatrix(A, pt2D, XPoints);
    ratio = NullspaceRatio(A, p);
  }
  if (bcheck) {
    if (ratio > 1e-5) //Assert that at least only one solution if found by SVD
    {
      Mat34 P = Map<Mat>(p.data(),4, 3).transpose();
      P = P * translationMatrix;
      P /= P(2,3);

      Mat3 K, R;
      Vec3 t;
      KRt_From_P(P, &K, &R, &t);

      //Assert point in front of the cam
      size_t cpt = 0;
      for (size_t i = 0; i < n; ++i) {
        cpt += (Depth(R, t, pt3d.col(i)) > 0) ? 1 : 0;
      }
      if (cpt == n) {
        Ps->push_back(P);
      }
    }
  }
  else  {
    Mat34 P = Map<Mat>(p.data(), 4, 3).transpose();
    P = P * translationMatrix;
    P /= P(2,3);
    Ps->push_back(P);
  }
}

double SixPointResectionSolver::Error
(
  const Mat34 & P,
  const Vec2 & pt2D,
  const Vec3 & pt3D
)
{
  return (pt2D - Project(P, pt3D)).norm();
}

}  // namespace kernel
}  // namespace resection
}  // namespace openMVG
