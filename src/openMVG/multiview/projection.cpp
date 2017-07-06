
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

#include "openMVG/multiview/projection.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {

/// Compute P = K[R|t]
void P_From_KRt(
  const Mat3 &K,  const Mat3 &R,  const Vec3 &t, Mat34 *P) {
  *P = K * HStack(R,t);
}

void KRt_From_P(const Mat34 &P, Mat3 *Kp, Mat3 *Rp, Vec3 *tp) {
  // Decompose using the RQ decomposition HZ A4.1.1 pag.579.
  Mat3 K = P.block(0, 0, 3, 3);

  Mat3 Q;
  Q.setIdentity();

  // Set K(2,1) to zero.
  if (K(2,1) != 0) {
    double c = -K(2,2);
    double s = K(2,1);
    double l = sqrt(c * c + s * s);
    c /= l; s /= l;
    Mat3 Qx;
    Qx << 1, 0, 0,
          0, c, -s,
          0, s, c;
    K = K * Qx;
    Q = Qx.transpose() * Q;
  }
  // Set K(2,0) to zero.
  if (K(2,0) != 0) {
    double c = K(2,2);
    double s = K(2,0);
    double l = sqrt(c * c + s * s);
    c /= l; s /= l;
    Mat3 Qy;
    Qy << c, 0, s,
          0, 1, 0,
         -s, 0, c;
    K = K * Qy;
    Q = Qy.transpose() * Q;
  }
  // Set K(1,0) to zero.
  if (K(1,0) != 0) {
    double c = -K(1,1);
    double s = K(1,0);
    double l = sqrt(c * c + s * s);
    c /= l; s /= l;
    Mat3 Qz;
    Qz << c,-s, 0,
          s, c, 0,
          0, 0, 1;
    K = K * Qz;
    Q = Qz.transpose() * Q;
  }

  Mat3 R = Q;

  //Mat3 H = P.block(0, 0, 3, 3);
  // RQ decomposition
  //Eigen::HouseholderQR<Mat3> qr(H);
  //Mat3 K = qr.matrixQR().triangularView<Eigen::Upper>();
  //Mat3 R = qr.householderQ();

  // Ensure that the diagonal is positive and R determinant == 1.
  if (K(2,2) < 0) {
    K = -K;
    R = -R;
  }
  if (K(1,1) < 0) {
    Mat3 S;
    S << 1, 0, 0,
         0,-1, 0,
         0, 0, 1;
    K = K * S;
    R = S * R;
  }
  if (K(0,0) < 0) {
    Mat3 S;
    S << -1, 0, 0,
          0, 1, 0,
          0, 0, 1;
    K = K * S;
    R = S * R;
  }

  // Compute translation.
  Eigen::PartialPivLU<Mat3> lu(K);
  Vec3 t = lu.solve(P.col(3));

  if (R.determinant()<0) {
    R = -R;
    t = -t;
  }

  // scale K so that K(2,2) = 1
  K = K / K(2,2);

  *Kp = K;
  *Rp = R;
  *tp = t;
}

Mat3 F_from_P(const Mat34 & P1, const Mat34 & P2)
{
  Mat3 F12;

  using Mat24 = Eigen::Matrix<double, 2, 4>;
  Mat24 X1 = P1.block<2, 4>(1, 0);
  Mat24 X2;  X2 << P1.row(2), P1.row(0);
  Mat24 X3 = P1.block<2, 4>(0, 0);
  Mat24 Y1 = P2.block<2, 4>(1, 0);
  Mat24 Y2;  Y2 << P2.row(2), P2.row(0);
  Mat24 Y3 = P2.block<2, 4>(0, 0);


  Mat4 X1Y1, X2Y1, X3Y1, X1Y2, X2Y2, X3Y2, X1Y3, X2Y3, X3Y3;
  X1Y1 << X1, Y1;  X2Y1 << X2, Y1;  X3Y1 << X3, Y1;
  X1Y2 << X1, Y2;  X2Y2 << X2, Y2;  X3Y2 << X3, Y2;
  X1Y3 << X1, Y3;  X2Y3 << X2, Y3;  X3Y3 << X3, Y3;


  F12 <<
    X1Y1.determinant(), X2Y1.determinant(), X3Y1.determinant(),
    X1Y2.determinant(), X2Y2.determinant(), X3Y2.determinant(),
    X1Y3.determinant(), X2Y3.determinant(), X3Y3.determinant();

  return F12;
}

Vec2 Project(const Mat34 &P, const Vec3 &X) {
  return Vec3(P * X.homogeneous()).hnormalized();
}

void Project(const Mat34 &P, const Mat3X &X, Mat2X *x) {
  x->resize(2, X.cols());
  for (size_t c = 0; c < static_cast<size_t>(X.cols()); ++c) {
    x->col(c) = Project(P, Vec3(X.col(c)));
  }
}

void Project(const Mat34 &P, const Mat4X &X, Mat2X *x) {
  x->resize(2, X.cols());
  for (Mat4X::Index c = 0; c < X.cols(); ++c) {
    const Vec3 hx = P * X.col(c);
    x->col(c) = hx.hnormalized();
  }
}

Mat2X Project(const Mat34 &P, const Mat3X &X) {
  Mat2X x(2, X.cols());
  Project(P, X, &x);
  return x;
}

Mat2X Project(const Mat34 &P, const Mat4X &X) {
  Mat2X x(2, X.cols());
  Project(P, X, &x);
  return x;
}

void HomogeneousToEuclidean(const Vec4 &H, Vec3 *X) {
  double w = H(3);
  *X << H(0) / w, H(1) / w, H(2) / w;
}

void EuclideanToHomogeneous(const Mat &X, Mat *H) {
  Mat::Index d = X.rows();
  Mat::Index n = X.cols();
  H->resize(d + 1, n);
  H->block(0, 0, d, n) = X;
  H->row(d).setOnes();
}

double Depth(const Mat3 &R, const Vec3 &t, const Vec3 &X) {
  return (R*X)[2] + t[2];
}

void HomogeneousToEuclidean(const Mat &H, Mat *X) {
  Mat::Index d = H.rows() - 1;
  Mat::Index n = H.cols();
  X->resize(d, n);
  for (Mat::Index i = 0; i < n; ++i) {
    double h = H(d, i);
    for (int j = 0; j < d; ++j) {
      (*X)(j, i) = H(j, i) / h;
    }
  }
}

Mat3X EuclideanToHomogeneous(const Mat2X &x) {
  Mat3X h(3, x.cols());
  h.block(0, 0, 2, x.cols()) = x;
  h.row(2).setOnes();
  return h;
}

void EuclideanToHomogeneous(const Mat2X &x, Mat3X *h) {
  h->resize(3, x.cols());
  h->block(0, 0, 2, x.cols()) = x;
  h->row(2).setOnes();
}

void HomogeneousToEuclidean(const Mat3X &h, Mat2X *e) {
  e->resize(2, h.cols());
  e->row(0) = h.row(0).array() / h.row(2).array();
  e->row(1) = h.row(1).array() / h.row(2).array();
}

void EuclideanToNormalizedCamera(const Mat2X &x, const Mat3 &K, Mat2X *n) {
  Mat3X x_image_h;
  EuclideanToHomogeneous(x, &x_image_h);
  Mat3X x_camera_h = K.inverse() * x_image_h;
  HomogeneousToEuclidean(x_camera_h, n);
}

void HomogeneousToNormalizedCamera(const Mat3X &x, const Mat3 &K, Mat2X *n) {
  Mat3X x_camera_h = K.inverse() * x;
  HomogeneousToEuclidean(x_camera_h, n);
}

/// Estimates the root mean square error (2D)
double RootMeanSquareError(const Mat2X &x_image,
  const Mat4X &X_world,
  const Mat34 &P) {
    const Mat2X::Index num_points = x_image.cols();
    const Mat2X dx = Project(P, X_world) - x_image;
    return std::sqrt(dx.squaredNorm() / num_points);
}

/// Estimates the root mean square error (2D)
double RootMeanSquareError(const Mat2X &x_image,
  const Mat3X &X_world,
  const Mat3 &K,
  const Mat3 &R,
  const Vec3 &t) {
    Mat34 P;
    P_From_KRt(K, R, t, &P);
    return RootMeanSquareError(x_image, X_world.colwise().homogeneous(), P);
}

} // namespace openMVG
