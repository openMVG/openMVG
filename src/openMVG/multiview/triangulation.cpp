// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {

// HZ 12.2 p.312
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

  const Eigen::JacobiSVD<Mat4> svd( design, Eigen::ComputeFullV );
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

bool TriangulateDLT
(
  const Mat3 &R0,
  const Vec3 &t0,
  const Vec3 &x0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &x1,
  Vec3 *X
)
{
  Mat34 P0, P1;
  P0.block<3,3>(0,0) = R0;
  P1.block<3,3>(0,0) = R1;
  P0.block<3,1>(0,3) = t0;
  P1.block<3,1>(0,3) = t1;
  TriangulateDLT(P0, x0, P1, x1, X);
  return x0.dot(R0 * (*X + R0.transpose() * t0)) > 0.0 &&
         x1.dot(R1 * (*X + R1.transpose() * t1)) > 0.0;
}

// Helper function
// Compute Relative motion between two absolute poses parameterized by Rt
// Rotate one bearing vector according to the relative motion
inline void AbsoluteToRelative(
  const Mat3 &R0,
  const Vec3 &t0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &x0,
  Mat3 &R,
  Vec3 &t,
  Vec3 &Rx0
)
{
  R = R1 * R0.transpose();
  t = t1 - R * t0;
  Rx0 = R * x0;
}

// Helper function
// Compute the 3D point from the outcome of an angular triangulation solver
// Return true if the point pass the cheirality test, false otherwise
// Lee & Civera Eq. (11) and  Table 1 - 4)
inline bool Compute3DPoint(
  const Vec3 &mprime0,
  const Vec3 &mprime1,
  const Vec3 &t,
  const Mat3 &R1,
  const Vec3 &t1,
  Vec3 *X)
{
  const Vec3 z = mprime1.cross(mprime0);
  const double z_squared = z.squaredNorm();
  const double lambda0 = z.dot(t.cross(mprime1)) / z_squared;
  const double lambda1 = z.dot(t.cross(mprime0)) / z_squared;
  const Vec3 xprime1 = t + lambda0 * mprime0;

  // x'1 is into the frame of camera1 convert it into the world frame in order to obtain the 3D point
  *X = R1.transpose() * (xprime1 - t1);

  // make and return the result of the cheirality test
  return lambda0 > 0.0 && lambda1 > 0.0;
}

bool TriangulateL1Angular
(
  const Mat3 &R0,
  const Vec3 &t0,
  const Vec3 &x0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &x1,
  Vec3 *X_euclidean
)
{
  // Table 1 - 1) we compute m0 (Rx0) and m1 (x1)
  Mat3 R;
  Vec3 t, Rx0;
  AbsoluteToRelative(R0, t0, R1, t1, x0, R, t, Rx0);

  // Table 1 - 2) obtain m'0 and m'1
  if(Rx0.normalized().cross(t).squaredNorm() <= x1.normalized().cross(t).squaredNorm())
  {
    const Vec3 n1 = x1.cross(t).normalized();
    // Eq. (12)
    const Vec3 mprime0 = Rx0 - Rx0.dot(n1) * n1;
    return Compute3DPoint(mprime0, x1, t, R1, t1, X_euclidean);
  }
  else
  {
    const Vec3 n0 = Rx0.cross(t).normalized();
    // Eq. (13)
    const Vec3 mprime1 = x1 - x1.dot(n0) * n0;
    return Compute3DPoint(Rx0, mprime1, t, R1, t1, X_euclidean);
  }
}

bool TriangulateLInfinityAngular
(
  const Mat3 &R0,
  const Vec3 &t0,
  const Vec3 &x0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &x1,
  Vec3 *X_euclidean
)
{
  // Table 1 - 1) we compute m0 (Rx0) and m1 (x1)
  Mat3 R;
  Vec3 t, Rx0;
  AbsoluteToRelative(R0, t0, R1, t1, x0, R, t, Rx0);

  // cf. 7. Lemma 2
  const Vec3 Rx0_norm = Rx0.normalized();
  const Vec3 x1_norm = x1.normalized();
  const Vec3 na = (Rx0_norm + x1_norm).cross(t);
  const Vec3 nb = (Rx0_norm - x1_norm).cross(t);

  const Vec3 nprime = na.squaredNorm() >= nb.squaredNorm() ? na.normalized() : nb.normalized();

  const Vec3 mprime0 = Rx0 - (Rx0.dot(nprime)) * nprime;
  const Vec3 mprime1 = x1 - (x1.dot(nprime)) * nprime;

  return Compute3DPoint(mprime0, mprime1, t, R1, t1, X_euclidean);
}

bool TriangulateIDWMidpoint(
  const Mat3 & R0,
  const Vec3 & t0,
  const Vec3 & x0,
  const Mat3 & R1,
  const Vec3 & t1,
  const Vec3 & x1,
  Vec3* X_euclidean
)
{
  // absolute to relative
  Mat3 R;
  Vec3 t, Rx0;
  AbsoluteToRelative(R0, t0, R1, t1, x0, R, t, Rx0);

  const double p_norm = Rx0.cross(x1).norm();
  const double q_norm = Rx0.cross(t).norm();
  const double r_norm = x1.cross(t).norm();

  // Eq. (10)
  const auto xprime1 = ( q_norm / (q_norm + r_norm) )
    * ( t + (r_norm / p_norm) * (Rx0 + x1) );

  // relative to absolute
  *X_euclidean = R1.transpose() * (xprime1 - t1);

  // Eq. (7)
  const Vec3 lambda0_Rx0 = (r_norm / p_norm) * Rx0;
  const Vec3 lambda1_x1 = (q_norm / p_norm) * x1;

  // Eq. (9) - test adequation
  return (t + lambda0_Rx0 - lambda1_x1).squaredNorm()
    <
    std::min(std::min(
      (t + lambda0_Rx0 + lambda1_x1).squaredNorm(),
      (t - lambda0_Rx0 - lambda1_x1).squaredNorm()),
      (t - lambda0_Rx0 + lambda1_x1).squaredNorm());
}

bool Triangulate2View
(
  const Mat3 &R0,
  const Vec3 &t0,
  const Vec3 &bearing0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &bearing1,
  Vec3 &X,
  ETriangulationMethod etri_method
)
{
  switch (etri_method)
  {
    case ETriangulationMethod::DIRECT_LINEAR_TRANSFORM:
      return TriangulateDLT(R0, t0, bearing0, R1, t1, bearing1, &X);
    break;
    case ETriangulationMethod::L1_ANGULAR:
      return TriangulateL1Angular(R0, t0, bearing0, R1, t1, bearing1, &X);
    break;
    case ETriangulationMethod::LINFINITY_ANGULAR:
      return TriangulateLInfinityAngular(R0, t0, bearing0, R1, t1, bearing1, &X);
    break;
    case ETriangulationMethod::INVERSE_DEPTH_WEIGHTED_MIDPOINT:
      return TriangulateIDWMidpoint(R0, t0, bearing0, R1, t1, bearing1, &X);
    break;
    default:
      return false;
  }
  return false;
}

}  // namespace openMVG
