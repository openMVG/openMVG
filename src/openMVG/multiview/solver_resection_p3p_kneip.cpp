// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2011, Laurent Kneip.
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_p3p_kneip.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/numeric/poly.h"

#include <cmath>

namespace openMVG {
namespace euclidean_resection {

/*
 *      Author: Laurent Kneip, adapted to openMVG by Pierre Moulon
 * Description: Compute the absolute pose of a camera using three 3D-to-2D correspondences
 *   Reference: [1] A Novel Parametrization of the P3P-Problem for a Direct Computation of
 *              Absolute Camera Position and Orientation
 *              Kneip, L.; Scaramuzza, D.; Siegwart, R.
 *              CVPR 2011
 *
 *       Input: featureVectors: 3x3 matrix with UNITARY feature vectors (each column is a vector)
 *              worldPoints: 3x3 matrix with corresponding 3D world points (each column is a point)
 *              solutions: 3x16 matrix that will contain the solutions
 *                         form: [ C1,R1, C2,R2 ... ]
 *                         the obtained orientation matrices are defined as transforming points from the cam to the world frame
 *      Output: bool: true if correct execution
 *                    false if world points aligned
 */

inline bool compute_P3P_Poses
(
  const Mat3 & featureVectors,
  const Mat3 & worldPoints,
  Mat & solutions
)
{
  solutions = Mat(3, 4*4);

  // Extraction of world points

  Vec3 P1 = worldPoints.col(0);
  Vec3 P2 = worldPoints.col(1);
  Vec3 P3 = worldPoints.col(2);

  // Verification that world points are not colinear

  if ( ((P2-P1).cross(P3-Vec3(1,1,1))).norm() == 0)
    return false;

  // Extraction of feature vectors

  Vec3 f1 = featureVectors.col(0);
  Vec3 f2 = featureVectors.col(1);
  Vec3 f3 = featureVectors.col(2);

  // Creation of intermediate camera frame

  Vec3 e1 = f1;
  Vec3 e3 = (f1.cross(f2)).normalized();
  Vec3 e2 = e3.cross(e1);

  Mat3 T;
  T.row(0) = e1;
  T.row(1) = e2;
  T.row(2) = e3;

  f3 = T * f3;

  // Reinforce that f3[2] > 0 for having theta in [0;pi]

  if (f3[2] > 0 )
  {
    f1 = featureVectors.col(1);
    f2 = featureVectors.col(0);
    f3 = featureVectors.col(2);

    e1 = f1;
    e3 = (f1.cross(f2)).normalized();
    e2 = e3.cross(e1);

    T.row(0) = e1;
    T.row(1) = e2;
    T.row(2) = e3;

    f3 = T * f3;

    P1 = worldPoints.col(1);
    P2 = worldPoints.col(0);
    P3 = worldPoints.col(2);
  }

  // Creation of intermediate world frame

  Vec3 n1 = (P2-P1).normalized();
  Vec3 n3 = (n1.cross(P3-P1)).normalized();
  Vec3 n2 = n3.cross(n1);

  Mat3 N;
  N.row(0) = n1;
  N.row(1) = n2;
  N.row(2) = n3;

  // Extraction of known parameters

  P3 = N * (P3 - P1);

  double d_12 = (P2 - P1).norm();
  double f_1 = f3[0] / f3[2];
  double f_2 = f3[1] / f3[2];
  double p_1 = P3[0];
  double p_2 = P3[1];

  double cos_beta = f1.transpose() * f2;
  double b = 1.0 / (1.0 - pow(cos_beta,2)) - 1.0;

  if (cos_beta < 0)
    b = -sqrt(b);
  else
    b = sqrt(b);

  // Definition of temporary variables for avoiding multiple computation

  double f_1_pw2 = pow(f_1,2);
  double f_2_pw2 = pow(f_2,2);
  double p_1_pw2 = pow(p_1,2);
  double p_1_pw3 = p_1_pw2 * p_1;
  double p_1_pw4 = p_1_pw3 * p_1;
  double p_2_pw2 = pow(p_2,2);
  double p_2_pw3 = p_2_pw2 * p_2;
  double p_2_pw4 = p_2_pw3 * p_2;
  double d_12_pw2 = pow(d_12,2);
  double b_pw2 = pow(b,2);

  // Computation of factors of 4th degree polynomial

  const std::array<double,5> factors = {
    {-f_2_pw2*p_2_pw4 - p_2_pw4*f_1_pw2 - p_2_pw4,

      2.*p_2_pw3*d_12*b +
      2.*f_2_pw2*p_2_pw3*d_12*b
      -2.*f_2*p_2_pw3*f_1*d_12,

      -f_2_pw2*p_2_pw2*p_1_pw2
      -f_2_pw2*p_2_pw2*d_12_pw2*b_pw2
      -f_2_pw2*p_2_pw2*d_12_pw2
      +f_2_pw2*p_2_pw4
      +p_2_pw4*f_1_pw2
      +2.*p_1*p_2_pw2*d_12
      +2.*f_1*f_2*p_1*p_2_pw2*d_12*b
      -p_2_pw2*p_1_pw2*f_1_pw2
      +2.*p_1*p_2_pw2*f_2_pw2*d_12
      -p_2_pw2*d_12_pw2*b_pw2
      -2.*p_1_pw2*p_2_pw2,

      2.*p_1_pw2*p_2*d_12*b
      +2.*f_2*p_2_pw3*f_1*d_12
      -2.*f_2_pw2*p_2_pw3*d_12*b
      -2.*p_1*p_2*d_12_pw2*b,

      -2.*f_2*p_2_pw2*f_1*p_1*d_12*b
      +f_2_pw2*p_2_pw2*d_12_pw2
      +2.*p_1_pw3*d_12
      -p_1_pw2*d_12_pw2
      +f_2_pw2*p_2_pw2*p_1_pw2
      -p_1_pw4
      -2.*f_2_pw2*p_2_pw2*p_1*d_12
      +p_2_pw2*f_1_pw2*p_1_pw2
      +f_2_pw2*p_2_pw2*d_12_pw2*b_pw2}};

  // Computation of roots
  std::array<double, 4> realRoots;
  solveQuarticPolynomial( factors, realRoots );

  // Backsubstitution of each solution

  for (int i=0; i<4; ++i)
  {
    double cot_alpha = (-f_1*p_1/f_2-realRoots[i]*p_2+d_12*b)/(-f_1*realRoots[i]*p_2/f_2+p_1-d_12);

    double cos_theta = realRoots[i];
    double sin_theta = sqrt(1-pow(realRoots[i],2));
    double sin_alpha = sqrt(1./(pow(cot_alpha,2)+1));
    double cos_alpha = sqrt(1.-pow(sin_alpha,2));

    if (cot_alpha < 0)
      cos_alpha = -cos_alpha;

    if (!is_finite(sin_theta))
      sin_theta = 0;

    Vec3 C (d_12*cos_alpha*(sin_alpha*b+cos_alpha),
        cos_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha),
        sin_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha));

    C = P1 + N.transpose()*C;

    Mat3 R;
    R<< -cos_alpha,    -sin_alpha*cos_theta,  -sin_alpha*sin_theta,
         sin_alpha,    -cos_alpha*cos_theta,  -cos_alpha*sin_theta,
                0,              -sin_theta,             cos_theta;

    R = N.transpose()*R.transpose()*T;

    solutions.col(i*4) = C;
    solutions.block<3,3>(0,i*4+1) = R.transpose();
  }
  return true;
}

void P3PSolver_Kneip::Solve
(
  const Mat & bearing_vectors,
  const Mat & pt3D,
  std::vector<Mat34> * models
)
{
  assert(3 == bearing_vectors.rows());
  assert(3 == pt3D.rows());
  assert(bearing_vectors.cols() == pt3D.cols());

  Mat solutions = Mat(3, 4*4);
  if (compute_P3P_Poses( bearing_vectors, pt3D, solutions))
  {
    Mat3 R;
    Vec3 t;
    Mat34 P;
    for (size_t i=0; i < 4; ++i)  {
      R = solutions.block<3,3>(0,i*4+1);
      t = -R * solutions.col(i*4);
      P_From_KRt(Mat3::Identity(), R, t, &P); // K = Id
      models->push_back(P);
    }
  }
}

double P3PSolver_Kneip::Error
(
  const Mat34 & P,
  const Vec3 & bearing_vector,
  const Vec3 & pt3D
)
{
  const auto new_bearing = (P * pt3D.homogeneous()).normalized();
  return 1.0 - (bearing_vector.dot(new_bearing));
}

}  // namespace euclidean_resection
}  // namespace openMVG
