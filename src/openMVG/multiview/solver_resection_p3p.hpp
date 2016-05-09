
/*
 * Copyright (c) 2011, Laurent Kneip, ETH Zurich
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ETH Zurich nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_RESECTION_P3P_H_
#define OPENMVG_MULTIVIEW_RESECTION_P3P_H_

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/projection.hpp"

#include <iostream>
#include <cmath>

namespace openMVG {
namespace euclidean_resection {

typedef Eigen::Matrix<double, 5, 1> Vec5;

inline void solveQuartic( const Vec5 & factors, Vec4 & realRoots)
{
  double A = factors[0];
  double B = factors[1];
  double C = factors[2];
  double D = factors[3];
  double E = factors[4];

  double A_pw2 = A*A;
  double B_pw2 = B*B;
  double A_pw3 = A_pw2*A;
  double B_pw3 = B_pw2*B;
  double A_pw4 = A_pw3*A;
  double B_pw4 = B_pw3*B;

  double alpha = -3*B_pw2/(8*A_pw2)+C/A;
  double beta = B_pw3/(8*A_pw3)-B*C/(2*A_pw2)+D/A;
  double gamma = -3*B_pw4/(256*A_pw4)+B_pw2*C/(16*A_pw3)-B*D/(4*A_pw2)+E/A;

  double alpha_pw2 = alpha*alpha;
  double alpha_pw3 = alpha_pw2*alpha;

  std::complex<double> P (-alpha_pw2/12-gamma,0);
  std::complex<double> Q (-alpha_pw3/108+alpha*gamma/3-pow(beta,2)/8,0);
  std::complex<double> R = -Q/2.0+sqrt(pow(Q,2.0)/4.0+pow(P,3.0)/27.0);

  std::complex<double> U = pow(R,(1.0/3.0));
  std::complex<double> y;

  if (U.real() == 0)
    y = -5.0*alpha/6.0-pow(Q,(1.0/3.0));
  else
    y = -5.0*alpha/6.0-P/(3.0*U)+U;

  std::complex<double> w = sqrt(alpha+2.0*y);

  std::complex<double> temp;

  temp = -B/(4.0*A) + 0.5*(w+sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
  realRoots[0] = temp.real();
  temp = -B/(4.0*A) + 0.5*(w-sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
  realRoots[1] = temp.real();
  temp = -B/(4.0*A) + 0.5*(-w+sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
  realRoots[2] = temp.real();
  temp = -B/(4.0*A) + 0.5*(-w-sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
  realRoots[3] = temp.real();
}

/*
 *      Author: Laurent Kneip, adapted to openMVG by Pierre Moulon
 * Description: Compute the absolute pose of a camera using three 3D-to-2D correspondences
 *   Reference: [1] A Novel Parametrization of the P3P-Problem for a Direct Computation of
 *              Absolute Camera Position and Orientation
 *              Kneip, L.; Scaramuzza, D. ; Siegwart, R.
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

inline bool compute_P3P_Poses( const Mat3 & featureVectors, const Mat3 & worldPoints, Mat & solutions )
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

  if( f3[2] > 0 )
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

  Vec5 factors;
  factors << -f_2_pw2*p_2_pw4 - p_2_pw4*f_1_pw2 - p_2_pw4,

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
    +f_2_pw2*p_2_pw2*d_12_pw2*b_pw2;

  // Computation of roots

  Vec4 realRoots;
  solveQuartic( factors, realRoots );

  // Backsubstitution of each solution

  for(int i=0; i<4; ++i)
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

struct P3PSolver {
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 4};
  // Solve the problem of camera pose.
  static void Solve(const Mat &pt2D, const Mat &pt3D, std::vector<Mat34> *models)
  {
    Mat3 R;
    Vec3 t;
    Mat34 P;
    assert(2 == pt2D.rows());
    assert(3 == pt3D.rows());
    assert(pt2D.cols() == pt3D.cols());
    Mat solutions = Mat(3, 4*4);
    Mat3 pt2D_3x3;
    pt2D_3x3.block<2,3>(0,0) = pt2D;
    pt2D_3x3.row(2).fill(1);
    pt2D_3x3.col(0).normalize();
    pt2D_3x3.col(1).normalize();
    pt2D_3x3.col(2).normalize();
    Mat3 pt3D_3x3 = pt3D;
    if (compute_P3P_Poses( pt2D_3x3, pt3D_3x3, solutions))
    {
      for (size_t i=0; i < 4; ++i)  {
        R = solutions.block<3,3>(0,i*4+1);
        t = -R * solutions.col(i*4);
        P_From_KRt(Mat3::Identity(), R, t, &P); // K = Id
        models->push_back(P);
      }
    }
  }

  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  static double Error(const Mat34 & P, const Vec2 & pt2D, const Vec3 & pt3D) {
    return (pt2D - Project(P, pt3D)).norm();
  }
};

class P3P_ResectionKernel_K {
 public:
  typedef Mat34 Model;
  enum { MINIMUM_SAMPLES = 3 };

  P3P_ResectionKernel_K(const Mat2X &x_camera, const Mat3X &X, const Mat3 &K = Mat3::Identity())
    :x_image_(x_camera), X_(X), K_(K)
  {
    assert(x_camera.cols() == X.cols());
    // Conversion from image coordinates to normalized camera coordinates
    Mat3X x_image_h;
    EuclideanToHomogeneous(x_image_, &x_image_h);
    x_camera_ = K_.inverse() * x_image_h;
    for(size_t i = 0; i < x_camera_.cols(); ++i)
      x_camera_.col(i).normalize();
  }

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const {
    const Mat3 pt2D_3x3 ( ExtractColumns(x_camera_, samples) );
    const Mat3 pt3D_3x3 ( ExtractColumns(X_, samples) );
    Mat solutions(3, 4*4);
    if (compute_P3P_Poses( pt2D_3x3, pt3D_3x3, solutions))
    {
      Mat34 P;
      Mat3 R;
      Vec3 t;
      for (size_t i=0; i < 4; ++i)  {
        R = solutions.block<3,3>(0,i*4+1);
        t = -R * solutions.col(i*4);
        P_From_KRt(K_, R, t, &P);
        models->push_back(P);
      }
    }
  }

  double Error(size_t sample, const Model &model) const {
    const Vec3 X = X_.col(sample);
    const Mat2X error = Project(model, X) - x_image_.col(sample);
    return error.col(0).norm();
  }

  size_t NumSamples() const {
    return static_cast<size_t>(x_camera_.cols());
  }

 private:
  Mat2X x_image_; // camera coordinates
  Mat3X x_camera_; // camera coordinates (normalized)
  Mat3X X_;        // 3D points
  Mat3 K_;
};

}  // namespace euclidean_resection
}  // namespace openMVG

#endif // OPENMVG_MULTIVIEW_RESECTION_P3P_H_

