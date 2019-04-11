// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Michael Persson
// Adapted to openMVG by Romain Janvier

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "openMVG/multiview/solver_resection_p3p_nordberg.hpp"
#include "openMVG/multiview/projection.hpp"

namespace openMVG
{
namespace euclidean_resection
{

bool computePosesNordberg(
    const Mat &bearing_vectors,
    const Mat &X,
    std::vector<std::tuple<Mat3, Vec3>> &rotation_translation_solutions)
{
  // Extraction of 3D points vectors
  Vec3 P1 = X.col(0);
  Vec3 P2 = X.col(1);
  Vec3 P3 = X.col(2);

  // Extraction of feature vectors
  Vec3 f1 = bearing_vectors.col(0);
  Vec3 f2 = bearing_vectors.col(1);
  Vec3 f3 = bearing_vectors.col(2);

  f1.normalize();
  f2.normalize();
  f3.normalize();

  double b12 = -2.0 * (f1.dot(f2));
  double b13 = -2.0 * (f1.dot(f3));
  double b23 = -2.0 * (f2.dot(f3));

  // implicit creation of Vec3, can be removed
  Vec3 d12 = P1 - P2;
  Vec3 d13 = P1 - P3;
  Vec3 d23 = P2 - P3;
  Vec3 d12xd13(d12.cross(d13));

  double a12 = d12.squaredNorm();
  double a13 = d13.squaredNorm();
  double a23 = d23.squaredNorm();

  //a*g^3 + b*g^2 + c*g + d = 0
  double c31 = -0.5 * b13;
  double c23 = -0.5 * b23;
  double c12 = -0.5 * b12;
  double blob = (c12 * c23 * c31 - 1.0);

  double s31_squared = 1.0 - c31 * c31;
  double s23_squared = 1.0 - c23 * c23;
  double s12_squared = 1.0 - c12 * c12;

  double p3 = a13 * (a23 * s31_squared - a13 * s23_squared);
  double p2 = 2.0 * blob * a23 * a13 + a13 * (2.0 * a12 + a13) * s23_squared + a23 * (a23 - a12) * s31_squared;
  double p1 = a23 * (a13 - a23) * s12_squared - a12 * a12 * s23_squared - 2.0 * a12 * (blob * a23 + a13 * s23_squared);
  double p0 = a12 * (a12 * s23_squared - a23 * s12_squared);

  double g = 0.0;

  // p3 is essentially det(D2) so it is definietly > 0 or it is degen
  //if (std::abs(p3) >= std::abs(p0) || true)
  //{
    p3 = 1.0 / p3;
    p2 *= p3;
    p1 *= p3;
    p0 *= p3;

    // get sharpest real root of above...
    g = cubick(p2, p1, p0);
  //}
  // else
  //{

    // lower numerical performance
    //g = 1.0 / (cubick(p1 / p0, p2 / p0, p3 / p0));
  //}

  //  cout<<"g: "<<g<<endl;

  // we can swap D1,D2 and the coeffs!
  // oki, Ds are:
  //D1=M12*XtX(2,2) - M23*XtX(1,1);
  //D2=M23*XtX(3,3) - M13*XtX(2,2);

  //[    a23 - a23*g,                 (a23*b12)/2,              -(a23*b13*g)/2]
  //[    (a23*b12)/2,           a23 - a12 + a13*g, (a13*b23*g)/2 - (a12*b23)/2]
  //[ -(a23*b13*g)/2, (a13*b23*g)/2 - (a12*b23)/2,         g*(a13 - a23) - a12]

  // gain 13 ns...
  double A00 = a23 * (1.0 - g);
  double A01 = (a23 * b12) * 0.5;
  double A02 = (a23 * b13 * g) * (-0.5);
  double A11 = a23 - a12 + a13 * g;
  double A12 = b23 * (a13 * g - a12) * 0.5;
  double A22 = g * (a13 - a23) - a12;

  Mat3 A;
  A << A00, A01, A02,
       A01, A11, A12,
       A02, A12, A22;

  // get sorted eigenvalues and eigenvectors given that one should be zero...
  Mat3 V;
  Vec3 L;

  eigwithknown0(A, V, L);

  double v = std::sqrt(std::max(0.0, -L(1) / L(0)));

  int valid = 0;
  std::array<Vec3, 4> Ls;

  // use the t=Vl with t2,st2,t3 and solve for t3 in t2
  { //+v
    double s = v;

    double w2 = 1.0 / (s * V(0,1) - V(0,0));
    double w0 = (V(1,0) - s * V(1,1)) * w2;
    double w1 = (V(2,0) - s * V(2,1)) * w2;

    double a = 1.0 / ((a13 - a12) * w1 * w1 - a12 * b13 * w1 - a12);
    double b = (a13 * b12 * w1 - a12 * b13 * w0 - 2.0 * w0 * w1 * (a12 - a13)) * a;
    double c = ((a13 - a12) * w0 * w0 + a13 * b12 * w0 + a13) * a;

    if (b * b - 4.0 * c >= 0.0)
    {
      double tau1, tau2;
      root2real(b, c, tau1, tau2);
      if (tau1 > 0)
      {
        double tau = tau1;
        double d = a23 / (tau * (b23 + tau) + 1.0);
        double l2 = std::sqrt(d);
        double l3 = tau * l2;

        double l1 = w0 * l2 + w1 * l3;
        if (l1 >= 0.0)
        {
          Ls[valid] = Vec3(l1, l2, l3);
          ++valid;
        }
      }
      if (tau2 > 0.0)
      {
        double tau = tau2;
        double d = a23 / (tau * (b23 + tau) + 1.0);
        double l2 = std::sqrt(d);
        double l3 = tau * l2;
        double l1 = w0 * l2 + w1 * l3;
        if (l1 >= 0.0)
        {
          Ls[valid] = Vec3(l1, l2, l3);
          ++valid;
        }
      }
    }
  }

  { //+v
    double s = -v;
    double w2 = 1.0 / (s * V(0, 1) - V(0, 0));
    double w0 = (V(1, 0) - s * V(1, 1)) * w2;
    double w1 = (V(2, 0) - s * V(2, 1)) * w2;

    double a = 1.0 / ((a13 - a12) * w1 * w1 - a12 * b13 * w1 - a12);
    double b = (a13 * b12 * w1 - a12 * b13 * w0 - 2.0 * w0 * w1 * (a12 - a13)) * a;
    double c = ((a13 - a12) * w0 * w0 + a13 * b12 * w0 + a13) * a;

    if (b * b - 4.0 * c >= 0)
    {
      double tau1, tau2;

      root2real(b, c, tau1, tau2);
      if (tau1 > 0)
      {
        double tau = tau1;
        double d = a23 / (tau * (b23 + tau) + 1.0);
        double l2 = std::sqrt(d);

        double l3 = tau * l2;

        double l1 = w0 * l2 + w1 * l3;
        if (l1 >= 0)
        {
          Ls[valid] = Vec3(l1, l2, l3);
          ++valid;
        }
      }
      if (tau2 > 0)
      {
        double tau = tau2;
        double d = a23 / (tau * (b23 + tau) + 1.0);
        double l2 = std::sqrt(d);

        double l3 = tau * l2;

        double l1 = w0 * l2 + w1 * l3;
        if (l1 >= 0)
        {
          Ls[valid] = Vec3(l1, l2, l3);
          ++valid;
        }
      }
    }
  }

  // if constexpr (refinement_iterations>0)
  for (int i = 0; i < valid; ++i)
  {
    gauss_newton_refineL(Ls[i], a12, a13, a23, b12, b13, b23);
  }

  Vec3 ry1, ry2, ry3;
  Vec3 yd1;
  Vec3 yd2;
  Vec3 yd1xd2;
  Mat3 Xmat;
  Xmat << d12(0), d13(0), d12xd13(0),
          d12(1), d13(1), d12xd13(1),
          d12(2), d13(2), d12xd13(2);
  
  Xmat = Xmat.inverse().eval();

  for (int i = 0; i < valid; ++i)
  {
    // compute the rotation:
    ry1 = f1 * Ls[i](0);
    ry2 = f2 * Ls[i](1);
    ry3 = f3 * Ls[i](2);

    yd1 = ry1 - ry2;
    yd2 = ry1 - ry3;
    yd1xd2 = yd1.cross(yd2);

    Mat3 Ymat;
    Ymat << yd1(0), yd2(0), yd1xd2(0),
            yd1(1), yd2(1), yd1xd2(1),
            yd1(2), yd2(2), yd1xd2(2);

    Mat3 Rs = Ymat * Xmat;
    rotation_translation_solutions.emplace_back(Rs, ry1 - Rs * P1);
  }
  return valid;
}

void P3PSolver_Nordberg::Solve(
    const Mat &bearing_vectors,
    const Mat &X, // 3D points
    std::vector<Mat34> *models)
{
  assert(3 == bearing_vectors.rows());
  assert(3 == X.rows());
  assert(bearing_vectors.cols() == X.cols());
  Mat34 P;
  std::vector<std::tuple<Mat3, Vec3>> rotation_translation_solutions;
  if (computePosesNordberg(bearing_vectors, X, rotation_translation_solutions))
  {
    for (const auto & rot_trans_it : rotation_translation_solutions) {
      Mat34 P;
      P_From_KRt(Mat3::Identity(),           // intrinsics
                  std::get<0>(rot_trans_it), // rotation
                  std::get<1>(rot_trans_it), // translation
                  &P);
      models->push_back(P);
    }
  }
};

double P3PSolver_Nordberg::Error
(
  const Mat34 & P,
  const Vec3 & bearing_vector,
  const Vec3 & pt3D
)
{
  const auto new_bearing = (P * pt3D.homogeneous()).normalized();
  return 1.0 - (bearing_vector.dot(new_bearing));
}

} // namespace euclidean_resection
} // namespace openMVG