// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Michael Persson
// Adapted to openMVG by Romain Janvier and Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_p3p_nordberg.hpp"
#include "openMVG/multiview/projection.hpp"

#include <array>

namespace openMVG
{
namespace euclidean_resection
{

/**
 * @brief Refine a valid solution with a Gauss-Newton Solver.
 * @param L Valid lambdas
 * @param a12 is the squared distance between X1 and X2
 * @param a13 is the squared distance between X1 and X3
 * @param a23 is the squared distance between X2 and X3
 * @param b12 is the cosine of the angle between bearing vector 1 and bearing vector 2
 * @param b13 is the cosine of the angle between bearing vector 1 and bearing vector 3
 * @param b23 is the cosine of the angle between bearing vector 2 and bearing vector 3
 * The paper note it rarely improve after two iterations. The original implementation use 5 iterations.
 */
static void gauss_newton_refineL(Vec3 &L,
                          const double & a12, const double & a13, const double & a23,
                          const double & b12, const double & b13, const double & b23)
{
  // const expr makes it easier for the compiler to unroll
  // TODO(RJ:) I have hardcoded the number of iterations here, it's a template parameter in the original implementation
  for (int i = 0; i < 5; ++i)
  {
    double l1 = L(0);
    double l2 = L(1);
    double l3 = L(2);
    double r1 = l1 * l1 + l2 * l2 + b12 * l1 * l2 - a12;
    double r2 = l1 * l1 + l3 * l3 + b13 * l1 * l3 - a13;
    double r3 = l2 * l2 + l3 * l3 + b23 * l2 * l3 - a23;

    if (std::abs(r1) + std::abs(r2) + std::abs(r3) < 1e-10)
      break;

    double dr1dl1 = 2.0 * l1 + b12 * l2;
    double dr1dl2 = 2.0 * l2 + b12 * l1;

    double dr2dl1 = 2.0 * l1 + b13 * l3;
    double dr2dl3 = 2.0 * l3 + b13 * l1;

    double dr3dl2 = 2.0 * l2 + b23 * l3;
    double dr3dl3 = 2.0 * l3 + b23 * l2;

    Vec3 r(r1, r2, r3);

    // or skip the inverse and make it explicit...
    {
      double v0 = dr1dl1;
      double v1 = dr1dl2;
      double v3 = dr2dl1;
      double v5 = dr2dl3;
      double v7 = dr3dl2;
      double v8 = dr3dl3;
      double det = 1.0 / (-v0 * v5 * v7 - v1 * v3 * v8);

      Mat3 Ji;
      Ji << -v5 * v7, -v1 * v8, v1 * v5,
            -v3 * v8, v0 * v8, -v0 * v5,
             v3 * v7, -v0 * v7, -v1 * v3;
      Vec3 L1 = Vec3(L) - det * (Ji * r);
      //%l=l - g*H\G;%inv(H)*G
      //L=L - g*J\r;
      //% works because the size is ok!
      {
        double l1 = L1(0);
        double l2 = L1(1);
        double l3 = L1(2);
        double r11 = l1 * l1 + l2 * l2 + b12 * l1 * l2 - a12;
        double r12 = l1 * l1 + l3 * l3 + b13 * l1 * l3 - a13;
        double r13 = l2 * l2 + l3 * l3 + b23 * l2 * l3 - a23;
        if (std::abs(r11) + std::abs(r12) + std::abs(r13) > std::abs(r1) + std::abs(r2) + std::abs(r3))
        {
          break;
        }
        else
          L = L1;
      }
    }
  }
};

static inline bool root2real(const double & b, const double & c, double & r1, double & r2)
{
  double v = b * b -4.0 * c;
  if (v < 0.0) {
      r1 = r2 = 0.5 * b;
      return false;
  }
  double y = std::sqrt(v);
  if (b < 0.0) {
      r1 = 0.5 * (-b + y);
      r2 = 0.5 * (-b - y);
  } else {
      r1 = 2.0 * c / (-b + y);
      r2 = 2.0 * c / (-b - y);
  }
  return true;
};

/**
 * @brief This function finds a single root of the cubic polynomial equation
 * @param b Coefficient of quadratic parameter
 * @param c Coefficient of linear parameter
 * @param d Coefficient of scalar parameter
 * @return the single root
 *
 * h(r) = r^3 + b*r^2 + c*r + d = 0
 *
 * The return root is as stable as possible in the sense that it has as high
 * derivative as possible.  The solution is found by simple Newton-Raphson
 * iterations, and the trick is to choose the intial solution r0 in a clever
 * way.
 *
 * The intial solution is found by considering 5 cases:
 *
 * Cases I and II: h has no stationary points. In this case its derivative
 * is positive.  The inital solution to the NR-iteration is r0 here h has
 * minimal derivative.
 *
 * Case III, IV, and V: has two stationary points, t1 < t2.  In this case,
 * h has negative derivative between t1 and t2.  In these cases, we can make
 * a second order approximation of h around each of t1 and t2, and choose r0
 * as the leftmost or rightmost root of these approximations, depending on
 * whether two, one, or both of h(t1) and h(t2) are > 0.
*/
static double cubick(const double &b, const double &c, const double &d)
{
  // Choose an initial solution
  double r0;
  // not monotonic
  if (b * b >= 3.0 * c)
  {
    // h has two stationary points, compute them
    // double t1 = t - std::sqrt(diff);
    double v = std::sqrt(b * b - 3.0 * c);
    double t1 = (-b - v) / (3.0);

    // Check if h(t1) > 0, in this case make a 2-order approx of h around t1
    double k = ((t1 + b) * t1 + c) * t1 + d;

    if (k > 0.0)
    {
      // Find leftmost root of 0.5*(r0 -t1)^2*(6*t1+2*b) +  k = 0
      r0 = t1 - std::sqrt(-k / (3.0 * t1 + b));
      // or use the linear comp too
      // r0 = t1 -
    }
    else
    {
      double t2 = (-b + v) / 3.0;
      k = ((t2 + b) * t2 + c) * t2 + d;
      // Find rightmost root of 0.5 * (r0 - t2)^2 * (6 * t2 +2 * b) + k1 = 0
      r0 = t2 + std::sqrt(-k / (3.0 * t2 + b));
    }
  }
  else
  {
    // r0=1.0/(cubick_inv(c/d,b/d,1.0/d));
    // about half work...
    // if(std::abs((((r0+b)*r0+c)*r0+d))>1e-10)
    r0 = -b / 3.0;
    if (std::abs(((3.0 * r0 + 2.0 * b) * r0 + c)) < 1e-4)
      r0 += 1;
    //else r0-=1;
    //double fx=(((r0+b)*r0+c)*r0+d); r0-=10; if(fx<0) r0+=20;
  }

  // Do ITER Newton-Raphson iterations
  // Break if position of root changes less than 1e-13
  // double starterr=std::abs(r0*(r0*(r0 + b) + c) + d);
  // TODO(RJ:) I have hardcoded the number of iteration here, it's a hardcoded in a define in the orginal implementation
  // according to the author, increasing it could lead to a better solution (more robust)
  for (unsigned int cnt = 0; cnt < 50; ++cnt)
  {
    double fx = (((r0 + b) * r0 + c) * r0 + d);

    if ((cnt < 7 || std::abs(fx) > 1e-13))
    {
      double fpx = ((3.0 * r0 + 2.0 * b) * r0 + c);
      r0 -= fx / fpx;
    }
    else
      break;
  }
  return r0;
};

/**
 * @brief eigwithknown0 eigen decomposition of a matrix which has a 0 eigen value
 * @param x the input matrix
 * @param E eigenvectors
 * @param L eigenvalues
 */
static void eigwithknown0(const Mat3 &x, Mat3 &E, Vec3 &L)
{
  // one eigenvalue is known to be 0.
  // the known one...
  L(2) = 0.0;

  Vec3 v3(x(3) * x(7) - x(6) * x(4),
          x(6) * x(1) - x(7) * x(0),
          x(4) * x(0) - x(3) * x(1));

  v3.normalize();

  double x01_squared = x(0, 1) * x(0, 1);
  // get the two other...
  double b = -x(0, 0) - x(1, 1) - x(2, 2);
  double c = -x01_squared - x(0, 2) * x(0, 2) - x(1, 2) * x(1, 2) +
             x(0, 0) * (x(1, 1) + x(2, 2)) + x(1, 1) * x(2, 2);
  double e1, e2;
  // roots(poly(x))
  root2real(b, c, e1, e2);

  if (std::abs(e1) < std::abs(e2))
    std::swap(e1, e2);
  L(0) = e1;
  L(1) = e2;

  double mx0011 = -x(0, 0) * x(1, 1);
  double prec_0 = x(0, 1) * x(1, 2) - x(0, 2) * x(1, 1);
  double prec_1 = x(0, 1) * x(0, 2) - x(0, 0) * x(1, 2);

  double e = e1;
  double tmp = 1.0 / (e * (x(0, 0) + x(1, 1)) + mx0011 - e * e + x01_squared);
  double a1 = -(e * x(0, 2) + prec_0) * tmp;
  double a2 = -(e * x(1, 2) + prec_1) * tmp;
  double rnorm = 1.0 / std::sqrt(a1 * a1 + a2 * a2 + 1.0);
  a1 *= rnorm;
  a2 *= rnorm;
  Vec3 v1(a1, a2, rnorm);

  // e = e2;
  double tmp2 = 1.0 / (e2 * (x(0, 0) + x(1, 1)) + mx0011 - e2 * e2 + x01_squared);
  double a21 = -(e2 * x(0, 2) + prec_0) * tmp2;
  double a22 = -(e2 * x(1, 2) + prec_1) * tmp2;
  double rnorm2 = 1.0 / std::sqrt(a21 * a21 + a22 * a22 + 1.0);
  a21 *= rnorm2;
  a22 *= rnorm2;
  Vec3 v2(a21, a22, rnorm2);

  // optionally remove axb from v1,v2
  // costly and makes a very small difference!
  // v1=(v1-v1.dot(v3)*v3);v1.normalize();
  // v2=(v2-v2.dot(v3)*v3);v2.normalize();
  // v2=(v2-v1.dot(v2)*v2);v2.normalize();
  E << v1(0), v2(0), v3(0),
      v1(1), v2(1), v3(1),
      v1(2), v2(2), v3(2);
};

/**
* @brief Compute the absolute pose of a camera using three 3D-to-2D correspondences.
*  Implementation of the paper "Lambda Twist: An Accurate Fast Robust Perspective Three Point (P3P) Solver". Persson, M. and Nordberg, K. ECCV 2018
*
* @authors Mikael Persson and Klas Nordberg
*
* @param[in] bearing_vectors 3x3 matrix with UNITARY feature vectors (each column is a vector)
* @param[in] X  3x3 matrix with corresponding 3D world points (each column is a point)
* @param[out] rotation_translation_solutions vector that will contain the solutions (up to 4 solutions)
*
* @return true if at least one solution is found, false if no solution was found
*
*/
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

  // double g = 0.0;

  //p3 is det(D2) so its definietly >0 or its a degenerate case
  //if (std::abs(p3) >= std::abs(p0) || true)
  //{
  p3 = 1.0 / p3;
  p2 *= p3;
  p1 *= p3;
  p0 *= p3;

  // get sharpest real root of above...
  double g = cubick(p2, p1, p0);
  //}
  // else
  //{

    // lower numerical performance
    //g = 1.0 / (cubick(p1 / p0, p2 / p0, p3 / p0));
  //}

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
      if (tau1 > 0.0)
      {
        double tau = tau1;
        double d = a23 / (tau * (b23 + tau) + 1.0);
        if(d > 0.0) {
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
      if (tau2 > 0.0)
      {
        double tau = tau2;
        double d = a23 / (tau * (b23 + tau) + 1.0);
        if(d > 0.0) {
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
  }

  { //-v
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
        if(d > 0.0) {
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
      if (tau2 > 0)
      {
        double tau = tau2;
        double d = a23 / (tau * (b23 + tau) + 1.0);
        if(d > 0.0) {
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

  return valid > 0;
}

void P3PSolver_Nordberg::Solve(
    const Mat &bearing_vectors,
    const Mat &X, // 3D points
    std::vector<Mat34> *models)
{
  assert(3 == bearing_vectors.rows());
  assert(3 == X.rows());
  assert(bearing_vectors.cols() == X.cols());
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

} // namespace euclidean_resection
} // namespace openMVG
