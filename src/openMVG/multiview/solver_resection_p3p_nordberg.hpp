// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Michael Persson
// Adapted to openMVG by Romain Janvier

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef OPENMVG_MULTIVIEW_RESECTION_P3P_PERSSON_HPP
#define OPENMVG_MULTIVIEW_RESECTION_P3P_PERSSON_HPP

#include "openMVG/multiview/two_view_kernel.hpp"

namespace openMVG
{
namespace euclidean_resection
{

struct P3PSolver_Persson
{
  enum
  {
    MINIMUM_SAMPLES = 3
  };
  enum
  {
    MAX_MODELS = 4
  };

  // Solve the absolute camera pose problem.
  // Use "Lambda Twist: An Accurate Fast Robust
  // Perspective Three Point (P3P) Solver
  // Persson, M.; Nordberg, K.
  static void Solve(
      const Mat &bearing_vectors,
      const Mat &X, // 3D points
      std::vector<Mat34> *models);

  // Compute the angular residual between the bearing vector and the 3d point projection vector
  static double Error(
      const Mat34 &P,
      const Vec3 &bearing_vector,
      const Vec3 &pt3D);
};

//-- Usable solver for robust estimation framework
using PoseResectionKernel_P3P_Persson =
    two_view::kernel::Kernel<
        P3PSolver_Persson, // Model estimator
        P3PSolver_Persson, // Error metric
        Mat34>;

static void gauss_newton_refineL(Vec3 &L,
                          const double & a12, const double & a13, const double & a23,
                          const double & b12, const double & b13, const double & b23)
{

  //TODO(RJ:) const expr makes it easier for the compiler to unroll
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
        double l1 = L(0);
        double l2 = L(1);
        double l3 = L(2);
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

static inline bool root2real(const double & b, const double & c, double & r1, double & r2){
    double v = b * b -4.0 * c;
    if(v < 0.0){
        r1 = r2 = 0.5 * b;
        return false;
    }
    double y = std::sqrt(v);
    if(b < 0.0){
        r1 = 0.5 * (-b +y);
        r2 = 0.5 * (-b -y);
    }else{
        r1= 2.0 * c / (-b + y);
        r2= 2.0 * c / (-b - y);
    }
    return true;
};

static double cubick(const double & b, const double & c, const double & d){
    // Choose an initial solution
    double r0;
    // not monotonic
    if (b * b >= 3.0 * c){
        // h has two stationary points, compute them
        // double t1 = t - std::sqrt(diff);
        double v = std::sqrt(b*b -3.0*c);
        double t1 = (-b - v)/(3.0);

        // Check if h(t1) > 0, in this case make a 2-order approx of h around t1
        double k = ((t1+b)*t1+c)*t1+d;

        if (k > 0.0) {
            // Find leftmost root of 0.5*(r0 -t1)^2*(6*t1+2*b) +  k = 0
            r0 = t1 - std::sqrt(-k/(3.0*t1 + b));
            // or use the linear comp too
            // r0 = t1 -
        } else {
            double t2 = (-b + v) / 3.0;
            k = ((t2 + b) * t2 + c) * t2 + d;
            // Find rightmost root of 0.5 * (r0 - t2)^2 * (6 * t2 +2 * b) + k1 = 0
            r0 = t2 + std::sqrt(-k/(3.0*t2 + b));
        }
    }
    else{
        // r0=1.0/(cubick_inv(c/d,b/d,1.0/d));
        // about half work...
        // if(std::abs((((r0+b)*r0+c)*r0+d))>1e-10)
        r0 = -b / 3.0;
        if(std::abs(((3.0 *r0 + 2.0 *b) * r0 + c)) < 1e-4) r0+=1;
        //else r0-=1;
        //double fx=(((r0+b)*r0+c)*r0+d); r0-=10; if(fx<0) r0+=20;

    }

    // Do ITER Newton-Raphson iterations
    // Break if position of root changes less than 1e-13
    // double starterr=std::abs(r0*(r0*(r0 + b) + c) + d);
    double fx,fpx;
    for (unsigned int cnt = 0; cnt < 50; ++cnt){ //TODO(RJ:) I have hardcoded the number of iter here
        fx = (((r0 + b) * r0 + c) * r0 + d);

        if((cnt < 7 || std::abs(fx) > 1e-13)  ){
            fpx = ((3.0 * r0 + 2.0 * b) * r0 + c);
            r0 -= fx / fpx;
        }
        else
            break;
    }
    return r0;
};

static void eigwithknown0(const Mat3& x, Mat3& E, Vec3 & L){
    // one eigenvalue is known to be 0.
    // the known one...
    L(2)=0.0;

    Vec3  v3(x(3) * x(7) - x(6) * x(4),
             x(6) * x(1) - x(7) * x(0),
             x(4) * x(0)- x(3) * x(1));

    v3.normalize();

    double x01_squared = x(0,1) * x(0,1);
    // get the two other...
    double b = -x(0,0) - x(1,1) - x(2,2);
    double c = -x01_squared - x(0,2) * x(0,2) - x(1,2) * x(1,2) +
                x(0,0) * (x(1,1) + x(2,2)) + x(1,1) * x(2,2);
    double e1, e2;
    //roots(poly(x))
    root2real(b,c,e1,e2);

    if(std::abs(e1) < std::abs(e2))
        std::swap(e1,e2);
    L(0) = e1;
    L(1) = e2;

    double mx0011 = -x(0,0) * x(1,1);
    double prec_0 = x(0,1) * x(1,2) - x(0,2) * x(1,1);
    double prec_1 = x(0,1) * x(0,2) - x(0,0) * x(1,2);

    double e = e1;
    double tmp = 1.0 / (e * (x(0,0) + x(1,1)) + mx0011 - e * e + x01_squared);
    double a1 = -(e * x(0,2) + prec_0) * tmp;
    double a2 = -(e * x(1,2) + prec_1) * tmp;
    double rnorm= 1.0 / std::sqrt(a1*a1 +a2*a2 + 1.0);
    a1 *= rnorm;
    a2 *= rnorm;
    Vec3 v1(a1,a2,rnorm);

    // e = e2;
    double tmp2 = 1.0 / (e2 * (x(0,0) + x(1,1)) + mx0011 - e2 * e2 + x01_squared);
    double a21 = -(e2 * x(0,2) + prec_0) * tmp2;
    double a22 = -(e2 * x (1,2) + prec_1) * tmp2;
    double rnorm2 = 1.0 / std::sqrt(a21*a21 +a22*a22 + 1.0);
    a21 *= rnorm2;
    a22 *= rnorm2;
    Vec3 v2(a21, a22, rnorm2);

    // optionally remove axb from v1,v2
    // costly and makes a very small difference!
    // v1=(v1-v1.dot(v3)*v3);v1.normalize();
    // v2=(v2-v2.dot(v3)*v3);v2.normalize();
    // v2=(v2-v1.dot(v2)*v2);v2.normalize();
    E << v1(0),v2(0),v3(0),  
         v1(1),v2(1),v3(1),  
         v1(2),v2(2),v3(2);
};

} // namespace euclidean_resection
} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_RESECTION_P3P_PERSSON_HPP
