// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLECAMERA_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLECAMERA_HPP

#include "openMVG/multiview/projection.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {
namespace cameras {

/// Pinhole camera P = K[R|t], t = -RC
struct PinholeCamera
{
  PinholeCamera(
    const Mat3 & K = Mat3::Identity(),
    const Mat3 & R = Mat3::Identity(),
    const Vec3 & t = Vec3::Zero())
    : _K(K), _R(R), _t(t)
  {
    _C = -R.transpose() * t;
    P_From_KRt(_K, _R, _t, &_P);
  }

  explicit PinholeCamera(const Mat34 & P)
  : _P(P)
  {
    KRt_From_P(_P, &_K, &_R, &_t);
    _C = -_R.transpose() * _t;
  }

  /// Projection matrix P = K[R|t]
  Mat34 _P;

  /// Intrinsic parameter (Focal, principal point)
  Mat3 _K;

  /// Extrinsic Rotation
  Mat3 _R;

  /// Extrinsic translation
  Vec3 _t;

  /// Camera center
  Vec3 _C;

  /// Projection of a 3D point into the camera plane
  static Vec2 Project(const Mat34 & P, const Vec3 & pt3D)
  {
    return openMVG::Project(P, pt3D);
  }

  /// Projection of a 3D point into the camera plane (member function)
  Vec2 Project(const Vec3 & pt3D) const
  {
    return openMVG::Project(_P, pt3D);
  }

  /// Return the residual value to the given 2d point
  static double Residual(
    const Mat34 & P,
    const Vec3 & pt3D,
    const Vec2 & ref) {
    return (ref - openMVG::Project(P, pt3D)).norm();
  }

  /// Return the residual value to the given 2d point
  double Residual(const Vec3 & pt3D, const Vec2 & ref) const  {
    return (ref - openMVG::Project(_P, pt3D)).norm();
  }

  double ResidualSquared(const Vec3 & pt3D, const Vec2 & ref) const  {
    return (ref - openMVG::Project(_P, pt3D)).squaredNorm();
  }

  // Compute the depth of the X point. R*X[2]+t[2].
  double Depth(const Vec3 &X) const{
    return openMVG::Depth(_R, _t, X);
  }

  /// Return the angle (degree) between two pinhole point rays
  static double AngleBetweenRay(
    const PinholeCamera & cam1,
    const PinholeCamera & cam2,
    const Vec2 & x1, const Vec2 & x2)
  {
    // x = (u, v, 1.0)  // image coordinates
    // X = R.t() * K.inv() * x + C // Camera world point
    // getting the ray:
    // ray = X - C = R.t() * K.inv() * x
    Vec3 ray1 = (cam1._R.transpose() *
      (cam1._K.inverse() * Vec3(x1(0), x1(1), 1.))).normalized();
    Vec3 ray2 = (cam2._R.transpose() *
      (cam2._K.inverse() * Vec3(x2(0), x2(1), 1.))).normalized();
    double mag = ray1.norm() * ray2.norm();
    double dotAngle = ray1.dot(ray2);
    return R2D(acos(clamp(dotAngle/mag, -1.0 + 1.e-8, 1.0 - 1.e-8)));
  }

};

} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLECAMERA_HPP
