// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/numeric/numeric.h"

namespace openMVG
{

// A pinhole camera with its associated rotation
// Used to sample the spherical image
struct PinholeCamera_R
{
  PinholeCamera_R
  (
    int focal,
    int w,
    int h,
    const Mat3 & R
  )
  : R_(R)
  {
    K_ <<
      focal, 0, w/2.0,
      0, focal, h/2.0,
      0, 0, 1;
  }

  Vec3 getLocalRay
  (
    double x,
    double y
  ) const
  {
    return (K_.inverse() * Vec3(x, y, 1.0)).normalized();
  }

  Vec3 getRay
  (
    double x,
    double y
  ) const
  {
    return R_ * getLocalRay(x, y);
  }

  /// Rotation matrix
  Mat3 R_;

  /// Intrinsic matrix
  Mat3 K_;
};

// Function to map 3D coordinates onto a 2D image according a spherical projection
class CsphericalMapping
{
public:

  static Vec2 Get2DPoint
  (
    const Vec3 & X,
    int wPano,
    int hPano
  )
  {
    const Vec3 polarcoord = Get3DPointPolar(X);

    const double phi   = polarcoord(0);
    const double theta = polarcoord(1);

    const double x = ((phi * wPano) / M_PI + wPano) / 2.0;  // between 0 and width
    const double y = theta * hPano / M_PI;                  // between 0 and height

    return Vec2(x, y);
  }

  static Vec3 Get3DPointPolar
  (
    const Vec3 & pos_3d
  )
  {
    const double
      x = pos_3d(0),
      y = pos_3d(1),
      z = pos_3d(2),
      theta = atan2(y, sqrt(Square(x) + Square(z))),
      phi = atan2(x, z);
    return Vec3 (phi, theta + M_PI/2.0, 1.0);
  }
};

// Function to map 3D coordinates onto a 2D image according a cylindrical projection
class CcylindricalMapping
{
public:
  static Vec2 Get2DPoint
  (
    const Vec3  & X,
    int wPano,
    int hPano
  )
  {
    //-- Project the vertex on polar coordinates
    const Vec3 polarcoord = CsphericalMapping::Get3DPointPolar(X);

    const double
      phi = polarcoord(0),
      theta = polarcoord(1),
      x = ((phi * wPano) / M_PI + wPano) / 2.0, // between 0 and width
      y = (tan(theta - M_PI/2.0)) * (double)(hPano / 2.0 / M_PI) + hPano / 2.0;

    return Vec2(x, y);
  }
};

/// Compute a rectilinear camera focal from an angular FoV
double focalFromPinholeHeight
(
  int h,
  double thetaMax = D2R(60) // Camera FoV
)
{
  float f = 1.f;
  while ( thetaMax < atan2( h / (2 * f) , 1))
  {
    f++;
  }
  return f;
}

} // namespace openMVG
