// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_FRUSTUM_HPP_
#define OPENMVG_GEOMETRY_FRUSTUM_HPP_

#include "openMVG/geometry/half_space_intersection.hpp"

namespace openMVG {
namespace geometry {

using namespace openMVG::geometry::halfPlane;

/// Define a camera Frustum:
///  - infinite Frustum (4 Half Spaces) (a pyramid)
///  - truncated Frustum (6 Half Spaces) (a truncated pyramid)
///  - This structure is used for testing frustum intersection (see if two cam can share visual content)
struct Frustum
{
  Vec3 cones[5]; // camera centre and the 4 points that define the image plane
  Half_planes planes; // Define infinite frustum planes + 2 optional Near and Far Half Space

  // Build a frustum from the image size, camera intrinsic and pose
  Frustum(const int w, const int h, const Mat3 & K, const Mat3 & R, const Vec3 & C)
  {
    const Mat3 Kinv = K.inverse();
    const Mat3 Rt = R.transpose();

    // Definition of the frustum with the supporting points
    cones[0] = C;
    cones[1] = Rt * ((Kinv * Vec3(0,0,1.0))) + C;
    cones[2] = Rt * ((Kinv * Vec3(w,0,1.0))) + C;
    cones[3] = Rt * ((Kinv * Vec3(w,h,1.0))) + C;
    cones[4] = Rt * ((Kinv * Vec3(0,h,1.0))) + C;

    // Definition of the supporting planes
    planes.push_back( Half_plane_p(cones[0], cones[4], cones[1]) );
    planes.push_back( Half_plane_p(cones[0], cones[1], cones[2]) );
    planes.push_back( Half_plane_p(cones[0], cones[2], cones[3]) );
    planes.push_back( Half_plane_p(cones[0], cones[3], cones[4]) );
  }

  Frustum(const int w, const int h, const Mat3 & K, const Mat3 & R, const Vec3 & C, const double zNear, const double zFar)
  {
    *this = Frustum(w,h,K,R,C);
    assert(zFar > zNear);

    // Add Znear and ZFar half plane using the cam looking direction
    const Vec3 camLookDirection_n = R.row(2).normalized();
    const double d_near = - zNear - camLookDirection_n.dot(C);
    planes.push_back( Half_plane(camLookDirection_n, d_near) );

    const double d_Far = zFar + camLookDirection_n.dot(C);
    planes.push_back( Half_plane(-camLookDirection_n, d_Far) );
  }

  /// Test if two frustums intersect or not
  bool intersect(const Frustum & f) const
  {
    // Concatenate the Half Planes and see if an intersection exists
    std::vector<Half_plane> vec_planes(planes.size() + f.planes.size());
    std::copy(&planes[0], &planes[0]+planes.size(), &vec_planes[0]);
    std::copy(&f.planes[0], &f.planes[0]+f.planes.size(), &vec_planes[planes.size()]);

    return halfPlane::isNotEmpty(vec_planes);
  }

  /// Return true if the Frustum is an infinite one
  bool isInfinite() const
  {
    return planes.size() == 4;
  }

  /// Return true if the Frustum is truncated
  bool isTruncated() const
  {
    return planes.size() == 6;
  }

}; // struct Frustum

} // namespace geometry
} // namespace openMVG

#endif // OPENMVG_GEOMETRY_FRUSTUM_HPP_
