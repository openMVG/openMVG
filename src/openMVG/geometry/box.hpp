// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_BOX_HPP_
#define OPENMVG_GEOMETRY_BOX_HPP_

#include "openMVG/geometry/half_space_intersection.hpp"

namespace openMVG
{

/**
* @brief Namespace holding various classes and methods for geometry manipulation
*/
namespace geometry
{

using namespace openMVG::geometry::halfPlane;

/**
* @brief Define a Bounding box:
*  - using 6 Half Spaces
*  - This structure is used for testing box-HalfPlaneObject intersection
*/
struct Box : public HalfPlaneObject
{
  /// Points that define the bounding box
  Vec3 points[8];
  
  Box
  (
    const Vec3 & center,
    const double radius, // Distance from center to corner (like sphere radius)
    const Mat3 R = Mat3::Identity()
  )
  {
    const Vec3 directions[8] =
    {
      // Top face
      (R * Vec3(-1,-1,1).normalized()) * radius + center,
      (R * Vec3(-1,1,1).normalized())  * radius + center,
      (R * Vec3(1,1,1).normalized())   * radius + center,
      (R * Vec3(1,-1,1).normalized())  * radius + center,
      // Bottom face
      (R * Vec3(-1,-1,-1).normalized()) * radius + center,
      (R * Vec3(-1,1,-1).normalized())  * radius + center,
      (R * Vec3(1,1,-1).normalized())   * radius + center,
      (R * Vec3(1,-1,-1).normalized())  * radius + center
    };
    std::copy(directions, directions + 8, points);
    
    // Defines the 6 half_planes that defines the supporting volume
    planes.push_back( Half_plane_p( directions[0], directions[1], directions[3] ) ); // Top
    planes.push_back( Half_plane_p( directions[4], directions[7], directions[5] ) ); // Bottom
    // Define remaining supporting planes (box borders)
    planes.push_back( Half_plane_p( directions[0], directions[4], directions[1] ) );
    planes.push_back( Half_plane_p( directions[0], directions[3], directions[4] ) );
    planes.push_back( Half_plane_p( directions[3], directions[2], directions[7] ) );
    planes.push_back( Half_plane_p( directions[1], directions[5], directions[2] ) );
  }
};

} // namespace geometry
} // namespace openMVG

#endif // OPENMVG_GEOMETRY_BOX_HPP_
