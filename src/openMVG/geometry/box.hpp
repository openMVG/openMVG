// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_BOX_HPP
#define OPENMVG_GEOMETRY_BOX_HPP

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

  Box() = default;

  /// Define a Square bounding box from a center position and radius (center to corner distance)
  Box
  (
    const Vec3 & center,
    const double radius, // Distance from center to corner (like sphere radius)
    const Mat3 R = Mat3::Identity()
  );

  /// Define a Rectangular parallelepiped bounding box from min/max axis aligned positions
  Box
  (
    const double x_min,
    const double y_min,
    const double z_min,
    const double x_max,
    const double y_max,
    const double z_max
  );

  /**
  * @brief Export the Box as a PLY file
  * @return true if the file can be saved on disk
  */
  static bool export_Ply
  (
    const Box & box,
    const std::string & filename,
    const Vec3 color = Vec3(0, 255, 0)
  );
};

} // namespace geometry
} // namespace openMVG

#endif // OPENMVG_GEOMETRY_BOX_HPP
