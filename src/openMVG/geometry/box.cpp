// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/box.hpp"

#include <fstream>

namespace openMVG{
namespace geometry{

Box::Box
(
  const Vec3 & center,
  const double radius, // Distance from center to corner (like sphere radius)
  const Mat3 R
)
{
  const Vec3 directions[8] =
  {
    // Top face
    (R * Vec3(-1,-1,1).normalized()) * radius + center,
    (R * Vec3(-1, 1,1).normalized()) * radius + center,
    (R * Vec3( 1, 1,1).normalized()) * radius + center,
    (R * Vec3( 1,-1,1).normalized()) * radius + center,
    // Bottom face
    (R * Vec3(-1,-1,-1).normalized()) * radius + center,
    (R * Vec3(-1, 1,-1).normalized()) * radius + center,
    (R * Vec3( 1, 1,-1).normalized()) * radius + center,
    (R * Vec3( 1,-1,-1).normalized()) * radius + center
  };
  std::copy(directions, directions + 8, points);

  // Defines the supporting volume thanks to 6 half_planes
  planes.push_back( Half_plane_p( points[0], points[1], points[3] ) ); // Top
  planes.push_back( Half_plane_p( points[4], points[7], points[5] ) ); // Bottom
  // Define remaining supporting planes (box borders)
  planes.push_back( Half_plane_p( points[0], points[4], points[1] ) );
  planes.push_back( Half_plane_p( points[0], points[3], points[4] ) );
  planes.push_back( Half_plane_p( points[3], points[2], points[7] ) );
  planes.push_back( Half_plane_p( points[1], points[5], points[2] ) );
}

Box::Box
(
  const double x_min,
  const double y_min,
  const double z_min,
  const double x_max,
  const double y_max,
  const double z_max
)
{
  const Vec3 points_[8] =
  {
    // Top face
    Vec3(x_min, y_min, z_max),
    Vec3(x_min, y_max, z_max),
    Vec3(x_max, y_max, z_max),
    Vec3(x_max, y_min, z_max),
    // Bottom face
    Vec3(x_min, y_min, z_min),
    Vec3(x_min, y_max, z_min),
    Vec3(x_max, y_max, z_min),
    Vec3(x_max, y_min, z_min),
  };
  std::copy(points_, points_ + 8, points);

  // Defines the supporting volume thanks to 6 half_planes
  planes.push_back( Half_plane_p( points[0], points[1], points[3] ) ); // Top
  planes.push_back( Half_plane_p( points[4], points[7], points[5] ) ); // Bottom
  // Define remaining supporting planes (box borders)
  planes.push_back( Half_plane_p( points[0], points[4], points[1] ) );
  planes.push_back( Half_plane_p( points[0], points[3], points[4] ) );
  planes.push_back( Half_plane_p( points[3], points[2], points[7] ) );
  planes.push_back( Half_plane_p( points[1], points[5], points[2] ) );
}

bool Box::export_Ply
(
  const Box & box,
  const std::string & filename,
  const Vec3 color
)
{
  std::ofstream of(filename.c_str());
  if (!of.is_open())
    return false;

  of << "ply" << '\n'
    << "format ascii 1.0" << '\n'
    << "element vertex 8\n"
    << "property float x" << '\n'
    << "property float y" << '\n'
    << "property float z" << '\n'
    << "property uchar red" << '\n'
    << "property uchar green" << '\n'
    << "property uchar blue" << '\n'
    << "element face 6 \n"
    << "property list uchar int vertex_index" << '\n'
    << "end_header" << '\n';

  // Export box points
  {
    for (int i = 0; i < 8; ++i)
      of << box.points[i].transpose() << " " << color.cast<int>() << '\n';
  }

  of
    // top & bottom planes
    << "3 0 1 3\n"
    << "3 4 7 5\n"
    // remaining planes
    << "3 0 4 1\n"
    << "3 0 3 4\n"
    << "3 3 2 7\n"
    << "3 1 5 2\n";

  of.flush();
  const bool bOk = of.good();
  of.close();
  return bOk;
}

} // namespace geometry
} // namespace openMVG
