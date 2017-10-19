// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_CONVEX_HULL_HPP
#define OPENMVG_GEOMETRY_CONVEX_HULL_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG
{
namespace geometry
{

  /// Define an alias in order to represent a Polygon as an array of 2d floating point coordinates.
  /// Polygon must be defined as CounterClockWise polygon.
  /// In order to check if your polygon is CounterClockWise you can use the Area Function.
  ///  - A positive area means CounterClockWise polygon.
  ///  - A negative area means a Clockwise polygon.
  using Polygon2d = std::vector<Eigen::Vector2d>;

///  @brief Test if a point is inside a convex polygon.
///         (Check if the point is on the same side of all the contiguous segment defined by the polygon points).
///
///  @param[in] p               Point that we want to test the polygon inclusion.
///  @param[in] convex_polygon  The convex polygon.
///
///  @return Returns true if the point is inside the polygon, else returns false.
bool IsIn(const Eigen::Vector2d & p, const Polygon2d & polygon);

///  @brief Compute the area of a convex polygon using the 'determinant method'
///         http://www.mathwords.com/a/area_convex_polygon.htm
///
///  @param[in]  polygon The convex polygon.
///  @param[out] are     The polygon area. It returns an empty area (0) if a polygon with no point or one point is used.
///
///  @return True if the area can be computed
///
bool ConvexPolygonArea(const Polygon2d & polygon, double & area);

///  @brief Compute the convex hull of a set of 2D points.
///         Note: the last point in the returned list is the same as the first one.
///         Implements Andrew's monotone chain algorithm. O(n log n) complexity.
///         https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
///
///  @param[in]  point_cloud  The point cloud for which the convex hull must be computed.
///  @param[out] convex_hull  The output convex hull.
///
///  @return True if the function was able to complete is job.
bool ComputeConvexHull(Polygon2d point_cloud, Polygon2d & convex_hull);


} // namespace geometry
} // namespace openMVG

#endif // OPENMVG_GEOMETRY_CONVEX_HULL_HPP
