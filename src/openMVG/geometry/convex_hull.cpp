// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/convex_hull.hpp"

#include <algorithm>

namespace openMVG
{
namespace geometry
{

namespace util
{
///  @brief Compute the 2D cross product of OA and OB vectors,
///  i.e. z-component of their 3D cross product.
///
///  @param[in] O The origin point.
///  @param[in] A First point.
///  @param[in] B Second point.
///
///  @return Returns a positive value, if OAB makes a counter-clockwise turn,
///          negative for clockwise turn, and zero if the points are collinear.
double CrossProductSign
(
  const Eigen::Vector2d &O,
  const Eigen::Vector2d &A,
  const Eigen::Vector2d &B
)
{
  return (A.x() - O.x()) * (B.y() - O.y()) - (A.y() - O.y()) * (B.x() - O.x());
}

} // namespace util

bool ConvexPolygonArea
(
  const Polygon2d & polygon,
  double & area
)
{
  area = 0.0;

  if (!polygon.empty())
  {
    if (polygon.front() != polygon.back())
      return false;
  }
  else
  {
    return false;
  }

  for (int i = 0; i < polygon.size() - 1 && polygon.size() >= 2; ++i)
    area += polygon[i].x() * polygon[i+1].y() - polygon[i].y() * polygon[i+1].x();
  area /= 2;

  return true;
}

bool IsIn
(
  const Eigen::Vector2d & p,
  const Polygon2d & polygon
)
{
  if (polygon.size() <= 2) // A point or line does not have any surface
    return false;

  for (int i = 0; i < polygon.size() - 1; ++i)
  {
    const auto & origin = polygon[i];
    const auto & edge = polygon[i+1];
    if (p != edge && p != origin)
    if (util::CrossProductSign(edge, origin, p) > 0)
      return false;
  }
  return true;
}

bool ComputeConvexHull
(
  Polygon2d point_cloud,
  Polygon2d & convex_hull
)
{
  convex_hull.clear();

  if (point_cloud.empty())
    return true;

  if (point_cloud.size() == 1)
  {
    convex_hull = point_cloud;
    return true;
  }

  // A. sort the point in lexicographic ordering (use x and then y)
  std::sort(point_cloud.begin(),
            point_cloud.end(),
            [](const Eigen::Vector2d & pa,
               const Eigen::Vector2d & pb) -> bool
            {
              return pa.x() < pb.x() || (pa.x() == pb.x() && pa.y() < pb.y());
            });

  convex_hull.resize(2 * point_cloud.size());

  int k = 0;
  // B. Build lower hull
  for (int i = 0; i < point_cloud.size(); ++i)
  {
    while (k >= 2 &&
           util::CrossProductSign(convex_hull[k-2],
                                  convex_hull[k-1],
                                  point_cloud[i]) <= 0)
    {
      --k;
    }
    convex_hull[k++] = point_cloud[i];
  }

  // C. Build upper hull
  for (int i = point_cloud.size() - 2, t = k+1; i >= 0; i--)
  {
    while (k >= t &&
           util::CrossProductSign(convex_hull[k-2],
                                  convex_hull[k-1],
                                  point_cloud[i]) <= 0)
    {
      --k;
    }
    convex_hull[k++] = point_cloud[i];
  }

  convex_hull.resize(k);
  return true;
}

} // namespace geometry
} // namespace openMVG
