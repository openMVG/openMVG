// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_HALF_SPACE_INTERSECTION_HPP
#define OPENMVG_GEOMETRY_HALF_SPACE_INTERSECTION_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST(Eigen::Hyperplane<double, 3>)

namespace openMVG
{
namespace geometry
{
namespace halfPlane
{

/// Define the Half_plane equation (abcd coefficients)
using Half_plane = Eigen::Hyperplane<double, 3>;

/// Define a collection of Half_plane
using Half_planes = std::vector<Half_plane>;

/**
* @brief Define a plane passing through the points (p, q, r).
* @param p First point
* @param q Second point
* @param r Third point
* @return Plane formed by p,q,r
* @note The plane is oriented such that p, q and r are oriented in a positive sense (that is counterclockwise).
*/
Half_plane
Half_plane_p
(
  const Vec3 & p,
  const Vec3 & q,
  const Vec3 & r
);

/**
* @brief Test if half-planes defines a non empty volume
* @param hplanes A list of half planes
* @retval true If half-planes define a non empty volume
* @retval false If half-planes define an empty volume
* @note Searched volume is defined by intersection of positive sides of half-planes
* @ref :
 [1] Paper: Finding the intersection of n half-spaces in time O(n log n).
 Author: F.P. Preparata, D.E. Muller
 Published in: Theoretical Computer Science, Volume 8, Issue 1, Pages 45-55
 Year: 1979
 More: ISSN 0304-3975, http://dx.doi.org/10.1016/0304-3975(79)90055-0.
*/
bool
isNotEmpty
(
  const Half_planes & hplanes
);

/**
* @brief Define a Volume 'object' thanks to a series of half_plane:
*  - This structure is used for testing generic HalfPlaneObject/HalfPlaneObject intersection
*/
struct HalfPlaneObject
{
  Half_planes planes;

  /**
  * @brief Test if two defined 'volume' intersects
  * @param rhs Another HalfPlaneObject
  * @retval true If an none empty intersection exists
  * @retval false If there's no intersection
  */
  bool intersect(const HalfPlaneObject & rhs) const;

  /**
  * @brief Test if a point is on the positive side of the HalfPlanes
  * @param rhs The 3D point to test
  * @retval true If The point is in the half plane defined 'volume'
  */
  bool contains(const Vec3 & rhs) const;
};

/**
* @brief Test if multiple defined 'volume' intersect
* @param hplanes A vector of HalfPlaneObject
* @retval true If a non-empty intersection exists
* @retval false If there's no intersection
*/
bool
intersect
(
  const std::vector<HalfPlaneObject> & hplanes
);

} // namespace halfPlane
} // namespace geometry
} // namespace openMVG

#endif // OPENMVG_GEOMETRY_HALF_SPACE_INTERSECTION_HPP
