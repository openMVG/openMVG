// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_HALF_SPACE_HPP_
#define OPENMVG_GEOMETRY_HALF_SPACE_HPP_

#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include <Eigen/Geometry>
#include <Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION( Eigen::Hyperplane<double, 3> )

namespace openMVG
{
namespace geometry
{
namespace halfPlane
{

/// Define the Half_plane equation (abcd coefficients)
typedef Eigen::Hyperplane<double, 3> Half_plane;

/// Define a collection of Half_plane
typedef std::vector<Half_plane> Half_planes;


/**
* @brief Define a plane passing through the points (p, q, r).
* @param p First point
* @param q Second point
* @param r Third point
* @return Plane formed by p,q,r
* @note The plane is oriented such that p, q and r are oriented in a positive sense (that is counterclockwise).
*/
inline Half_plane
Half_plane_p
(
  const Vec3 & p,
  const Vec3 & q,
  const Vec3 & r
)
{
  return Half_plane::Through(r,q,p);
}

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
inline bool
isNotEmpty
(
  const Half_planes & hplanes
)
{
  // Check if it exists a point on all positive side of the half plane thanks to a Linear Program formulation [1].
  // => If a point exists: there is a common subspace defined and so intersections.
  // The LP formulation consists in set the Half_plane as constraint and check if a point can fit the equations.

  using namespace openMVG;
  using namespace openMVG::linearProgramming;

  LP_Constraints cstraint;
  {
    cstraint.nbParams_ = 3; // {X,Y,Z}
    cstraint.vec_bounds_.resize( cstraint.nbParams_ );
    std::fill( cstraint.vec_bounds_.begin(), cstraint.vec_bounds_.end(),
               std::make_pair( ( double ) - 1e+30, ( double )1e+30 ) ); // [X,Y,Z] => -inf, +inf
    cstraint.bminimize_ = true;

    // Configure constraints
    const size_t nbConstraints = hplanes.size();
    cstraint.constraint_mat_.resize( nbConstraints, 3 );
    cstraint.vec_sign_.resize( nbConstraints );
    cstraint.constraint_objective_ = Vec( nbConstraints );

    // Fill the constrains (half-space equations)
    for ( unsigned char i = 0; i < hplanes.size(); ++i )
    {
      const Vec & half_plane_coeff = hplanes[i].coeffs();
      // add the half plane equation to the system
      cstraint.constraint_mat_.row( i ) =
        Vec3( half_plane_coeff( 0 ),
              half_plane_coeff( 1 ),
              half_plane_coeff( 2 ) );
      cstraint.vec_sign_[i] = LP_Constraints::LP_GREATER_OR_EQUAL;
      cstraint.constraint_objective_( i ) = - half_plane_coeff( 3 );
    }
  }

  // Solve in order to see if a point exists within the half spaces positive side?
  OSI_CLP_SolverWrapper solver( cstraint.nbParams_ );
  solver.setup( cstraint );
  const bool bIntersect = solver.solve(); // Status of the solver tell if there is an intersection or not
  return bIntersect;
}

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
  bool intersect(const HalfPlaneObject & rhs) const
  {
    // Concatenate the Half Planes and see if an intersection exists
    std::vector<Half_plane> vec_planes(planes.size() + rhs.planes.size());
    std::copy(&planes[0], &planes[0] + planes.size(), &vec_planes[0]);
    std::copy(&rhs.planes[0], &rhs.planes[0] + rhs.planes.size(), &vec_planes[planes.size()]);

    return halfPlane::isNotEmpty(vec_planes);
  }

  /**
  * @brief Test if a point is on the positive side of the HalfPlanes
  * @param rhs The 3D point to test
  * @retval true If The point is in the half plane defined 'volume'
  */
  bool contains(const Vec3 & rhs) const
  {
    unsigned int count_positive = 0;
    for (const Half_plane & hp : planes)
    {
      count_positive += (hp.signedDistance(rhs) > 0) ? 1 : 0;
    }
    return (count_positive == planes.size()) && !planes.empty();
  }
};

} // namespace halfPlane
} // namespace geometry
} // namespace openMVG

#endif // OPENMVG_GEOMETRY_HALF_SPACE_HPP_
