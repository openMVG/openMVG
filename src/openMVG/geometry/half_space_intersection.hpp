// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_HALF_SPACE_HPP_
#define OPENMVG_GEOMETRY_HALF_SPACE_HPP_

#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include <Eigen/Geometry>

namespace openMVG {
namespace geometry {
namespace halfPlane {

/// Define the Half_plane equation (abcd coefficients)
typedef Eigen::Hyperplane<double,3> Half_plane;
/// Define a collection of Half_plane
typedef std::vector<Half_plane> Half_planes;

// Define a plane passing through the points (p, q, r).
// The plane is oriented such that p, q and r are oriented in a positive sense (that is counterclockwise).
Half_plane Half_plane_p(const Vec3 & p, const Vec3 & q, const Vec3 & r)
{
  const Vec3 abc = (p-r).cross(q-r);
  const double d = - abc.dot(r);
  Half_plane hp;
  hp.coeffs() << abc(0),abc(1),abc(2),d;
  return hp;
}

// [1] Paper: Finding the intersection of n half-spaces in time O(n log n).
// Author: F.P. Preparata, D.E. Muller
// Published in: Theoretical Computer Science, Volume 8, Issue 1, Pages 45-55
// Year: 1979
// More: ISSN 0304-3975, http://dx.doi.org/10.1016/0304-3975(79)90055-0.


/// Return true if the half_planes define a not empty volume (an intersection exists)
static bool isNotEmpty(const Half_planes & hplanes)
{
  // Check if it exists a point on all positive side of the half plane thanks to a Linear Program formulation [1].
  // => If a point exists: there is a common subspace defined and so intersections.
  // The LP formulation consists in set the Half_plane as constraint and check if a point can fit the equations.

  using namespace openMVG;
  using namespace openMVG::linearProgramming;

  LP_Constraints cstraint;
  {
    cstraint._nbParams = 3; // {X,Y,Z}
    cstraint._vec_bounds.resize(cstraint._nbParams);
    std::fill(cstraint._vec_bounds.begin(),cstraint._vec_bounds.end(),
      std::make_pair((double)-1e+30, (double)1e+30)); // [X,Y,Z] => -inf, +inf
    cstraint._bminimize = true;

    // Configure constraints
    const size_t nbConstraints = hplanes.size();
    cstraint._constraintMat = Mat(nbConstraints,3);
    cstraint._vec_sign.resize(nbConstraints);
    cstraint._Cst_objective = Vec(nbConstraints);

    // Fill the constrains (half-space equations)
    for (unsigned char i= 0; i < hplanes.size(); ++i)
    {
      const Vec & half_plane_coeff = hplanes[i].coeffs();
      // add the half plane equation to the system
      cstraint._constraintMat.row(i) =
        Vec3(half_plane_coeff(0),
          half_plane_coeff(1),
          half_plane_coeff(2));
      cstraint._vec_sign[i] = LP_Constraints::LP_GREATER_OR_EQUAL;
      cstraint._Cst_objective(i) = - half_plane_coeff(3);
    }
  }

  // Solve in order to see if a point exists within the half spaces positive side?
  OSI_CLP_SolverWrapper solver(cstraint._nbParams);
  solver.setup(cstraint);
  const bool bIntersect = solver.solve(); // Status of the solver tell if there is an intersection or not
  return bIntersect;
}

} // namespace geometry
} // namespace openMVG
} // namespace halfPlane

#endif // OPENMVG_GEOMETRY_HALF_SPACE_HPP_
