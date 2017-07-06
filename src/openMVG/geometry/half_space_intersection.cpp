// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/geometry/half_space_intersection.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"

namespace openMVG{
namespace geometry{
namespace halfPlane{

Half_plane Half_plane_p
(
  const Vec3 & p,
  const Vec3 & q,
  const Vec3 & r
)
{
  return Half_plane::Through(r,q,p);
}

bool isNotEmpty
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
    cstraint.vec_bounds_.assign( cstraint.nbParams_,
      {std::numeric_limits<double>::lowest(),
      std::numeric_limits<double>::max()}); // [X,Y,Z] => -inf, +inf
    cstraint.bminimize_ = true;

    // Configure constraints
    const size_t nbConstraints = hplanes.size();
    cstraint.constraint_mat_.resize( nbConstraints, 3 );
    cstraint.vec_sign_.resize( nbConstraints );
    cstraint.constraint_objective_.resize( nbConstraints );

    // Fill the constrains (half-space equations)
    for ( unsigned char i = 0; i < hplanes.size(); ++i )
    {
      const Vec & half_plane_coeff = hplanes[i].coeffs();
      // add the half plane equation to the system
      cstraint.constraint_mat_.row( i ) <<
        half_plane_coeff( 0 ),
        half_plane_coeff( 1 ),
        half_plane_coeff( 2 );
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

bool HalfPlaneObject::intersect(const HalfPlaneObject & rhs) const
{
  // Concatenate the Half Planes and see if an intersection exists
  Half_planes vec_planes(planes.size() + rhs.planes.size());
  std::copy(&planes[0], &planes[0] + planes.size(), &vec_planes[0]);
  std::copy(&rhs.planes[0], &rhs.planes[0] + rhs.planes.size(), &vec_planes[planes.size()]);

  return halfPlane::isNotEmpty(vec_planes);
}

bool HalfPlaneObject::contains(const Vec3 & rhs) const
{
  unsigned int count_positive = 0;
  for (const Half_plane & hp : planes)
  {
    count_positive += (hp.signedDistance(rhs) > 0) ? 1 : 0;
  }
  return (count_positive == planes.size()) && !planes.empty();
}

bool intersect
(
  const std::vector<HalfPlaneObject> & hplanes
)
{
  // Compute the total number of half-planes
  std::size_t s = 0;
  std::for_each(hplanes.cbegin(), hplanes.cend(),
                [&s](const HalfPlaneObject & hp) { s += hp.planes.size(); });

  // Concatenate the half-planes and see if an intersection exists
  Half_planes vec_planes;
  vec_planes.reserve(s);
  for (const HalfPlaneObject & hp : hplanes)
  {
    vec_planes.insert(vec_planes.end(), hp.planes.cbegin(), hp.planes.cend());
  }

  return halfPlane::isNotEmpty(vec_planes);
}

} // namespace halfPlane
} // namespace geometry
} // namespace openMVG
