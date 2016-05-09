// Copyright (c) 2013,2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/half_space_intersection.hpp"
#include "openMVG/geometry/frustum.hpp"
#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"
#include <iostream>

using namespace openMVG;
using namespace openMVG::geometry::halfPlane;
using namespace std;

TEST(HALF_PLANE, ExistingSubspace) {

  std::vector<Half_plane> vec_hplanes;

  Vec3 a,b,c;
  a << 0,0,0;
  b << 1,0,0;
  c << 0,1,0;

  const Vec3 offset(0,0,2);
  vec_hplanes.push_back(Half_plane_p(a,b,c));
  vec_hplanes.push_back(Half_plane_p(a+offset,b+offset,c+offset));

  //    /\
  // ___|____ z = 2
  //
  //    /\
  // ___|____ z = 0

  EXPECT_TRUE( isNotEmpty(vec_hplanes) );
}

TEST(HALF_PLANE, EmptyIntersection) {

  std::vector<Half_plane> vec_hplanes;

  Vec3 a,b,c;
  a << 0,0,0;
  b << 1,0,0;
  c << 0,1,0;

  const Vec3 offset(0,0,2);
  vec_hplanes.push_back(Half_plane_p(a,b,c));
  vec_hplanes.push_back(Half_plane_p(a+offset,b+offset,c+offset));
  vec_hplanes[1].normal() *= -1; //invert the side of the half plane

  //    /\
  // ___|____ z = 0
  //
  //
  // _______ z = -2
  //    |
  //   \/

  EXPECT_FALSE( isNotEmpty(vec_hplanes) );
}

TEST(HALF_PLANE, Side)
{
  HalfPlaneObject half_planes_obj;

  Vec3 a,b,c;
  a << 0,0,0;
  b << 1,0,0;
  c << 0,1,0;

  const Vec3 offset(0,0,2);
  half_planes_obj.planes.push_back(Half_plane_p(a,b,c));
  half_planes_obj.planes.push_back(Half_plane_p(a+offset,b+offset,c+offset));

  //    /\
  // ___|____ z = 2
  //
  //    /\
  // ___|____ z = 0


  // Test with a point that is visible by the two half plane
  EXPECT_TRUE( half_planes_obj.contains(Vec3(0,0,3)) );

  // Test with a point that is visible by only one of the half plane
  EXPECT_FALSE( half_planes_obj.contains(Vec3(0,0,1)) );

  // Test with a point that is invisible by the two half plane
  EXPECT_FALSE( half_planes_obj.contains(Vec3(0,0,-1)) );
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
