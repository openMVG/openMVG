// Copyright (c) 2013,2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/half_space_intersection.hpp"
#include "openMVG/geometry/frustum.hpp"
#include "openMVG/geometry/box.hpp"

#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/multiview/projection.hpp"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"
#include <iostream>
#include <fstream>
#include <Eigen/Geometry>

using namespace openMVG;
using namespace openMVG::geometry;
using namespace openMVG::geometry::halfPlane;
using namespace std;

//--
// Box/Camera frustum intersection unit test
//--

TEST(box_point, intersection)
{
  const Box box(Vec3::Zero(), sqrt(1.));

  // Test with a point that is inside the defined volume
  EXPECT_TRUE( box.contains(Vec3(0,0,0)) );

  // Test with a point that is outside the defined volume
  EXPECT_FALSE( box.contains(Vec3(1,1,1)) );
}

TEST(box_frustum, intersection)
{
  const int focal = 1000;
  const int principal_Point = 500;
  //-- Setup a circular camera rig or "cardioid".
  const int iNviews = 4;
  const int iNbPoints = 6;
  const NViewDataSet d =
    NRealisticCamerasRing(
    iNviews, iNbPoints,
    nViewDatasetConfigurator(focal, focal, principal_Point, principal_Point, 5, 0));

  const Box box(Vec3::Zero(), sqrt(1.));
  {
    std::ostringstream os;
    os << "box.ply";
    Box::export_Ply(box, os.str());
  }

  // Test with infinite Frustum for each camera
  {
    std::vector<Frustum> vec_frustum;
    for (int i=0; i < iNviews; ++i)
    {
      const Frustum f (principal_Point*2, principal_Point*2, d._K[i], d._R[i], d.C[i]);
      EXPECT_TRUE(f.intersect(box));
      EXPECT_TRUE(box.intersect(f));

      std::ostringstream os;
      os << i << "frust.ply";
      Frustum::export_Ply(f, os.str());
    }
  }

  // Test with truncated frustum
  {
    // Build frustum with near and far plane defined by min/max depth per camera
    for (int i=0; i < iNviews; ++i)
    {
      double minDepth = std::numeric_limits<double>::max();
      double maxDepth = std::numeric_limits<double>::min();
      for (int j=0; j < iNbPoints; ++j)
      {
        const double depth = Depth(d._R[i], d._t[i], d.X.col(j));
        if (depth < minDepth)
          minDepth = depth;
        if (depth > maxDepth)
          maxDepth = depth;
      }
      const Frustum f(principal_Point*2, principal_Point*2,
          d._K[i], d._R[i], d.C[i], minDepth, maxDepth);

      EXPECT_TRUE(f.intersect(box));
      EXPECT_TRUE(box.intersect(f));
    }
  }
}

TEST(box_frustum, no_intersection)
{
  const int focal = 1000;
  const int principal_Point = 500;
  //-- Setup a circular camera rig or "cardioid".
  const int iNviews = 4;
  const int iNbPoints = 6;
  const NViewDataSet d =
    NRealisticCamerasRing(
      iNviews, iNbPoints,
      nViewDatasetConfigurator(focal, focal, principal_Point, principal_Point, 5, 0));

  // Put the box out of field of the camera
  // (since camera are Y up, we move the box along Y axis)
  const Vec3 position(0, -4, 0);
  const Box box(position, sqrt(1.));
  {
    std::ostringstream os;
    os << "box.ply";
    Box::export_Ply(box, os.str());
  }

  // Test with infinite Frustum for each camera
  {
    std::vector<Frustum> vec_frustum;
    for (int i = 0; i < iNviews; ++i)
    {
      const Frustum f(principal_Point * 2, principal_Point * 2, d._K[i], d._R[i], d.C[i]);
      EXPECT_FALSE(f.intersect(box));
      EXPECT_FALSE(box.intersect(f));

      std::ostringstream os;
      os << i << "frust.ply";
      Frustum::export_Ply(f, os.str());
    }
  }

  // Test with truncated frustum
  {
    // Build frustum with near and far plane defined by min/max depth per camera
    for (int i = 0; i < iNviews; ++i)
    {
      double minDepth = std::numeric_limits<double>::max();
      double maxDepth = std::numeric_limits<double>::min();
      for (int j = 0; j < iNbPoints; ++j)
      {
        const double depth = Depth(d._R[i], d._t[i], d.X.col(j));
        if (depth < minDepth)
          minDepth = depth;
        if (depth > maxDepth)
          maxDepth = depth;
      }
      const Frustum f(principal_Point * 2, principal_Point * 2,
        d._K[i], d._R[i], d.C[i], minDepth, maxDepth);

      EXPECT_FALSE(f.intersect(box));
      EXPECT_FALSE(box.intersect(f));
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
