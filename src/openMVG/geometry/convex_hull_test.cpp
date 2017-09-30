// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/convex_hull.hpp"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <limits>
#include <random>

using namespace openMVG;
using namespace openMVG::geometry;

static const double KEps = std::numeric_limits<double>::epsilon();

// Test corner case (an empty point cloud)
TEST(ConvexHull, Empty) {
  const Polygon2d point_cloud;

  Polygon2d convex_hull;
  EXPECT_TRUE(ComputeConvexHull(point_cloud, convex_hull));
  EXPECT_TRUE(convex_hull.empty());
}

TEST(ConvexHull, Line) {
  Polygon2d point_cloud;
  point_cloud.emplace_back(0., 0.);
  point_cloud.emplace_back(0., 2.);
  // Since a line is not a polygon (no area) we expect false for the IsIn test
  EXPECT_FALSE(IsIn(point_cloud[1], point_cloud));
  EXPECT_FALSE(IsIn(point_cloud[0], point_cloud));
}

// Test corner case (point cloud with one point)
TEST(ConvexHull, OnePoint) {
  const Polygon2d point_cloud(1, {0., 0.});
  // A one point polygon have no area
  double area = -1.0;
  EXPECT_TRUE(ConvexPolygonArea(point_cloud, area));
  EXPECT_NEAR(0.0, area, KEps);

  Polygon2d convex_hull;
  EXPECT_TRUE(ComputeConvexHull(point_cloud, convex_hull));
  EXPECT_EQ(1, convex_hull.size());  // Equal to the existing point

  // An empty polygon should have no area
  EXPECT_TRUE(ConvexPolygonArea(convex_hull, area));
  EXPECT_NEAR(0.0, area, KEps);
}

TEST(ConvexHull, Random_Point_ConvexHull) {
  // a. Generate some random point in a given area
  // b. Compute its convex hull
  // c. Check that all the points are inside the defined convex hull

  std::random_device rd;
  std::mt19937 gen(std::mt19937::default_seed);
  std::uniform_real_distribution<double> dis_x(100, 200), dis_y(300, 400);

  const int kIteration = 10;
  for (int i = 0; i < kIteration; ++i) {
    const int kNbPoints = 5;
    Polygon2d point_cloud(kNbPoints);
    for (auto & point : point_cloud) {
      point << dis_x(gen), dis_y(gen);
    }

    Polygon2d convex_hull;
    EXPECT_TRUE(ComputeConvexHull(point_cloud, convex_hull));
    EXPECT_TRUE(!convex_hull.empty());

    // Test that all the points are inside the convex_hull
    for (const auto & it_pt : point_cloud)
      EXPECT_TRUE(IsIn(it_pt, convex_hull));
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
