// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/numeric/numeric.h"

#include "testing/testing.h"

using namespace openMVG;
using namespace std;

TEST(Triangulation, TriangulateDLT) {

  const NViewDataSet d = NRealisticCamerasRing(2, 12);

  for (int i = 0; i < d._X.cols(); ++i) {
    const Vec3 X_gt = d._X.col(i);
    Vec3 X_estimated = Vec3::Random();
    const Vec2 x1 = d._x[0].col(i);
    const Vec2 x2 = d._x[1].col(i);
    TriangulateDLT(d.P(0), x1.homogeneous(), d.P(1), x2.homogeneous(), &X_estimated);
    EXPECT_NEAR(0, DistanceLInfinity(X_estimated, X_gt), 1e-8);
    X_estimated = Vec3::Random();
    TriangulateDLT(d._R[0], d._t[0], d._K[0].inverse() * x1.homogeneous(), d._R[1], d._t[1], d._K[1].inverse() * x2.homogeneous(), &X_estimated);
    EXPECT_NEAR(0, DistanceLInfinity(X_estimated, X_gt), 1e-8);
  }
}

TEST(Triangulation, TriangulateL1Angular)
{
  const NViewDataSet d = NRealisticCamerasRing(2, 12);
  for (int i = 0; i < d._X.cols(); ++i) {
    const Vec3 X_gt = d._X.col(i);
    Vec3 X_estimated;
    const Vec2 x1 = d._x[0].col(i);
    const Vec2 x2 = d._x[1].col(i);
    const bool cheir_ok = TriangulateL1Angular(d._R[0], d._t[0], d._K[0].inverse() * x1.homogeneous(), d._R[1], d._t[1], d._K[1].inverse() * x2.homogeneous(), &X_estimated);
    EXPECT_TRUE(cheir_ok);
    EXPECT_NEAR(0, DistanceLInfinity(X_estimated, X_gt), 1e-8);
  }
}

  TEST(Triangulation, TriangulateLInfinityAngular)
{
  const NViewDataSet d = NRealisticCamerasRing(2, 12);
  for (int i = 0; i < d._X.cols(); ++i) {
    const Vec3 X_gt = d._X.col(i);
    Vec3 X_estimated;
    const Vec2 x1 = d._x[0].col(i);
    const Vec2 x2 = d._x[1].col(i);
    const bool cheir_ok = TriangulateLInfinityAngular(d._R[0], d._t[0], d._K[0].inverse() * x1.homogeneous(), d._R[1], d._t[1], d._K[1].inverse() * x2.homogeneous(), &X_estimated);
    EXPECT_TRUE(cheir_ok);
    EXPECT_NEAR(0, DistanceLInfinity(X_estimated, X_gt), 1e-8);
  }
}

TEST(Triangulation, TriangulateIDWMidpoint)
{
  const NViewDataSet d = NRealisticCamerasRing(2, 12);
  for (int i = 0; i < d._X.cols(); ++i) {
    const Vec3 X_gt = d._X.col(i);
    Vec3 X_estimated;
    const Vec2 x1 = d._x[0].col(i);
    const Vec2 x2 = d._x[1].col(i);
    const bool cheir_ok = TriangulateIDWMidpoint(d._R[0], d._t[0], d._K[0].inverse() * x1.homogeneous(), d._R[1], d._t[1], d._K[1].inverse() * x2.homogeneous(), &X_estimated);
    EXPECT_TRUE(cheir_ok);
    EXPECT_NEAR(0, DistanceLInfinity(X_estimated, X_gt), 1e-8);
  }
}

TEST(Triangulation, Factoring)
{
  const NViewDataSet d = NRealisticCamerasRing(2, 12);
  const std::vector<ETriangulationMethod> triangulation_methods =
    {
      ETriangulationMethod::DIRECT_LINEAR_TRANSFORM,
      ETriangulationMethod::L1_ANGULAR,
      ETriangulationMethod::LINFINITY_ANGULAR,
      ETriangulationMethod::INVERSE_DEPTH_WEIGHTED_MIDPOINT
    };
  for (const auto& method_it : triangulation_methods)
  {
    for (int i = 0; i < d._X.cols(); ++i)
    {
      const Vec3 X_gt = d._X.col(i);
      Vec3 X_estimated = Vec3::Random();
      const Vec2 x1 = d._x[0].col(i);
      const Vec2 x2 = d._x[1].col(i);
      const bool cheir_ok =
        Triangulate2View(
          d._R[0], d._t[0], d._K[0].inverse() * x1.homogeneous(),
          d._R[1], d._t[1], d._K[1].inverse() * x2.homogeneous(),
          X_estimated,
          method_it);
      EXPECT_TRUE(cheir_ok);
      EXPECT_NEAR(0, DistanceLInfinity(X_estimated, X_gt), 1e-8);
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
