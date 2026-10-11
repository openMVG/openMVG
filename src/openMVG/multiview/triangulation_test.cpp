// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Spherical.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/numeric/numeric.h"

#include "testing/testing.h"

using namespace openMVG;

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

TEST(Triangulation, TriangulateDLTSphericalCheirality)
{
  const cameras::Intrinsic_Spherical camera(1000, 500);
  const Mat3 R0 = RotationAroundY(0.2) * RotationAroundX(-0.1);
  const Vec3 t0(0.25, -0.1, 0.15);
  const Vec3 t1(-0.75, 0.15, -0.1);
  const Mat3 rotations[] = {
    RotationAroundY(-0.15) * RotationAroundZ(0.1),
    RotationAroundY(M_PI) * R0
  };
  const Vec3 points[] = {Vec3(0.4, 0.3, 4.0), Vec3(-0.3, 0.2, -4.0)};

  // Cover front, rear, and opposite camera hemispheres with nontrivial poses.
  for (const Mat3 &R1 : rotations)
  {
    for (const Vec3 &X_gt : points)
    {
      const Vec3 Y0 = R0 * X_gt + t0;
      const Vec3 Y1 = R1 * X_gt + t1;
      const Vec3 b0 = camera(camera.project(Y0)).col(0);
      const Vec3 b1 = camera(camera.project(Y1)).col(0);
      EXPECT_MATRIX_NEAR(Y0.normalized(), b0, 1e-12);
      EXPECT_MATRIX_NEAR(Y1.normalized(), b1, 1e-12);

      Vec3 X_estimated;
      EXPECT_EQ(true, TriangulateDLT(R0, t0, b0, R1, t1, b1, &X_estimated));
      EXPECT_NEAR(0.0, DistanceLInfinity(X_estimated, X_gt), 1e-8);
      EXPECT_EQ(true, Triangulate2View(
        R0, t0, b0, R1, t1, b1, X_estimated,
        ETriangulationMethod::DIRECT_LINEAR_TRANSFORM));
      EXPECT_NEAR(0.0, DistanceLInfinity(X_estimated, X_gt), 1e-8);

      // Positive scaling changes neither the rays nor their cheirality.
      EXPECT_EQ(true, TriangulateDLT(
        R0, t0, 3.0 * b0, R1, t1, 0.25 * b1, &X_estimated));
      EXPECT_NEAR(0.0, DistanceLInfinity(X_estimated, X_gt), 1e-8);

      // Antipodal spherical pixels have the same DLT projective constraints,
      // but the reconstructed point lies behind at least one observed ray.
      const Vec3 opposite0 = camera(camera.project(-Y0)).col(0);
      const Vec3 opposite1 = camera(camera.project(-Y1)).col(0);
      EXPECT_MATRIX_NEAR((-b0), opposite0, 1e-12);
      EXPECT_MATRIX_NEAR((-b1), opposite1, 1e-12);
      for (int antipodal = 1; antipodal < 4; ++antipodal)
      {
        EXPECT_EQ(false, Triangulate2View(
          R0, t0, (antipodal & 1) ? opposite0 : b0,
          R1, t1, (antipodal & 2) ? opposite1 : b1,
          X_estimated, ETriangulationMethod::DIRECT_LINEAR_TRANSFORM));
        EXPECT_NEAR(0.0, DistanceLInfinity(X_estimated, X_gt), 1e-8);
      }

      // Small image noise must preserve the positive ray depths.
      const Vec2 pixel0 = camera.project(Y0) + Vec2(0.05, -0.025);
      const Vec2 pixel1 = camera.project(Y1) + Vec2(-0.025, 0.05);
      const Vec3 noisy0 = camera(pixel0).col(0);
      const Vec3 noisy1 = camera(pixel1).col(0);
      EXPECT_TRUE(noisy0.dot(Y0) > 0.0);
      EXPECT_TRUE(noisy1.dot(Y1) > 0.0);
      EXPECT_EQ(true, Triangulate2View(
        R0, t0, noisy0, R1, t1, noisy1, X_estimated,
        ETriangulationMethod::DIRECT_LINEAR_TRANSFORM));
      EXPECT_NEAR(0.0, DistanceLInfinity(X_estimated, X_gt), 0.05);
      EXPECT_NEAR(0.0, (camera.project(R0 * X_estimated + t0) - pixel0).norm(), 0.5);
      EXPECT_NEAR(0.0, (camera.project(R1 * X_estimated + t1) - pixel1).norm(), 0.5);
    }
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
