// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/cameras.hpp"
using namespace openMVG;
using namespace openMVG::cameras;

#include "testing/testing.h"
#include "openMVG/cameras/Camera_Unit_Test.inl"

TEST(Cameras_Spherical, disto_undisto) {

  const Intrinsic_Spherical cam(2000, 1000);

  Test_camera(cam);
}

TEST(Cameras_Spherical, ResidualUsesShortestLongitudeDistance) {
  const Intrinsic_Spherical cam(2000, 1000);
  const IntrinsicBase & intrinsic = cam;
  EXPECT_MATRIX_NEAR(Vec2(2, 3),
    intrinsic.residual(cam(Vec2(1999, 500)), Vec2(1, 503)), 1e-10);
  EXPECT_MATRIX_NEAR(Vec2(-2, -3),
    intrinsic.residual(cam(Vec2(1, 500)), Vec2(1999, 497)), 1e-10);
  EXPECT_MATRIX_NEAR(Vec2(2, 3),
    intrinsic.residual(cam(Vec2(999, 500)), Vec2(1001, 503)), 1e-10);

  // Pinhole image boundaries are not periodic.
  const Pinhole_Intrinsic pinhole(2000, 1000, 1000, 1000, 500);
  EXPECT_MATRIX_NEAR(Vec2(-1998, 3),
    pinhole.residual(pinhole(Vec2(1999, 500)), Vec2(1, 503)), 1e-10);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
