// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Sida Li, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole_Brown.hpp"
using namespace openMVG;
using namespace openMVG::cameras;

#include "testing/testing.h"
#include "openMVG/cameras/Camera_Unit_Test.inl"

TEST(Cameras_Brown, disto_undisto_T2) {

  const Pinhole_Intrinsic_Brown_T2 cam(1000, 1000, 1000, 500, 500,
    // K1, K2, K3, T1, T2
    -0.054, 0.014, 0.006, 0.001, -0.001);

  Test_camera(cam);
}

TEST(Cameras_Brown, InvertNonContractiveDistortion) {
  // Fixed-point subtraction oscillates at (1, 1) for k1 = 0.2,
  // although the radial mapping is monotonic and has a unique inverse.
  for (const auto & coefficients : std::vector<std::vector<double>>{
    {0.2, 0, 0, 0, 0}, {0.2, 0.03, 0.004, 0.01, -0.02}})
  {
    const Pinhole_Intrinsic_Brown_T2 cam(2000, 2000, 500, 1000, 1000,
      coefficients[0], coefficients[1], coefficients[2],
      coefficients[3], coefficients[4]);
    for (int y = -10; y <= 10; ++y)
    {
      for (int x = -10; x <= 10; ++x)
      {
        const Vec2 point(x / 10.0, y / 10.0);
        const Vec2 distorted = cam.add_disto(point);
        const Vec2 recovered = cam.remove_disto(distorted);
        EXPECT_TRUE(recovered.allFinite());
        EXPECT_MATRIX_NEAR(point, recovered, 1e-9);
        EXPECT_MATRIX_NEAR(distorted, cam.add_disto(recovered), 1e-9);
        EXPECT_MATRIX_NEAR(cam.cam2ima(point),
          cam.get_ud_pixel(cam.cam2ima(distorted)), 1e-6);
      }
    }
  }
}

TEST(Cameras_Brown, ZeroDistortionPreservesPoints) {
  const Pinhole_Intrinsic_Brown_T2 cam(2000, 2000, 500, 1000, 1000);
  const Vec2 point(1.2, -0.8);
  EXPECT_MATRIX_NEAR(point, cam.remove_disto(point), 0.0);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
