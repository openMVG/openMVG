// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Romain Janvier <romain.janvier~AT~univ-orleans.fr>

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole_Fisheye.hpp"
#include <array>
using namespace openMVG;
using namespace openMVG::cameras;

#include "testing/testing.h"
#include "openMVG/cameras/Camera_Unit_Test.inl"

TEST(Cameras_Fisheye, disto_undisto_Fisheye) {

  const Pinhole_Intrinsic_Fisheye cam(1000, 1000, 1000, 500, 500,
                                      -0.054, 0.014, 0.006, 0.011); // K1, K2, K3, K4

  Test_camera(cam);
}

TEST(Cameras_Fisheye, StrongDistortionPixelRoundTrip) {
  const std::array<std::array<double, 4>, 7> coefficients = {{
    {{0, 0, 0, 0}}, {{0.1, 0, 0, 0}}, {{0.2, 0, 0, 0}},
    {{0.5, 0, 0, 0}}, {{-0.054, 0.014, 0.006, 0.011}},
    {{0.5, -0.1, 0, 0}}, {{0, 0, 0, 0.1}}
  }};
  for (const auto & k : coefficients) {
    const Pinhole_Intrinsic_Fisheye cam(20000, 20000, 500, 10000, 10000,
                                       k[0], k[1], k[2], k[3]);
    for (const Vec2 & point : {Vec2(0, 0), Vec2(1, 0), Vec2(-1, 0),
                              Vec2(2, 0), Vec2(0, -2), Vec2(1, 1),
                              Vec2(3, 2), Vec2(5, 0), Vec2(0, 10)}) {
      const Vec2 pixel = cam.cam2ima(point);
      const Vec2 distorted = cam.get_d_pixel(pixel);
      EXPECT_MATRIX_NEAR(pixel, cam.get_ud_pixel(distorted), 1e-7);
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
