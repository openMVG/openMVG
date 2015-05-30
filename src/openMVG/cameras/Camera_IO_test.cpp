// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_IO.hpp"
using namespace openMVG;
using namespace openMVG::cameras;

#include "testing/testing.h"

using std::string;

TEST(Camera_IO, PinholeSaveRead) {

  const Mat3 R = Mat3
    (Eigen::AngleAxisd(rand(), Vec3::UnitX())
    * Eigen::AngleAxisd(rand(), Vec3::UnitY())
    * Eigen::AngleAxisd(rand(), Vec3::UnitZ()));
  const PinholeCamera camGT( Mat3::Identity(), R, Vec3(0,1,2));

  EXPECT_TRUE( save( "pinholeCam.bin", camGT));
  EXPECT_FALSE( save( "pinholeCam.txt", camGT)); // extension must be .bin

  PinholeCamera cam;
  EXPECT_TRUE( load( "pinholeCam.bin", cam));
  EXPECT_MATRIX_NEAR(camGT._P, cam._P, 1e-3);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
