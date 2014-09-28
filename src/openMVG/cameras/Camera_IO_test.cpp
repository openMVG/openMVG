// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_IO.hpp"
using namespace openMVG;

#include "testing/testing.h"

using std::string;

TEST(Camera_IO, PinholeSaveRead) {

  Mat3 R = Mat3
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

TEST(Camera_IO, BrownPinholeSaveRead) {

  const BrownPinholeCamera camGT(
    1000, 500, 500,
    Mat3(Eigen::AngleAxisd(rand(), Vec3::UnitX())
      * Eigen::AngleAxisd(rand(), Vec3::UnitY())
      * Eigen::AngleAxisd(rand(), Vec3::UnitZ())),
    Vec3(0,1,2), 0.2, 0.1, 0.05);

  EXPECT_TRUE( save( "BrownPinholeCam.txt", camGT));
  EXPECT_FALSE( save( "BrownPinholeCam.bin", camGT)); // extension must be .bin

  BrownPinholeCamera cam;
  EXPECT_TRUE( load( "BrownPinholeCam.txt", cam));
  EXPECT_MATRIX_NEAR(camGT._P, cam._P, 1e-3);
  EXPECT_EQ( camGT._f, cam._f);
  EXPECT_EQ( camGT._ppx, cam._ppx);
  EXPECT_EQ( camGT._ppy, cam._ppy);
  EXPECT_EQ( camGT._k1, cam._k1);
  EXPECT_EQ( camGT._k2, cam._k2);
  EXPECT_EQ( camGT._k3, cam._k3);
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
