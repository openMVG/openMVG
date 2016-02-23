// Copyright (c) 2015 Romain Janvier <romain.janvier~AT~univ-orleans.fr>

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/cameras.hpp"
using namespace openMVG;
using namespace openMVG::cameras;

#include "testing/testing.h"

//-----------------
// Test summary:
//-----------------
// - Create a Pinhole_Intrinsic_Fisheye
// - Generate random point inside the image domain
// - Add and remove distortion and assert we found back the generated point
// - Check the last point in the camera & image domain
// - Assert that the tested distortion is not null (in order to ensure validity of the test)
//-----------------
TEST(Cameras_Fisheye, disto_undisto_Fisheye) {

const Pinhole_Intrinsic_Fisheye cam(1000, 1000, 1000, 500, 500,
                                    -0.054, 0.014, 0.006, 0.011); // K1, K2, K3, K4

  const double epsilon = 1e-4;
  for (int i = 0; i < 10; ++i)
  {
    // generate random point inside the image domain (last random to avoid 0,0)
    const Vec2 ptImage = (Vec2::Random() * 800./2.) + Vec2(500,500) + Vec2::Random();
    const Vec2 ptCamera = cam.ima2cam(ptImage);
    // Check that adding and removing distortion allow to recover the provided point
    EXPECT_MATRIX_NEAR( ptCamera, cam.remove_disto(cam.add_disto(ptCamera)), epsilon);
    EXPECT_MATRIX_NEAR( ptImage, cam.cam2ima(cam.remove_disto(cam.add_disto(ptCamera))), epsilon);

    // Assert that distortion field is not null and it has moved the initial provided point
    EXPECT_FALSE( (cam.add_disto(ptCamera) == cam.remove_disto(cam.add_disto(ptCamera))) ) ;
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */