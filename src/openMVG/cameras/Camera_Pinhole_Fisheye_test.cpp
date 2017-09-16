// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Romain Janvier <romain.janvier~AT~univ-orleans.fr>

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole_Fisheye.hpp"
using namespace openMVG;
using namespace openMVG::cameras;

#include "testing/testing.h"
#include "openMVG/cameras/Camera_Unit_Test.inl"

TEST(Cameras_Fisheye, disto_undisto_Fisheye) {

  const Pinhole_Intrinsic_Fisheye cam(1000, 1000, 1000, 500, 500,
                                      -0.054, 0.014, 0.006, 0.011); // K1, K2, K3, K4

  Test_camera(cam);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
