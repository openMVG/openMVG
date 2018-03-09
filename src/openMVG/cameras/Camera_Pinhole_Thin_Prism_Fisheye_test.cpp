// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong,Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole_Thin_Prism_Fisheye.hpp"
using namespace openMVG;
using namespace openMVG::cameras;

#include "testing/testing.h"
#include "openMVG/cameras/Camera_Unit_Test.inl"

TEST(Cameras_Fisheye, disto_undisto_Fisheye) {

  const Pinhole_Intrinsic_Thin_Prism_Fisheye cam(1000, 1000, 1000, 500, 500,
                                      0.21047,0.21102,-5.36231e-06,0.00051541,-0.158023,0.406856,-8.46499e-05,0.000861313); // K1, K2, p1, p2, K3, K4, sx1, sx2

  Test_camera(cam);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
