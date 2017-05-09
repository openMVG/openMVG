// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
using namespace openMVG;
using namespace openMVG::cameras;

#include "testing/testing.h"
#include "openMVG/cameras/Camera_Unit_Test.inl"

TEST(Cameras_Radial, disto_undisto_K1) {

  const Pinhole_Intrinsic cam(1000, 1000, 1000, 500, 500);

  Test_camera(cam);
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
