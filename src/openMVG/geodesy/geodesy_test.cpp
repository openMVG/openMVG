// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geodesy/geodesy.hpp"
using namespace openMVG::geodesy;

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

TEST(GEODESY, LLA_TO_ECEF_TO_LLA)
{
  const double lat = 10, lon = 20, alt = 30;
  const openMVG::Vec3 ecef_xyz = lla_to_ecef(lat, lon, alt);
  const openMVG::Vec3 lla = ecef_to_lla(ecef_xyz(0), ecef_xyz(1), ecef_xyz(2));
  EXPECT_NEAR(lat, lla(0), 1e-6);
  EXPECT_NEAR(lon, lla(1), 1e-6);
  EXPECT_NEAR(alt, lla(2), 1e-6);
  EXPECT_TRUE(ecef_xyz(0) != lla(0));
  EXPECT_TRUE(ecef_xyz(1) != lla(1));
  EXPECT_TRUE(ecef_xyz(2) != lla(2));
}


TEST(GEODESY, LLA_TO_UTM)
{
  const double lat = 10, lon = 20, alt = 30;
  const openMVG::Vec3 utm = lla_to_utm(lat, lon, alt);
  EXPECT_NEAR(390399.227485604, utm(0), 1e-6);
  EXPECT_NEAR(1105578.58919181, utm(1), 1e-6);
  EXPECT_NEAR(alt, utm(2), 1e-6);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
