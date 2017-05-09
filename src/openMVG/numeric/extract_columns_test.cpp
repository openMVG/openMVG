// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/numeric/extract_columns.hpp"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

using namespace openMVG;

TEST(Numeric, ExtractColumns) {
  Mat2X A(2, 5);
  A << 1, 2, 3, 4, 5,
       6, 7, 8, 9, 10;
  Vec2i columns; columns << 0, 2;
  Mat2X extracted = ExtractColumns(A, columns);
  EXPECT_NEAR(1, extracted(0,0), 1e-15);
  EXPECT_NEAR(3, extracted(0,1), 1e-15);
  EXPECT_NEAR(6, extracted(1,0), 1e-15);
  EXPECT_NEAR(8, extracted(1,1), 1e-15);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
