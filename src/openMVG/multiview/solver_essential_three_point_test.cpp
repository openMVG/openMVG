// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romain JANVIER.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/solver_essential_three_point.hpp"
#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/numeric/numeric.h"

#include "testing/testing.h"

#include <iostream>

using namespace openMVG;

TEST(ThreePointsRelativePose, reprod) {
  Mat2X x1(2,3), x2(2,3);
  x1 << 0.023397,  0.0028923, 0.022274,
        -0.094774, -0.0434725, -0.047761;
  x2 << 0.0057636, -2.9155e-04, 0.021184,
        -0.0927769, -4.5890e-02, -0.044213;
  // Functional testing (only assert that some model are found)
  std::vector<Mat3> models;
  ThreePointsRelativePose(x1, x2, &models);
  EXPECT_EQ(2, models.size());
  std::cerr
    << models[0] << std::endl
    << "=============" << std::endl
    << models[1] << std::endl;
}




/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
