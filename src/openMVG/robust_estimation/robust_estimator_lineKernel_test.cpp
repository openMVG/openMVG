// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/robust_estimation/robust_estimator_lineKernel_test.hpp"

#include "testing/testing.h"

#include <numeric>
#include <vector>

using namespace openMVG;
using namespace openMVG::robust;

// Since the line fitter isn't so simple, test it in isolation.
TEST(LineFitter, ItWorks) {

  Mat2X xy(2, 5);
  // y = 2x + 1
  xy << 1, 2, 3, 4,  5,
        3, 5, 7, 9, 11;
  std::vector<Vec2> models;
  LineKernel kernel(xy);
  std::vector<uint32_t> samples(static_cast<size_t>(xy.cols()));
  std::iota(samples.begin(), samples.end(), 0);
  kernel.Fit(samples, &models);
  CHECK_EQUAL(1, models.size());
  EXPECT_NEAR(2.0, models[0][1], 1e-9);
  EXPECT_NEAR(1.0, models[0][0], 1e-9);
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
