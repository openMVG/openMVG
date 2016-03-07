
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/robust_estimation/rand_sampling.hpp"
#include "testing/testing.h"
#include <set>

using namespace openMVG;
using namespace openMVG::robust;

// Assert that each time exactly N random number are picked (no repetition)
TEST(UniformSampleTest, NoRepetions) {

  std::vector<size_t> samples;
  for (size_t total = 1; total < 500; total *= 2) { //Size of the data set
    for (size_t num_samples = 1; num_samples <= total; num_samples *= 2) { //Size of the consensus set
      UniformSample(num_samples, total, &samples);
      std::set<size_t> myset;
      for (size_t i = 0; i < num_samples; ++i) {
        myset.insert(samples[i]);
      }
      CHECK_EQUAL(num_samples, myset.size());
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
