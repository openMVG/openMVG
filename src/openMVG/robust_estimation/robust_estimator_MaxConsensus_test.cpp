// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/robust_estimation/robust_estimator_lineKernel_test.hpp"
#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::robust;

// Test without outlier
TEST(MaxConsensusLineFitter, OutlierFree) {

  Mat2X xy(2, 5);
  // y = 2x + 1
  xy << 1, 2, 3, 4,  5,
        3, 5, 7, 9, 11;

  // The base estimator
  LineKernel kernel(xy);

  // Check the best model that fit the most of the data
  //  in a robust framework (Max-consensus).
  std::vector<uint32_t> vec_inliers;
  const Vec2 model = MaxConsensus(kernel,
    ScorerEvaluator<LineKernel>(0.3), &vec_inliers);
  EXPECT_NEAR(2.0, model[1], 1e-9);
  EXPECT_NEAR(1.0, model[0], 1e-9);
  CHECK_EQUAL(5, vec_inliers.size());
}

// Test without getting back the model
TEST(MaxConsensusLineFitter, OutlierFree_DoNotGetBackModel) {

  Mat2X xy(2, 5);
  // y = 2x + 1
  xy << 1, 2, 3, 4,  5,
        3, 5, 7, 9, 11;

  LineKernel kernel(xy);
  std::vector<uint32_t> vec_inliers;
  const Vec2 model = MaxConsensus(kernel,
    ScorerEvaluator<LineKernel>(0.3), &vec_inliers);
  CHECK_EQUAL(5, vec_inliers.size());
  EXPECT_NEAR(2.0, model[1], 1e-9);
  EXPECT_NEAR(1.0, model[0], 1e-9);
}

// Test efficiency of MaxConsensus to find (inlier/outlier) in contamined data
TEST(MaxConsensusLineFitter, OneOutlier) {

  Mat2X xy(2, 6);
  // y = 2x + 1 with an outlier
  xy << 1, 2, 3, 4,  5, 100, // outlier!
        3, 5, 7, 9, 11, -123; // outlier!

  LineKernel kernel(xy);

  std::vector<uint32_t> vec_inliers;
  const Vec2 model = MaxConsensus(kernel,
    ScorerEvaluator<LineKernel>(0.3), &vec_inliers);
  EXPECT_NEAR(2.0, model[1], 1e-9);
  EXPECT_NEAR(1.0, model[0], 1e-9);
  CHECK_EQUAL(5, vec_inliers.size());
}

// Critical test:
// Test if the robust estimator do not return inlier if too few point
// was given for an estimation.
TEST(MaxConsensusLineFitter, TooFewPoints) {

  Mat2X xy(2, 1);
  xy << 1,
        3;   // y = 2x + 1 with x = 1
  LineKernel kernel(xy);
  std::vector<uint32_t> vec_inliers;
  MaxConsensus(kernel, ScorerEvaluator<LineKernel>(0.3), &vec_inliers);
  CHECK_EQUAL(0, vec_inliers.size());
}

// From a GT model :
//  Compute a list of point that fit the model.
//  Add white noise to given amount of points in this list.
//  Check that the number of inliers and the model are correct.
TEST(MaxConsensusLineFitter, RealisticCase) {

  constexpr int NbPoints = 30;
  constexpr double inlierRatio = 30.0 / 100.0; // works with 40
  Mat2X xy(2, NbPoints);

  Vec2 GTModel; // y = 2x + 1
  GTModel <<  -2.0, 6.3;

  //-- Build the point list according the given model
  for (int i = 0; i < NbPoints; ++i)  {
    xy.col(i) << i, static_cast<double>(i)*GTModel[1] + GTModel[0];
  }

  //-- Add some noise (for the asked percentage amount)
  constexpr auto nbPtToNoise = static_cast<uint32_t>(NbPoints*inlierRatio);
  std::vector<uint32_t> vec_samples; // fit with unique random index
  std::mt19937 random_generator(std::mt19937::default_seed);
  UniformSample(nbPtToNoise, NbPoints, random_generator, &vec_samples);

  std::uniform_int_distribution<int> d0(-3, 2);
  std::uniform_int_distribution<int> d1(-6, 8);
  for (const auto index : vec_samples) {
    // additive random noise
    xy.col(index) << xy.col(index)(0) + static_cast<double>(d0(random_generator)),
                     xy.col(index)(1) + static_cast<double>(d1(random_generator));
  }

  LineKernel kernel(xy);
  std::vector<uint32_t> vec_inliers;
  const Vec2 model = MaxConsensus(kernel,
    ScorerEvaluator<LineKernel>(0.3), &vec_inliers);
  CHECK_EQUAL(NbPoints-nbPtToNoise, vec_inliers.size());
  EXPECT_NEAR(-2.0, model[0], 1e-9);
  EXPECT_NEAR( 6.3, model[1], 1e-9);
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
