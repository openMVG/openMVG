
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/robust_estimation/robust_estimator_lineKernel_test.hpp"
#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

#include "openMVG/numeric/numeric.h"

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
  std::vector<size_t> vec_inliers;
  Vec2 model = MaxConsensus(kernel,
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
  std::vector<size_t> vec_inliers;
  Vec2 model = MaxConsensus(kernel,
    ScorerEvaluator<LineKernel>(0.3), &vec_inliers);
  CHECK_EQUAL(5, vec_inliers.size());
}

// Test efficiency of MaxConsensus to find (inlier/outlier) in contamined data
TEST(MaxConsensusLineFitter, OneOutlier) {

  Mat2X xy(2, 6);
  // y = 2x + 1 with an outlier
  xy << 1, 2, 3, 4,  5, 100, // outlier!
        3, 5, 7, 9, 11, -123; // outlier!

  LineKernel kernel(xy);

  std::vector<size_t> vec_inliers;
  Vec2 model = MaxConsensus(kernel,
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
  std::vector<size_t> vec_inliers;
  Vec2 model = MaxConsensus(kernel,
    ScorerEvaluator<LineKernel>(0.3), &vec_inliers);
  CHECK_EQUAL(0, vec_inliers.size());
}

// From a GT model :
//  Compute a list of point that fit the model.
//  Add white noise to given amount of points in this list.
//  Check that the number of inliers and the model are correct.
TEST(MaxConsensusLineFitter, RealisticCase) {

  const int NbPoints = 30;
  const int inlierPourcentAmount = 30; //works with 40
  Mat2X xy(2, NbPoints);

  Vec2 GTModel; // y = 2x + 1
  GTModel <<  -2.0, 6.3;

  //-- Build the point list according the given model
  for(int i = 0; i < NbPoints; ++i)  {
    xy.col(i) << i, (double)i*GTModel[1] + GTModel[0];
  }

  //-- Add some noise (for the asked percentage amount)
  int nbPtToNoise = (int) NbPoints*inlierPourcentAmount/100.0;
  std::vector<size_t> vec_samples; // Fit with unique random index
  UniformSample(nbPtToNoise, NbPoints, &vec_samples);
  for(size_t i = 0; i <vec_samples.size(); ++i)
  {
    const size_t randomIndex = vec_samples[i];
    //Additive random noise
    xy.col(randomIndex) << xy.col(randomIndex)(0)+rand()%2-3,
                           xy.col(randomIndex)(1)+rand()%8-6;
  }

  LineKernel kernel(xy);
  std::vector<size_t> vec_inliers;
  Vec2 model = MaxConsensus(kernel,
    ScorerEvaluator<LineKernel>(0.3), &vec_inliers);
  CHECK_EQUAL(NbPoints-nbPtToNoise, vec_inliers.size());
  EXPECT_NEAR(-2.0, model[0], 1e-9);
  EXPECT_NEAR( 6.3, model[1], 1e-9);
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
