// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/features/feature.hpp"
#include "openMVG/features/feature_container.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"

#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::matching;

TEST(DeduplicateLeftOrRight, Empty)
{
  matching::IndMatches matches;
  EXPECT_FALSE(DeduplicateLeftAndRight({}, {}, matches));
  EXPECT_TRUE(matches.empty());
}

TEST(DeduplicateLeftOrRight, NoDuplicates)
{
  const features::PointFeatures left_features = {{0,1}, {1,2}};
  const features::PointFeatures right_features = {{2,3}, {4,5}};
  matching::IndMatches input_matches = {{0,0}, {1,1}};
  EXPECT_FALSE(DeduplicateLeftAndRight(left_features,
                                       right_features,
                                       input_matches));
  EXPECT_EQ(2, input_matches.size());
}

TEST(DeduplicateLeftAndRight, Duplicates)
{
  const features::PointFeatures left_features = {{1,1}, {1,1}};
  const features::PointFeatures right_features = {{0,0}, {0,0}};
  {
    matching::IndMatches input_matches = {{0,0}, {1,1}};
    EXPECT_TRUE(DeduplicateLeftAndRight(left_features,
                                        right_features,
                                        input_matches));
    EXPECT_EQ(1, input_matches.size());
  }

  {
    matching::IndMatches input_matches = {{0,0}, {0,0}};
    EXPECT_TRUE(DeduplicateLeftAndRight(left_features,
                                        right_features,
                                        input_matches));
    EXPECT_EQ(1, input_matches.size());
  }
}

TEST(indMatchDecoratorXY, DeduplicateLeftOrRight)
{
  const features::PointFeature x1(1,1), x2(2,2), x3(3,3);

  // Corner case (empty arrays)
  {
    matching::IndMatches matches = {};
    EXPECT_FALSE(
      DeduplicateLeftOrRight(
        {}, {}, matches));
  }
  // No deduplication check
  {
    matching::IndMatches matches = {{0,0}, {1,1}};
    EXPECT_FALSE(
      DeduplicateLeftOrRight(
        {x1, x2}, {x1, x3}, matches)
    );
  }

  // Deduplication check
  {
    matching::IndMatches matches = {{0,0}, {1,0}};
    EXPECT_TRUE(
      DeduplicateLeftOrRight(
        {x1, x2}, {x1}, matches));
  }
  {
    matching::IndMatches matches = {{0,0}, {1,2}};
    EXPECT_TRUE(
      DeduplicateLeftOrRight(
        {x1, x2}, {x1, x2, x1}, matches));
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
