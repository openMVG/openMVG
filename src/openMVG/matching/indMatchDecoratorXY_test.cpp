// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"

#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::matching;

TEST(IndMatchDecorator, Empty)
{
  const std::vector<features::PointFeature> features;
  matching::IndMatches matches;
  const IndMatchDecorator<float> decorator(matches, features, features);
  EXPECT_FALSE(decorator.getDeduplicated(matches));
}

TEST(IndMatchDecorator, NoDuplicates)
{
  const std::vector<features::PointFeature> left_features = {{0,1}, {1,2}};
  const std::vector<features::PointFeature> right_features = {{2,3}, {4,5}};
  const matching::IndMatches input_matches = {{0,0}, {1,1}};
  const IndMatchDecorator<float> decorator(input_matches, left_features, right_features);
  matching::IndMatches output_matches;
  EXPECT_FALSE(decorator.getDeduplicated(output_matches));
  EXPECT_EQ(2, output_matches.size());
}

TEST(IndMatchDecorator, RightFeaturesAreTheSame)
{
  const std::vector<features::PointFeature> left_features = {{1,1}, {1,1}};
  const std::vector<features::PointFeature> right_features = {{0,0}};
  const matching::IndMatches input_matches = {{0,0}, {1,0}};
  const IndMatchDecorator<float> decorator(input_matches, left_features, right_features);
  matching::IndMatches output_matches;
  EXPECT_TRUE(decorator.getDeduplicated(output_matches));
  EXPECT_EQ(1, output_matches.size());
}

TEST(IndMatchDecorator, LeftFeaturesAreTheSame)
{
  const std::vector<features::PointFeature> left_features = {{0,0}};
  const std::vector<features::PointFeature> right_features = {{1,1}, {1,0}};
  const matching::IndMatches input_matches = {{0,0}, {0,1}};
  const IndMatchDecorator<float> decorator(input_matches, left_features, right_features);
  matching::IndMatches output_matches;
  EXPECT_TRUE(decorator.getDeduplicated(output_matches));
  EXPECT_EQ(1, output_matches.size());
}

TEST(IndMatchDecorator, PartialDeletion)
{
  const std::vector<features::PointFeature> left_features = {{0,0}, {0,0}, {5,15}};
  const std::vector<features::PointFeature> right_features = {{1,1}, {2,2}};
  const matching::IndMatches input_matches = {{0,0}, {1,0}, {2,1}};
  const IndMatchDecorator<float> decorator(input_matches, left_features, right_features);
  matching::IndMatches output_matches;
  EXPECT_TRUE(decorator.getDeduplicated(output_matches));
  EXPECT_EQ(2, output_matches.size());
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
