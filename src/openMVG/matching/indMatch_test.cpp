// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "testing/testing.h"
#include "openMVG/matching/indMatch.hpp"

using namespace openMVG;
using namespace matching;

TEST(IndMatch, DuplicateRemoval_NoRemoval)
{
  std::vector<IndMatch> vec_indMatch;

  vec_indMatch.push_back(IndMatch(2,3)); // 0
  vec_indMatch.push_back(IndMatch(0,1)); // 1

  // Check no removal
  EXPECT_FALSE(IndMatch::getDeduplicated(vec_indMatch));

  // Check lexigraphical sorting
  // Due to lexigraphical sorting (0,1) must appears first
  EXPECT_EQ(IndMatch(0,1), vec_indMatch[0]);
  EXPECT_EQ(IndMatch(2,3), vec_indMatch[1]);
}

TEST(IndMatch, DuplicateRemoval_Simple)
{
  std::vector<IndMatch> vec_indMatch;

  vec_indMatch.push_back(IndMatch(0,1)); // 0
  vec_indMatch.push_back(IndMatch(0,1)); // 1: error with addition 0

  vec_indMatch.push_back(IndMatch(1,2)); // 2
  vec_indMatch.push_back(IndMatch(1,2)); // 3: error with addition 2

  EXPECT_TRUE(IndMatch::getDeduplicated(vec_indMatch));
  // Two matches must remain (line 0 and 2)
  EXPECT_EQ(2, vec_indMatch.size());
}

TEST(IndMatch, DuplicateRemoval)
{
  std::vector<IndMatch> vec_indMatch;

  vec_indMatch.push_back(IndMatch(0,1));
  vec_indMatch.push_back(IndMatch(0,1)); // Error defined before

  // Add some other matches
  vec_indMatch.push_back(IndMatch(0,2));
  vec_indMatch.push_back(IndMatch(1,1));
  vec_indMatch.push_back(IndMatch(2,3));
  vec_indMatch.push_back(IndMatch(3,3));

  EXPECT_TRUE(IndMatch::getDeduplicated(vec_indMatch));
  // Five matches must remain (one (0,1) must disappear)
  EXPECT_EQ(5, vec_indMatch.size());

  EXPECT_EQ(IndMatch(0,1), vec_indMatch[0]);
  EXPECT_EQ(IndMatch(0,2), vec_indMatch[1]);
  EXPECT_EQ(IndMatch(1,1), vec_indMatch[2]);
  EXPECT_EQ(IndMatch(2,3), vec_indMatch[3]);
  EXPECT_EQ(IndMatch(3,3), vec_indMatch[4]);
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
