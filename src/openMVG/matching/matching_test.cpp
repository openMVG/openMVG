// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "testing/testing.h"
#include "openMVG/numeric/numeric.h"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matcher_kdtree_flann.hpp"
#include "openMVG/matching/matcher_cascade_hashing.hpp"
#include <iostream>
using namespace std;

using namespace openMVG;
using namespace matching;

TEST(Matching, ArrayMatcherBruteForce_Simple_Dim1)
{
  const float array[] = {0, 1, 2, 3, 4};
  ArrayMatcherBruteForce<float> matcher;
  EXPECT_TRUE( matcher.Build(array, 5, 1) );

  const float query[] = {2};
  int nIndice = -1;
  float fDistance = -1.0f;
  EXPECT_TRUE( matcher.SearchNeighbour( query, &nIndice, &fDistance) );

  EXPECT_EQ( 2, nIndice); // index of the found nearest neighbor
  EXPECT_NEAR( 0.0f, fDistance, 1e-8); //distance
}

TEST(Matching, ArrayMatcherBruteForce_NN)
{
  const float array[] = {0, 1, 2, 5, 6};
  // no 3, because it involve the same dist as 1,1
  ArrayMatcherBruteForce<float> matcher;
  EXPECT_TRUE( matcher.Build(array, 5, 1) );

  const float query[] = {2};
  IndMatches vec_nIndice;
  vector<float> vec_fDistance;
  EXPECT_TRUE( matcher.SearchNeighbours(query,1, &vec_nIndice, &vec_fDistance, 5) );

  EXPECT_EQ( 5, vec_nIndice.size());
  EXPECT_EQ( 5, vec_fDistance.size());

  // Check distances:
  EXPECT_NEAR( vec_fDistance[0], Square(2.0f-2.0f), 1e-6);
  EXPECT_NEAR( vec_fDistance[1], Square(1.0f-2.0f), 1e-6);
  EXPECT_NEAR( vec_fDistance[2], Square(0.0f-2.0f), 1e-6);
  EXPECT_NEAR( vec_fDistance[3], Square(5.0f-2.0f), 1e-6);
  EXPECT_NEAR( vec_fDistance[4], Square(6.0f-2.0f), 1e-6);

  // Check indexes:
  EXPECT_EQ(IndMatch(0,2), vec_nIndice[0]);
  EXPECT_EQ(IndMatch(0,1), vec_nIndice[1]);
  EXPECT_EQ(IndMatch(0,0), vec_nIndice[2]);
  EXPECT_EQ(IndMatch(0,3), vec_nIndice[3]);
  EXPECT_EQ(IndMatch(0,4), vec_nIndice[4]);
}

TEST(Matching, ArrayMatcherBruteForce_Simple_Dim4)
{
  const float array[] = {
    0, 1, 2, 3,
    4, 5, 6, 7,
    8, 9, 10, 11};
  ArrayMatcherBruteForce<float> matcher;
  EXPECT_TRUE( matcher.Build(array, 3, 4) );

  const float query[] = {4, 5, 6, 7};
  int nIndice = -1;
  float fDistance = -1.0f;
  EXPECT_TRUE( matcher.SearchNeighbour( query, &nIndice, &fDistance) );

  EXPECT_EQ( 1, nIndice); // index of the found nearest neighbor
  EXPECT_NEAR( 0.0f, fDistance, 1e-8); //distance
}

TEST(Matching, ArrayMatcher_Kdtree_Flann_Simple__NN)
{
  const float array[] = {0, 1, 2, 5, 6};
  // no 3, because it involve the same dist as 1,1

  ArrayMatcher_Kdtree_Flann<float> matcher;
  EXPECT_TRUE( matcher.Build(array, 5, 1) );

  const float query[] = {2};
  IndMatches vec_nIndice;
  vector<float> vec_fDistance;
  const int NN = 5;
  EXPECT_TRUE( matcher.SearchNeighbours(query, 1, &vec_nIndice, &vec_fDistance, NN) );

  EXPECT_EQ( 5, vec_nIndice.size());
  EXPECT_EQ( 5, vec_fDistance.size());

  // Check distances:
  EXPECT_NEAR( vec_fDistance[0], Square(2.0f-2.0f), 1e-6);
  EXPECT_NEAR( vec_fDistance[1], Square(1.0f-2.0f), 1e-6);
  EXPECT_NEAR( vec_fDistance[2], Square(0.0f-2.0f), 1e-6);
  EXPECT_NEAR( vec_fDistance[3], Square(5.0f-2.0f), 1e-6);
  EXPECT_NEAR( vec_fDistance[4], Square(6.0f-2.0f), 1e-6);

  // Check indexes:
  EXPECT_EQ(IndMatch(0,2), vec_nIndice[0]);
  EXPECT_EQ(IndMatch(0,1), vec_nIndice[1]);
  EXPECT_EQ(IndMatch(0,0), vec_nIndice[2]);
  EXPECT_EQ(IndMatch(0,3), vec_nIndice[3]);
  EXPECT_EQ(IndMatch(0,4), vec_nIndice[4]);
}

//-- Test LIMIT case (empty arrays)

TEST(Matching, ArrayMatcherBruteForce_Simple_EmptyArrays)
{
  std::vector<float> array;
  ArrayMatcherBruteForce<float> matcher;
  EXPECT_FALSE( matcher.Build(&array[0], 0, 4) );

  int nIndice = -1;
  float fDistance = -1.0f;
  EXPECT_FALSE( matcher.SearchNeighbour( &array[0], &nIndice, &fDistance) );
}

TEST(Matching, ArrayMatcher_Kdtree_Flann_Simple_EmptyArrays)
{
  std::vector<float> array;
  ArrayMatcher_Kdtree_Flann<float> matcher;
  EXPECT_FALSE( matcher.Build(&array[0], 0, 4) );

  int nIndice = -1;
  float fDistance = -1.0f;
  EXPECT_FALSE( matcher.SearchNeighbour( &array[0], &nIndice, &fDistance) );
}

TEST(Matching, Cascade_Hashing_Simple_EmptyArrays)
{
  std::vector<float> array;
  ArrayMatcherCascadeHashing<float> matcher;
  EXPECT_FALSE( matcher.Build(&array[0], 0, 4) );

  int nIndice = -1;
  float fDistance = -1.0f;
  EXPECT_FALSE( matcher.SearchNeighbour( &array[0], &nIndice, &fDistance) );
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
