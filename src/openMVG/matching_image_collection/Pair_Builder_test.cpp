// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "testing/testing.h"

#include <iostream>
using namespace std;

using namespace openMVG;

// Check pairs follow a weak ordering pair.first < pair.second
bool checkPairOrder(const PairsT & pairs)
{
  for (PairsT::const_iterator iterP = pairs.begin(); iterP != pairs.end();
    ++iterP)
  {
    if (iterP->first >= iterP->second)
      return false;
  }
  return true;
}

TEST(matching_image_collection, exhaustivePairs)
{
  PairsT pairSet = exhaustivePairs(0);
  EXPECT_EQ( 0, pairSet.size());

  pairSet = exhaustivePairs(4);
  EXPECT_TRUE( checkPairOrder(pairSet) );
  EXPECT_EQ( 6, pairSet.size());
  EXPECT_TRUE( pairSet.find(std::make_pair(0,1)) != pairSet.end() );
  EXPECT_TRUE( pairSet.find(std::make_pair(0,2)) != pairSet.end() );
  EXPECT_TRUE( pairSet.find(std::make_pair(0,3)) != pairSet.end() );
  EXPECT_TRUE( pairSet.find(std::make_pair(1,2)) != pairSet.end() );
  EXPECT_TRUE( pairSet.find(std::make_pair(1,3)) != pairSet.end() );
  EXPECT_TRUE( pairSet.find(std::make_pair(2,3)) != pairSet.end() );
}

TEST(matching_image_collection, contiguousWithOverlap)
{
  PairsT pairSet = exhaustivePairs(0);
  EXPECT_EQ( 0, pairSet.size());

  pairSet = contiguousWithOverlap(4,1);
  EXPECT_TRUE( checkPairOrder(pairSet) );
  EXPECT_EQ( 3, pairSet.size());
  EXPECT_TRUE( pairSet.find(std::make_pair(0,1)) != pairSet.end() );
  EXPECT_TRUE( pairSet.find(std::make_pair(1,2)) != pairSet.end() );
  EXPECT_TRUE( pairSet.find(std::make_pair(2,3)) != pairSet.end() );
}

TEST(matching_image_collection, IO)
{
  PairsT pairSetGT;
  pairSetGT.insert( std::make_pair(0,2) );
  pairSetGT.insert( std::make_pair(2,4) );

  EXPECT_TRUE( savePairs("pairsT_IO.txt", pairSetGT));

  PairsT loaded_Pairs;
  EXPECT_TRUE( loadPairs("pairsT_IO.txt", loaded_Pairs));
  EXPECT_TRUE( loaded_Pairs == pairSetGT);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
