// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "testing/testing.h"

#include <iostream>
#include <algorithm>
using namespace std;

using namespace openMVG;

// Check pairs follow a weak ordering pair.first < pair.second
template<typename IterablePairs>
bool checkPairOrder(const IterablePairs & pairs)
{
  for (typename IterablePairs::const_iterator iterP = pairs.begin(); iterP != pairs.end();
    ++iterP)
  {
    if (iterP->first >= iterP->second)
      return false;
  }
  return true;
}

TEST(matching_image_collection, exhaustivePairs)
{
  Pair_Set pairSet = exhaustivePairs(0);
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
  Pair_Set pairSet = exhaustivePairs(0);
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
  Pair_Set pairSetGT;
  pairSetGT.insert( std::make_pair(0,1) );
  pairSetGT.insert( std::make_pair(1,2) );
  pairSetGT.insert( std::make_pair(2,0) );

  Pair_Set pairSetGTsorted;
  pairSetGTsorted.insert( std::make_pair(0,1) );
  pairSetGTsorted.insert( std::make_pair(0,2) );
  pairSetGTsorted.insert( std::make_pair(1,2) );

  EXPECT_TRUE( savePairs("pairsT_IO.txt", pairSetGT));

  Pair_Set loaded_Pairs;
  EXPECT_TRUE( loadPairs(3, "pairsT_IO.txt", loaded_Pairs));
  EXPECT_TRUE( std::equal(loaded_Pairs.begin(), loaded_Pairs.end(), pairSetGTsorted.begin()) );
}

TEST(matching_image_collection, IO_InvalidInput)
{
  // A pair with index superior to the expected picture count
  Pair_Set pairSetGT;
  pairSetGT.insert( std::make_pair(0,1) );
  pairSetGT.insert( std::make_pair(10,20) );

  EXPECT_TRUE( savePairs("pairsT_IO_InvalidInput.txt", pairSetGT));

  Pair_Set loaded_Pairs;
  const int expectedPicCount = 2;
  EXPECT_FALSE( loadPairs(expectedPicCount, "pairsT_IO_InvalidInput.txt", loaded_Pairs));

  // A pair with equal index
  pairSetGT.clear();
  pairSetGT.insert( std::make_pair(0,1) );
  pairSetGT.insert( std::make_pair(0,0) );
  EXPECT_FALSE( loadPairs(expectedPicCount, "pairsT_IO_InvalidInput.txt", loaded_Pairs));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
