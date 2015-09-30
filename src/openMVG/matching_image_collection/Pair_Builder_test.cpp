// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "testing/testing.h"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_view.hpp"

#include <iostream>
#include <algorithm>
#include <memory>

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
  sfm::Views views;
  {
    // Empty
    Pair_Set pairSet = exhaustivePairs(views);
    EXPECT_EQ( 0, pairSet.size());
  }
  {
    std::vector<IndexT> indexes = {{ 12, 54, 89, 65 }};
    for( IndexT i: indexes )
    {
      views[i] = std::make_shared<sfm::View>("filepath", i);
    }


    Pair_Set pairSet = exhaustivePairs(views);
    EXPECT_TRUE( checkPairOrder(pairSet) );
    EXPECT_EQ( 6, pairSet.size());
    EXPECT_TRUE( pairSet.find(std::make_pair(12,54)) != pairSet.end() );
    EXPECT_TRUE( pairSet.find(std::make_pair(12,89)) != pairSet.end() );
    EXPECT_TRUE( pairSet.find(std::make_pair(12,65)) != pairSet.end() );
    EXPECT_TRUE( pairSet.find(std::make_pair(54,89)) != pairSet.end() );
    EXPECT_TRUE( pairSet.find(std::make_pair(54,65)) != pairSet.end() );
    EXPECT_TRUE( pairSet.find(std::make_pair(65,89)) != pairSet.end() );
  }
}

TEST(matching_image_collection, contiguousWithOverlap)
{
  sfm::Views views;
  {
    // Empty
    Pair_Set pairSet = exhaustivePairs(views);
    EXPECT_EQ( 0, pairSet.size());
  }
  {
    std::vector<IndexT> indexes = {{ 12, 54, 89, 65 }};
    for( IndexT i: indexes )
    {
      views[i] = std::make_shared<sfm::View>("filepath", i);
    }

    Pair_Set pairSet = contiguousWithOverlap(views, 1);
    EXPECT_TRUE( checkPairOrder(pairSet) );
    EXPECT_EQ( 3, pairSet.size());
    EXPECT_TRUE( pairSet.find(std::make_pair(12,54)) != pairSet.end() );
    EXPECT_TRUE( pairSet.find(std::make_pair(54,65)) != pairSet.end() );
    EXPECT_TRUE( pairSet.find(std::make_pair(65,89)) != pairSet.end() );
  }
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
  EXPECT_TRUE( loadPairs("pairsT_IO.txt", loaded_Pairs));
  EXPECT_TRUE( std::equal(loaded_Pairs.begin(), loaded_Pairs.end(), pairSetGTsorted.begin()) );
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
