// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015, 2017 WhuAegeanSea, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_graph_utils.hpp"

#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;


TEST(SFM_DATA_GRAPH, PairsToConnectedComponents)
{
  Pair_Set pair;

  // connected graph follows the biEdge condition: 0 -> 1 -> 2 ->3 ->4 -> 0
  pair.emplace(std::make_pair(0, 1));
  pair.emplace(std::make_pair(0, 2));
  pair.emplace(std::make_pair(0, 3));
  pair.emplace(std::make_pair(1, 2));
  pair.emplace(std::make_pair(1, 3));
  pair.emplace(std::make_pair(2, 3));

  // connected graph: 4 -> 5 -> 6 -> 7 -> 8
  pair.emplace(std::make_pair(4, 5));
  pair.emplace(std::make_pair(5, 6));
  pair.emplace(std::make_pair(6, 7));
  pair.emplace(std::make_pair(7, 8));

  // connected graph: 9 -> 10
  pair.emplace(std::make_pair(9, 10));

  std::map<IndexT, std::set<IndexT>> subgraphs_ids;

  // test for GlobalSFM with the biEdge condition
  const bool flag1 = PairsToConnectedComponents(pair, true, 3, subgraphs_ids);
  EXPECT_TRUE(flag1);
  EXPECT_EQ(1, subgraphs_ids.size());
  if (subgraphs_ids.size() == 1)
  {
    const auto & iter_begin = subgraphs_ids.begin();
    const std::set<IndexT> & pairs0 = iter_begin->second;
    EXPECT_EQ(4, pairs0.size());
    if (pairs0.size() == 4)
    {
      EXPECT_TRUE(pairs0.find(0) != pairs0.end());
      EXPECT_TRUE(pairs0.find(1) != pairs0.end());
      EXPECT_TRUE(pairs0.find(2) != pairs0.end());
      EXPECT_TRUE(pairs0.find(3) != pairs0.end());
    }
  }

  // test for IncrementalSFM with no biEdge condition
  const bool flag2 = PairsToConnectedComponents(pair, false, 3, subgraphs_ids);
  EXPECT_TRUE(flag2);
  EXPECT_EQ(2, subgraphs_ids.size());
  if (subgraphs_ids.size() == 2)
  {
    auto & iter_first = subgraphs_ids.begin();
    const std::set<IndexT> & pairs0 = iter_first->second;
    EXPECT_EQ(5, pairs0.size());
    if (pairs0.size() == 5)
    {
      EXPECT_TRUE(pairs0.find(4) != pairs0.end());
      EXPECT_TRUE(pairs0.find(5) != pairs0.end());
      EXPECT_TRUE(pairs0.find(6) != pairs0.end());
      EXPECT_TRUE(pairs0.find(7) != pairs0.end());
      EXPECT_TRUE(pairs0.find(8) != pairs0.end());
    }

    auto & iter_second = ++iter_first;
    const std::set<IndexT> & pairs1 = iter_second->second;
    EXPECT_EQ(4, pairs1.size());
    if (pairs1.size() == 4)
    {
      EXPECT_TRUE(pairs1.find(0) != pairs1.end());
      EXPECT_TRUE(pairs1.find(1) != pairs1.end());
      EXPECT_TRUE(pairs1.find(2) != pairs1.end());
      EXPECT_TRUE(pairs1.find(3) != pairs1.end());
    }
  }
}

/* ************************************************************************* */
int main() {
  TestResult tr; return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */
