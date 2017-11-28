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

TEST(SFM_DATA_GRAPH, PairsToConnectedComponents_EmptyInput)
{
  std::map<IndexT, std::set<IndexT>> subgraphs_ids;
  EXPECT_FALSE(PairsToConnectedComponents({}, true, 1, subgraphs_ids));
}

TEST(SFM_DATA_GRAPH, PairsToConnectedComponents)
{
  const Pair_Set pairs =
  {
    // the first connected graph follows the biEdge condition: 0 -> 1 -> 2 ->3 ->4 -> 0
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3},
    // the second connected graph: 4 -> 5 -> 6 -> 7 -> 8
    {4, 5}, {5, 6}, {6, 7}, {7, 8},
    // the third connected graph: 9 -> 10
    {9, 10}
  };

  std::map<IndexT, std::set<IndexT>> subgraphs_ids;

  // test for GlobalSFM with the biEdge condition
  EXPECT_TRUE(PairsToConnectedComponents(pairs, true, 3, subgraphs_ids));
  EXPECT_EQ(1, subgraphs_ids.size());
  if (subgraphs_ids.size() == 1)
  {
    const auto & iter_begin = subgraphs_ids.cbegin();
    const std::set<IndexT> & kept_idx = iter_begin->second;
    EXPECT_EQ(4, kept_idx.size());
    {
      EXPECT_TRUE(kept_idx.find(0) != kept_idx.cend());
      EXPECT_TRUE(kept_idx.find(1) != kept_idx.cend());
      EXPECT_TRUE(kept_idx.find(2) != kept_idx.cend());
      EXPECT_TRUE(kept_idx.find(3) != kept_idx.cend());
    }
  }

  // test for IncrementalSFM with no biEdge condition
  // The smallest graph must disappear since it is too small
  EXPECT_TRUE(PairsToConnectedComponents(pairs, false, 3, subgraphs_ids));
  EXPECT_EQ(2, subgraphs_ids.size());
  if (subgraphs_ids.size() == 2)
  {
    const auto iter_first = subgraphs_ids.cbegin();
    const std::set<IndexT> & kept_idx = iter_first->second;
    EXPECT_EQ(5, kept_idx.size());
    {
      EXPECT_TRUE(kept_idx.find(4) != kept_idx.cend());
      EXPECT_TRUE(kept_idx.find(5) != kept_idx.cend());
      EXPECT_TRUE(kept_idx.find(6) != kept_idx.cend());
      EXPECT_TRUE(kept_idx.find(7) != kept_idx.cend());
      EXPECT_TRUE(kept_idx.find(8) != kept_idx.cend());
    }

    auto iter_second = subgraphs_ids.begin();
    std::advance(iter_second, 1);
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
