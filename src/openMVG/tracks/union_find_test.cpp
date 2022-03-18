// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/tracks/union_find.hpp"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <set>

using namespace openMVG;

TEST(Tracks, union_find_logic) {

  // Initialize the union_find tree the track four element: {0, 1, 2, 3}
  UnionFind uf_tree;
  uf_tree.InitSets(4); // ids from 0 to 4
  EXPECT_EQ(4, uf_tree.GetNumNodes());

  // Link 0 and 1
  uf_tree.Union(0, 1);
  EXPECT_EQ(uf_tree.Find(0), uf_tree.Find(1));

  // The component id of 1,2,3 should not be linked
  CHECK(uf_tree.Find(1) != uf_tree.Find(2));
  CHECK(uf_tree.Find(2) != uf_tree.Find(3));

  // Link 0 to self (should not crash)
  uf_tree.Union(0, 0);

  // Link 2 to 3
  uf_tree.Union(2, 3);
  EXPECT_EQ(uf_tree.Find(2), uf_tree.Find(3));
  CHECK(uf_tree.Find(1) != uf_tree.Find(2));

  // Link 1 and 2
  //  It creates the link 1->2 and link the two previously existing component (0->1; 2->3)
  uf_tree.Union(1, 2);
  EXPECT_EQ(uf_tree.Find(1), uf_tree.Find(2));
  EXPECT_EQ(uf_tree.Find(0), uf_tree.Find(3));
}

TEST(Tracks, one_cc) {

  UnionFind uf_tree;
  uf_tree.InitSets(5);

  // Create the connections: 0-1-2-3-4
  uf_tree.Union(0, 1);
  uf_tree.Union(1, 2);
  uf_tree.Union(3, 4);
  uf_tree.Union(1, 3);

  // Run path compression to identify all the point belonging to the CC
  for (unsigned int i = 0; i < uf_tree.GetNumNodes(); ++i)
    uf_tree.Find(i);

  // Count the number of CC
  const std::set<unsigned int> parent_id(uf_tree.m_cc_parent.cbegin(), uf_tree.m_cc_parent.cend());
  EXPECT_EQ(1, parent_id.size());
}

TEST(Tracks, two_cc) {

  UnionFind uf_tree;
  uf_tree.InitSets(4);

  // Create the following connections
  // 0-1
  // 2-3
  uf_tree.Union(0, 1);
  uf_tree.Union(2, 3);

  // Run path compression to identify all the point belonging to the CC
  for (unsigned int i = 0; i < uf_tree.GetNumNodes(); ++i)
    uf_tree.Find(i);

  // Count the number of CC
  const std::set<unsigned int> parent_id(uf_tree.m_cc_parent.cbegin(), uf_tree.m_cc_parent.cend());
  EXPECT_EQ(2, parent_id.size());
}

TEST(Tracks, four_cc) {

  UnionFind uf_tree;
  uf_tree.InitSets(6);

  // Create the following connections
  // 0-1
  // 2-3
  // 4
  // 5
  uf_tree.Union(0, 1);
  uf_tree.Union(2, 3);

  // Run path compression to identify all the point belonging to the CC
  for (unsigned int i = 0; i < uf_tree.GetNumNodes(); ++i)
    uf_tree.Find(i);

  // Count the number of CC
  const std::set<unsigned int> parent_id(uf_tree.m_cc_parent.cbegin(), uf_tree.m_cc_parent.cend());
  EXPECT_EQ(4, parent_id.size());
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
