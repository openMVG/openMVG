// Copyright (c) 2012, 2013 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <vector>

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include "openMVG/graph/graph.hpp"
using namespace lemon;

TEST(connectedComponents, Empty) {
  lemon::ListGraph graph;

  const int connectedComponentCount = lemon::countConnectedComponents(graph);
  EXPECT_EQ(0, connectedComponentCount);
}

TEST(connectedComponents, OneCC) {
  lemon::ListGraph graph;
  lemon::ListGraph::Node a = graph.addNode(), b = graph.addNode();
  graph.addEdge(a,b);
  const int connectedComponentCount = lemon::countConnectedComponents(graph);
  EXPECT_EQ(1, connectedComponentCount);
}

TEST(connectedComponents, TwoCC) {
  lemon::ListGraph graph;

  lemon::ListGraph::Node a = graph.addNode(), b = graph.addNode();
  graph.addEdge(a,b);
  lemon::ListGraph::Node a2 = graph.addNode(), b2 = graph.addNode();
  graph.addEdge(a2,b2);
  const int connectedComponentCount = lemon::countConnectedComponents(graph);
  EXPECT_EQ(2, connectedComponentCount);
}


TEST(connectedComponents, TwoCC_Parsing) {
  lemon::ListGraph graph;

  lemon::ListGraph::Node a = graph.addNode(), b = graph.addNode();
  graph.addEdge(a,b);
  lemon::ListGraph::Node a2 = graph.addNode(), b2 = graph.addNode();
  graph.addEdge(a2,b2);

  typedef ListGraph::NodeMap<size_t> IndexMap;
  IndexMap connectedNodeMap(graph);
  const int connectedComponentCount =  lemon::connectedComponents(graph, connectedNodeMap);
  EXPECT_EQ(2, connectedComponentCount);
  for (IndexMap::MapIt it(connectedNodeMap); it != INVALID; ++it)
  {
    std::cout << *it << "\t" << graph.id(it) << "\n";
  }
}


/// Test to get back node id of each CC
// a
//
// b-c
//
// d-g
// | |
// e-f
//
// h-i-j-k
//   |/
//   l
TEST(exportGraphToMapSubgraphs, CC_Subgraph) {
  lemon::ListGraph graph;

  // single
  lemon::ListGraph::Node a = graph.addNode();

  // two
  lemon::ListGraph::Node b = graph.addNode(), c = graph.addNode();
  graph.addEdge(b,c);

  // four
  lemon::ListGraph::Node d = graph.addNode(), e = graph.addNode(),
    f = graph.addNode(), g = graph.addNode();
  graph.addEdge(d,e);
  graph.addEdge(e,f);
  graph.addEdge(f,g);
  graph.addEdge(g,d);

  // five
  lemon::ListGraph::Node h = graph.addNode(), i = graph.addNode(),
    j = graph.addNode(), k = graph.addNode(),l = graph.addNode();
  graph.addEdge(h,i);
  graph.addEdge(i,j);
  graph.addEdge(j,k);
  graph.addEdge(i,l);
  graph.addEdge(j,l);

  const std::map<size_t, std::set<lemon::ListGraph::Node> > map_subgraphs =
    openMVG::graph::exportGraphToMapSubgraphs<lemon::ListGraph, size_t>(graph);

  EXPECT_EQ(4, map_subgraphs.size());
  EXPECT_EQ(5, map_subgraphs.at(0).size());
  EXPECT_EQ(4, map_subgraphs.at(1).size());
  EXPECT_EQ(2, map_subgraphs.at(2).size());
  EXPECT_EQ(1, map_subgraphs.at(3).size());
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
