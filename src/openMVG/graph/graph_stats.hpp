// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GRAPH_GRAPH_STATS_HPP
#define OPENMVG_GRAPH_GRAPH_STATS_HPP

#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/graph/graph_builder.hpp"

#include "openMVG/numeric/numeric.h"

namespace openMVG
{
namespace graph
{

/// Export graph statistics to the console:
/// - number of node
/// - number of cc
/// - node degree statistics (min, max, mean, median)
template <typename IterablePairs>
void getGraphStatistics
(
  const int nb_nodes,
  const IterablePairs & pairs
)
{
  const indexedGraph graph(pairs);

  if (lemon::countNodes(graph.g) > nb_nodes)
    return;

  // Compute node degree for each node
  std::vector<int> nodes_degree;
  for (indexedGraph::GraphT::NodeIt n(graph.g); n != lemon::INVALID; ++n)
  {
    nodes_degree.emplace_back(lemon::countOutArcs(graph.g, n));
  }

  int min, max, mean, median;
  const bool stats = minMaxMeanMedian(
    nodes_degree.cbegin(), nodes_degree.cend(), min, max, mean, median);

  // Compute the number of connected component
  const auto connected_components =
    exportGraphToMapSubgraphs<indexedGraph::GraphT, IndexT>(graph.g);

  std::cout
    << "Graph statistics:\n"
    << "\t#nodes: " << nb_nodes<< "\n"
    << "\t#cc: " << connected_components.size() + (nb_nodes - lemon::countNodes(graph.g)) << "\n"
    << "\t#singleton: " << nb_nodes - lemon::countNodes(graph.g)<< "\n"
    << "\tNode degree statistics:\t"
    << "min: " << min << ", "
    << "max: " << max << ", "
    << "mean: " << mean << ", "
    << "median: " << median << "\n";
}

}
}


#endif // OPENMVG_GRAPH_GRAPH_STATS_HPP
