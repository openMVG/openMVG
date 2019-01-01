// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015, 2017 WhuAegeanSea, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_graph_utils.hpp"

#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/graph/graph_builder.hpp"
#include "openMVG/types.hpp"


namespace openMVG {
namespace sfm {

static void KeepOnlyReferencedElement
(
  const std::set<IndexT> & set_remainingIds,
  const matching::PairWiseMatches & matches_in,
  matching::PairWiseMatches & matches_out
)
{
  matching::PairWiseMatches kept_matches;
  for (const auto & iter_pairwise_matches : matches_in)
  {
    if (set_remainingIds.count(iter_pairwise_matches.first.first) &&
        set_remainingIds.count(iter_pairwise_matches.first.second))
    {
      kept_matches.insert(iter_pairwise_matches);
    }
  }
  matches_out.swap(kept_matches);
}

bool PairsToConnectedComponents
(
  const Pair_Set & pairs,
  bool is_biedge,
  int min_nodes,
  std::map<IndexT, std::set<IndexT>>& subgraphs_ids
)
{
  subgraphs_ids.clear();

  using Graph = graph::indexedGraph::GraphT;
  graph::indexedGraph putativeGraph(pairs);

  // For global SFM, firstly remove the not bi-edge element
  if (is_biedge)
  {
    using EdgeMapAlias = Graph::EdgeMap<bool>;
    EdgeMapAlias cutMap(putativeGraph.g);
    if (lemon::biEdgeConnectedCutEdges(putativeGraph.g, cutMap) > 0)
    {
      // Some edges must be removed because they don't follow the biEdge condition.
      using EdgeIterator = Graph::EdgeIt;
      EdgeIterator itEdge(putativeGraph.g);
      for (EdgeMapAlias::MapIt it(cutMap); it != lemon::INVALID; ++it, ++itEdge)
      {
        if (*it)
        {
          putativeGraph.g.erase(itEdge);  // remove the not bi-edge element
        }
      }
    }
  }

  // Compute all subgraphs in the putative graph
  const int connectedComponentCount = lemon::countConnectedComponents(putativeGraph.g);
  if (connectedComponentCount >= 1)
  {
    const auto map_subgraphs = graph::exportGraphToMapSubgraphs<Graph, IndexT>(putativeGraph.g);
    for (const auto & iter_map_subgraphs : map_subgraphs)
    {
      if (iter_map_subgraphs.second.size() > min_nodes)
      {
        std::set<IndexT> subgraphNodes;
        const std::set<lemon::ListGraph::Node> & ccSet = iter_map_subgraphs.second;
        for (const auto & iter2 : ccSet)
        {
          const IndexT Id = (*putativeGraph.node_map_id)[iter2];
          subgraphNodes.insert(Id);
        }
        subgraphs_ids.emplace(iter_map_subgraphs.first, subgraphNodes);
      }
    }
  }
  return !subgraphs_ids.empty();
}

bool SplitMatchesIntoSubgraphMatches
(
  const Pair_Set & pairs,
  const matching::PairWiseMatches & matches,
  bool is_biedge,
  int min_nodes,
  std::vector<matching::PairWiseMatches> & subgraphs_matches
)
{
  if (pairs.size() == 0 || matches.size() == 0)
  {
    return false;
  }
  subgraphs_matches.clear();
  std::map<IndexT, std::set<IndexT>> subgraphs_ids;
  if (PairsToConnectedComponents(pairs, is_biedge, min_nodes, subgraphs_ids))
  {
    for (const auto & iter_subgraphs_ids : subgraphs_ids)
    {
      matching::PairWiseMatches component_matches;
      KeepOnlyReferencedElement(iter_subgraphs_ids.second, matches, component_matches);
      subgraphs_matches.emplace_back(component_matches);
    }
  }
  return !subgraphs_matches.empty();
}

} // namespace sfm
} // namespace openMVG
