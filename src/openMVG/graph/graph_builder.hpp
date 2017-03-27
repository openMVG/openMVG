// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GRAPH_GRAPH_BUILDER_HPP
#define OPENMVG_GRAPH_GRAPH_BUILDER_HPP

#include <map>
#include <memory>
#include <set>

#include "openMVG/types.hpp"

#include "lemon/list_graph.h"

namespace openMVG
{
namespace graph
{


/**
* @brief Structure used to keep information of an image graph:
  - Build a graph (add nodes and connection between nodes)
*/
struct indexedGraph
{
  /// Type of graph
  using GraphT = lemon::ListGraph;

  /// Type of index of nodes
  using map_NodeMapIndex = GraphT::NodeMap<IndexT>;

  /// The graph
  GraphT g;

  /// Association of data to graph Node
  std::unique_ptr<map_NodeMapIndex> node_map_id;


  /**
  * @brief Build a graph from a list of pair
  * @param pairs List of pairs
  */
  template <typename IterablePairs>
  indexedGraph( const IterablePairs & pairs )
  {
    node_map_id.reset( new map_NodeMapIndex( g ) );

    //A-- Compute the number of node we need
    std::set<IndexT> nodes;
    for (const auto & pair_it : pairs)
    {
      nodes.insert( pair_it.first );
      nodes.insert( pair_it.second );
    }

    //B-- Create a node graph for each element of the set
    std::map<IndexT, GraphT::Node> id_to_node;
    for ( const auto & node_it : nodes)
    {
      const GraphT::Node n = g.addNode();
      ( *node_map_id ) [n] = node_it;
      id_to_node[node_it] = n;
    }

    //C-- Add weighted edges from the pairs object
    for (const auto & pair_it : pairs)
    {
      g.addEdge( id_to_node[pair_it.first], id_to_node[pair_it.second] );
    }
  }

  /**
  * @brief Create a graph from node index and pairs (edges)
  * @param nodes List of nodes
  * @param pairs List of pairs
  * @note pairs must contains valid nodes indexes
  */
  template <typename IterableNodes, typename IterablePairs>
  indexedGraph( const IterableNodes & nodes, const IterablePairs & pairs )
  {
    node_map_id.reset( new map_NodeMapIndex( g ) );

    std::set<IndexT> node_set (nodes.begin(), nodes.end());
    //A-- Create a node graph for each element of the set
    std::map<IndexT, GraphT::Node> id_to_node;
    for ( const auto & node_it : node_set)
    {
      const GraphT::Node n = g.addNode();
      ( *node_map_id ) [n] = node_it;
      id_to_node[node_it] = n;
    }

    //B-- Add weighted edges from the pairs object
    for (const auto & pair_it : pairs)
    {
      g.addEdge( id_to_node[pair_it.first], id_to_node[pair_it.second] );
    }
  }
};

} // namespace graph
} // namespace openMVG

#endif // OPENMVG_GRAPH_GRAPH_BUILDER_HPP
