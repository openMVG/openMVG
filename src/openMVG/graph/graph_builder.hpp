// Copyright (c) 2012, 2013 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GRAPH_BUILDER_H_
#define OPENMVG_GRAPH_BUILDER_H_

#include <memory>
#include <set>

namespace openMVG {
namespace graph  {

// Structure used to keep information of an image graph:
//  - Build a graph (add nodes and connection between nodes)
//
struct indexedGraph
{
  typedef lemon::ListGraph GraphT;
  typedef std::map<IndexT, GraphT::Node> map_Size_t_Node;
  typedef GraphT::NodeMap<IndexT> map_NodeMapIndex;

  GraphT g;
  map_Size_t_Node map_size_t_to_node; // Original image index to graph node
  std::unique_ptr<map_NodeMapIndex> map_nodeMapIndex; // Association of data to graph Node

  template <typename IterablePairs>
  indexedGraph(const IterablePairs & pairs)
  {
    map_nodeMapIndex.reset( new map_NodeMapIndex(g) );

    //A-- Compute the number of node we need
    std::set<IndexT> setNodes;
    for (typename IterablePairs::const_iterator iter = pairs.begin();
      iter != pairs.end();
      ++iter)
    {
      setNodes.insert(iter->first);
      setNodes.insert(iter->second);
    }

    //B-- Create a node graph for each element of the set
    for (std::set<IndexT>::const_iterator iter = setNodes.begin();
      iter != setNodes.end();
      ++iter)
    {
      map_size_t_to_node[*iter] = g.addNode();
      (*map_nodeMapIndex) [map_size_t_to_node.at(*iter)] = *iter;
    }

    //C-- Add weighted edges from the pairs object
    for (typename IterablePairs::const_iterator iter = pairs.begin();
      iter != pairs.end();
      ++iter)
    {
      const IndexT i = iter->first;
      const IndexT j = iter->second;
      g.addEdge(map_size_t_to_node.at(i), map_size_t_to_node.at(j));
    }
  }

  /// Create a graph from node index and pairs (edges)
  // /!\ pairs must contains valid nodes indexes
  template <typename IterableNodes, typename IterablePairs>
  indexedGraph(const IterableNodes & nodes, const IterablePairs & pairs)
  {
    map_nodeMapIndex.reset( new map_NodeMapIndex(g) );

    //A-- Create a node graph for each element of the set
    for (typename IterableNodes::const_iterator iter = nodes.begin();
      iter != nodes.end();
      ++iter)
    {
      map_size_t_to_node[*iter] = g.addNode();
      (*map_nodeMapIndex) [map_size_t_to_node.at(*iter)] = *iter;
    }

    //B-- Add weighted edges from the pairs object
    for (typename IterablePairs::const_iterator iter = pairs.begin();
      iter != pairs.end();
      ++iter)
    {
      const IndexT i = iter->first;
      const IndexT j = iter->second;
      g.addEdge(map_size_t_to_node.at(i), map_size_t_to_node.at(j));
    }
  }
};

} // namespace graph
} // namespace openMVG

#endif // OPENMVG_GRAPH_BUILDER__H_
