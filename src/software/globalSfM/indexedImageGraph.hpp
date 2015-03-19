
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_INDEXED_IMAGE_GRAPH_H_
#define OPENMVG_INDEXED_IMAGE_GRAPH_H_

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <memory>

#include "lemon/list_graph.h"
using namespace lemon;
#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::matching;

namespace openMVG  {
namespace imageGraph  {

// Structure used to keep information of an image graph:
//  - A graph (connection between nodes)
//  - Node => linked to String (the name of the image)
//  - EdgeMap => Number of point connection between the source and target
//
struct indexedImageGraph
{
  typedef lemon::ListGraph GraphT;
  typedef std::map<size_t, GraphT::Node> map_IndexT_Node;
  typedef GraphT::NodeMap<size_t> map_NodeMapIndex;
  typedef GraphT::NodeMap<std::string> map_NodeMapName;
  typedef GraphT::EdgeMap<size_t> map_EdgeMap;

  GraphT g;
  map_IndexT_Node map_IndexT_to_node; // Original image index to graph node
  std::auto_ptr<map_NodeMapIndex> map_nodeMapIndex; // Association of data to graph Node
  std::auto_ptr<map_NodeMapName> map_codeMapName; // Association of data to graph Node
  std::auto_ptr<map_EdgeMap> map_edgeMap; // Number of point matches between the source and the target

  // Build from pairwise matches and corresponding view filename
  indexedImageGraph( const PairWiseMatches & map_indexedMatches,
    const std::vector<std::string> &vec_fileNames)
  {
    map_nodeMapIndex = std::auto_ptr<map_NodeMapIndex>(new map_NodeMapIndex(g));
    map_codeMapName = std::auto_ptr<map_NodeMapName>(new map_NodeMapName(g));
    map_edgeMap = std::auto_ptr<map_EdgeMap>(new map_EdgeMap(g));

    //A-- Compute the number of node we need
    std::set<IndexT> setNodes;
    for (PairWiseMatches::const_iterator iter = map_indexedMatches.begin();
      iter != map_indexedMatches.end();
      ++iter)
    {
      setNodes.insert(iter->first.first);
      setNodes.insert(iter->first.second);
    }

    //B-- Create a node graph for each element of the set
    for (std::set<IndexT>::const_iterator iter = setNodes.begin();
      iter != setNodes.end();
      ++iter)
    {
      map_IndexT_to_node[*iter] = g.addNode();
      (*map_nodeMapIndex) [map_IndexT_to_node[*iter]] = *iter;
      (*map_codeMapName) [map_IndexT_to_node[*iter]] = vec_fileNames[*iter];
    }

    //C-- Add weighted edges from the "map_indexedMatches" object
    for (PairWiseMatches::const_iterator iter = map_indexedMatches.begin();
      iter != map_indexedMatches.end();
      ++iter)
    {
      const std::vector<IndMatch> & vec_FilteredMatches = iter->second;
      if (vec_FilteredMatches.size() > 0)
      {
        const IndexT i = iter->first.first;
        const IndexT j = iter->first.second;
        GraphT::Edge edge =  g.addEdge(map_IndexT_to_node[i], map_IndexT_to_node[j]);
        (*map_edgeMap)[ edge ] = vec_FilteredMatches.size();
      }
    }
  }

  // Build the view graph from a sequence of pair
  template <typename IterablePairSequence>
  indexedImageGraph(const IterablePairSequence & pairs)
  {
    map_nodeMapIndex = std::auto_ptr<map_NodeMapIndex>( new map_NodeMapIndex(g) );
    map_edgeMap = std::auto_ptr<map_EdgeMap>( new map_EdgeMap(g) );

    //A-- Compute the number of node we need
    std::set<IndexT> setNodes;
    for (typename IterablePairSequence::const_iterator iter = pairs.begin();
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
      map_IndexT_to_node[*iter] = g.addNode();
      (*map_nodeMapIndex) [map_IndexT_to_node[*iter]] = *iter;
    }

    //C-- Add weighted edges from the pairs object
    for (typename IterablePairSequence::const_iterator iter = pairs.begin();
      iter != pairs.end();
      ++iter)
    {
      const IndexT i = iter->first;
      const IndexT j = iter->second;
      GraphT::Edge edge =  g.addEdge(map_IndexT_to_node[i], map_IndexT_to_node[j]);
      (*map_edgeMap)[ edge ] = 1;
    }
  }
};

} // namespace imageGraph
} // namespace openMVG

#endif // OPENMVG_INDEXED_IMAGE_GRAPH_H_
