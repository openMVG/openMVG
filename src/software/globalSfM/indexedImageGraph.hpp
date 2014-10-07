
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

#include "lemon/list_graph.h"
using namespace lemon;
#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::matching;

namespace openMVG  {
namespace imageGraph  {
  using namespace std;

// Structure used to keep information of an image graph :
//  - A graph (connection between nodes
//  - Node => linked to String (the name of the image)
//  - EdgeMap => Number of point connection between the source and target
//
struct indexedImageGraph
{
  typedef lemon::ListGraph GraphT;
  typedef std::map<size_t, GraphT::Node> map_Size_t_Node;
  typedef GraphT::NodeMap<size_t> map_NodeMapIndex;
  typedef GraphT::NodeMap<std::string> map_NodeMapName;
  typedef GraphT::EdgeMap<size_t> map_EdgeMap;

  GraphT g;
  map_Size_t_Node map_size_t_to_node; // Original image index to graph node
  auto_ptr<map_NodeMapIndex> map_nodeMapIndex; // Association of data to graph Node
  auto_ptr<map_NodeMapName> map_codeMapName; // Association of data to graph Node
  auto_ptr<map_EdgeMap> map_edgeMap; // Number of point matches between the source and the target

  indexedImageGraph( const PairWiseMatches & map_indexedMatches,
    const std::vector<string> &vec_fileNames)
  {
    map_nodeMapIndex =  auto_ptr<map_NodeMapIndex>( new map_NodeMapIndex(g) );
    map_codeMapName =  auto_ptr<map_NodeMapName>( new map_NodeMapName(g) );
    map_edgeMap = auto_ptr<map_EdgeMap>( new map_EdgeMap(g) );

    //A-- Compute the number of node we need
    set<size_t> setNodes;
    for (PairWiseMatches::const_iterator iter = map_indexedMatches.begin();
      iter != map_indexedMatches.end();
      ++iter)
    {
      setNodes.insert(iter->first.first);
      setNodes.insert(iter->first.second);
    }

    //B-- Create a node graph for each element of the set
    for (set<size_t>::const_iterator iter = setNodes.begin();
      iter != setNodes.end();
      ++iter)
    {
      map_size_t_to_node[*iter] = g.addNode();
      (*map_nodeMapIndex) [map_size_t_to_node[*iter]] = *iter;
      (*map_codeMapName) [map_size_t_to_node[*iter]] = vec_fileNames[*iter];
    }

    //C-- Add weighted edges from the "map_indexedMatches" object
    for (PairWiseMatches::const_iterator iter = map_indexedMatches.begin();
      iter != map_indexedMatches.end();
      ++iter)
    {
      const std::vector<IndMatch> & vec_FilteredMatches = iter->second;
      if (vec_FilteredMatches.size() > 0)
      {
        const size_t i = iter->first.first;
        const size_t j = iter->first.second;
        GraphT::Edge edge =  g.addEdge(map_size_t_to_node[i], map_size_t_to_node[j]);
        (*map_edgeMap)[ edge ] = vec_FilteredMatches.size();
      }
    }
  }

  indexedImageGraph( const std::vector<std::pair<size_t, size_t> > & map_pairs,
    const std::vector<string> &vec_fileNames)
  {
    typedef std::vector<std::pair<size_t, size_t> > Pairs_T;
    map_nodeMapIndex =  auto_ptr<map_NodeMapIndex>( new map_NodeMapIndex(g) );
    map_codeMapName =  auto_ptr<map_NodeMapName>( new map_NodeMapName(g) );
    map_edgeMap = auto_ptr<map_EdgeMap>( new map_EdgeMap(g) );

    //A-- Compute the number of node we need
    set<size_t> setNodes;
    for (Pairs_T::const_iterator iter = map_pairs.begin();
      iter != map_pairs.end();
      ++iter)
    {
      setNodes.insert(iter->first);
      setNodes.insert(iter->second);
    }

    //B-- Create a node graph for each element of the set
    for (set<size_t>::const_iterator iter = setNodes.begin();
      iter != setNodes.end();
      ++iter)
    {
      map_size_t_to_node[*iter] = g.addNode();
      (*map_nodeMapIndex) [map_size_t_to_node[*iter]] = *iter;
      (*map_codeMapName) [map_size_t_to_node[*iter]] = vec_fileNames[*iter];
    }

    //C-- Add weighted edges from the pairs object
    for (Pairs_T::const_iterator iter = map_pairs.begin();
      iter != map_pairs.end();
      ++iter)
    {
      const size_t i = iter->first;
      const size_t j = iter->second;
      GraphT::Edge edge =  g.addEdge(map_size_t_to_node[i], map_size_t_to_node[j]);
      (*map_edgeMap)[ edge ] = 1;
    }
  }
};

} // namespace imageGraph
} // namespace openMVG

#endif // OPENMVG_INDEXED_IMAGE_GRAPH_H_
