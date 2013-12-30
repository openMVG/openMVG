// Copyright (c) 2012, 2013 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GRAPH_CONNECTED_COMPONENT_H_
#define OPENMVG_GRAPH_CONNECTED_COMPONENT_H_

#include <lemon/connectivity.h>
#include <set>

namespace openMVG
{
namespace graphUtils
{

/// Export node of each CC (Connected Component) in a map
template <typename GraphT>
std::map<size_t, std::set<lemon::ListGraph::Node> >  exportGraphToMapSubgraphs(
  const GraphT & g)
{
  typedef lemon::ListGraph::NodeMap<size_t> IndexMap;
  IndexMap connectedNodeMap(g);
  int connectedComponentCount =  lemon::connectedComponents(g, connectedNodeMap);

  std::map<size_t, std::set<lemon::ListGraph::Node> > map_subgraphs;

  // Create subgraphs' map
  typedef lemon::ListGraph::NodeIt NodeIterator;
  NodeIterator itNode(g);
  for (IndexMap::MapIt it(connectedNodeMap);
    it != lemon::INVALID; ++it, ++itNode) {
    map_subgraphs[*it].insert(itNode);
  }
  return map_subgraphs;
}

} // namespace graphUtils
} // namespace openMVG

#endif // OPENMVG_GRAPH_CONNECTED_COMPONENT_H_
