// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GRAPH_TRIPLET_FINDER_H
#define OPENMVG_GRAPH_TRIPLET_FINDER_H

#include "lemon/list_graph.h"
using namespace lemon;

#include <algorithm>
#include <vector>

namespace openMVG{
namespace graphUtils{

/// Simple container for tuple of three value
/// It is used to store the node id of triplets of a graph.
struct Triplet
{
  Triplet(size_t ii, size_t jj, size_t kk)
    :i(ii), j(jj), k(kk)
  { }
  size_t i,j,k;

  bool contain(const std::pair<size_t,size_t> & edge) const
  {
    size_t It = edge.first;
    size_t Jt = edge.second;
    if ( (It == i || It == j || It == k ) &&
         (Jt == i || Jt == j || Jt == k ) && It != Jt)
      return true;
    else
      return false;
  }

  friend bool operator==(const Triplet& m1, const Triplet& m2)  {
    return m1.contain(std::make_pair(m2.i, m2.j))
     && m1.contain(std::make_pair(m2.i, m2.k));
  }

  friend bool operator!=(const Triplet& m1, const Triplet& m2)  {
    return !(m1 == m2);
  }

  friend std::ostream & operator<<(std::ostream & os, const Triplet & t)
  {
    os << t.i << " " << t.j << " " << t.k << std::endl;
    return os;
  }
};

/// Function that return all the triplets found in a graph
/// vec_triplets must be empty.
template<typename GraphT>
bool List_Triplets(const GraphT & g, std::vector< Triplet > & vec_triplets)
{
  // Algorithm
  //
  //-- For each node
  //    - list the outgoing not visited edge
  //    -  for each tuple of edge
  //       - if their end are connected
  //          Detected cyle of length 3
  //          Mark first edge as visited

  typedef typename GraphT::OutArcIt OutArcIt;
  typedef typename GraphT::NodeIt NodeIterator;
  typedef typename GraphT::template EdgeMap<bool> BoolEdgeMap;

  BoolEdgeMap map_edge(g, false); // Visited edge map

  // For each nodes
  for (NodeIterator itNode(g); itNode != INVALID; ++itNode)
  {
    // For each edges (list the not visited outgoing edges)
    std::vector<OutArcIt> vec_edges;
    for (OutArcIt e(g, itNode); e!=INVALID; ++e)
    {
      if (!map_edge[e]) // If not visited
        vec_edges.push_back(e);
    }

    // For all tuples look of ends of edges are linked
    while(vec_edges.size()>1)
    {
      OutArcIt itPrev = vec_edges[0]; // For all tuple (0,Inth)
      for(size_t i=1; i < vec_edges.size(); ++i)
      {
        // Check if the extremity is linked
        typename GraphT::Arc cycleEdge = findArc(g, g.target(itPrev), g.target(vec_edges[i]));
        if (cycleEdge!= INVALID && !map_edge[cycleEdge])
        {
          // Elementary cycle found (make value follow a monotonic ascending serie)
          double triplet[3] = {
            g.id(itNode),
            g.id(g.target(itPrev)),
            g.id(g.target(vec_edges[i]))};
          std::sort(&triplet[0], &triplet[3]);
          vec_triplets.push_back(Triplet(triplet[0],triplet[1],triplet[2]));
        }
      }
      // Mark the current ref edge as visited
      map_edge[itPrev] = true;
      // remove head to list remaining tuples
      vec_edges.erase(vec_edges.begin());
    }
  }
  return (!vec_triplets.empty());
}

} // namespace graphUtils
} // namespace openMVG

#endif // OPENMVG_GRAPH_TRIPLET_FINDER_H
