// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GRAPH_GRAPH_TRIPLET_FINDER_HPP
#define OPENMVG_GRAPH_GRAPH_TRIPLET_FINDER_HPP

#include <algorithm>
#include <utility>
#include <vector>

#include "openMVG/graph/graph.hpp"
#include "openMVG/types.hpp"

#include "lemon/list_graph.h"

namespace openMVG
{
namespace graph
{

/**
* @brief Simple container for tuple of three value
* @note It is used to store the node id of triplets of a graph.
* @todo Why not using std::tuple ?
*/
struct Triplet
{

  /**
  * @brief Constructor
  * @param ii First element of the triplet
  * @param jj Second element of the triplet
  * @param kk Third element of the triplet
  */
  Triplet( IndexT ii, IndexT jj, IndexT kk )
    : i( ii ), j( jj ), k( kk )
  { }

  /**
  * @brief Indicate if an edge contains one of the element of the triplet
  * @param edge Edge to test
  * @retval true if edge contains at least one of index of the triplet
  * @retval false if edge contains none of the index of the triplet
  */
  bool contain( const std::pair<IndexT, IndexT> & edge ) const
  {
    const IndexT It = edge.first;
    const IndexT Jt = edge.second;
    if ( ( It == i || It == j || It == k ) &&
         ( Jt == i || Jt == j || Jt == k ) && It != Jt )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  /**
  * @brief Indicate if triplet are the same
  * @param m1 First triplet
  * @param m2 Second triplet
  * @retval true if triplets are the same
  * @retval false if triplets are differents
  */
  friend bool operator==( const Triplet& m1, const Triplet& m2 )
  {
    return m1.contain( {m2.i, m2.j} )
           && m1.contain( {m2.i, m2.k} );
  }

  /**
  * @brief Indicate if triplet are differents
  * @param m1 First triplet
  * @param m2 Second triplet
  * @retval true if triplets are differents
  * @retval false if triplets are the same
  */
  friend bool operator!=( const Triplet& m1, const Triplet& m2 )
  {
    return !( m1 == m2 );
  }

  /**
  * @brief Stream operation
  * @param os Stream in which triplet is written
  * @param t Triplet to write
  * @return stream after write operation
  */
  friend std::ostream & operator<<( std::ostream & os, const Triplet & t )
  {
    os << t.i << " " << t.j << " " << t.k << std::endl;
    return os;
  }

  /// the three triplet index id
  IndexT i, j, k;
};


/**
* @brief Function that return all the triplets found in a graph
* @param g Input graph
* @param[out] vec_triplets Output list of triplet found in the graph
* @retval true If at least one triplet is found
* @retval false If no triplet is found
* @note vec_triplets must be empty.
*/
template<typename GraphT>
bool List_Triplets( const GraphT & g, std::vector< Triplet > & vec_triplets )
{
  // Algorithm
  //
  //-- For each node
  //    - list the outgoing not visited edge
  //    -  for each tuple of edge
  //       - if their end are connected
  //          Detected cycle of length 3
  //          Mark first edge as visited

  /// Type of graph iterator
  using OutArcIt = typename GraphT::OutArcIt;

  /// Type of node iterator
  using NodeIterator = typename GraphT::NodeIt;

  /// Type for edge maps
  using BoolEdgeMap = typename GraphT::template EdgeMap<bool>;

  /// List of visited edge map
  BoolEdgeMap map_edge( g, false ); // Visited edge map

  // For each nodes
  for ( NodeIterator itNode( g ); itNode != lemon::INVALID; ++itNode )
  {

    // For each edges (list the not visited outgoing edges)
    std::vector<OutArcIt> vec_edges;
    for ( OutArcIt e( g, itNode ); e != lemon::INVALID; ++e )
    {
      if ( !map_edge[e] ) // If not visited
      {
        vec_edges.push_back( e );
      }
    }

    // For all tuples look of ends of edges are linked
    while (vec_edges.size() > 1)
    {
      OutArcIt itPrev = vec_edges[0]; // For all tuple (0,Inth)
      for (size_t i = 1; i < vec_edges.size(); ++i )
      {
        // Check if the extremity is linked
        typename GraphT::Arc cycleEdge = findArc( g, g.target( itPrev ), g.target( vec_edges[i] ) );
        if ( cycleEdge != lemon::INVALID && !map_edge[cycleEdge] )
        {
          // Elementary cycle found (make value follow a monotonic ascending serie)
          int triplet[3] =
          {
            g.id( itNode ),
            g.id( g.target( itPrev ) ),
            g.id( g.target( vec_edges[i] ) )
          };
          std::sort( &triplet[0], &triplet[3] );
          vec_triplets.emplace_back( triplet[0], triplet[1], triplet[2] );
        }
      }
      // Mark the current ref edge as visited
      map_edge[itPrev] = true;
      // remove head to list remaining tuples
      vec_edges.erase( vec_edges.begin() );
    }
  }
  return ( !vec_triplets.empty() );
}


/**
* @brief Return triplets contained in the graph build from IterablePairs
* @param pairs Graph pairs
* @return List of triplet found in graph
*/
template <typename IterablePairs>
static std::vector< graph::Triplet > tripletListing
(
  const IterablePairs & pairs
)
{
  indexedGraph putativeGraph( pairs );
  std::vector< graph::Triplet > vec_triplets;
  graph::List_Triplets<indexedGraph::GraphT>( putativeGraph.g, vec_triplets );

  //Change triplets to ImageIds
  for ( auto & triplet : vec_triplets )
  {
    const IndexT I = ( *putativeGraph.node_map_id )[putativeGraph.g.nodeFromId( triplet.i )];
    const IndexT J = ( *putativeGraph.node_map_id )[putativeGraph.g.nodeFromId( triplet.j )];
    const IndexT K = ( *putativeGraph.node_map_id )[putativeGraph.g.nodeFromId( triplet.k )];
    IndexT triplet_[3] = { I, J, K };
    std::sort( &triplet_[0], &triplet_[3] );
    triplet = graph::Triplet( triplet_[0], triplet_[1], triplet_[2] );
  }
  return vec_triplets;
}

} // namespace graph
} // namespace openMVG

#endif // OPENMVG_GRAPH_GRAPH_TRIPLET_FINDER_HPP
