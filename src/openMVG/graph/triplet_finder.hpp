// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013. 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GRAPH_GRAPH_TRIPLET_FINDER_HPP
#define OPENMVG_GRAPH_GRAPH_TRIPLET_FINDER_HPP

#include <algorithm>
#include <array>
#include <ostream>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "openMVG/types.hpp"

namespace openMVG
{
namespace graph
{

/**
* @brief Simple container for a tuple of three value
* @note It is used to store the node id of triplets of a graph.
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
    return ( ( It == i || It == j || It == k ) &&
             ( Jt == i || Jt == j || Jt == k ) && It != Jt );
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
    return m1.contain( { m2.i, m2.j } ) &&
           m1.contain( { m2.i, m2.k } );
  }

  /**
  * @brief Indicate if triplet are differents
  * @param m1 First triplet
  * @param m2 Second triplet
  * @retval true if triplets are differents
  * @retval false if triplets are the same
  */
  friend bool operator!=( const Triplet & m1, const Triplet & m2 )
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
    os << t.i << ' ' << t.j << ' ' << t.k << '\n';
    return os;
  }

  /// the three triplet index id
  IndexT i, j, k;
};

/**
* @brief Return triplets contained in the graph build from IterablePairs
* @param[in] pairs A list of pairs
* @param[out] triplets List of triplet found in graph
* @return boolean return true if some triplet are found
**/
template <typename IterablePairs, class TTripletContainer>
bool ListTriplets
(
  const IterablePairs & pairs,
  TTripletContainer & triplets
)
{
  triplets.clear();

  // Build an adjacency list corresponding to the edge list
  std::unordered_map<IndexT, std::set<IndexT>> adjacency_list;
  for (const auto & edge : pairs)
  {
    adjacency_list[edge.first].insert(edge.second);
    adjacency_list[edge.second].insert(edge.first);
  }

  std::vector<IndexT> node_candidate_for_triplet;

  // List the pair and find all triplets thanks to the adjacency list
  for (const auto & edge_it : pairs)
  {
    // Find any targeting edge that contains the first and the second node index
    const auto & node1_edges = adjacency_list.find(edge_it.first)->second;
    const auto & node2_edges = adjacency_list.find(edge_it.second)->second;

    // Compute the intersection between the two adjacency lists to find
    //  triplets (it will list the nodes that are connected to the first and
    //  second)
    node_candidate_for_triplet.clear();
    std::set_intersection(node1_edges.cbegin(),
                          node1_edges.cend(),
                          node2_edges.cbegin(),
                          node2_edges.cend(),
                          std::back_inserter(node_candidate_for_triplet));
    // Add a triplet
    for (const auto & node_index_it : node_candidate_for_triplet)
    {
      std::array<IndexT, 3> triplet_indexes {{
        static_cast<IndexT>(edge_it.first),
        static_cast<IndexT>(edge_it.second),
        node_index_it}};
      // sort the triplet indexes as i<j<k (monotonic ascending sorting)
      std::sort(triplet_indexes.begin(), triplet_indexes.end());
      triplets.emplace_back(triplet_indexes[0],
                            triplet_indexes[1],
                            triplet_indexes[2]);
    }

    // Since we have already listed all the triplets than contain this edge.
    // We can now remove the edge from the adjacency graph to reduce the
    // node array size for the next iterations.
    adjacency_list[edge_it.first].erase(edge_it.second);
    adjacency_list[edge_it.second].erase(edge_it.first);
  }
  return ( !triplets.empty() );
}

/**
* @brief Return triplets contained in the graph build from IterablePairs
* @param pairs Graph pairs
* @return List of triplet found in graph
*/
template <typename IterablePairs>
static std::vector<graph::Triplet> TripletListing
(
  const IterablePairs & pairs
)
{
  std::vector<graph::Triplet> triplets;
  graph::ListTriplets( pairs, triplets );
  return triplets;
}

} // namespace graph
} // namespace openMVG

#endif // OPENMVG_GRAPH_GRAPH_TRIPLET_FINDER_HPP
