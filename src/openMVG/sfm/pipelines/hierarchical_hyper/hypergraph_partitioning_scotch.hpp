// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/**
 * This file contains interfaces to the SCOTCH hypergraph partitioning library
 * that we use for hypersfm.
 */

#pragma once

#include "openMVG/types.hpp"
#include <fstream>
#include <iostream>

// scotch library for hypergraph partitioning
#include "scotch.h"

/**
 * @note WARNING : this function USED TO create a memory leak (the incriminated part is SCOTCH_meshGraph() in
 * the scotch library)
 * you have to use our version of the scotch library which corrects this memory leak (in third_party/scotch_6.0.4).
 * @brief Use the scotch library to partition a hypergraph made of view_id nodes and tracks set hyperedges
 * into two sets of nodes (view ids in that case)
 * @param hyper_graph The hypergraph to be partitioned, defined as a map of <view_ids,tracks_ids>.
 * @param view_id_partitions The returned partition, in the shape of a pair of view_ids sets.
 * @retval true If partitioning succeeded
 * @retval false If partitioning failed
 */
inline bool ScotchPartitionHyperGraph(
    const std::map<std::set<openMVG::IndexT>, std::set<size_t>> & hyper_graph,
    std::pair<std::set<openMVG::IndexT>, std::set<openMVG::IndexT>> & view_id_partitions)
{
  // check hyper graph has been computed
  if (hyper_graph.size() == 0)
  {
    std::cerr << "empty hyper graph ! " << std::endl;
    return false;
  }

  // regroup all view indices of the hypergraph
  std::set<openMVG::IndexT> view_ids_all;
  for (const auto & h_edge : hyper_graph)
    view_ids_all.insert(h_edge.first.begin(), h_edge.first.end());

  // find total number of edges
  int edge_number_total(0);
  for (const auto & h_edge : hyper_graph)
    edge_number_total += h_edge.first.size();

  // the mesh data container
  SCOTCH_Mesh meshdat;

  // note : SCOTCH_Num is just an int
  // note2 : these obscure variable names come from the scotch documentation...
  SCOTCH_Num velmbas = 0;
  SCOTCH_Num vnodbas = velmbas + hyper_graph.size();
  SCOTCH_Num velmnbr = hyper_graph.size();
  SCOTCH_Num vnodnbr = view_ids_all.size();
  SCOTCH_Num * verttab = new SCOTCH_Num[velmnbr + vnodnbr + 1]; // array of start indices in edgetab of vertex
  SCOTCH_Num * vendtab = verttab + 1; // this points to verttab + 1
  SCOTCH_Num * velotab = new SCOTCH_Num[velmnbr]; // element weights
  //SCOTCH_Num vnlotab[vnodnbr]; // node weights NOTE : not used because no weights on the nodes in our case
  //SCOTCH_Num vlbltab[velmnbr + vnodnbr]; // vertex label array ... we don't use this
  SCOTCH_Num edgenbr = 2*edge_number_total; // should be 2*number of edges
  SCOTCH_Num * edgetab = new SCOTCH_Num[edgenbr]; // adjacency array

  // map hyper edges with element ids
  int el_id(velmbas);
  std::map<int, std::set<openMVG::IndexT>> map_el_id_h_edge;
  for (const auto & h_edge : hyper_graph)
  {
    map_el_id_h_edge[el_id] = h_edge.first;
    el_id++;
  }

  // map view ids with node ids
  int node_id(vnodbas);
  std::map<openMVG::IndexT, int> map_view_id_node_id;
  std::map<int, openMVG::IndexT> map_node_id_view_id;
  for ( const auto & view_id : view_ids_all)
  {
    map_view_id_node_id[view_id] = node_id;
    map_node_id_view_id[node_id] = view_id;
    node_id++;
  }

  // each hyper edge is set of view ids. the weight of it is the size of
  // its tracks ids.
  int i(0), j(0);
  for (const auto & h_edge : hyper_graph)
  {
    const auto & edge_view_ids = h_edge.first;
    const auto & weight = h_edge.second.size();

    verttab[j] = i;
    velotab[j] = weight;
    for (const auto & view_id : edge_view_ids)
    {
      edgetab[i] = map_view_id_node_id.at(view_id);
      i++;
    }
    j++;
  }
  
  for (const auto & view_id : view_ids_all)
  {
    const int current_node_id = map_view_id_node_id[view_id];

    verttab[current_node_id] = i;
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
    for (const auto & el_id_h_edge_pair : map_el_id_h_edge)
    {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
      {
        const int element_id = el_id_h_edge_pair.first;
        const std::set<openMVG::IndexT> & h_edge = el_id_h_edge_pair.second;

        // if view is part of the hyper edge
        if (h_edge.find(view_id) != h_edge.end())
        {
#ifdef OPENMVG_USE_OPENMP
    #pragma omp critical
#endif
          {
            edgetab[i] = element_id;
            i++;
          }
        }
      }// omp section
    }
  }
  // last index for verttab
  verttab[velmnbr + vnodnbr] = i;

  // build the mesh (hypergraph equivalent in scotch)
  SCOTCH_meshInit(&meshdat);

  if (SCOTCH_meshBuild(&meshdat, velmbas, vnodbas, velmnbr, vnodnbr, verttab,
      vendtab, velotab, NULL, NULL, edgenbr, edgetab) == 0)
    std::cout << "built mesh (hypergraph)!" << std::endl;

  std::cout << "number of edges : " << edgenbr << std::endl
    << "number of nodes : " << vnodnbr << std::endl
    << "number of elements : " << velmnbr << std::endl;

  SCOTCH_Graph grafdat;
  SCOTCH_graphInit(&grafdat);

  // convert mesh to graph
  if (SCOTCH_meshGraph(&meshdat, &grafdat) == 0)
    std::cout << "mesh converted to graph !" << std::endl;
  std::cout << "check grafdat  : " << SCOTCH_graphCheck(&grafdat) << std::endl;

  SCOTCH_Num nvert, nedges;
  SCOTCH_graphSize(&grafdat, &nvert, &nedges);
  std::cout << "graph size : " << std::endl
    << "vertices : " << nvert << std::endl
    << "edges : " << nedges << std::endl;

  // define default strategy for partitioning
  SCOTCH_Strat * stratptr = SCOTCH_stratAlloc();
  SCOTCH_stratInit(stratptr);

  // array for partitioned data
  SCOTCH_Num parttab[vnodnbr];
  for (auto & part_element : parttab)
    part_element = -1;

  // partition in 2
  if (SCOTCH_graphPart(&grafdat, 2, stratptr, parttab) == 0)
    std::cout << "graph partitioned !" << std::endl;

  // free memory
  SCOTCH_graphSize(&grafdat, &nvert, &nedges);
  std::cout << "graph size : " << std::endl
    << "vertices : " << nvert << std::endl
    << "edges : " << nedges << std::endl;
  SCOTCH_graphFree(&grafdat);
  SCOTCH_graphExit(&grafdat);
  SCOTCH_graphSize(&grafdat, &nvert, &nedges);
  std::cout << "graph size : " << std::endl
    << "vertices : " << nvert << std::endl
    << "edges : " << nedges << std::endl;
  SCOTCH_meshExit(&meshdat);
  SCOTCH_stratExit(stratptr);
  SCOTCH_memFree(stratptr);

  // cannot free those before since they are used by meshdat
  delete [] edgetab;
  delete [] verttab;
  delete [] velotab;

  // partition data the view ids
  int current_node_id(vnodbas);
  for (const auto & part_element : parttab)
  {
    if (part_element == 0)
      view_id_partitions.first.insert(map_node_id_view_id.at(current_node_id));
    if (part_element == 1)
      view_id_partitions.second.insert(map_node_id_view_id.at(current_node_id));
    current_node_id++;
  }


  return true;
}
