/**
 * This file contains interfaces to the hypergraph partitioning libraries
 * that we use for hypersfm. 
 * contains : 
 *  - interface to scotch lib
 *  - interface to hmetis lib
 */

#pragma once

// scotch library for hypergraph partitioning
#include <scotch.h>
#include <fstream>

// get all permutations between the values of a set, as a set of pairs <x1,x2>, where x2 > x1
std::set<std::pair<unsigned,unsigned>> getPermutations(const std::set<unsigned> & id_set)
{
  std::set<std::pair<unsigned,unsigned>> permutations_set;
  for (std::set<unsigned>::const_iterator id = id_set.begin(); id != id_set.end(); id++)
  {
    for (std::set<unsigned>::const_iterator id_2 = std::next(id,1); id_2 != id_set.end(); id_2++)
      permutations_set.insert(std::make_pair(*id, *id_2));
  }
  return permutations_set;
}

/**
 * @note Alternative to ScotchPartitionHyperGraph function which cause no memory leak. Here we don't
 * build a mesh (hypergraph) but directly a graph for scotch to partition SHIT RESULTS THO
 * @brief Use the scotch library to partition a hypergraph made of view_id nodes and tracks set hyperedges 
 * into two sets of nodes (view ids in that case)
 * @param hyper_graph The hypergraph to be partitioned, defined as a map of <view_ids,tracks_ids>. 
 * @param view_id_partitions The returned partition, in the shape of a pair of view_ids sets.
 * @retval true If partitioning succeeded
 * @retval false If partitioning failed
 */
inline bool ScotchPartitionHyperGraphNoMesh(
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
  {
    view_ids_all.insert(h_edge.first.begin(), h_edge.first.end());
  }

  std::set<std::pair<unsigned,unsigned>> all_edges;
  for (const auto & h_edge : hyper_graph)
  {
    std::set<std::pair<unsigned,unsigned>> new_perms = getPermutations(h_edge.first);
    all_edges.insert(new_perms.begin(), new_perms.end());
  }
  std::cout << " LKAJDSLKSAJF " << all_edges.size() << std::endl;


  SCOTCH_Graph grafdat;
  SCOTCH_graphInit(&grafdat);

  SCOTCH_Num baseval = 0; // array base number
  SCOTCH_Num vertnbr = view_ids_all.size(); // number of vertices
  SCOTCH_Num verttab[vertnbr + 1]; // adjacency index array
  SCOTCH_Num * vendtab = verttab + 1;// end thing
  SCOTCH_Num velotab[vertnbr];// vertex load array (TODO we don't use this ?)
  SCOTCH_Num edgenbr = 2*all_edges.size();//TODO compute edge_tot // number of arcs
  SCOTCH_Num edgetab[edgenbr]; // adjacency array
  SCOTCH_Num edlotab[edgenbr];// arc load array

  // fill load arrays with zeros
  std::fill(edlotab, edlotab + edgenbr, 0);
  std::fill(velotab, velotab + vertnbr, 0);

  int node_id(baseval);
  std::map<openMVG::IndexT, int> map_view_id_node_id;
  std::map<int, openMVG::IndexT> map_node_id_view_id;
  for ( const auto & view_id : view_ids_all)
  {
    map_view_id_node_id[view_id] = node_id;
    map_node_id_view_id[node_id] = view_id;
    node_id++;
  }

  int k(0), l(0);
  for (const auto & view_id : view_ids_all)
  {
    verttab[l] = k;
    std::set<int> paired_indices;
    velotab[l] += 1;
    std::map<int, int> edge_loads;
    for (const auto & h_edge : hyper_graph)
    {
      // if view is not in hyper edge-> try next
      if (h_edge.first.count(view_id) == 0)
        continue;

      for (const auto & id : h_edge.first)
      {
        if (id != view_id)
        {
          const auto & n_id = map_view_id_node_id.at(id);
          paired_indices.insert(n_id);
          if (edge_loads.find(n_id) == edge_loads.end())
            edge_loads[n_id] = 1;
          else
            edge_loads[n_id] += 1;
        }
      }
    }
    for (const auto & ind : paired_indices)
    {
      edgetab[k] = ind;
      //edlotab[k] += 1;
      edlotab[k] = edge_loads.at(ind);
      k++;
    }
    l++;
  }
  verttab[l] = edgenbr;
  velotab[l] += 1;

  //SCOTCH_graphBuild(&grafdat, baseval, vertnbr, verttab, vendtab, velotab, vlbltab, edgenbr, edgetab, edlotab);
  //SCOTCH_graphBuild(&grafdat, baseval, vertnbr, verttab, vendtab, velotab, NULL, edgenbr, edgetab, edlotab);
  SCOTCH_graphBuild(&grafdat, baseval, vertnbr, verttab, vendtab, velotab, NULL, edgenbr, edgetab, NULL);
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
  SCOTCH_Num parttab[edgenbr]; 
  for (auto & part_element : parttab)
    part_element = -1;

  // partition in 2
  if (SCOTCH_graphPart(&grafdat, 2, stratptr, parttab) == 0)
    std::cout << "graph partitioned !" << std::endl;

  // partition data the view ids
  int current_node_id(baseval);
  for (const auto & part_element : parttab)
  {
    if (part_element == 0)
      view_id_partitions.first.insert(map_node_id_view_id.at(current_node_id));
    if (part_element == 1)
      view_id_partitions.second.insert(map_node_id_view_id.at(current_node_id));
    current_node_id++;
  }
  
  // free memory
  SCOTCH_graphExit(&grafdat);
  SCOTCH_stratExit(stratptr);
  SCOTCH_memFree(stratptr);

  return true;
}

/**
 * @note WARNING : this function USED TO create a memory leak (the incriminated part is SCOTCH_meshGraph())
 * you have to use our version of the scotch library which corrects this memory leak (in third_party/scotch_6.0.4).
 * otherwise, use the ScotchPartitionHyperGraphNoMesh function 
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
  //SCOTCH_Num vlbltab[velmnbr + vnodnbr]; // vertex label array ... do we use this ?
  SCOTCH_Num edgenbr = 2*edge_number_total; // to fill, should be 2*number of edges (multiple edge per h_edge) TODO : << THIS !!  
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
  delete edgetab;
  delete verttab;
  delete velotab;

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
