// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef _OPENMVG_CLUSTERING_K_MEANS_HPP_
#define _OPENMVG_CLUSTERING_K_MEANS_HPP_

#include "openMVG/clustering/kmeans_trait.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

#include <limits>
#include <random>
#include <vector>

namespace openMVG
{
namespace clustering
{

/**
* @brief Kind of initialization of the kmeans centers
*/
enum KMeansInitType
{
  KMEANS_INIT_RANDOM , /* Standard Llyod algoritm */
  KMEANS_INIT_PP ,     /* Kmeans++ initialisation */
} ;

/**
* @brief Compute minimum distance to any center
* @param pts Input points
* @param centers Centers
* @param[out] dists computed (minimum) distence to any center
*/
template< typename DataType >
void MinimumDistanceToAnyCenter( const std::vector< DataType > & pts ,
                                 const std::vector< DataType > & centers ,
                                 std::vector< typename KMeansVectorDataTrait<DataType>::scalar_type > & dists )
{
  using trait = KMeansVectorDataTrait<DataType>;

  dists.resize( pts.size() , std::numeric_limits<typename trait::scalar_type>::max() ) ;

  for( size_t id_pt = 0 ; id_pt < pts.size() ; ++id_pt )
  {
    const auto & pt = pts[ id_pt ] ;
    for( const auto & c : centers )
    {
      const typename trait::scalar_type cur_d = trait::L2( pt , c ) ;
      dists[ id_pt ] = std::min( dists[ id_pt ] , cur_d ) ;
    }
  }
}

/**
* @brief Compute simple kmeans clustering on specified data
* @param source_data Input data
* @param[out] cluster_assignment index for each point in the input set to a specified cluster
* @param[out] centers Centers of the clusters
* @param nb_cluster requested number of cluster in the output
* @param max_nb_iteration maximum number of iteration to do for clustering
* @note This is the standard llyod algorithm
*/
template< typename DataType >
void KMeans( const std::vector< DataType > & source_data ,
             std::vector< uint32_t > & cluster_assignment ,
             std::vector< DataType > & centers ,
             const uint32_t nb_cluster ,
             const uint32_t max_nb_iteration = std::numeric_limits<uint32_t>::max() ,
             const KMeansInitType init_type = KMEANS_INIT_PP )
{
  if( source_data.size() == 0 )
  {
    return ;
  }

  using trait = KMeansVectorDataTrait<DataType>;

  // 1 - Compute range of the input
  DataType min , max ;
  trait::minMax( source_data , min , max ) ;

  std::mt19937_64 rng ;

  // 2 - init center of mass
  if( init_type == KMEANS_INIT_PP )
  {
    // Kmeans++ init :
    // first one is a random one
    // the others based on the importance probability (Di / \sum_i Di) where :
    // Di is the minimum distance to any created centers already created
    std::uniform_int_distribution<size_t> distrib_first( 0 , source_data.size() - 1 ) ;
    centers.emplace_back( source_data[ distrib_first( rng ) ] );

    std::vector< typename trait::scalar_type > dists ;

    for( size_t id_center = 1 ; id_center < nb_cluster ; ++id_center )
    {
      // Compute Di / \sum Di pdf
      MinimumDistanceToAnyCenter( source_data , centers , dists ) ;
      std::discrete_distribution<size_t> distrib_c( dists.begin() , dists.end() );

      // Sample a point from this distribution
      centers.emplace_back( source_data[distrib_c( rng )] ) ;
    }
  }
  else
  {
    // Standard Llyod init
    for( uint32_t id_center = 0 ; id_center < nb_cluster ; ++id_center )
    {
      centers.emplace_back( trait::random( min , max , rng ) ) ;
    }
  }

  // Assign all element to the first one
  cluster_assignment.resize( source_data.size() , nb_cluster ) ;


  // 3 - Perform kmeans
  bool changed ;
  uint32_t id_iteration = 0 ;

  do
  {
    changed = false ;

    // 3.1 affect center to each points
    for( size_t id_pt = 0 ; id_pt < source_data.size() ; ++id_pt )
    {
      const DataType & cur_pt = source_data[id_pt] ;
      typename trait::scalar_type min_dist = std::numeric_limits<typename trait::scalar_type>::max() ;
      uint32_t nearest_center = nb_cluster ;

      for( uint32_t cur_center = 0 ; cur_center < nb_cluster ; ++cur_center )
      {
        const typename trait::scalar_type cur_dist = trait::L2( cur_pt , centers[ cur_center ] ) ;
        if( cur_dist < min_dist )
        {
          min_dist = cur_dist ;
          nearest_center = cur_center ;
        }
      }

      // Assign id to the nearest center
      if( cluster_assignment[id_pt] != nearest_center )
      {
        cluster_assignment[id_pt] = nearest_center ;
        changed = true ;
      }
    }

    // 3.2 Compute new centers of mass
    std::vector< DataType > new_centers( nb_cluster , trait::null( centers[0] ) ) ;
    std::vector< uint32_t > nb_per_center( nb_cluster , 0 ) ;

    for( size_t id_pt = 0 ; id_pt < source_data.size() ; ++id_pt )
    {
      const uint32_t id_center = cluster_assignment[id_pt] ;
      trait::accumulate( new_centers[ id_center ] , source_data[id_pt] ) ;
      ++nb_per_center[id_center] ;
    }
    for( uint32_t id_center = 0 ; id_center < nb_cluster ; ++id_center )
    {
      if( nb_per_center[id_center] == 0 )
      {
        new_centers[id_center] = trait::random( min , max , rng ) ;
      }
      else
      {
        trait::divide( new_centers[id_center] , nb_per_center[id_center] ) ;
      }
    }

    centers = new_centers ;
    ++id_iteration ;
  }
  while( changed && id_iteration < max_nb_iteration ) ;
}

} // namespace clustering
} // namespace openMVG

#endif