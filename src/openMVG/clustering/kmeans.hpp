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
enum class KMeansInitType
{
  KMEANS_INIT_RANDOM, /* Standard Llyod algoritm */
  KMEANS_INIT_PP,     /* Kmeans++ initialization */
};

/**
* @brief Compute minimum distance to any center
* @param pts Input points
* @param centers Centers
* @param[out] dists computed (minimum) distance to any center
*/
template< typename DataType >
void MinimumDistanceToAnyCenter( const std::vector< DataType > & pts,
                                 const std::vector< DataType > & centers,
                                 std::vector< typename KMeansVectorDataTrait<DataType>::scalar_type > & dists )
{
  using trait = KMeansVectorDataTrait<DataType>;

  dists.resize( pts.size(), std::numeric_limits<typename trait::scalar_type>::max() );

  #pragma omp parallel for
  for( int id_pt = 0; id_pt < static_cast<int>(pts.size()); ++id_pt )
  {
    const auto & pt = pts[ id_pt ];
    for( const auto & c : centers )
    {
      const typename trait::scalar_type cur_d = trait::L2( pt, c );
      dists[ id_pt ] = std::min( dists[ id_pt ], cur_d );
    }
  }
}

/**
* @brief Compute Nearest center Id of a given point
* @param pt Query point
* @param centers list of test centers
* @return id of the nearest center (0-based)
*/
template< typename DataType >
uint32_t NearestCenterID( const DataType & pt,
                          const std::vector< DataType > & centers )
{
  using trait = KMeansVectorDataTrait<DataType>;
  const uint32_t nb_cluster = static_cast<uint32_t>( centers.size() );

  typename trait::scalar_type min_dist = std::numeric_limits<typename trait::scalar_type>::max();
  uint32_t nearest_center = nb_cluster;

  for( uint32_t cur_center = 0; cur_center < nb_cluster; ++cur_center )
  {
    const typename trait::scalar_type cur_dist = trait::L2( pt, centers[ cur_center ] );
    if( cur_dist < min_dist )
    {
      min_dist = cur_dist;
      nearest_center = cur_center;
    }
  }
  return nearest_center;
}

/**
* @brief Compute center of mass of a set a points
* @param pts List of points
* @param assigned_center Id of the center to be affected to a given point
* @param nb_center Number of center of mass in the result
* @return New centers of mass
*/
template< typename DataType >
std::vector< DataType > ComputeCenterOfMass( const std::vector< DataType > & pts,
    const std::vector< uint32_t > & assigned_center,
    const uint32_t nb_center )
{
  using trait = KMeansVectorDataTrait<DataType>;

  std::vector< DataType > new_centers( nb_center, trait::null( pts[0] ) );
  std::vector< uint32_t > nb_per_center( nb_center, 0 );

  // Affect points to centers
  for( size_t id_pt = 0; id_pt < pts.size(); ++id_pt )
  {
    const uint32_t id_center = assigned_center[id_pt];
    trait::accumulate( new_centers[ id_center ], pts[id_pt] );
    ++nb_per_center[id_center];
  }

  // Compute mean of centers based on the number of points affected to each centers
  #pragma omp parallel for
  for( int id_center = 0; id_center < static_cast<int>(nb_center); ++id_center )
  {
    trait::divide( new_centers[id_center], nb_per_center[id_center] );
  }

  return new_centers;
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
void KMeans( const std::vector< DataType > & source_data,
             std::vector< uint32_t > & cluster_assignment,
             std::vector< DataType > & centers,
             const uint32_t nb_cluster,
             const uint32_t max_nb_iteration = std::numeric_limits<uint32_t>::max(),
             const KMeansInitType init_type = KMeansInitType::KMEANS_INIT_PP )
{
  if( source_data.size() == 0 )
  {
    return;
  }

  using trait = KMeansVectorDataTrait<DataType>;

  std::mt19937_64 rng(std::mt19937_64::default_seed);


  // 1 - init center of mass
  if( init_type == KMeansInitType::KMEANS_INIT_PP )
  {
    // Kmeans++ init:
    // first one is a random one
    // the others based on the importance probability (Di / \sum_i Di) where:
    // Di is the minimum distance to any created centers already created
    std::uniform_int_distribution<size_t> distrib_first( 0, source_data.size() - 1 );
    centers.reserve( nb_cluster );
    centers.emplace_back( source_data[ distrib_first( rng ) ] );

    std::vector< typename trait::scalar_type > dists;

    for( uint32_t id_center = 1; id_center < nb_cluster; ++id_center )
    {
      // Compute Di / \sum Di pdf
      MinimumDistanceToAnyCenter( source_data, centers, dists );
      std::discrete_distribution<size_t> distrib_c( dists.cbegin(), dists.cend() );

      // Sample a point from this distribution
      centers.emplace_back( source_data[distrib_c( rng )] );
    }
  }
  else if (init_type == KMeansInitType::KMEANS_INIT_RANDOM)
  {
    DataType min, max;
    trait::minMax( source_data, min, max );

    // Standard Llyod init
    centers.resize( nb_cluster );
    std::uniform_int_distribution<size_t> distrib( 0, source_data.size() - 1 );
    for( auto & cur_center : centers  )
    {
      cur_center = source_data[distrib( rng )];
    }
  }
  else // Invalid Kmeans initialization type
  {
    return;
  }

  // Assign all element to the first center
  cluster_assignment.resize( source_data.size(), nb_cluster );

  bool changed;
  uint32_t id_iteration = 0;

  // 2 - Perform kmeans
  do
  {
    changed = false;

    // 2.1 affect center to each points
    #pragma omp parallel for shared(changed)
    for( int id_pt = 0; id_pt < static_cast<int>(source_data.size()); ++id_pt )
    {
      const DataType & cur_pt = source_data[id_pt];
      // Compute nearest center of this point
      const uint32_t nearest_center = NearestCenterID( cur_pt, centers );
      if( cluster_assignment[id_pt] != nearest_center )
      {
        cluster_assignment[id_pt] = nearest_center;
        changed = true;
      }
    }

    // 2.2 Compute new centers of mass
    centers = ComputeCenterOfMass( source_data, cluster_assignment, nb_cluster );

    ++id_iteration;
  }
  while( changed && id_iteration < max_nb_iteration );
}

} // namespace clustering
} // namespace openMVG

#endif
