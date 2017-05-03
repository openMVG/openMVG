// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_REGISTRATION_ICP_HPP
#define OPENMVG_GEOMETRY_REGISTRATION_ICP_HPP

#include "openMVG/geometry/kd_tree_3d.hpp"
#include "openMVG/geometry/registration/rigid_motion_3d_3d_estimation.hpp"
#include "openMVG/numeric/numeric.h"

#include <flann/flann.hpp>

#include <memory>
#include <numeric>
#include <random>

namespace openMVG
{
namespace geometry
{
namespace registration
{

/**
* @brief Apply a rigid transformation on a given set
* @param[in,out] data a data point
* @param q Quaternion
* @param t Translation
*/
template <typename Scalar>
void Transform( Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &data, const Mat4 &tra )
{
  for ( int id_pt = 0; id_pt < data.rows(); ++id_pt )
  {
    const Scalar x = data( id_pt, 0 );
    const Scalar y = data( id_pt, 1 );
    const Scalar z = data( id_pt, 2 );

    Vec3 pt_tra( x * tra( 0, 0 ) + y * tra( 0, 1 ) + z * tra( 0, 2 ) + tra( 0, 3 ),
                 x * tra( 1, 0 ) + y * tra( 1, 1 ) + z * tra( 1, 2 ) + tra( 1, 3 ),
                 x * tra( 2, 0 ) + y * tra( 2, 1 ) + z * tra( 2, 2 ) + tra( 2, 3 ) );
    //        const Scalar w = x * tra( 3, 0 ) + y * tra( 3, 1 ) + z * tra( 3, 2 ) + tra( 3, 3 );

    data( id_pt, 0 ) = pt_tra[ 0 ]; // / w;
    data( id_pt, 1 ) = pt_tra[ 1 ]; // / w;
    data( id_pt, 2 ) = pt_tra[ 2 ]; // / w;
  }
}

/**
* @brief Apply a rigid transformation on a given set
* @param[in,out] data a data point
* @param q Quaternion
* @param t Translation
*/
template<typename Scalar>
static inline void Transform( Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> & data ,
                              const Eigen::Quaternion<Scalar> & q ,
                              const Eigen::Matrix<Scalar, 3, 1> & t )
{
  for( int id_pt = 0 ; id_pt < data.rows() ; ++id_pt )
  {
    //    const Eigen::Matrix<Scalar, 3, 1> cur_pt = data.row( id_pt ).transpose() ;
    //    cur_pt << data( id_pt , 0 ) , data( id_pt , 1 ) , data( id_pt , 2 ) ;
    data.row( id_pt ) = ( ( q * data.row( id_pt ).transpose() ) + t ).transpose() ;
  }
}

/**
* @brief Compute mean square error between two set of points
* @param target target point list
* @param data source point list
* @param corresp pair (index, corresp[index]) for error computation
* @retval mse if number of pair is > 0
* @retval +inf if number of pair == 0
*/
template< typename Scalar>
static inline Scalar ComputeMSE( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &target,
                                 const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &data ,
                                 const std::vector<int> & corresp )
{
  int nb_valid = 0 ;
  Scalar mse = 0 ;
  for ( int id_pt = 0; id_pt < data.rows(); ++id_pt )
  {
    const int corr = corresp[ id_pt ];
    if ( corr >= 0 )
    {
      // Add distance between the two points
      mse += ( data.row( id_pt ) - target.row( corr ) ) .squaredNorm() ;
      ++nb_valid;
    }
  }
  if( nb_valid > 0 )
  {
    return mse / static_cast<Scalar>( nb_valid );
  }
  else
  {
    return std::numeric_limits<Scalar>::max() ;
  }
}

/**
* @brief Compute a random subset of a list
* @param highest_value Maximum index for the output
* @param percentage Percentage of value to keep in the set [0;highest_value]
* @param[out] samples Output samples
* @param rng Random number generator
* @todo : move to numeric ?
*/
static inline void RandomSubset( const size_t highest_value , // included
                                 const double percentage , // between 0 and 1
                                 std::vector< int > & samples , // selected samples
                                 std::mt19937_64 & rng ) // random generator
{
  std::vector< int > res( highest_value ) ;
  std::iota( res.begin() , res.end() , 0 ) ;
  std::shuffle( res.begin() , res.end() , rng ) ;

  const int nb_values = std::min( static_cast<int>( highest_value ) , static_cast<int>( std::ceil( static_cast<double>( highest_value ) * percentage ) ) ) ;
  samples.resize( nb_values ) ;
  std::copy( res.begin() , res.begin() + nb_values , samples.begin() ) ;
}

/**
* @brief Compute standard deviation of a given vector
* @param v a vector
* @return standard deviation of the vector
* @todo : move to numeric ?
* @note : this assume that at least one value is in the vector
*/
template< typename Scalar >
static inline Scalar StdDev( const std::vector< Scalar > & v )
{
  const double sum = std::accumulate( v.begin(), v.end(), Scalar( 0 ) );
  const double mean = sum / v.size();

  std::vector<double> diff( v.size() );
  std::transform( v.begin(), v.end(), diff.begin(),
                  std::bind2nd( std::minus<double>(), mean ) );
  const double sq_sum = std::inner_product( diff.begin(), diff.end(), diff.begin(), Scalar( 0 ) );
  return std::sqrt( sq_sum / v.size() );
}


/**
* @brief estimate transformation between two set of points
* @param target Target point
* @param data data point
* @param corresp list of correspondonce (index, corresp[index]) map data to target
* @param use_ceres indicate if ceres will be used for estimation
*/
template< typename Scalar >
Mat4 EstimateTransformation( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &target ,
                             const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &data ,
                             const std::vector<int> &corresp ,
                             const bool use_ceres = false )
{
  if( use_ceres )
  {
    RigidMotion3d3dEstimation<Scalar> estimator ;
    std::pair< Eigen::Quaternion<Scalar> , Eigen::Matrix<Scalar, 3, 1> > tra = estimator( target , data , corresp ) ;

    Mat4 tmp = Mat4::Identity() ;
    tmp.block( 0 , 0 , 3 , 3 ) = tra.first.toRotationMatrix().template cast<double>() ;
    tmp.block( 0 , 3 , 3 , 1 ) = tra.second.template cast<double>() ;

    return tmp ;
  }
  else
  {
    Eigen::Matrix<Scalar, 3 , Eigen::Dynamic, Eigen::RowMajor> tmp_target ;
    Eigen::Matrix<Scalar, 3 , Eigen::Dynamic, Eigen::RowMajor> tmp_data ;

    int nb_valid = 0 ;
    for( const auto & cur_corresp : corresp )
    {
      if( cur_corresp >= 0 )
      {
        ++nb_valid ;
      }
    }

    if( nb_valid == 0 )
    {
      return Mat4::Identity() ;
    }

    tmp_target.resize( 3 , nb_valid ) ;
    tmp_data.resize( 3 , nb_valid ) ;

    int id_out = 0 ;
    for( int id_pt = 0 ; id_pt < corresp.size() ; ++id_pt )
    {
      const int cur_corresp = corresp[ id_pt ] ;
      if( cur_corresp >= 0 )
      {
        tmp_target.col( id_out ) = target.row( cur_corresp ).transpose() ;
        tmp_data.col( id_out ) = data.row( id_pt ).transpose() ;
        ++id_out ;
      }
    }
    
    const Eigen::Matrix<Scalar,4,4> res = Eigen::umeyama( tmp_data , tmp_target , false ) ;
    return res.template cast<double>() ; 
  }
}


/**
* @brief Given two sets of points: target and data
* This function computes rigid transformation that maps model on target minimizing MSE distance
* @param target Target shape : a matrix of 3d points (one point per row)
* @param data data shape : a matrix of 3d points (one point per row)
* @param nb_iteration Maximum number of iteration
* @param mse_threshold Threshold use to stop computation
* @param[out] t Translation transformation (3d vector)
* @param[out] R Rotation transformation (3x3 rotation matrix)
*/
template <typename Scalar>
void ICP( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &target,
          const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &data_,
          const unsigned long max_nb_iteration,
          const Scalar mse_threshold,
          openMVG::Vec3 &t,
          openMVG::Mat3 &R )
{

  // Build Kd-Tree
  KDTree3d<Scalar> tree( target );

  unsigned long id_iteration = 0;
  Scalar cur_mse             = std::numeric_limits<Scalar>::max();

  // Working sample
  Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> data = data_;

  Mat4 final_tra = Mat4::Identity() ;

  std::vector<int> corresp( data.rows() );
  std::vector<Scalar> distance( data.rows() );

  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> indices;
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dists;

  std::vector< int > subset_indice ;
  const double percentage = 0.10 ; // only sample 10% of the subset
  std::mt19937_64 rng ;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> subset_data ;

  while ( id_iteration < max_nb_iteration && cur_mse > mse_threshold )
  {
    std::fill( corresp.begin(), corresp.end(), -1 );
    std::fill( distance.begin(), distance.end(), Scalar( -1.0 ) );

    // 1 - Pick a random subset
    RandomSubset( data.rows() - 1 , percentage , subset_indice , rng ) ;
    const int nb_subset_elts = std::min( static_cast<int>( data.rows() ) - 1 , static_cast<int>( subset_indice.size() ) ) ;
    subset_data.resize( nb_subset_elts , 3 ) ;
    for( size_t id_sample = 0 ; id_sample < subset_indice.size() ; ++id_sample )
    {
      const int indice = subset_indice[ id_sample ] ;
      subset_data.row( id_sample ) = data.row( indice ) ;
    }

    // 2 - Establish pairs based on nearest neighbor search
    tree.Search( subset_data , 1 , indices, dists );
    std::vector< Scalar > compact_dist( subset_data.rows() ) ;
    int id_dist = 0 ;
    for ( int id_pt = 0; id_pt < subset_data.rows(); ++id_pt )
    {
      const int id_point = indices( id_pt, 0 );

      if ( ( id_point >= 0 ) &&
           ( id_point < target.rows() ) )
      {
        distance[ id_pt ] = dists( id_pt, 0 );
        corresp[ id_pt ]  = id_point;
        compact_dist[ id_dist ] = dists( id_pt, 0 );

        ++id_dist ;
      }
    }
    compact_dist.resize( id_dist ) ;

    // 3 - Filter points based on 3 * stddev
    const Scalar t = 3.0 * StdDev( compact_dist ) ;
    for ( int id_pt = 0; id_pt < subset_data.rows(); ++id_pt )
    {
      if ( corresp[ id_pt ] >= 0 )
      {
        if ( distance[ id_pt ] > t )
        {
          corresp[ id_pt ]  = -1;
          distance[ id_pt ] = -1.0;
        }
      }
    }
    const double mse_before = ComputeMSE( target , subset_data , corresp ) ;
    if( id_iteration == 0 )
    {
      cur_mse = mse_before ;
    }

    // 4 - Compute best rigid transformation based on pairs
    const Mat4 tra = EstimateTransformation( target , subset_data , corresp ) ;

    // 5 - Update data points and final transformation
    Transform( subset_data , tra ) ; //.first , tra.second );
    const double mse_after = ComputeMSE( target , subset_data , corresp ) ;
    if( mse_after < mse_before )
    {
      // Update the whole set
      Transform( data , tra ) ; //.first , tra.second );

      // Update global transformation
      final_tra = tra * final_tra ;
      // Update mse
      cur_mse = mse_after ;
    }

    ++id_iteration;
  }

  // Compute final transformation
  R = final_tra.block( 0 , 0 , 3 , 3 ) ;
  t = final_tra.block( 0 , 3 , 3 , 1 ) ;
}


/**
* @brief Given two sets of points: target and data
* This function computes rigid transformation that maps model on target minimizing MSE distance
* This use the point to normal distance for computation
* @param target Target shape : a matrix of 3d points (one point per row)
* @param target_n Target shape normals : a matrix of 3d vectors (one vector per row)
* @param data data shape : a matrix of 3d points (one point per row)
* @param nb_iteration Maximum number of iteration
* @param mse_threshold Threshold use to stop computation
* @param[out] t Translation transformation (3d vector)
* @param[out] R Rotation transformation (3x3 rotation matrix)
*/
template <typename Scalar>
void ICP( const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &target,
          const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &target_n,
          const Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> &data_,
          const unsigned long max_nb_iteration,
          const Scalar mse_threshold,
          openMVG::Vec3 &t,
          openMVG::Mat3 &R )
{
  // Build Kd-Tree
  KDTree3d<Scalar> tree( target );

  unsigned long id_iteration = 0;
  Scalar cur_mse             = std::numeric_limits<Scalar>::max();

  // Working sample
  Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> data = data_;

  Mat4 final_tra = Mat4::Identity() ;

  std::vector<int> corresp( data.rows() );
  std::vector<Scalar> distance( data.rows() );

  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> indices;
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dists;


  std::vector< int > subset_indice ;
  const double percentage = 0.10 ; // only sample 10% of the subset
  std::mt19937_64 rng ;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> subset_data ;

  while ( id_iteration < max_nb_iteration && cur_mse > mse_threshold )
  {
    std::fill( corresp.begin(), corresp.end(), -1 );
    std::fill( distance.begin(), distance.end(), Scalar( -1.0 ) );

    // 1 - Pick a random subset
    RandomSubset( data.rows() - 1 , percentage , subset_indice , rng ) ;
    const int nb_subset_elts = std::min( static_cast<int>( data.rows() ) - 1 , static_cast<int>( subset_indice.size() ) ) ;
    subset_data.resize( nb_subset_elts , 3 ) ;
    for( size_t id_sample = 0 ; id_sample < subset_indice.size() ; ++id_sample )
    {
      const int indice = subset_indice[ id_sample ] ;
      subset_data.row( id_sample ) = data.row( indice ) ;
    }

    // 2 - Establish pairs based on nearest neighbor search
    tree.Search( subset_data , 1 , indices, dists );
    std::vector< Scalar > compact_dist( subset_data.rows() ) ;
    int id_dist = 0 ;
    for ( int id_pt = 0; id_pt < subset_data.rows(); ++id_pt )
    {
      const int id_point = indices( id_pt, 0 );

      if ( ( id_point >= 0 ) &&
           ( id_point < target.rows() ) )
      {
        distance[ id_pt ] = dists( id_pt, 0 );
        corresp[ id_pt ]  = id_point;
        compact_dist[ id_dist ] = dists( id_pt, 0 );

        ++id_dist ;
      }
    }
    compact_dist.resize( id_dist ) ;

    // 3 - Filter points based on 3 * stddev
    const Scalar t = 3.0 * StdDev( compact_dist ) ;
    for ( int id_pt = 0; id_pt < subset_data.rows(); ++id_pt )
    {
      if ( corresp[ id_pt ] >= 0 )
      {
        if ( distance[ id_pt ] > t )
        {
          corresp[ id_pt ]  = -1;
          distance[ id_pt ] = -1.0;
        }
      }
    }
    const double mse_before = ComputeMSE( target , subset_data , corresp ) ;
    if( id_iteration == 0 )
    {
      cur_mse = mse_before ;
    }

    // 4 - Compute best rigid transformation based on pairs
    RigidMotion3d3dEstimation<Scalar> estimator ;
    std::pair< Eigen::Quaternion<Scalar> , Eigen::Matrix<Scalar, 3, 1> > tra = estimator( target , target_n , subset_data , corresp ) ;

    // 5 - Update data points and final transformation
    Transform( subset_data , tra.first , tra.second );
    const double mse_after = ComputeMSE( target , subset_data , corresp ) ;
    if( mse_after < mse_before )
    {
      // Update whole set
      Transform( data , tra.first , tra.second );

      // Update global transformation
      Mat4 tmp = Mat4::Identity() ;
      tmp.block( 0 , 0 , 3 , 3 ) = tra.first.toRotationMatrix() ;
      tmp.block( 0 , 3 , 3 , 1 ) = tra.second ;
      final_tra = tmp * final_tra ;

      // Update mse
      cur_mse = mse_after ;
    }

    ++id_iteration;
  }

  // Compute final transformation
  R = final_tra.block( 0 , 0 , 3 , 3 ) ;
  t = final_tra.block( 0 , 3 , 3 , 1 ) ;
}



} // namespace registration
} // namespace geometry
} // namespace openMVG

#endif