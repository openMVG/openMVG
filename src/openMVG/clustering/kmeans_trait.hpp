// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef _OPENMVG_KMEANS_TRAIT_HPP_
#define _OPENMVG_KMEANS_TRAIT_HPP_

#include "openMVG/numeric/eigen_alias_definition.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <random>
#include <vector>

namespace openMVG
{
namespace clustering
{

/**
* @brief Class used to detail parts of a vector in a generic way
* @note this is only tested with floating point type
*/
template< typename VectorType >
class KMeansVectorDataTrait
{
  public:
    /// base type
    typedef VectorType type;

    /// Type of a scalar element
    typedef VectorType scalar_type;

    /**
    * @brief number of element in the vector
    * @param aVector a Vector
    * @return number of scalar element in the vector
    */
    static size_t size( const type & aVector );

    /**
    * @brief Square euclidean distance between two vectors
    * @param aVec1 first vector
    * @param aVec2 second vector
    * @return square euclidean distance
    */
    static scalar_type L2( const type & aVec1, const type & aVec2 );

    /**
    * @brief Draw a random vector in a range
    * @param min minimum bound of the sampling range
    * @param max maximum bound of the sampling range
    * @param rng A c++11 random generator
    * @return a Random vector in the given range
    */
    template < typename RngType >
    static type random( const type & min, const type & max, RngType & rng );

    /**
    * @brief Compute minimum and maximum value of a set of points
    * @params elts list of points
    * @param[out] min minimum (component-wise) of the points
    * @param[out] max maximum (component-wise) of the points
    */
    static void minMax( const std::vector<type> & elts, type & min, type & max );

    /**
    * @brief get a zero valued vector data
    * @param dummy a dummy vector
    * @return a null vector
    */
    static type null( const type & dummy );

    /**
    * @brief Accumulate value inside a vector
    * @param self vector used for accumulation
    * @param data vector to add to the self vector
    * @note this perform self += data (component-wise)
    */
    static void accumulate( type & self, const type & data );

    /**
    * @brief Scalar division
    * @param self Vector to divide
    * @param val scalar divisor
    * @note this perform self /= data (component-wise)
    */
    static void divide( type & self, const size_t val );
};

/**
* @brief overloading for std::array
*/
template< typename T, size_t N >
class KMeansVectorDataTrait<std::array<T, N>>
{
  public:
    // base type
    typedef std::array<T, N> type;

    typedef T scalar_type;

    /**
    * @brief number of element in the vector
    * @param aVector a Vector
    * @return number of scalar element in the vector
    */
    static size_t size( const type & aVector )
    {
      return N;
    }

    /**
    * @brief Square euclidean distance between two vectors
    * @param aVec1 first vector
    * @param aVec2 second vector
    * @return square euclidean distance
    */
    static scalar_type L2( const type & aVec1, const type & aVec2 )
    {
      typedef Eigen::Matrix<scalar_type, 1, Eigen::Dynamic> VecType;
      typedef Eigen::Map<const VecType> VecTypeConst;
      VecTypeConst map_aVec1(aVec1.data(), aVec1.size());
      VecTypeConst map_aVec2(aVec2.data(), aVec2.size());
      return ( map_aVec1 - map_aVec2 ).squaredNorm();
    }


    /**
    * @brief Draw a random vector in a range
    * @param min minimum bound of the sampling range
    * @param max maximum bound of the sampling range
    * @param rng A c++11 random generator
    * @return a Random vector in the given range
    */
    template < typename RngType >
    static type random( const type & min, const type & max, RngType & rng )
    {
      type res;
      for( size_t id_dim = 0; id_dim < N; ++id_dim )
      {
        std::uniform_real_distribution<scalar_type> distrib( min[id_dim], max[id_dim] );
        res[ id_dim ] = distrib( rng );
      }
      return res;
    }



    /**
    * @brief Compute minimum and maximum value of a set of points
    * @params elts list of points
    * @param[out] min minimum (component-wise) of the points
    * @param[out] max maximum (component-wise) of the points
    */
    static void minMax( const std::vector<type> & elts, type & min, type & max )
    {
      if( elts.size() == 0 )
      {
        return;
      }

      // Init
      std::fill( min.begin(), min.end(), std::numeric_limits<scalar_type>::max() );
      std::fill( max.begin(), max.end(), std::numeric_limits<scalar_type>::lowest() );

      // min/max search
      for( size_t id_pt = 0; id_pt < elts.size(); ++id_pt )
      {
        const type & cur_elt = elts[ id_pt ];
        for( size_t id_dim = 0; id_dim < cur_elt.size(); ++id_dim )
        {
          // Get min_max for ith dim
          min[ id_dim ] = std::min( min[id_dim], cur_elt[ id_dim ] );
          max[ id_dim ] = std::max( max[id_dim], cur_elt[ id_dim ] );
        }
      }
    }

    /**
    * @brief get a zero valued vector data
    * @param dummy a dummy vector
    * @return a null vector
    */
    static type null( const type & dummy )
    {
      type res;
      res.fill( scalar_type( 0 ) );
      return res;
    }


    /**
    * @brief Accumulate value inside a vector
    * @param self vector used for accumulation
    * @param data vector to add to the self vector
    * @note this perform self += data (component-wise)
    */
    static void accumulate( type & self, const type & data )
    {
      for( size_t id_dim = 0; id_dim < N; ++id_dim )
      {
        self[id_dim] += data[id_dim];
      }
    }

    /**
    * @brief Scalar division
    * @param self Vector to divide
    * @param val scalar divisor
    * @note this perform self /= data (component-wise)
    */
    static void divide( type & self, const size_t val )
    {
      for( size_t id_dim = 0; id_dim < N; ++id_dim )
      {
        self[id_dim] /= static_cast<scalar_type>( val );
      }
    }

};

/**
* @brief Overloading for std::vector
*/
template< typename T>
class KMeansVectorDataTrait<std::vector<T>>
{
  public:
    // base type
    typedef std::vector<T> type;

    typedef T scalar_type;

    /**
    * @brief number of element in the vector
    * @param aVector a Vector
    * @return number of scalar element in the vector
    */
    static size_t size( const type & aVector )
    {
      return aVector.size();
    }

    /**
    * @brief Square euclidean distance between two vectors
    * @param aVec1 first vector
    * @param aVec2 second vector
    * @return square euclidean distance
    */
    static scalar_type L2( const type & aVec1, const type & aVec2 )
    {
      typedef Eigen::Matrix<scalar_type, 1, Eigen::Dynamic> VecType;
      typedef Eigen::Map<const VecType> VecTypeConst;
      VecTypeConst map_aVec1(aVec1.data(), aVec1.size());
      VecTypeConst map_aVec2(aVec2.data(), aVec2.size());
      return ( map_aVec1 - map_aVec2 ).squaredNorm();
    }


    /**
    * @brief Draw a random vector in a range
    * @param min minimum bound of the sampling range
    * @param max maximum bound of the sampling range
    * @param rng A c++11 random generator
    * @return a Random vector in the given range
    */
    template < typename RngType >
    static type random( const type & min, const type & max, RngType & rng )
    {
      type res( min );
      for( size_t id_dim = 0; id_dim < res.size(); ++id_dim )
      {
        std::uniform_real_distribution<scalar_type> distrib( min[id_dim], max[id_dim] );
        res[ id_dim ] = distrib( rng );
      }
      return res;
    }


    /**
    * @brief Compute minimum and maximum value of a set of points
    * @params elts list of points
    * @param[out] min minimum (component-wise) of the points
    * @param[out] max maximum (component-wise) of the points
    */
    static void minMax( const std::vector<type> & elts, type & min, type & max )
    {
      if( elts.size() == 0 )
      {
        return;
      }

      // Init
      min.resize( elts[0].size(), std::numeric_limits<scalar_type>::max() );
      max.resize( elts[0].size(), std::numeric_limits<scalar_type>::lowest() );

      // min/max search
      for( size_t id_pt = 0; id_pt < elts.size(); ++id_pt )
      {
        const type & cur_elt = elts[ id_pt ];
        for( size_t id_dim = 0; id_dim < cur_elt.size(); ++id_dim )
        {
          // Get min_max for ith dim
          min[ id_dim ] = std::min( min[id_dim], cur_elt[ id_dim ] );
          max[ id_dim ] = std::max( max[id_dim], cur_elt[ id_dim ] );
        }
      }
    }

    /**
    * @brief get a zero valued vector data
    * @param dummy a dummy vector
    * @return a null vector
    */
    static type null( const type & dummy )
    {
      type res( dummy.size(), scalar_type( 0 ) );
      return res;
    }


    /**
    * @brief Accumulate value inside a vector
    * @param self vector used for accumulation
    * @param data vector to add to the self vector
    * @note this perform self += data (component-wise)
    */
    static void accumulate( type & self, const type & data )
    {
      for( size_t id_dim = 0; id_dim < self.size(); ++id_dim )
      {
        self[id_dim] += data[id_dim];
      }
    }


    /**
    * @brief Scalar division
    * @param self Vector to divide
    * @param val scalar divisor
    * @note this perform self /= data (component-wise)
    */
    static void divide( type & self, const size_t val )
    {
      for( size_t id_dim = 0; id_dim < self.size(); ++id_dim )
      {
        self[id_dim] /= static_cast<scalar_type>( val );
      }
    }
};


/**
* @brief Overloading for Vec2
*/
template<>
class KMeansVectorDataTrait<Vec2>
{
  public:
    // base type
    typedef Vec2 type;

    typedef double scalar_type;

    /**
    * @brief number of element in the vector
    * @param aVector a Vector
    * @return number of scalar element in the vector
    */
    static size_t size( const type & aVector )
    {
      return 2;
    }

    /**
    * @brief Square euclidean distance between two vectors
    * @param aVec1 first vector
    * @param aVec2 second vector
    * @return square euclidean distance
    */
    static scalar_type L2( const type & aVec1, const type & aVec2 )
    {
      return ( aVec1 - aVec2 ).squaredNorm();
    }


    /**
    * @brief Draw a random vector in a range
    * @param min minimum bound of the sampling range
    * @param max maximum bound of the sampling range
    * @param rng A c++11 random generator
    * @return a Random vector in the given range
    */
    template < typename RngType >
    static type random( const type & min, const type & max, RngType & rng )
    {
      std::uniform_real_distribution<scalar_type> distrib_x( min[0], max[0] );
      std::uniform_real_distribution<scalar_type> distrib_y( min[1], max[1] );

      return type( distrib_x( rng ), distrib_y( rng ) );
    }



    /**
    * @brief Compute minimum and maximum value of a set of points
    * @params elts list of points
    * @param[out] min minimum (component-wise) of the points
    * @param[out] max maximum (component-wise) of the points
    */
    static void minMax( const std::vector<type> & elts, type & min, type & max )
    {
      if( elts.size() == 0 )
      {
        return;
      }

      // Init
      min.Constant( std::numeric_limits<scalar_type>::max() );
      max.Constant( std::numeric_limits<scalar_type>::lowest() );

      // min/max search
      for( size_t id_pt = 0; id_pt < elts.size(); ++id_pt )
      {
        const type & cur_elt = elts[ id_pt ];
        for( size_t id_dim = 0; id_dim < 2; ++id_dim )
        {
          // Get min_max for ith dim
          min[ id_dim ] = std::min( min[id_dim], cur_elt[ id_dim ] );
          max[ id_dim ] = std::max( max[id_dim], cur_elt[ id_dim ] );
        }
      }
    }


    /**
    * @brief get a zero valued vector data
    * @param dummy a dummy vector
    * @return a null vector
    */
    static type null( const type & dummy )
    {
      return type( 0, 0 );
    }

    /**
    * @brief Accumulate value inside a vector
    * @param self vector used for accumulation
    * @param data vector to add to the self vector
    * @note this perform self += data (component-wise)
    */
    static void accumulate( type & self, const type & data )
    {
      self += data;
    }


    /**
    * @brief Scalar division
    * @param self Vector to divide
    * @param val scalar divisor
    * @note this perform self /= data (component-wise)
    */
    static void divide( type & self, const size_t val )
    {
      self /= static_cast<scalar_type>( val );
    }

};

/**
* @brief Overloading for Vec3
*/
template<>
class KMeansVectorDataTrait<Vec3>
{
  public:
    // base type
    typedef Vec3 type;

    typedef double scalar_type;

    /**
    * @brief number of element in the vector
    * @param aVector a Vector
    * @return number of scalar element in the vector
    */
    static size_t size( const type & aVector )
    {
      return 3;
    }

    /**
    * @brief Square euclidean distance between two vectors
    * @param aVec1 first vector
    * @param aVec2 second vector
    * @return square euclidean distance
    */
    static scalar_type L2( const type & aVec1, const type & aVec2 )
    {
      return ( aVec1 - aVec2 ).squaredNorm();
    }


    /**
    * @brief Draw a random vector in a range
    * @param min minimum bound of the sampling range
    * @param max maximum bound of the sampling range
    * @param rng A c++11 random generator
    * @return a Random vector in the given range
    */
    template < typename RngType >
    static type random( const type & min, const type & max, RngType & rng )
    {
      std::uniform_real_distribution<scalar_type> distrib_x( min[0], max[0] );
      std::uniform_real_distribution<scalar_type> distrib_y( min[1], max[1] );
      std::uniform_real_distribution<scalar_type> distrib_z( min[2], max[2] );

      return type( distrib_x( rng ), distrib_y( rng ), distrib_z( rng ) );
    }

    /**
    * @brief Compute minimum and maximum value of a set of points
    * @params elts list of points
    * @param[out] min minimum (component-wise) of the points
    * @param[out] max maximum (component-wise) of the points
    */
    static void minMax( const std::vector<type> & elts, type & min, type & max )
    {
      if( elts.size() == 0 )
      {
        return;
      }

      // Init
      min.Constant( std::numeric_limits<scalar_type>::max() );
      max.Constant( std::numeric_limits<scalar_type>::lowest() );

      // min/max search
      for( size_t id_pt = 0; id_pt < elts.size(); ++id_pt )
      {
        const type & cur_elt = elts[ id_pt ];
        for( size_t id_dim = 0; id_dim < 3; ++id_dim )
        {
          // Get min_max for ith dim
          min[ id_dim ] = std::min( min[id_dim], cur_elt[ id_dim ] );
          max[ id_dim ] = std::max( max[id_dim], cur_elt[ id_dim ] );
        }
      }
    }


    /**
    * @brief get a zero valued vector data
    * @param dummy a dummy vector
    * @return a null vector
    */
    static type null( const type & dummy )
    {
      return type( 0, 0, 0 );
    }


    /**
    * @brief Accumulate value inside a vector
    * @param self vector used for accumulation
    * @param data vector to add to the self vector
    * @note this perform self += data (component-wise)
    */
    static void accumulate( type & self, const type & data )
    {
      self += data;
    }


    /**
    * @brief Scalar division
    * @param self Vector to divide
    * @param val scalar divisor
    * @note this perform self /= data (component-wise)
    */
    static void divide( type & self, const size_t val )
    {
      self /= static_cast<scalar_type>( val );
    }
};

/**
* @brief Overloading for Vec4
*/
template<>
class KMeansVectorDataTrait<Vec4>
{
  public:
    // base type
    typedef Vec4 type;

    typedef double scalar_type;

    /**
    * @brief number of element in the vector
    * @param aVector a Vector
    * @return number of scalar element in the vector
    */
    static size_t size( const type & aVector )
    {
      return 4;
    }

    /**
    * @brief Square euclidean distance between two vectors
    * @param aVec1 first vector
    * @param aVec2 second vector
    * @return square euclidean distance
    */
    static scalar_type L2( const type & aVec1, const type & aVec2 )
    {
      return ( aVec1 - aVec2 ).squaredNorm();
    }


    /**
    * @brief Draw a random vector in a range
    * @param min minimum bound of the sampling range
    * @param max maximum bound of the sampling range
    * @param rng A c++11 random generator
    * @return a Random vector in the given range
    */
    template < typename RngType >
    static type random( const type & min, const type & max, RngType & rng )
    {
      std::uniform_real_distribution<scalar_type> distrib_x( min[0], max[0] );
      std::uniform_real_distribution<scalar_type> distrib_y( min[1], max[1] );
      std::uniform_real_distribution<scalar_type> distrib_z( min[2], max[2] );
      std::uniform_real_distribution<scalar_type> distrib_w( min[3], max[3] );

      return type( distrib_x( rng ), distrib_y( rng ), distrib_z( rng ), distrib_w( rng ) );
    }

    /**
    * @brief Compute minimum and maximum value of a set of points
    * @params elts list of points
    * @param[out] min minimum (component-wise) of the points
    * @param[out] max maximum (component-wise) of the points
    */
    static void minMax( const std::vector<type> & elts, type & min, type & max )
    {
      if( elts.size() == 0 )
      {
        return;
      }

      // Init
      min.Constant( std::numeric_limits<scalar_type>::max() );
      max.Constant( std::numeric_limits<scalar_type>::lowest() );

      // min/max search
      for( size_t id_pt = 0; id_pt < elts.size(); ++id_pt )
      {
        const type & cur_elt = elts[ id_pt ];
        for( size_t id_dim = 0; id_dim < 4; ++id_dim )
        {
          // Get min_max for ith dim
          min[ id_dim ] = std::min( min[id_dim], cur_elt[ id_dim ] );
          max[ id_dim ] = std::max( max[id_dim], cur_elt[ id_dim ] );
        }
      }
    }

    /**
    * @brief get a zero valued vector data
    * @param dummy a dummy vector
    * @return a null vector
    */
    static type null( const type & dummy )
    {
      return type( 0, 0, 0, 0 );
    }


    /**
    * @brief Accumulate value inside a vector
    * @param self vector used for accumulation
    * @param data vector to add to the self vector
    * @note this perform self += data (component-wise)
    */
    static void accumulate( type & self, const type & data )
    {
      self += data;
    }


    /**
    * @brief Scalar division
    * @param self Vector to divide
    * @param val scalar divisor
    * @note this perform self /= data (component-wise)
    */
    static void divide( type & self, const size_t val )
    {
      self /= static_cast<scalar_type>( val );
    }
};


/**
* @brief Overloading for Vec2
*/
template<>
class KMeansVectorDataTrait<Vec2f>
{
  public:
    // base type
    typedef Vec2f type;

    typedef float scalar_type;

    /**
    * @brief number of element in the vector
    * @param aVector a Vector
    * @return number of scalar element in the vector
    */
    static size_t size( const type & aVector )
    {
      return 2;
    }

    /**
    * @brief Square euclidean distance between two vectors
    * @param aVec1 first vector
    * @param aVec2 second vector
    * @return square euclidean distance
    */
    static scalar_type L2( const type & aVec1, const type & aVec2 )
    {
      return ( aVec1 - aVec2 ).squaredNorm();
    }


    /**
    * @brief Draw a random vector in a range
    * @param min minimum bound of the sampling range
    * @param max maximum bound of the sampling range
    * @param rng A c++11 random generator
    * @return a Random vector in the given range
    */
    template < typename RngType >
    static type random( const type & min, const type & max, RngType & rng )
    {
      std::uniform_real_distribution<scalar_type> distrib_x( min[0], max[0] );
      std::uniform_real_distribution<scalar_type> distrib_y( min[1], max[1] );

      return type( distrib_x( rng ), distrib_y( rng ) );
    }


    /**
    * @brief Compute minimum and maximum value of a set of points
    * @params elts list of points
    * @param[out] min minimum (component-wise) of the points
    * @param[out] max maximum (component-wise) of the points
    */
    static void minMax( const std::vector<type> & elts, type & min, type & max )
    {
      if( elts.size() == 0 )
      {
        return;
      }

      // Init
      min.Constant( std::numeric_limits<scalar_type>::max() );
      max.Constant( std::numeric_limits<scalar_type>::lowest() );

      // min/max search
      for( size_t id_pt = 0; id_pt < elts.size(); ++id_pt )
      {
        const type & cur_elt = elts[ id_pt ];
        for( size_t id_dim = 0; id_dim < 2; ++id_dim )
        {
          // Get min_max for ith dim
          min[ id_dim ] = std::min( min[id_dim], cur_elt[ id_dim ] );
          max[ id_dim ] = std::max( max[id_dim], cur_elt[ id_dim ] );
        }
      }
    }

    /**
    * @brief get a zero valued vector data
    * @param dummy a dummy vector
    * @return a null vector
    */
    static type null( const type & dummy )
    {
      return type( 0, 0 );
    }


    /**
    * @brief Accumulate value inside a vector
    * @param self vector used for accumulation
    * @param data vector to add to the self vector
    * @note this perform self += data (component-wise)
    */
    static void accumulate( type & self, const type & data )
    {
      self += data;
    }


    /**
    * @brief Scalar division
    * @param self Vector to divide
    * @param val scalar divisor
    * @note this perform self /= data (component-wise)
    */
    static void divide( type & self, const size_t val )
    {
      self /= static_cast<scalar_type>( val );
    }
};

/**
* @brief Overloading for Vec3
*/
template<>
class KMeansVectorDataTrait<Vec3f>
{
  public:
    // base type
    typedef Vec3f type;

    typedef float scalar_type;

    /**
    * @brief number of element in the vector
    * @param aVector a Vector
    * @return number of scalar element in the vector
    */
    static size_t size( const type & aVector )
    {
      return 3;
    }

    /**
    * @brief Square euclidean distance between two vectors
    * @param aVec1 first vector
    * @param aVec2 second vector
    * @return square euclidean distance
    */
    static scalar_type L2( const type & aVec1, const type & aVec2 )
    {
      return ( aVec1 - aVec2 ).squaredNorm();
    }


    /**
    * @brief Draw a random vector in a range
    * @param min minimum bound of the sampling range
    * @param max maximum bound of the sampling range
    * @param rng A c++11 random generator
    * @return a Random vector in the given range
    */
    template < typename RngType >
    static type random( const type & min, const type & max, RngType & rng )
    {
      std::uniform_real_distribution<scalar_type> distrib_x( min[0], max[0] );
      std::uniform_real_distribution<scalar_type> distrib_y( min[1], max[1] );
      std::uniform_real_distribution<scalar_type> distrib_z( min[2], max[2] );

      return type( distrib_x( rng ), distrib_y( rng ), distrib_z( rng ) );
    }

    /**
    * @brief Compute minimum and maximum value of a set of points
    * @params elts list of points
    * @param[out] min minimum (component-wise) of the points
    * @param[out] max maximum (component-wise) of the points
    */
    static void minMax( const std::vector<type> & elts, type & min, type & max )
    {
      if( elts.size() == 0 )
      {
        return;
      }

      // Init
      min.Constant( std::numeric_limits<scalar_type>::max() );
      max.Constant( std::numeric_limits<scalar_type>::lowest() );

      // min/max search
      for( size_t id_pt = 0; id_pt < elts.size(); ++id_pt )
      {
        const type & cur_elt = elts[ id_pt ];
        for( size_t id_dim = 0; id_dim < 3; ++id_dim )
        {
          // Get min_max for ith dim
          min[ id_dim ] = std::min( min[id_dim], cur_elt[ id_dim ] );
          max[ id_dim ] = std::max( max[id_dim], cur_elt[ id_dim ] );
        }
      }
    }


    /**
    * @brief get a zero valued vector data
    * @param dummy a dummy vector
    * @return a null vector
    */
    static type null( const type & dummy )
    {
      return type( 0, 0, 0 );
    }


    /**
    * @brief Accumulate value inside a vector
    * @param self vector used for accumulation
    * @param data vector to add to the self vector
    * @note this perform self += data (component-wise)
    */
    static void accumulate( type & self, const type & data )
    {
      self += data;
    }


    /**
    * @brief Scalar division
    * @param self Vector to divide
    * @param val scalar divisor
    * @note this perform self /= data (component-wise)
    */
    static void divide( type & self, const size_t val )
    {
      self /= static_cast<scalar_type>( val );
    }
};


/**
* @brief Specialization for Vec
*/
/**
* @brief Overloading for Vec3
*/
template<>
class KMeansVectorDataTrait<Vec>
{
  public:
    // base type
    typedef Vec type;

    typedef double scalar_type;

    /**
    * @brief number of element in the vector
    * @param aVector a Vector
    * @return number of scalar element in the vector
    */
    static size_t size( const type & aVector )
    {
      return aVector.size();
    }

    /**
    * @brief Square euclidean distance between two vectors
    * @param aVec1 first vector
    * @param aVec2 second vector
    * @return square euclidean distance
    */
    static scalar_type L2( const type & aVec1, const type & aVec2 )
    {
      return ( aVec1 - aVec2 ).squaredNorm();
    }


    /**
    * @brief Draw a random vector in a range
    * @param min minimum bound of the sampling range
    * @param max maximum bound of the sampling range
    * @param rng A c++11 random generator
    * @return a Random vector in the given range
    */
    template < typename RngType >
    static type random( const type & min, const type & max, RngType & rng )
    {
      Vec res( min.size() );
      for( size_t id_dim = 0; id_dim < res.size(); ++id_dim )
      {
        std::uniform_real_distribution<scalar_type> distrib( min[0], max[0] );
        res[id_dim] = distrib( rng );
      }
      return res;
    }

    /**
    * @brief Compute minimum and maximum value of a set of points
    * @params elts list of points
    * @param[out] min minimum (component-wise) of the points
    * @param[out] max maximum (component-wise) of the points
    */
    static void minMax( const std::vector<type> & elts, type & min, type & max )
    {
      if( elts.size() == 0 )
      {
        return;
      }

      // Init
      min = Vec( elts[0].size() );
      min.fill( std::numeric_limits<scalar_type>::max() );
      max = Vec( elts[0].size() );
      max.fill( std::numeric_limits<scalar_type>::lowest() );

      // min/max search
      for( size_t id_pt = 0; id_pt < elts.size(); ++id_pt )
      {
        const type & cur_elt = elts[ id_pt ];
        min = min.cwiseMin( cur_elt );
        max = max.cwiseMax( cur_elt );
      }
    }

    /**
    * @brief get a zero valued vector data
    * @param dummy a dummy vector
    * @return a null vector
    */
    static type null( const type & dummy )
    {
      Vec res( dummy.size() );
      res.fill( scalar_type( 0 ) );
      return res;
    }


    /**
    * @brief Accumulate value inside a vector
    * @param self vector used for accumulation
    * @param data vector to add to the self vector
    * @note this perform self += data (component-wise)
    */
    static void accumulate( type & self, const type & data )
    {
      self += data;
    }

    /**
    * @brief Scalar division
    * @param self Vector to divide
    * @param val scalar divisor
    * @note this  perform self /= data (component-wise)
    */
    static void divide( type & self, const size_t val )
    {
      self /= static_cast<scalar_type>( val );
    }
};


} // namespace clustering
} // namespace openMVG

#endif
