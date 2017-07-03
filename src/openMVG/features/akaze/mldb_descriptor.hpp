// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Romuald PERROT, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_DESCRIPTION_MLDB_DESCRIPTOR_HPP
#define OPENMVG_IMAGE_DESCRIPTION_MLDB_DESCRIPTOR_HPP

#include <algorithm>
#include <cmath>

#include "openMVG/features/descriptor.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/image/image_container.hpp"

namespace openMVG {
namespace features {

  /**
  ** @brief Compute mean values (Li,Lx,Ly) in each subdivisions
  ** @param samples_Li input values on Li
  ** @param samples_Lx input values on Lx
  ** @param samples_Ly input values on Ly
  ** @param nb_subdiv number of subdivision (on each 2d axis)
  ** @param subdiv_size size of a subdivision (on each 2d axis)
  ** @param pattern_size source size (on each 2d axis)
  ** @param c cosinus of main orientation (for computing new derivatives values)
  ** @param s sinus of main orientation (for computing new derivatives values)
  ** @param mean_Li mean of Li in each subdivision
  ** @param mean_Lx mean of Lx in each subdivision
  ** @param mean_Ly mean of Ly in each subdivision
  **/
  template< typename Real>
  static inline void ComputeMeanValuesInSubdivisions(
      const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & samples_Li ,
      const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & samples_Lx ,
      const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & samples_Ly ,
      const int nb_subdiv ,
      const int subdiv_size ,
      const int pattern_size ,
      const Real c , // cos( theta )
      const Real s , // sin( theta )
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & mean_Li ,
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & mean_Lx ,
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & mean_Ly )
  {
    mean_Li.resize( nb_subdiv , nb_subdiv );
    mean_Lx.resize( nb_subdiv , nb_subdiv );
    mean_Ly.resize( nb_subdiv , nb_subdiv );

    const int max_w = samples_Li.cols();
    const int max_h = samples_Li.rows();

    for (int i = 0; i < nb_subdiv; ++i )
    {
      for (int j = 0; j < nb_subdiv; ++j )
      {
        // Compute subdivision extend
        const int min_x = j * subdiv_size;
        const int min_y = i * subdiv_size;
        const int max_x = std::min( ( j + 1 ) * subdiv_size , max_w );
        const int max_y = std::min( ( i + 1 ) * subdiv_size , max_h );

        // Sum every elements of this subdivision
        mean_Li( i , j ) = mean_Lx( i , j ) = mean_Ly( i , j ) = 0;

        size_t nb_elt = 0;
        for (int ii = min_y; ii < max_y; ++ii )
        {
          for (int jj = min_x; jj < max_x; ++jj )
          {
            mean_Li( i , j ) += samples_Li( ii , jj );

            const Real dx = samples_Lx( ii , jj );
            const Real dy = samples_Ly( ii , jj );

            // Rotate derivatives
            // a is original angle, b is keypoint angle
            // Cos( a - b ) = cosA cosB + sinA sinB
            //              = dx * c + dy * s
            // Sin( a - b ) = sinA cosB - cosA sinB
            //              = dy * c - dx * s
            mean_Ly( i , j ) +=   dx * c + dy * s;
            mean_Lx( i , j ) +=   dy * c - dx * s;

            ++nb_elt;
          }
        }

        mean_Li( i , j ) /= static_cast<Real>( nb_elt );
        mean_Lx( i , j ) /= static_cast<Real>( nb_elt );
        mean_Ly( i , j ) /= static_cast<Real>( nb_elt );
      }
    }
  }

  /**
  ** @brief Compute binary description
  ** @param mean_Li input mean values of Li values (mean per subdivision)
  ** @param mean_Lx input mean values of Lx values (mean per subdivision)
  ** @param mean_Ly input mean values of Ly values (mean per subdivision)
  ** @param nb_subdiv Number of subdivision (in 2d so it's a nb_subdivxnb_subdiv pattern)
  ** @param outIndex input/ouput index to store description
  ** @param desc ouput vector (idealy a std::bitset) containing binary description of theses regions
  **/
  template< typename DescriptorType , typename Real>
  static inline void ComputeBinaryValues(
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & mean_Li ,
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & mean_Lx ,
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & mean_Ly ,
    const int nb_subdiv ,
    size_t & outIndex ,
    DescriptorType & desc )
  {
    // Binary comparisons (ie (0,0) with (0,1), (O,0) with (0,2), ... )
    for (int i = 0; i < ( nb_subdiv * nb_subdiv ); ++i )
    {
      for (int j = i + 1; j < ( nb_subdiv * nb_subdiv ); ++j )
      {
        const int src_i = i / nb_subdiv;
        const int src_j = i % nb_subdiv;

        const int dst_i = j / nb_subdiv;
        const int dst_j = j % nb_subdiv;

        // Compare (src_i,src_j) with (dst_i,dst_j) on the three values
        desc[ outIndex++ ] = mean_Li( src_i , src_j ) > mean_Li( dst_i , dst_j );
        desc[ outIndex++ ] = mean_Lx( src_i , src_j ) > mean_Lx( dst_i , dst_j );
        desc[ outIndex++ ] = mean_Ly( src_i , src_j ) > mean_Ly( dst_i , dst_j );
      }
    }
  }

  /**
    ** @brief Compute final keypoint (ie interest point + description) for a given interest point
    ** @param Li Input Octave slice
    ** @param Lx Input X-derivative
    ** @param Ly Input Y-derivative
    ** @param id_octave Id of current octave
    ** @param ipt Input interest point
    ** @param desc output descriptor (binary descriptor)
    **/
  template< typename Real>
  void ComputeMLDBDescriptor(
    const image::Image<Real> & Li,
    const image::Image<Real> &Lx,
    const image::Image<Real> &Ly,
    const int id_octave ,
    const SIOPointFeature & ipt ,
    Descriptor<bool, 486> & desc )
  {
    // // Note : in KAZE description we compute descriptor of previous slice and never the current one

    // See if it's necessary to change this value (pass it into a param ?)
    const int pattern_size = 10;

    // Sampling size according to the scale value
    const Real inv_octave_scale = static_cast<Real>( 1 ) / static_cast<Real>( 1 << id_octave );
    const Real sigma_scale = std::round( ipt.scale() * inv_octave_scale );

    // Get every samples inside 2pattern x 2pattern square region
    // Memory efficient (get samples then work in aligned)
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      samples_Li( 2 * pattern_size + 1 , 2 * pattern_size + 1 ),
      samples_Lx( 2 * pattern_size + 1 , 2 * pattern_size + 1 ),
      samples_Ly( 2 * pattern_size + 1 , 2 * pattern_size + 1 );

    // Compute cos and sin values for this point orientation
    const Real c = std::cos( ipt.orientation() );
    const Real s = std::sin( ipt.orientation() );

    const Real cur_x = ipt.x() * inv_octave_scale;
    const Real cur_y = ipt.y() * inv_octave_scale;

    // Retrieve samples in oriented pattern
    for (int i = - pattern_size; i <= pattern_size; ++i )
    {
      for (int j = - pattern_size; j <= pattern_size; ++j )
      {
        // sample_y = yf + (l*scale*co + k*scale*si);
        // sample_x = xf + (-l*scale*si + k*scale*co);

        // Need to do a function for that rotate of angle -theta ( minus theta because image frame is Y down)
        const Real delta_y =   static_cast<Real>( j ) * c + static_cast<Real>( i ) * s;
        const Real delta_x =  -static_cast<Real>( j ) * s + static_cast<Real>( i ) * c;

        // Compute new real position
        const Real dx = cur_x + sigma_scale * delta_x;
        const Real dy = cur_y + sigma_scale * delta_y;

        // Compute new discrete position
        const int y = std::round( dy );
        const int x = std::round( dx );

        samples_Li( i + pattern_size , j + pattern_size ) = Li( y , x );
        samples_Lx( i + pattern_size , j + pattern_size ) = Lx( y , x );
        samples_Ly( i + pattern_size , j + pattern_size ) = Ly( y , x );
      }
    }

    size_t outIndex = 0; // Index to store next binary value

    // Grid 1 : 2x2 subdivision
    int subdiv_size = pattern_size;
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> sumLi , sumLx , sumLy;
    ComputeMeanValuesInSubdivisions( samples_Li , samples_Lx , samples_Ly , 2 , subdiv_size , pattern_size , c , s , sumLi , sumLx , sumLy );
    ComputeBinaryValues( sumLi , sumLx , sumLy , 2 , outIndex , desc );

    // Grid 2 : 3x3 subdivision
    subdiv_size = static_cast<int>( std::ceil( static_cast<Real>( 2 * pattern_size ) / static_cast<Real>( 3 ) ) );
    ComputeMeanValuesInSubdivisions( samples_Li , samples_Lx , samples_Ly , 3 , subdiv_size , pattern_size , c , s , sumLi , sumLx , sumLy );
    ComputeBinaryValues( sumLi , sumLx , sumLy , 3 , outIndex , desc );

    // Grid 3 : 4x4 subdivision
    subdiv_size = pattern_size / 2;
    ComputeMeanValuesInSubdivisions( samples_Li , samples_Lx , samples_Ly , 4 , subdiv_size , pattern_size , c , s , sumLi , sumLx , sumLy );
    ComputeBinaryValues( sumLi , sumLx , sumLy , 4 , outIndex , desc );

    assert( outIndex == 486 ); // Just to be sure (and we are sure ! completly sure !)
  }

} // namespace features
} // namespace openMVG

#endif // OPENMVG_IMAGE_DESCRIPTION_MLDB_DESCRIPTOR_HPP
