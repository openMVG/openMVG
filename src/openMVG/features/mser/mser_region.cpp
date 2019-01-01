// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "mser_region.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

namespace openMVG
{
  namespace features
  {
    namespace MSER
    {
      /**
      * @brief region_level Region value (grayscale value)
      * @param base_x Base pixel of the region in x-coordinate
      * @param base_y Base pixel of the region in y-coordinate
      */
      MSERRegion::MSERRegion( const int region_level , const int base_x , const int base_y )
        : m_parent( nullptr ) ,
          m_child( nullptr ) ,
          m_sister( nullptr ) ,
          m_is_stable( false ) ,
          m_variation( std::numeric_limits<double>::infinity() ) ,
          m_level( region_level ) ,
          m_base_x( base_x ) ,
          m_base_y( base_y ) ,
          m_area( 0 ) ,
          m_moment_x( 0 ) ,
          m_moment_x2( 0 ) ,
          m_moment_y( 0 ) ,
          m_moment_y2( 0 ) ,
          m_moment_xy( 0 )
      {

      }

      /**
      * Compute ellipse position of the given region
      * @param[out] (x,y) -> Center of the ellipse
      */
      void MSERRegion::FitEllipse( double & x , double & y ) const
      {
        /*
        * Compute center of the ellipse
        * x_mean = int_{region} X / area
        * y_mean = int_{region} Y / area
        */
        x = m_moment_x / static_cast<double>( m_area );
        y = m_moment_y / static_cast<double>( m_area );
      }

      /**
      * Compute covariance matrix corresponding to the ellipse
      *
      * cov = | a b | = | dx * dx    dx * dy |
              | b c |   | dx * dy    dy * dy |
      * @param[out] (a,b,c) -> Value of the covariant matrix
      */
      void MSERRegion::FitEllipse( double & a , double & b , double & c) const
      {
        /*
        * Compute center of the ellipse
        * x_mean = int_{region} X / area
        * y_mean = int_{region} Y / area
        */
        const double ell_x = m_moment_x / static_cast<double>( m_area );
        const double ell_y = m_moment_y / static_cast<double>( m_area );

        /*
        * Compute covariance matrix corresponding to the ellipse
        */
        const double i20 = m_moment_x2 - static_cast<double>( m_area ) * ell_x * ell_x;
        const double i02 = m_moment_y2 - static_cast<double>( m_area ) * ell_y * ell_y;
        const double i11 = m_moment_xy - static_cast<double>( m_area ) * ell_x * ell_y;
        const double n = i20*i02 - i11*i11;

        if (n != 0)
        {
          a = i02/n * (m_area-1)/4.;
          b = -i11/n * (m_area-1)/4.;
          c = i20/n * (m_area-1)/4.;
        }
       }

      /**
      * @brief Fit ellipse for the given region
      * @param[out] (ell_x,ell_y) -> Center of the ellipse
      * @param[out] (ell_major_x , ell_major_y) -> major axis base vector (normalized)
      * @param[out] (ell_minor_x , ell_minor_y) -> minor axis base vector (normalized)
      * @param[out] ell_major_length -> (half) Length of the major axis
      * @param[out] ell_minor_length -> (half) Length of the minor axis
      */
      void MSERRegion::FitEllipse( double & ell_x , double & ell_y ,
                                   double & ell_major_x , double & ell_major_y ,
                                   double & ell_minor_x , double & ell_minor_y ,
                                   double & ell_major_length , double & ell_minor_length ) const
      {
        /*
        * Compute center of the ellipse
        * x_mean = int_{region} X / area
        * y_mean = int_{region} Y / area
        */
        ell_x = m_moment_x / static_cast<double>( m_area );
        ell_y = m_moment_y / static_cast<double>( m_area );

        /*
        * Compute covariance matrix corresponding to the ellipse
        *
        * cov = | a b | = | dx * dx    dx * dy |
                | b c |   | dx * dy    dy * dy |
        * dx = x - x_mean
        * dy = x - y_mean
        */
        const double a = m_moment_x2 / static_cast<double>( m_area ) - ell_x * ell_x;
        const double b = m_moment_xy / static_cast<double>( m_area ) - ell_x * ell_y;
        const double c = m_moment_y2 / static_cast<double>( m_area ) - ell_y * ell_y;

        /**
        * Extract major and minor axis of the ellipse.
        * axis are given by the eigen vector of the ellipse
        * axis length are eigen values
        */
        // http://math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
        const double d  = a + c; // Trace
        const double e  = a - c; // Simplified det (todo : need to check if it's correct since b could be 0)
        const double f  = std::sqrt( 4.0 * b * b + e * e );

        // Eigenvalues (square of the axis lengths)
        const double e0 = ( d + f ) * 0.5;
        const double e1 = ( d - f ) * 0.5;

        ell_major_length = std::sqrt( e0 );
        ell_minor_length = std::sqrt( e1 );

        // Check if ellipse is axis aligned or not
        if (std::abs( b ) < 0.0001 )
        {
          // Axis aligned
          ell_major_x = 1.0;
          ell_major_y = 0.0;

          ell_minor_x = 0.0;
          ell_minor_y = 1.0;
        }
        else
        {
          // Not axis aligned
          ell_major_x = e0 - c;
          ell_major_y = b;

          ell_minor_x = e1 - c;
          ell_minor_y = b;

          // Normalize axis
          const double inv_norm_major = 1.0 / std::hypot( ell_major_x, ell_major_y );
          const double inv_norm_minor = 1.0 / std::hypot( ell_minor_x, ell_minor_y );

          ell_major_x *= inv_norm_major;
          ell_major_y *= inv_norm_major;

          ell_minor_x *= inv_norm_minor;
          ell_minor_y *= inv_norm_minor;
        }

        // Ensure correct ordering on axis
        if (ell_minor_length > ell_major_length )
        {
          std::swap( ell_minor_length , ell_major_length );
          std::swap( ell_major_x , ell_minor_x );
          std::swap( ell_major_y , ell_minor_y );
        }
      }

      /**
       * @brief Add a new pixel in the region
       * @param x X-Coord of the pixel to append
       * @param y Y-Coord of the pixel to append
       */
      void MSERRegion::AppendPixel( const int x , const int y )
      {
        // One more pixel -> area is growing (what a deduction !)
        ++m_area;

        // update geometric moments
        m_moment_x += static_cast<double>( x );
        m_moment_y += static_cast<double>( y );

        m_moment_x2 += static_cast<double>( x * x );
        m_moment_y2 += static_cast<double>( y * y );

        m_moment_xy += static_cast<double>( x * y );
      }

      /**
      * @brief Add a new child to this region (append infos and make connectivity)
      * @param child Region to append
      */
      void MSERRegion::MergeRegion( MSERRegion * child )
      {
        // Update area
        m_area += child->m_area;

        // Update moments
        m_moment_x += child->m_moment_x;
        m_moment_y += child->m_moment_y;
        m_moment_x2 += child->m_moment_x2;
        m_moment_y2 += child->m_moment_y2;
        m_moment_xy += child->m_moment_xy;

        // Enque child in the linked child list
        child->m_sister = m_child;
        m_child = child;

        // Update child parent
        child->m_parent = this;
      }


      /**
      * @brief Compute MSER stability , compress parent path, and check if the region (and its children) is stable
      * @param delta distance (in level) to the other region to compute stability (stability is not made up to delta level)
      * @param minimumArea Minimum area to be accepted as stable
      * @param maximumArea Maximum area to be accepted as stable
      * @param maxVariation Maximum variation (difference between parent area and region area) to be accepted as stable
      */
      void MSERRegion::ComputeStability( const int delta , const int minimumArea , const int maximumArea , const double maxVariation )
      {
        // Found the upmost region with region no more than delta (ie: reference for stability)
        MSERRegion * ref = this;

        while (ref->m_parent && ( ref->m_parent->m_level <= ( m_level + delta ) ))
        {
          ref = ref->m_parent; // Climb hierarchy
        }

        // Compute variation between the reference region and the current region
        m_variation = static_cast<double>( ref->m_area - m_area ) / static_cast<double>( m_area );

        // Ensure the region could be qualified as stable
        const bool stable =
          // Hierarchy MSER criterion
          ( ( !ref ) || m_variation <= ref->m_variation ) && // Top level region or less variation than its parent
          // Basic MSER criterion
          ( ( m_variation <= maxVariation ) && ( m_area >= minimumArea ) && ( m_area <= maximumArea ) );

        // Process the children and check at least one child is stable
        MSERRegion * cur_child = m_child;
        if (! cur_child )
        {
          if (stable )
          {
            m_is_stable = true;
          }
        }
        else
        {
          // Process child list
          while (cur_child != nullptr)
          {
            cur_child->ComputeStability( delta , minimumArea , maximumArea , maxVariation );

            // At least one child is stable -> the region is stable
            if (stable && ( m_variation < cur_child->m_variation ))
            {
              m_is_stable = true;
            }

            cur_child = cur_child->m_sister;
          }
        }
      }

      /**
      * @brief check criteria on the whole region
      * @param variation Minimum variation to validate the criteria
      * @param area Maximum area to validate the criteria
      * @return true if it as minimum variation and maximum area (on the whole region, ie: region and its child)
      */
      bool MSERRegion::CheckCriteria( const double variation , const int area ) const
      {
        // Has minimum area -> that's enough to decide
        if (m_area <= area )
        {
          return true;
        }

        // It's stable and with enough variation
        if (m_is_stable && ( m_variation < variation ) )
        {
          return false;
        }

        MSERRegion * cur_child = m_child;
        while (cur_child != nullptr)
        {
          if (! cur_child->CheckCriteria( variation , area ) )
          {
            return false;
          }

          cur_child = cur_child->m_sister;
        }

        // Every child pass the criteria
        return true;
      }

      /**
      * Cut all reference to any region in the hierarchy
      */
      void MSERRegion::Isolate( )
      {
        m_parent = nullptr;
        m_child = nullptr;
        m_sister = nullptr;
      }


      /**
      * @brief Export region and it's children if it validate the MSER diversity criterion
      * @param minDiversity Minimum diversity to validate the region
      * @param[out] regions Region that validate the criterion
      */
      void MSERRegion::Export( const double minDiversity , std::vector<MSERRegion> & regions )
      {
        // Only export if it's stable
        if (m_is_stable )
        {
          // But don't export only if it's stable !
          // Need to check the hierarchy if parent is more stable than this region

          // 1 - Parent must have minimum size
          const int minimum_parent_area = static_cast<int>( static_cast<double>( m_area ) / ( 1.0 - minDiversity ) + 0.5 );

          MSERRegion * cur_parent = this;
          while (cur_parent->m_parent != nullptr && ( cur_parent->m_parent->m_area < minimum_parent_area ))
          {
            cur_parent = cur_parent->m_parent;

            // 2 - Parent must have move diversity than its children
            // If not it's more stable than this region
            if (cur_parent->m_is_stable && ( cur_parent->m_variation <= m_variation ) )
            {
              m_is_stable = false;
              break;
            }
          }

          // If it's still stable, check children area
          if (m_is_stable )
          {
            const int max_children_area = static_cast<int>( static_cast<double>( m_area ) * ( 1.0 - minDiversity ) + 0.5 );
            if (! CheckCriteria( m_variation , max_children_area ) )
            {
              // Does'nt pass stability criteria
              m_is_stable = false;
            }
          }

          // If it successfully pass all the test, add the region
          if (m_is_stable )
          {
            // Make a copy
            regions.push_back( *this );
            regions.back().Isolate();
          }
        }

        // Try to export the child regions
        MSERRegion * cur_child = m_child;
        while (cur_child != nullptr)
        {
          cur_child->Export( minDiversity , regions );
          cur_child = cur_child->m_sister;
        }

      }

      /**
      * @brief Compute valid MSER regions from a given region
      * @param delta distance (in level) to the other region to compute stability (stability is not made up to delta level)
      * @param minimumArea Minimum area to be accepted as stable
      * @param maximumArea Maximum area to be accepted as stable
      * @param maxVariation Maximum variation (difference between parent area and region area) to be accepted as stable
      * @param minDiversity Minimum diversity to validate the region
      * @param[out] regions Region that validate the criterion
      */
      void MSERRegion::ComputeMSER( const int delta , const int minArea , const int maxArea , const double maxVariation , const double minDiversity , std::vector<MSERRegion> & regions )
      {
        // Compute MSER stability stats on all the hierarchy
        ComputeStability( delta , minArea , maxArea , maxVariation );
        // Export regions that pass all MSER stability checks (region and its direct hierarchy)
        Export( minDiversity , regions );
      }

    }

  }

}
