// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MSER_MSER_REGION_HPP
#define OPENMVG_MSER_MSER_REGION_HPP

#include <vector>

namespace openMVG
{
namespace features
{
namespace MSER
{

    /**
     * @brief Class holding a MSER region
     * @note User should not deal with construction of the class but only use MSERExtractor to build regions
     * User should only use the FitEllipse method
     */
    class MSERRegion
    {
      friend class MSERExtractor;

      public:
      /**
        * @brief region_level Region value (grayscale value)
        * @param base_x Base pixel of the region in x-coordinate
        * @param base_y Base pixel of the region in y-coordinate
        */
      MSERRegion( const int region_level = 256, const int base_x = 0, const int base_y = 0 );

      /**
        * @brief Fit ellipse for the given region
        * @param[out] (ell_x,ell_y) -> Center of the ellipse
        * @param[out] (ell_major_x , ell_major_y) -> major axis base vector (normalized)
        * @param[out] (ell_minor_x , ell_minor_y) -> minor axis base vector (normalized)
        * @param[out] ell_major_length -> (half) Length of the major axis
        * @param[out] ell_minor_length -> (half) Length of the minor axis
        */
      void FitEllipse( double &ell_x, double &ell_y,
                       double &ell_major_x, double &ell_major_y,
                       double &ell_minor_x, double &ell_minor_y,
                       double &ell_major_length, double &ell_minor_length ) const;

      /**
        * Compute covariance matrix corresponding to the ellipse
        *
        * cov = | a b | = | dx * dx    dx * dy |
                | b c |
        * @param[out] (a,b,c) -> Value of the covariant matrix
        */
      void FitEllipse( double &a, double &b, double &c ) const;

      /**
        * Compute ellipse position of the given region
        * @param[out] (x,y) -> Center of the ellipse
        */
      void FitEllipse( double &x, double &y ) const;

      private:

      // Connectivity
      MSERRegion *m_parent; // Pointer to the parent region
      MSERRegion *m_child;  // Pointer to the first child of the region
      MSERRegion *m_sister; // Pointer to the next region (at the same level)

      // MSER statistics
      bool   m_is_stable;
      double m_variation; // Variation between parent area and this area

      // Region info
      int m_level;  // Level (grayscale value) of the region
      int m_base_x; // X-coord of the initial pixel of the region
      int m_base_y; // Y-coord of the initial pixel of the region
      int m_area;   // Area of the region (ie: nb pixels)

      // Geometric Moments : used for computation of the underlaying ellipse of the region
      double m_moment_x;  // int_{Region} X
      double m_moment_x2; // int_{Region} X*X
      double m_moment_y;  // int_{Region} Y
      double m_moment_y2; // int_{Region} Y * Y
      double m_moment_xy; // int_{Region} X * Y

      /**
        * Cut all reference to any region in the hierarchy
        */
      void Isolate( void );

      /**
        * @brief Add a new pixel in the region
        * @param x X-Coord of the pixel to append
        * @param y Y-Coord of the pixel to append
        */
      void AppendPixel( const int x, const int y );

      /**
        * @brief Add a new child to this region (append infos and make connectivity)
        * @param child Region to append
        */
      void MergeRegion( MSERRegion *child );

      /**
        * @brief Compute MSER stability , compress parent path, and check if the region (and it's children) is stable
        * @param delta distance (in level) to the other region to compute stability (stability is not made up to delta level)
        * @param minimumArea Minimum area to be accepted as stable
        * @param maximumArea Maximum area to be accepted as stable
        * @param maxVariation Maximum variation (difference between parent area and region area) to be accepted as stable
        */
      void ComputeStability( const int delta, const int minimumArea, const int maximumArea, const double maxVariation );

      /**
        * @brief check criteria on the whole region
        * @param variation Minimum variation to validate the criteria
        * @param area Maximum area to validate the criteria
        * @return true if it as minimum variation and maximum area (on the whole region, ie: region and it's child)
        */
      bool CheckCriteria( const double variation, const int area ) const;

      /**
        * @brief Export region and it's children if it validate the MSER diversity criterion
        * @param minDiversity Minimum diversity to validate the region
        * @param[out] regions Region that validate the criterion
        */
      void Export( const double minDiversity, std::vector<MSERRegion> &regions );

      /**
        * @brief Compute valid MSER regions from a given region
        * @param delta distance (in level) to the other region to compute stability (stability is not made up to delta level)
        * @param minimumArea Minimum area to be accepted as stable
        * @param maximumArea Maximum area to be accepted as stable
        * @param maxVariation Maximum variation (difference between parent area and region area) to be accepted as stable
        * @param minDiversity Minimum diversity to validate the region
        * @param[out] regions Region that validate the criterion
        */
      void ComputeMSER( const int delta, const int minArea, const int maxArea, const double maxVariation, const double minDiversity, std::vector<MSERRegion> &regions );
    };

} // namespace MSER
} // namespace features
} // namespace openMVG

#endif // OPENMVG_MSER_MSER_REGION_HPP
