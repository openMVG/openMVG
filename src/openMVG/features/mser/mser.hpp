// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_MSER_MSER_HPP
#define OPENMVG_FEATURES_MSER_MSER_HPP

#include <vector>
namespace openMVG { namespace image { template <typename T> class Image; } }

/**
* @note Implementation based of the D. Nister method : "Linear Time Maximally Stable Extremal Regions" [1]
*/
//------------------
//-- Bibliography --
//------------------
//- [1] "Linear Time Maximally Stable Extremal Regions"
//- Authors: David Nistér, Henrik Stewénius
//- Date: 2008.
//- Conference: ECCV.

namespace openMVG
{
  namespace features
  {
    namespace MSER
    {

      class MSERRegion;

      /**
      * @brief Function used to extract MSER regions of an image
      * @note This should used in a two step approach :
      * MSERExtractor extr8( ... , MSER_8_CONNECTIVITY );
      * MSERExtractor extr4( ... , MSER_4_CONNECTIVITY );
      * Image img;
      * std::vector< MSERRegion > regs;
      * extr8( img , regs );
      * extr4( ~img , regs ); // ~is the inverse color image
      *
      * We do this because MSER is extracted from the lowest grayscale level to the upper and so is extracts only the black regions
      */
      class MSERExtractor
      {
      public:

        /*
        * Connectivity used to find neighboring of a pixel
        *
        * 4-connect :
                 |   | 1 |   |
                 | 2 | X | 0 |
                 |   | 3 |   |
        * 8-connect :
                 | 3 | 2 | 1 |
                 | 4 | X | 0 |
                 | 5 | 6 | 7 |
        */
        static const int MSER_4_CONNECTIVITY;
        static const int MSER_8_CONNECTIVITY;


        /**
        * @brief MSER extractor feature
        * @param delta Distance between levels used to check stability (area(i) - area(i - delta) ) / area( i ) )
        * @param min_area Minimum area of the extracted region (relative to the image area)
        * @param max_area Maximum area of the extracted region (relative to the image area)
        * @param max_variation Maximum difference (in stability score) between the region
        * @param min_diversity Threshold to keep only the more stable region (when one region is inside another)
        * @param connectivity Kind of connectivity to process region (MSER_4_CONNECTIVITY or MSER_8_CONNECTIVITY)
        */
        MSERExtractor( const int delta = 2 ,
                       const double min_area = 0.0001 , const double max_area = 0.5 ,
                       const double max_variation = 0.5 ,
                       const double min_diversity = 1.0 / 3.0 ,
                       const int connectivity = MSER_4_CONNECTIVITY );

        /**
        * @brief Extract MSER regions
        * @param img Input image
        * @param[out] regions Output regions
        */
        void Extract( const image::Image<unsigned char> & img , std::vector< MSERRegion > & regions ) const;

      private:

        /**
        * Try to merge region with same level
        * @param nextLevel Level used to limit search
        * @param pixel_x X-coord of the base of the merged region
        * @param pixel_y Y-coord of the base of the merged region
        */
        void ProcessStack( const int nextLevel , const int pixel_x , const int pixel_y , std::vector< MSERRegion * > & regionStack ) const;

        int m_delta; // Maximum level distance to check stability
        double m_minimum_area; // Minimum area (relative to the image) of the valid regions
        double m_maximum_area; // Maximum area (relative to the image) of the valid regions
        double m_max_variation; // Maximum variation between two regions
        double m_min_diversity; // Stability distance between two region in the same hierarchy
        int m_connectivity;
      };
    } // namespace MSER
  } // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_MSER_MSER_HPP
