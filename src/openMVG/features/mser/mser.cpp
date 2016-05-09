// Copyright (c) 2015 Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "mser.hpp"

#include <iostream>
#include <stack>

namespace openMVG
{
  namespace features
  {
    namespace MSER
    {

      const int MSERExtractor::MSER_4_CONNECTIVITY = 0 ;
      const int MSERExtractor::MSER_8_CONNECTIVITY = 1 ;


      // Ordering of neighboring
      /*
          3     2     1
          4           0
          5     6     7
      */
      enum PixelNeighborsDirection
      {
        PIXEL_RIGHT = 0 ,
        PIXEL_TOP_RIGHT    = 1 ,
        PIXEL_TOP = 2 ,
        PIXEL_TOP_LEFT = 3 ,
        PIXEL_LEFT = 4 ,
        PIXEL_BOTTOM_LEFT = 5 ,
        PIXEL_BOTTOM = 6 ,
        PIXEL_BOTTOM_RIGHT = 7
      } ;


      /**
      * @brief Given preceding direction used to search a pixel, get the next pixel direction
      * @param dir Preceding direction
      * @param neighboring_type Kind of connectivity used (MSER_4_CONNECTIVITY or MSER_8_CONNECTIVITY)
      */
      static inline PixelNeighborsDirection NextDirection( const PixelNeighborsDirection dir , const int neighboring_type )
      {
        if( neighboring_type == MSERExtractor::MSER_4_CONNECTIVITY )
        {
          // Assuming dir is even
          return PixelNeighborsDirection( ( static_cast<int>( dir ) + 2 ) ) ;
        }
        else // neighboring_type == MSERExtractor::MSER_8_CONNECTIVITY
        {
          return PixelNeighborsDirection( ( static_cast<int>( dir ) + 1 ) ) ;
        }
      }

      /**
      * @brief Indicate if pixel is inside the image range
      * @param x X-coord of the tested pixel
      * @param y Y-coord of the tested pixel
      * @param img_width Width of the image
      * @param img_height Height of the image
      * @retval true if the pixel is inside the image
      * @retval false if the pixel is outside the image
      */
      static inline bool ValidPixel( const int x , const int y ,
                                     const int img_width , const int img_height )
      {
        return x >= 0 && y >= 0 && x < img_width && y < img_height ;
      }

      /*
      * Given a pixel and image dimension, retrive it's neighbor
      * If neighbor is out of the image range, return (-1,-1) coord
      * @param x Input x pixel coord
      * @param y Input y pixel coord
      * @param dir Neighbor direction to retrieve
      * @param img_width Width of the image
      * @param img_height Height of the image
      * @param[out] nx Neighbor x pixel coord
      * @param[out] ny Neighbor y pixel coord
      * @return true if input params are all good | false if dir is incorrect
      */
      static inline bool GetNeighbor( const int x , const int y ,
                                      const PixelNeighborsDirection dir ,
                                      const int img_width , const int img_height ,
                                      int & nx , int & ny )
      {
        nx = x ;
        ny = y ;

        switch( dir )
        {
        case PIXEL_TOP_LEFT :
        {
          nx -= 1 ;
          ny -= 1 ;
          break ;
        }
        case PIXEL_TOP :
        {
          ny -= 1 ;
          break ;
        }
        case PIXEL_TOP_RIGHT :
        {
          nx += 1 ;
          ny -= 1 ;
          break ;
        }
        case PIXEL_LEFT :
        {
          nx -= 1 ;
          break ;
        }
        case PIXEL_RIGHT:
        {
          nx += 1 ;
          break ;
        }
        case PIXEL_BOTTOM_LEFT:
        {
          nx -= 1 ;
          ny += 1 ;
          break ;
        }
        case PIXEL_BOTTOM :
        {
          ny += 1 ;
          break ;
        }
        case PIXEL_BOTTOM_RIGHT :
        {
          nx += 1 ;
          ny += 1 ;
          break ;
        }
        default :
        {
          std::cerr << "Unhandled pixel direction" << std::endl ;
          return false ;
        }
        }

        if( ! ValidPixel( nx , ny , img_width , img_height ) )
        {
          nx = -1 ;
          ny = -1 ;
        }
        return true ;
      }

      /**
      * @brief Pixel info (for easy manipulation)
      */
      struct PixelStackElt
      {
        // Position of the processed pixel
        int pix_x ;
        int pix_y ;
        // Level of the processed pixel
        int pix_level ;
        // Current neighhobor to process
        PixelNeighborsDirection edge_index ;
      } ;


      /**
      * @brief MSER extractor feature
      * @param delta Distance between levels used to check stability (area(i) - area(i - delta) ) / area( i ) )
      * @param min_area Minimum area of the extracted region (relative to the image area)
      * @param max_area Maximum area of the extracted region (relative to the image area)
      * @param max_variation Maximum difference (in stability score) between the region
      * @param min_diversity Threshold to keep only the more stable region (when one region is inside another)
      */
      MSERExtractor::MSERExtractor( const int delta ,
                                    const double min_area , const double max_area ,
                                    const double max_variation ,
                                    const double min_diversity ,
                                    const int connectivity )
        :  m_delta( delta ) ,
           m_minimum_area( min_area ) ,
           m_maximum_area( max_area ) ,
           m_max_variation( max_variation ) ,
           m_min_diversity( min_diversity ) ,
           m_connectivity( connectivity )
      {

      }

      /**
      * @brief Extract MSER regions
      * @param img Input image
      * @param[out] regions Output regions
      */
      void MSERExtractor::Extract( const image::Image<unsigned char> & img , std::vector< MSERRegion > & regions ) const
      {
        // Compute minimum and maximum region area relative to this image
        const int minRegArea = img.Width() * img.Height() * m_minimum_area ;
        const int maxRegArea = img.Width() * img.Height() * m_maximum_area ;

        // List of processed pixels (maybe we can use a more efficient structure)
        std::vector< std::vector< bool > > processed ;
        processed.resize( img.Width() ) ;
        for( int i = 0 ; i < img.Width() ; ++i )
        {
          processed[ i ].resize( img.Height() ) ;
          std::fill( processed[ i ].begin() , processed[ i ].end() , false ) ;
        }

        // Holds the boundary of given grayscale value (boundary[0] -> pixels in the boundary with 0 grayscale value)
        std::vector< PixelStackElt > boundary[ 256 ] ;

        // List of regions computed so far (not only valid MSER regions)
        std::vector< MSERRegion * > regionStack ;

        // Push en empty region
        regionStack.push_back( new MSERRegion ) ;

        // Start processing from top left pixel
        PixelStackElt cur_pix ;
        cur_pix.pix_x = 0 ;
        cur_pix.pix_y = 0 ;
        cur_pix.pix_level = img( 0 , 0 ) ;
        cur_pix.edge_index = PIXEL_RIGHT ;

        processed[ cur_pix.pix_x ][ cur_pix.pix_y ] = true ;

        regionStack.push_back( new MSERRegion( cur_pix.pix_level , cur_pix.pix_x , cur_pix.pix_y ) ) ;

        int priority = 256 ;

        // Start process
        while( 1 )
        {
          bool restart = false ;

          // Process neighboring to see if there's something to search with lower grayscale level
          for(  PixelNeighborsDirection curDir = cur_pix.edge_index ;
                curDir <= PIXEL_BOTTOM_RIGHT ;
                curDir = NextDirection( curDir , m_connectivity ) )
          {
            int nx , ny ;
            GetNeighbor( cur_pix.pix_x , cur_pix.pix_y , curDir , img.Width() , img.Height() , nx , ny ) ;

            // Pixel was not processed before
            if( ValidPixel( nx , ny , img.Width() , img.Height() ) && ! processed[ nx ][ ny ] )
            {
              const int nLevel = img( ny , nx ) ;
              processed[ nx ][ ny ] = true ;

              // Info of the neighboring pixel
              PixelStackElt n_elt ;
              n_elt.pix_x = nx ;
              n_elt.pix_y = ny ;
              n_elt.pix_level = nLevel ;
              n_elt.edge_index = PIXEL_RIGHT ;

              // Now look from which pixel do we have to continue
              if( nLevel >= cur_pix.pix_level )
              {
                // Continue from the same pixel
                boundary[ nLevel ].push_back( n_elt ) ;

                // Store the lowest value so far
                priority = std::min( nLevel , priority ) ;
              }
              else
              {
                // Go on with the neighboring pixel (go down)
                cur_pix.edge_index = NextDirection( curDir , m_connectivity ) ; // Next time we have to process the next boundary pixel
                boundary[ cur_pix.pix_level ].push_back( cur_pix ) ;

                // Store the lowest value so far
                priority = std::min( cur_pix.pix_level , priority ) ;

                // Push the next pixel to process
                cur_pix = n_elt ;
                restart = true ;
                break ;
              }
            }
          }
          // Do we have to restart from a new pixel ?
          if( restart )
          {
            // If so it's that because we found a lower grayscale value so let's start a new region
            regionStack.push_back( new MSERRegion( cur_pix.pix_level , cur_pix.pix_x , cur_pix.pix_y ) ) ;
            continue ;
          }

          // We have process all the neighboring pixels, current pixel is the lowest we have found so far
          // now process the current pixel
          regionStack.back()->AppendPixel( cur_pix.pix_x , cur_pix.pix_y ) ;

          // End of the process : we have no boundary region, compute MSER from graph
          if( priority == 256 )
          {
            regionStack.back()->ComputeMSER( m_delta , minRegArea , maxRegArea , m_max_variation , m_min_diversity , regions ) ;
            break ;
          }

          PixelStackElt next_pix = boundary[ priority ].back() ;
          boundary[ priority ].pop_back() ;

          // Get the next pixel level
          while( boundary[ priority ].empty() && ( priority < 256 ) )
          {
            ++priority ;
          }

          // Clear the stack
          int newLevel = next_pix.pix_level ;

          // Process the current stack of pixels if the next processing pixel is not at the same curent level
          if( newLevel != cur_pix.pix_level )
          {
            // Try to merge the regions to fomr a tree
            ProcessStack( newLevel , next_pix.pix_x , next_pix.pix_y , regionStack ) ;
          }

          // Update next pixel for processing
          cur_pix = next_pix ;
        }

        // Clear region stack created so far
        for( size_t i = 0 ; i < regionStack.size() ; ++i )
        {
          delete regionStack[ i ] ;
        }
      }

      /**
      * Try to merge region with same level
      * @param nextLevel Level used to limit search
      * @param pixel_x X-coord of the base of the merged region
      * @param pixel_y Y-coord of the base of the merged region
      */
      void MSERExtractor::ProcessStack( const int nextLevel , const int pixel_x , const int pixel_y , std::vector< MSERRegion * > & regionStack ) const
      {
        do
        {
          MSERRegion * top = regionStack.back() ;
          regionStack.pop_back() ;

          if( nextLevel < regionStack.back()->m_level )
          {
            regionStack.push_back( new MSERRegion( nextLevel , pixel_x , pixel_y ) ) ;

            regionStack.back()->MergeRegion( top ) ;
            return ;
          }

          regionStack.back()->MergeRegion( top ) ;


        }
        while( nextLevel > regionStack.back()->m_level ) ;
      }


    } // namespace MSER
  } // namespace features
} // namespace openMVG
