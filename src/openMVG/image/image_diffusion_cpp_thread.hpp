// Copyright (c) 2014 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_DIFFUSION_CPP_THREAD_HPP
#define OPENMVG_IMAGE_IMAGE_DIFFUSION_CPP_THREAD_HPP

#include <thread>
#include <vector>

namespace openMVG
{
  /**
   ** Apply Fast Explicit Diffusion to an Image (on central part)
   ** @param src input image
   ** @param diff diffusion coefficient image
   ** @param half_t Half diffusion time
   ** @param out Output image
   ** @param row_start Row range beginning (range is [row_start ; row_end [ )
   ** @param row_end Row range end (range is [row_start ; row_end [ )
   **/
  template< typename Image >
  void ImageFEDCentral( const Image & src , const Image & diff , const typename Image::Tpixel half_t , Image & out , 
    const int row_start , const int row_end )
  {
    typedef typename Image::Tpixel Real ;

    const int width  = src.Width() ;

    Real n_diff[4] ;
    Real n_src[4] ;

    // Compute FED step on general range
    for( int i = row_start ; i < row_end ; ++i )
    {
      for( int j = 1 ; j < width - 1 ; ++j )
      {
        // Retrieve neighbors : TODO check if we need a cache efficient version ?
        n_diff[0] = diff( i     , j + 1 ) ;
        n_diff[1] = diff( i - 1 , j ) ;
        n_diff[2] = diff( i     , j - 1 ) ;
        n_diff[3] = diff( i + 1 , j ) ;

        n_src[0] = src( i     , j + 1 ) ;
        n_src[1] = src( i - 1 , j ) ;
        n_src[2] = src( i     , j - 1 ) ;
        n_src[3] = src( i + 1 , j ) ;

        // Compute diffusion factor for given pixel
        const Real cur_src  = src( i , j ) ;
        const Real cur_diff = diff( i , j ) ;

        const Real a = ( cur_diff + n_diff[0] ) * ( n_src[0] - cur_src ) ;
        const Real b = ( cur_diff + n_diff[1] ) * ( cur_src - n_src[1] ) ;
        const Real c = ( cur_diff + n_diff[2] ) * ( cur_src - n_src[2] ) ;
        const Real d = ( cur_diff + n_diff[3] ) * ( n_src[3] - cur_src ) ;

        const Real value = half_t * ( a - c + d - b ) ;
        out( i , j ) = value ;
      }
    }
  }


  /**
   ** Apply Fast Explicit Diffusion to an Image (on central part)
   ** @param src input image
   ** @param diff diffusion coefficient image
   ** @param half_t Half diffusion time
   ** @param out Output image
   **/
  template< typename Image >
  void ImageFEDCentralCPPThread( const Image & src , const Image & diff , const typename Image::Tpixel half_t , Image & out )
  {
    const int nb_thread = std::thread::hardware_concurrency() ;

    // Compute ranges
    std::vector< int > range;
    SplitRange( 1 , (int) ( src.rows() - 1 ) , nb_thread , range ) ;

    std::vector< std::thread > threads ;

    // Force use of a function (std::thread use auto and cannot deduce overloading)
    void ( *fct )( const Image & , const Image & , const typename Image::Tpixel , Image & , const int , const int ) = ImageFEDCentral<Image> ;

    for( size_t i = 1 ; i < range.size() ; ++i )  {
      threads.push_back( std::thread( fct , std::cref( src ) , std::cref( diff ) , half_t , std::ref( out ) , range[i-1] , range[i] ) ) ;
    }

    for( size_t i = 0 ; i < threads.size() ; ++i )  {
      threads[i].join() ;
    }
  }

}  // namespace openMVG

#endif // OPENMVG_IMAGE_IMAGE_DIFFUSION_CPP_THREAD_HPP
