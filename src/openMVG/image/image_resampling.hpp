// Copyright (c) 2014 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_RESAMPLING_HPP_
#define OPENMVG_IMAGE_IMAGE_RESAMPLING_HPP_


namespace openMVG
{
	/** 
   ** Half sample an image (ie reduce it's size by a factor 2) using bilinear interpolation
   ** @param src input image 
   ** @param out output image 
   **/
  template < typename Image > 
  void ImageHalfSample( const Image & src , Image & out )
  {
    const int new_width  = src.Width() / 2 ;
    const int new_height = src.Height() / 2 ;

    out.resize( new_width , new_height ) ;

    for( int i = 0 ; i < new_height ; ++i )
    {
      for( int j = 0 ; j < new_width ; ++j )
      {
        // Use .5f offset to ensure mid pixel and correct bilinear sampling
        out( i , j ) =
          SampleLinear( src, 2.f * (i+.5f), 2.f * (j+.5f) );
      }
    }
  }

  // TODO : provide a Double size resampling image 

}

#endif