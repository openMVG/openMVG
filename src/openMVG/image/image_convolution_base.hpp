// Copyright (c) 2014 Romuald Perrot, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_CONVOLUTION_BASE_HPP
#define OPENMVG_IMAGE_IMAGE_CONVOLUTION_BASE_HPP

#ifdef USE_SSE
#include "openMVG/image/image_convolution_sse.hpp"
#endif


namespace openMVG
{

  /**
   ** Type of border management for convolutions
   **/
  enum EBorderManagement
  {
    BORDER_COPY , // Copy border pixels
    BORDER_CROP   // Crop final image to avoid border mgmt (result is less large than source)
  } ;


  /**
   ** Filter an extended row [halfKernelSize][row][halfKernelSize]
   ** @param buffer data to filter
   ** @param kernel kernel array
   ** @param rsize buffer length
   ** @param ksize kernel length
  **/ 
  template<class T1, class T2> inline
  void conv_buffer_( T1* buffer, const T2* kernel, int rsize, int ksize )
  {
    for ( size_t i = 0; i < rsize; ++i )
    {
      T2 sum( 0 );
      for ( size_t j = 0; j < ksize; ++j )
      {
        sum += buffer[i + j] * kernel[j];
      }
      buffer[i] = sum;
    }
  }

  // float specialization
  template<> inline
  void conv_buffer_( float* buffer, const float * kernel, int rsize, int ksize )
  {
#ifdef USE_SSE
   convolve_sse_fast( buffer , buffer , rsize , kernel , ksize ) ; 
#else
    for ( size_t i = 0; i < rsize; ++i )
    {
      float sum( 0.0f );
      for ( size_t j = 0; j < ksize; ++j )
      {
        sum += buffer[i + j] * kernel[j];
      }
      buffer[i] = sum;
    }
#endif
  }
} // namespace openMVG

#endif // OPENMVG_IMAGE_IMAGE_CONVOLUTION_BASE_HPP
