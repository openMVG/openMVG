// Copyright (c) 2014 Pierre Moulon, Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_CONVOLUTION_SSE_HPP
#define OPENMVG_IMAGE_IMAGE_CONVOLUTION_SSE_HPP

#include <Eigen/StdVector>

namespace openMVG
{
  /**
   ** @brief Compute 1d convolution using SSE
   ** @param in Input array
   ** @param out Output array
   ** @param length Length of input arrays
   ** @param kernel Convolution kernel
   ** @param kernel_length Length of kernel
   **/
  static inline void convolve_sse_fast( const float * const in , float * const out , const int length , const float * kernel , const int kernel_length )
  {
    /*
     * Computes 8 values at a time using two accumulators
     */
    int i = 0 ; 
    for( ; i < length - 8 ; i += 8 )
    {
      // Src = start of pixel
      const float * src = in + i ; 

      // Kernel values
      __m128 f ;

      // Accumulators
      __m128 s0 = _mm_setzero_ps() , s1 = s0 ;

      for( int k = 0 ; k < kernel_length ; ++k , ++src )
      {
        // Load kernel value
        // Faster than _mm_set1_ps()
        f = _mm_load_ss( kernel + k ) ;
        f = _mm_shuffle_ps( f , f , 0 ) ; 

        s0 = _mm_add_ps( s0 , _mm_mul_ps( _mm_loadu_ps( src ) , f ) ) ;
        s1 = _mm_add_ps( s1 , _mm_mul_ps( _mm_loadu_ps( src + 4 ) , f ) ) ;  
      } 
      _mm_store_ps( out + i , s0 ) ;
      _mm_store_ps( out + i + 4 , s1 ) ; 
    }

    // Handle the last values
    for( ; i < length ; ++i ) {
      float sum = 0.f ;
      for( int k = 0; k < kernel_length; ++k )  {
        sum += in[i + k] * kernel[k];
      }
      out[i] = sum ;
    }
  }

} // namespace openMVG

#endif // OPENMVG_IMAGE_IMAGE_CONVOLUTION_SSE_HPP

