
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_METRIC_H
#define OPENMVG_MATCHING_METRIC_H

#include "openMVG/matching/metric_hamming.hpp"

namespace openMVG {
namespace matching {

template<typename T>
struct Accumulator { typedef T Type; };
template<>
struct Accumulator<unsigned char>  { typedef float Type; };
template<>
struct Accumulator<unsigned short> { typedef float Type; };
template<>
struct Accumulator<unsigned int> { typedef float Type; };
template<>
struct Accumulator<char>   { typedef float Type; };
template<>
struct Accumulator<short>  { typedef float Type; };
template<>
struct Accumulator<int> { typedef float Type; };

/// Squared Euclidean distance functor.
template<class T>
struct L2_Simple
{
    typedef T ElementType;
    typedef typename Accumulator<T>::Type ResultType;

    template <typename Iterator1, typename Iterator2>
    inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
    {
        ResultType result = ResultType();
        ResultType diff;
        for(size_t i = 0; i < size; ++i ) {
            diff = *a++ - *b++;
            result += diff*diff;
        }
        return result;
    }
};

/// Squared Euclidean distance functor (vectorized version)
template<class T>
struct L2_Vectorized
{
    typedef T ElementType;
    typedef typename Accumulator<T>::Type ResultType;

    template <typename Iterator1, typename Iterator2>
    inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
    {
      ResultType result = ResultType();
      ResultType diff0, diff1, diff2, diff3;
      Iterator1 last = a + size;
      Iterator1 lastgroup = last - 3;

      // Process 4 items with each loop for efficiency.
      while (a < lastgroup) {
          diff0 = a[0] - b[0];
          diff1 = a[1] - b[1];
          diff2 = a[2] - b[2];
          diff3 = a[3] - b[3];
          result += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
          a += 4;
          b += 4;
      }
      // Process last 0-3 pixels.  Not needed for standard vector lengths.
      while (a < last) {
          diff0 = *a++ - *b++;
          result += diff0 * diff0;
      }
      return result;
    }
};

#ifdef USE_SSE

namespace optim_ss2{

#include <xmmintrin.h>

  /// Union to switch between SSE and float array
  union sseRegisterHelper
  {
      __m128 m;
      float f[4];
  };

  // Euclidean distance (SSE method) (squared result)
  inline float l2_sse(float * b1, float * b2, int size)
  {
    float* b1Pt = (float*)b1;
    float* b2Pt = (float*)b2;
    if(size%4 == 0)
    {
      __m128 srcA, srcB, temp, cumSum;
      float zeros[4] = {0.f,0.f,0.f,0.f};
      cumSum = _mm_load_ps( zeros );
      for(int i = 0 ; i < size; i+=4)
      {
        srcA = _mm_load_ps(b1Pt+i);
        srcB = _mm_load_ps(b2Pt+i);
        //-- Subtract
        temp = _mm_sub_ps( srcA, srcB );
        //-- Multiply
        temp =  _mm_mul_ps( temp, temp );
        //-- sum
        cumSum = _mm_add_ps( cumSum, temp );
      }
      sseRegisterHelper res;
      res.m = cumSum;
      return (res.f[0]+res.f[1]+res.f[2]+res.f[3]);
    }
    else
    {
      std::cout <<"\n/!\\ size is not modulus 4,"
        << " distance cannot be performed in SSE"<< std::endl;
      return 0.0f;
    }
  }
} // namespace optim_ss2

// Template specification to run SSE L2 squared distance
//  on float vector
template<>
struct L2_Vectorized<float>
{
    typedef float ElementType;
    typedef Accumulator<float>::Type ResultType;

    template <typename Iterator1, typename Iterator2>
    inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
    {
      return optim_ss2::l2_sse(a,b,size);
    }
};

#endif // USE_SSE

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_METRIC_H
