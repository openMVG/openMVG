
// Copyright (c) 2014-2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_METRIC_HAMMING_H
#define OPENMVG_MATCHING_METRIC_HAMMING_H

#include "openMVG/matching/metric.hpp"
#include <bitset>

#ifdef _MSC_VER
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
#include <intrin.h>
#else
#include <stdint.h>
#endif

#ifdef __ARM_NEON__
#include "arm_neon.h"
#endif

// Brief:
// Hamming distance count the number of bits in common between descriptors
//  by using a XOR operation + a count.
// For maximal performance SSE4 must be enable for builtin popcount activation.

namespace openMVG {
namespace matching {

#undef PLATFORM_64_BIT
#undef PLATFORM_32_BIT
#if __amd64__ || __x86_64__ || _WIN64 || _M_X64
#define PLATFORM_64_BIT
#else
#define PLATFORM_32_BIT
#endif

/// Hamming distance:
///  Working for STL fixed size BITSET and boost DYNAMIC_BITSET
template<typename TBitset>
struct HammingBitSet
{
  typedef TBitset ElementType;
  typedef size_t ResultType;

  // Returns the Hamming Distance between two binary descriptors
  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    return (*a ^ *b).count();
  }
};

// https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable
// Lookup table to count the number of common 1 bits on unsigned char values
static const unsigned char pop_count_LUT[256] =
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};

// Hamming distance to work on raw memory
//  like unsigned char *
template<typename T>
struct Hamming
{
  typedef T ElementType;
  typedef unsigned int ResultType;

  /** This is popcount_3() from:
   * http://en.wikipedia.org/wiki/Hamming_weight */
  static inline unsigned int popcnt32(uint32_t n)
  {
#ifdef _MSC_VER
    return __popcnt(n);
#else
#if (defined __GNUC__ || defined __clang__)
    return __builtin_popcountl(n);
#endif
    n -= ((n >> 1) & 0x55555555);
    n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
    return (((n + (n >> 4))& 0xF0F0F0F)* 0x1010101) >> 24;
#endif
  }

  static inline unsigned int popcnt64(uint64_t n)
  {
#if defined _MSC_VER && defined PLATFORM_64_BIT
    return __popcnt64(n);
#else
#if (defined __GNUC__ || defined __clang__)
    return __builtin_popcountll(n);
#endif
    n -= ((n >> 1) & 0x5555555555555555LL);
    n = (n & 0x3333333333333333LL) + ((n >> 2) & 0x3333333333333333LL);
    return (((n + (n >> 4))& 0x0f0f0f0f0f0f0f0fLL)* 0x0101010101010101LL) >> 56;
#endif
  }

  // Size must be equal to number of ElementType
  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    ResultType result = 0;
// Windows & generic platforms:

#ifdef PLATFORM_64_BIT
    if(size%sizeof(uint64_t) == 0)
    {
      const uint64_t* pa = reinterpret_cast<const uint64_t*>(a);
      const uint64_t* pb = reinterpret_cast<const uint64_t*>(b);
      size /= (sizeof(uint64_t)/sizeof(unsigned char));
      for(size_t i = 0; i < size; ++i, ++pa, ++pb ) {
        result += popcnt64(*pa ^ *pb);
      }
    }
    else if(size%sizeof(uint32_t) == 0)
    {
      const uint32_t* pa = reinterpret_cast<const uint32_t*>(a);
      const uint32_t* pb = reinterpret_cast<const uint32_t*>(b);
      size /= (sizeof(uint32_t)/sizeof(unsigned char));
      for(size_t i = 0; i < size; ++i, ++pa, ++pb ) {
        result += popcnt32(*pa ^ *pb);
      }
    }
    else
    {
      const ElementType * a2 = reinterpret_cast<const ElementType*> (a);
      const ElementType * b2 = reinterpret_cast<const ElementType*> (b);
      for (size_t i = 0;
           i < size / (sizeof(unsigned char)); ++i) {
        result += pop_count_LUT[a2[i] ^ b2[i]];
      }
    }
#else // PLATFORM_64_BIT
    if(size%sizeof(uint32_t) == 0)
    {
      const uint32_t* pa = reinterpret_cast<const uint32_t*>(a);
      const uint32_t* pb = reinterpret_cast<const uint32_t*>(b);
      size /= (sizeof(uint32_t)/sizeof(unsigned char));
      for(size_t i = 0; i < size; ++i, ++pa, ++pb ) {
        result += popcnt32(*pa ^ *pb);
      }
    }
    else
    {
      const ElementType * a2 = reinterpret_cast<const ElementType*> (a);
      const ElementType * b2 = reinterpret_cast<const ElementType*> (b);
      for (size_t i = 0;
           i < size / (sizeof(unsigned char)); ++i) {
        result += pop_count_LUT[a2[i] ^ b2[i]];
      }
    }
#endif // PLATFORM_64_BIT
    return result;
  }
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_METRIC_HAMMING_H
