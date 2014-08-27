
// Copyright (c) 2014 Pierre MOULON.

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

// Hamming distance to work on raw memory
//  like unsigned char *
template<typename T>
struct Hamming
{
  typedef T ElementType;
  typedef unsigned int ResultType;

  /** This is popcount_3() from:
   * http://en.wikipedia.org/wiki/Hamming_weight */
  inline unsigned int popcnt32(uint32_t n) const
  {
#ifdef _MSC_VER
    return __popcnt(n);
#else
    n -= ((n >> 1) & 0x55555555);
    n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
    return (((n + (n >> 4))& 0xF0F0F0F)* 0x1010101) >> 24;
#endif
  }

  inline unsigned int popcnt64(uint64_t n) const
  {
#ifdef _MSC_VER
    return __popcnt64(n);
#else
    n -= ((n >> 1) & 0x5555555555555555LL);
    n = (n & 0x3333333333333333LL) + ((n >> 2) & 0x3333333333333333LL);
    return (((n + (n >> 4))& 0x0f0f0f0f0f0f0f0fLL)* 0x0101010101010101LL) >> 56;
#endif     
  }

  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    ResultType result = 0;
#if (defined __GNUC__ || defined __clang__) && defined USE_SSE 
#ifdef __ARM_NEON__
    {
      uint32x4_t bits = vmovq_n_u32(0);
      for (size_t i = 0; i < size; i += 16) {
        uint8x16_t A_vec = vld1q_u8 (a + i);
        uint8x16_t B_vec = vld1q_u8 (b + i);
        uint8x16_t AxorB = veorq_u8 (A_vec, B_vec);
        uint8x16_t bitsSet = vcntq_u8 (AxorB);
        uint16x8_t bitSet8 = vpaddlq_u8 (bitsSet);
        uint32x4_t bitSet4 = vpaddlq_u16 (bitSet8);
        bits = vaddq_u32(bits, bitSet4);
      }
      uint64x2_t bitSet2 = vpaddlq_u32 (bits);
      result = vgetq_lane_s32 (vreinterpretq_s32_u64(bitSet2),0);
      result += vgetq_lane_s32 (vreinterpretq_s32_u64(bitSet2),2);
    }
#else
    {
      //for portability just use unsigned long -- and use the __builtin_popcountll (see docs for __builtin_popcountll)
      typedef unsigned long long pop_t;
      const size_t modulo = size % sizeof(pop_t);
      const pop_t* a2 = reinterpret_cast<const pop_t*> (a);
      const pop_t* b2 = reinterpret_cast<const pop_t*> (b);
      const pop_t* a2_end = a2 + (size / sizeof(pop_t));

      for (; a2 != a2_end; ++a2, ++b2) result += __builtin_popcountll((*a2) ^ (*b2));

      if (modulo) {
        //in the case where size is not dividable by sizeof(size_t)
        //need to mask off the bits at the end
        pop_t a_final = 0, b_final = 0;
        memcpy(&a_final, a2, modulo);
        memcpy(&b_final, b2, modulo);
        result += __builtin_popcountll(a_final ^ b_final);
      }
    }
#endif //NEON
    return result;
#endif
#ifdef PLATFORM_64_BIT
    if(size%64 == 0)
    {
      const uint64_t* pa = reinterpret_cast<const uint64_t*>(a);
      const uint64_t* pb = reinterpret_cast<const uint64_t*>(b);
      size /= (sizeof(uint64_t)/sizeof(unsigned char));
      for(size_t i = 0; i < size; ++i ) {
        result += popcnt64(*pa ^ *pb);
        ++pa;
        ++pb;
      }
    }
    else
    {
      const uint32_t* pa = reinterpret_cast<const uint32_t*>(a);
      const uint32_t* pb = reinterpret_cast<const uint32_t*>(b);
      size /= (sizeof(uint32_t)/sizeof(unsigned char));
      for(size_t i = 0; i < size; ++i ) {
        result += popcnt32(*pa ^ *pb);
        ++pa;
        ++pb;
      }
    }    
#else
    const uint32_t* pa = reinterpret_cast<const uint32_t*>(a);
    const uint32_t* pb = reinterpret_cast<const uint32_t*>(b);
    size /= (sizeof(uint32_t)/sizeof(unsigned char));
    for(size_t i = 0; i < size; ++i ) {
      result += popcnt32(*pa ^ *pb);
      ++pa;
      ++pb;
    }
#endif
    return result;
  }
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_METRIC_HAMMING_H
