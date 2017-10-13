// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014-2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_METRIC_HAMMING_HPP
#define OPENMVG_MATCHING_METRIC_HAMMING_HPP

#include "openMVG/matching/metric.hpp"

#include <bitset>
#include <cstdint>
#include <type_traits>

#ifdef _MSC_VER
#include "nmmintrin.h"
#endif

// Brief:
// Hamming distance count the number of bits in common between descriptors
//  by using a XOR operation + a count.
// For maximal performance SSE4 must be enable for builtin popcount activation.

namespace openMVG {
namespace matching {

/// Hamming distance:
///  Working for STL fixed size BITSET and boost DYNAMIC_BITSET
template<typename TBitset>
struct HammingBitSet
{
  using ElementType = TBitset;
  using ResultType = size_t;

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
  using ElementType = T;
  using ResultType = unsigned int;

  template<typename U>
  static inline std::size_t constexpr popcnt(const U & rhs)
  {
    static_assert(std::is_integral<U>::value && std::is_unsigned<T>::value,
      "U must be an unsigned integral type.");
    return std::bitset<sizeof(U) * 8>(rhs).count();
  }

#ifdef _MSC_VER
  inline std::size_t constexpr popcnt(const uint32_t & rhs)
  {
    return _mm_popcnt_u32(rhs);
  }

  inline std::size_t constexpr popcnt(const uint64_t & rhs)
  {
#if __amd64__ || __x86_64__ || _WIN64 || _M_X64
    return _mm_popcnt_u64(rhs);
#else
    // Process low and high bits
    return popcnt(static_cast<std::uint32_t>(rhs)) +
           popcnt(static_cast<std::uint32_t>(rhs >> 32));
#endif
  }
#endif

  template <typename U>
  inline ResultType popcntLoop(const U * pa, const U * pb, size_t size) const
  {
    ResultType result = 0;
    size /= (sizeof(U) / sizeof(uint8_t));
    for (size_t i = 0; i < size; ++i) {
      result += popcnt(pa[i] ^ pb[i]);
    }
    return result;
  }

  // Size must be equal to number of ElementType
  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    if (size % sizeof(uint64_t) == 0)
    {
      const uint64_t* pa = reinterpret_cast<const uint64_t*>(a);
      const uint64_t* pb = reinterpret_cast<const uint64_t*>(b);
      return popcntLoop(pa, pb, size);
    }
    else if (size % sizeof(uint32_t) == 0)
    {
      const uint32_t* pa = reinterpret_cast<const uint32_t*>(a);
      const uint32_t* pb = reinterpret_cast<const uint32_t*>(b);
      return popcntLoop(pa, pb, size);
    }
    else
    {
      const ElementType * pa = reinterpret_cast<const ElementType*> (a);
      const ElementType * pb = reinterpret_cast<const ElementType*> (b);
      return popcntLoop(pa, pb, size);
    }
    return 0;
  }
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_METRIC_HAMMING_HPP
