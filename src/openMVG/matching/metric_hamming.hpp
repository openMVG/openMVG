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
    static_assert(std::is_integral<U>::value, "U must be an integral type.");
    return std::bitset<sizeof(U) * 8>(rhs).count();
  }

  // Size must be equal to number of ElementType
  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    ResultType result = 0;

    if (size % sizeof(uint64_t) == 0)
    {
      const uint64_t* pa = reinterpret_cast<const uint64_t*>(a);
      const uint64_t* pb = reinterpret_cast<const uint64_t*>(b);
      size /= (sizeof(uint64_t)/sizeof(unsigned char));
      for (size_t i = 0; i < size; ++i, ++pa, ++pb ) {
        result += popcnt(*pa ^ *pb);
      }
    }
    else if (size % sizeof(uint32_t) == 0)
    {
      const uint32_t* pa = reinterpret_cast<const uint32_t*>(a);
      const uint32_t* pb = reinterpret_cast<const uint32_t*>(b);
      size /= (sizeof(uint32_t)/sizeof(unsigned char));
      for (size_t i = 0; i < size; ++i, ++pa, ++pb ) {
        result += popcnt(*pa ^ *pb);
      }
    }
    else
    {
      const ElementType * a2 = reinterpret_cast<const ElementType*> (a);
      const ElementType * b2 = reinterpret_cast<const ElementType*> (b);
      for (size_t i = 0; i < size / (sizeof(ElementType)); ++i) {
        result += popcnt(a2[i] ^ b2[i]);
      }
    }
    return result;
  }
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_METRIC_HAMMING_HPP
