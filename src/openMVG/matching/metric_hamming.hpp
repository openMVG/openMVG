
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_METRIC_HAMMING_H
#define OPENMVG_MATCHING_METRIC_HAMMING_H

#include "openMVG/matching/metric.hpp"
#include <bitset>

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
  typedef TBitset ElementType;
  typedef size_t ResultType;

  // Returns the Hamming Distance between two binary descriptors
  template <typename Iterator1, typename Iterator2>
  ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    return (*a ^ *b).count();
  }
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_METRIC_HAMMING_H
