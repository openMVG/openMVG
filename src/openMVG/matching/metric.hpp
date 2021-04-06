// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_METRIC_HPP
#define OPENMVG_MATCHING_METRIC_HPP

#include "openMVG/matching/metric_simd.hpp"
#include "openMVG/matching/metric_hamming.hpp"
#include "openMVG/numeric/accumulator_trait.hpp"
#include <cstdint>

namespace openMVG {
namespace matching {

/// Squared Euclidean distance functor
template<class T>
struct L2
{
  using ElementType = T;
  using ResultType = typename Accumulator<ElementType>::Type;

  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    ResultType result = ResultType();
    ResultType diff0, diff1, diff2, diff3;
    Iterator1 last = a + size;
    Iterator1 lastgroup = last - 3;

    // Process 4 items for each loop for efficiency.
    while (a < lastgroup) {
      diff0 = a[0] - b[0];
      diff1 = a[1] - b[1];
      diff2 = a[2] - b[2];
      diff3 = a[3] - b[3];
      result += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
      a += 4;
      b += 4;
    }
    // Process last 0-3 elements.  Not needed for standard vector lengths.
    while (a < last) {
      diff0 = *a++ - *b++;
      result += diff0 * diff0;
    }
    return result;
  }
};

// Template specialization for the uint8_t type
template<>
struct L2<uint8_t>
{
  using ElementType = uint8_t;
  using ResultType = int;

  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    #ifdef OPENMVG_USE_AVX2
    if (size == 128)
    {
      return L2_AVX2(a, b, size);
    }
    #endif

    ResultType result = ResultType();
    ResultType diff0, diff1, diff2, diff3;
    Iterator1 last = a + size;
    Iterator1 lastgroup = last - 3;

    // Process 4 items for each loop for efficiency.
    while (a < lastgroup) {
      diff0 = a[0] - b[0];
      diff1 = a[1] - b[1];
      diff2 = a[2] - b[2];
      diff3 = a[3] - b[3];
      result += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
      a += 4;
      b += 4;
    }
    // Process last 0-3 elements.  Not needed for standard vector lengths.
    while (a < last) {
      diff0 = *a++ - *b++;
      result += diff0 * diff0;
    }
    return result;
  }
};

// Template specialization for the float type
template<>
struct L2<float>
{
  using ElementType = float;
  using ResultType = typename Accumulator<ElementType>::Type;

  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    #ifdef OPENMVG_USE_AVX
    if (size > 4 && size % 8)
    {
      return L2_AVX(a, b, size);
    }
    #endif

    ResultType result = ResultType();
    ResultType diff0, diff1, diff2, diff3;
    Iterator1 last = a + size;
    Iterator1 lastgroup = last - 3;

    // Process 4 items for each loop for efficiency.
    while (a < lastgroup) {
      diff0 = a[0] - b[0];
      diff1 = a[1] - b[1];
      diff2 = a[2] - b[2];
      diff3 = a[3] - b[3];
      result += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
      a += 4;
      b += 4;
    }
    // Process last 0-3 elements.  Not needed for standard vector lengths.
    while (a < last) {
      diff0 = *a++ - *b++;
      result += diff0 * diff0;
    }
    return result;
  }
};

template <class T>
struct L1
{
  using ElementType = T;
  using ResultType = typename Accumulator<ElementType>::Type;

  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    ResultType result = ResultType();
    Iterator1 last = a + size;
    Iterator1 lastgroup = last - 3;

    // Process 4 items for each loop for efficiency.
    while (a < lastgroup) {
      result += std::abs(a[0] - b[0]);
      result += std::abs(a[1] - b[1]);
      result += std::abs(a[2] - b[2]);
      result += std::abs(a[3] - b[3]);
      a += 4;
      b += 4;
    }
    // Process last 0-3 elements.  Not needed for standard vector lengths.
    while (a < last) {
      result += std::abs(*a++ - *b++);
    }
    return result;
  }
};

template<>
struct L1<uint8_t>
{
  using ElementType = uint8_t;
  using ResultType = int;

  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    ResultType result = ResultType();
    Iterator1 last = a + size;
    Iterator1 lastgroup = last - 3;

    // Process 4 items for each loop for efficiency.
    while (a < lastgroup) {
      result += std::abs(a[0] - b[0]);
      result += std::abs(a[1] - b[1]);
      result += std::abs(a[2] - b[2]);
      result += std::abs(a[3] - b[3]);
      a += 4;
      b += 4;
    }
    // Process last 0-3 elements.  Not needed for standard vector lengths.
    while (a < last) {
      result += std::abs(*a++ - *b++);
    }
    return result;
  }
};

/// Inner product similarity. It's equivalent to the cosine similarity for normalized vector.
/// The inner product is multiplied by -1 for sorting purpose w.r.t. other metrics.
template<class T>
struct LInner
{
  using ElementType = T;
  using ResultType = typename Accumulator<ElementType>::Type;

  template <typename Iterator1, typename Iterator2>
  inline ResultType operator()(Iterator1 a, Iterator2 b, size_t size) const
  {
    ResultType result = ResultType();
    ResultType prod0, prod1, prod2, prod3;
    Iterator1 last = a + size;
    Iterator1 lastgroup = last - 3;

    // Process 4 items for each loop for efficiency.
    while (a < lastgroup) {
      prod0 = a[0] * b[0];
      prod1 = a[1] * b[1];
      prod2 = a[2] * b[2];
      prod3 = a[3] * b[3];
      result += prod0 + prod1 + prod2 + prod3;
      a += 4;
      b += 4;
    }
    // Process last 0-3 elements.  Not needed for standard vector lengths.
    while (a < last) {
      prod0 = *a++ * *b++;
      result += prod0;
    }
    return ElementType(-1) * result; // TODO: inner product is a similarity metric, not a distance
  }
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_METRIC_HPP
