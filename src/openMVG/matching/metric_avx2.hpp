// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/*
*
* Define fast AVX2 squared euclidean distance computation for SIFT array
*/

#ifndef OPENMVG_MATCHING_METRIC_AVX2_HPP
#define OPENMVG_MATCHING_METRIC_AVX2_HPP

#include <array>
#include <numeric>

#ifdef OPENMVG_USE_AVX2
#include <cstdint>
#include <immintrin.h>

namespace openMVG {
namespace matching {

#ifdef _MSC_VER
#define ALIGNED32 __declspec(align(32))
#else
#define ALIGNED32 __attribute__((aligned(32)))
#endif

inline int L2_AVX2
(
  const uint8_t * a,
  const uint8_t * b,
  size_t size
)
{
  // Accumulator
  __m256i acc (_mm256_setzero_si256());

  const ALIGNED32 __m256i* ad = reinterpret_cast<const ALIGNED32 __m256i*>(a);
  const ALIGNED32 __m256i* bd = reinterpret_cast<const ALIGNED32 __m256i*>(b);

  // Compute (A-B) * (A-B) on 32 components per iteration
  for (int i = 0; i < 4; ++i) {
    // In order to avoid overflow, process low and high order value
    const __m256i min = _mm256_min_epu8(ad[i], bd[i]);
    const __m256i max = _mm256_max_epu8(ad[i], bd[i]);
    const __m256i d = _mm256_sub_epi8(max, min);

    // Squared elements in range [0,15]
    __m256i dl = _mm256_unpacklo_epi8(d, _mm256_setzero_si256());
    dl = _mm256_madd_epi16(dl, dl);
    // Squared elements in range [15,31]
    __m256i dh = _mm256_unpackhi_epi8(d, _mm256_setzero_si256());
    dh = _mm256_madd_epi16(dh, dh);
    acc = _mm256_add_epi32(acc, _mm256_add_epi32(dl, dh));
  }
  // Compute the sum in the accumulator
  __m128i l = _mm256_extracti128_si256(acc, 0);
  __m128i h = _mm256_extracti128_si256(acc, 1);
  __m128i r = _mm_hadd_epi32(_mm_add_epi32(h, l), _mm_setzero_si128());
  return _mm_extract_epi32(r, 0) + _mm_extract_epi32(r, 1);
}

inline float L2_AVX2
(
  const float * a,
  const float * b,
  size_t size
)
{
  // Accumulator
  __m256 acc (_mm256_setzero_ps());

  // Compute (A-B) * (A-B) on 16 components per iteration
  for (int j = 0; j < size; j += 8)
  {
    const __m256 t0 = _mm256_sub_ps(_mm256_loadu_ps(a + j), _mm256_loadu_ps(b + j));
    acc = _mm256_add_ps(acc, _mm256_mul_ps(t0, t0));
  }
  float ALIGNED32 acc_float[8];
  _mm256_store_ps(acc_float, acc);
  return std::accumulate(acc_float, acc_float + 8, 0.f);
}

}  // namespace matching
}  // namespace openMVG
#endif

#endif // OPENMVG_MATCHING_METRIC_AVX2_HPP
