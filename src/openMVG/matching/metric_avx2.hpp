
// Copyright (c) 2017 Kareem Omar, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/*
*
* Define fast AVX2 & SSEs squared euclidean distance computation for SIFT array
*
* Author: Kareem Omar
* kareem.omar@uah.edu
* https://github.com/komrad36/KfNN
*/

#ifndef OPENMVG_MATCHING_METRIC_AVX2_HPP
#define OPENMVG_MATCHING_METRIC_AVX2_HPP

#ifdef OPENMVG_USE_AVX2
#include <cstdint>
#include <immintrin.h>

namespace openMVG {
namespace matching {

inline int L2_128_uint8_t_AVX2
(
  const uint8_t * a,
  const uint8_t * b,
  size_t size
)
{
  const uint8_t* __restrict tset = reinterpret_cast<const uint8_t* __restrict>(a);
  const uint8_t* __restrict qset = reinterpret_cast<const uint8_t* __restrict>(b);
  __m256i q1 = _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(qset)));
  __m256i q2 = _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(qset + 16)));
  __m256i q3 = _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(qset + 32)));
  __m256i q4 = _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(qset + 48)));
  __m256i q5 = _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(qset + 64)));
  __m256i q6 = _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(qset + 80)));
  __m256i q7 = _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(qset + 96)));
  __m256i q8 = _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(qset + 112)));

  q1 = _mm256_sub_epi16(q1, _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(tset))));
  q1 = _mm256_madd_epi16(q1, q1);
  q2 = _mm256_sub_epi16(q2, _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(tset + 16))));
  q2 = _mm256_madd_epi16(q2, q2);
  q3 = _mm256_sub_epi16(q3, _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(tset + 32))));
  q3 = _mm256_madd_epi16(q3, q3);
  q4 = _mm256_sub_epi16(q4, _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(tset + 48))));
  q4 = _mm256_madd_epi16(q4, q4);
  q5 = _mm256_sub_epi16(q5, _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(tset + 64))));
  q5 = _mm256_madd_epi16(q5, q5);
  q6 = _mm256_sub_epi16(q6, _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(tset + 80))));
  q6 = _mm256_madd_epi16(q6, q6);
  q7 = _mm256_sub_epi16(q7, _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(tset + 96))));
  q7 = _mm256_madd_epi16(q7, q7);
  q8 = _mm256_sub_epi16(q8, _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i*)(tset + 112))));
  q8 = _mm256_madd_epi16(q8, q8);

  q1 = _mm256_add_epi32(q1, q2);
  q3 = _mm256_add_epi32(q3, q4);
  q5 = _mm256_add_epi32(q5, q6);
  q7 = _mm256_add_epi32(q7, q8);

  q1 = _mm256_add_epi32(q1, q3);
  q5 = _mm256_add_epi32(q5, q7);

  q1 = _mm256_add_epi32(q1, q5);

  __m128i x = _mm_add_epi32(_mm256_extractf128_si256(q1, 1), _mm256_castsi256_si128(q1));
  x = _mm_add_epi32(x, _mm_shuffle_epi32(x, 14));
  x = _mm_add_epi32(x, _mm_shuffle_epi32(x, 1));
  return _mm_cvtsi128_si32(x);
}

inline float L2_128_float_AVX2
(
  const float * a,
  const float * b,
  size_t size
)
{
  const float* __restrict tset = reinterpret_cast<const float* __restrict>(a);
  const float* __restrict qset = reinterpret_cast<const float* __restrict>(b);
  __m256 q1 = _mm256_load_ps(qset);
  __m256 q2 = _mm256_load_ps(qset + 8);
  __m256 q3 = _mm256_load_ps(qset + 16);
  __m256 q4 = _mm256_load_ps(qset + 24);
  __m256 q5 = _mm256_load_ps(qset + 32);
  __m256 q6 = _mm256_load_ps(qset + 40);
  __m256 q7 = _mm256_load_ps(qset + 48);
  __m256 q8 = _mm256_load_ps(qset + 56);
  q1 = _mm256_sub_ps(q1, _mm256_load_ps(tset));
  q1 = _mm256_mul_ps(q1, q1);
  q2 = _mm256_sub_ps(q2, _mm256_load_ps(tset + 8));
  q2 = _mm256_mul_ps(q2, q2);
  q1 = _mm256_add_ps(q1, q2);
  q3 = _mm256_sub_ps(q3, _mm256_load_ps(tset + 16));
  q3 = _mm256_mul_ps(q3, q3);
  q4 = _mm256_sub_ps(q4, _mm256_load_ps(tset + 24));
  q4 = _mm256_mul_ps(q4, q4);
  q3 = _mm256_add_ps(q3, q4);
  q1 = _mm256_add_ps(q1, q3);
  q5 = _mm256_sub_ps(q5, _mm256_load_ps(tset + 32));
  q5 = _mm256_mul_ps(q5, q5);
  q6 = _mm256_sub_ps(q6, _mm256_load_ps(tset + 40));
  q6 = _mm256_mul_ps(q6, q6);
  q5 = _mm256_add_ps(q5, q6);
  q7 = _mm256_sub_ps(q7, _mm256_load_ps(tset + 48));
  q7 = _mm256_mul_ps(q7, q7);
  q8 = _mm256_sub_ps(q8, _mm256_load_ps(tset + 56));
  q8 = _mm256_mul_ps(q8, q8);
  q7 = _mm256_add_ps(q7, q8);
  q5 = _mm256_add_ps(q5, q7);
  __m256 res1 = _mm256_add_ps(q1, q5);

  q1 = _mm256_load_ps(qset + 64);
  q2 = _mm256_load_ps(qset + 72);
  q3 = _mm256_load_ps(qset + 80);
  q4 = _mm256_load_ps(qset + 88);
  q5 = _mm256_load_ps(qset + 96);
  q6 = _mm256_load_ps(qset + 104);
  q7 = _mm256_load_ps(qset + 112);
  q8 = _mm256_load_ps(qset + 120);
  q1 = _mm256_sub_ps(q1, _mm256_load_ps(tset + 64));
  q1 = _mm256_mul_ps(q1, q1);
  q2 = _mm256_sub_ps(q2, _mm256_load_ps(tset + 72));
  q2 = _mm256_mul_ps(q2, q2);
  q1 = _mm256_add_ps(q1, q2);
  q3 = _mm256_sub_ps(q3, _mm256_load_ps(tset + 80));
  q3 = _mm256_mul_ps(q3, q3);
  q4 = _mm256_sub_ps(q4, _mm256_load_ps(tset + 88));
  q4 = _mm256_mul_ps(q4, q4);
  q3 = _mm256_add_ps(q3, q4);
  q1 = _mm256_add_ps(q1, q3);
  q5 = _mm256_sub_ps(q5, _mm256_load_ps(tset + 96));
  q5 = _mm256_mul_ps(q5, q5);
  q6 = _mm256_sub_ps(q6, _mm256_load_ps(tset + 104));
  q6 = _mm256_mul_ps(q6, q6);
  q5 = _mm256_add_ps(q5, q6);
  q7 = _mm256_sub_ps(q7, _mm256_load_ps(tset + 112));
  q7 = _mm256_mul_ps(q7, q7);
  q8 = _mm256_sub_ps(q8, _mm256_load_ps(tset + 120));
  q8 = _mm256_mul_ps(q8, q8);
  q7 = _mm256_add_ps(q7, q8);
  q5 = _mm256_add_ps(q5, q7);
  res1 = _mm256_add_ps(res1, _mm256_add_ps(q1, q5));

  __m128 x = _mm_add_ps(_mm256_extractf128_ps(res1, 1), _mm256_castps256_ps128(res1));
  __m128 shuf = _mm_movehdup_ps(x);
  __m128 sums = _mm_add_ps(x, shuf);
  shuf = _mm_movehl_ps(shuf, sums);
  return _mm_cvtss_f32(_mm_add_ss(sums, shuf));
}

}  // namespace matching
}  // namespace openMVG
#endif

#endif // OPENMVG_MATCHING_METRIC_AVX2_HPP

