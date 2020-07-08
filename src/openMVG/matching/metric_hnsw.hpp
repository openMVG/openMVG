// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Romain JANVIER

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_METRIC_HNSW_HPP
#define OPENMVG_MATCHING_METRIC_HNSW_HPP

#include "openMVG/matching/metric_hamming.hpp"

#include "third_party/hnswlib/hnswlib.h"

/*
* This file define specialized HNSW kernels for differents metrics/spaces
*/
namespace openMVG {
namespace matching {
namespace custom_hnsw{

template <typename U>
static unsigned int HammingKernel(const void *__restrict pVect1, const void *__restrict pVect2, const void *__restrict qty_ptr) 
{
  constexpr openMVG::matching::Hamming<U> hamming;
  const U *a = reinterpret_cast<const U *>(pVect1);
  const U *b = reinterpret_cast<const U *>(pVect2);
  return hamming(a, b,*(reinterpret_cast<const size_t*>(qty_ptr)));
}

template <typename U>
class HammingSpace : public hnswlib::SpaceInterface<unsigned int>
{
  hnswlib::DISTFUNC<unsigned int> fstdistfunc_;
  size_t data_size_;
  size_t dim_;

public:
  explicit HammingSpace(size_t dim)
  {
    fstdistfunc_ = HammingKernel<U>;
    dim_ = dim;
    data_size_ = dim_ * sizeof(U);
  }

  ~HammingSpace() {}

  size_t get_data_size() override
  {
    return data_size_;
  }

  hnswlib::DISTFUNC<unsigned int> get_dist_func() override
  {
    return fstdistfunc_;
  }

  void *get_dist_func_param() override
  {
    return &dim_;
  }
};

static int L1Kernel(const void *__restrict pVect1, const void *__restrict pVect2, const void *__restrict qty_ptr) 
{
  size_t size = *(reinterpret_cast<const size_t*>(qty_ptr));
  size = size >> 2;
  int result = 0;
  uint8_t *a = (uint8_t *)(pVect1); // discard const
  uint8_t *b = (uint8_t *)(pVect2); // discard const

  // Process 4 items for each loop for efficiency.
  for (size_t i = 0; i < size; i++) {
    result += std::abs(*a - *b);
    a++;b++;
    result += std::abs(*a - *b);
    a++;b++;
    result += std::abs(*a - *b);
    a++;b++;
    result += std::abs(*a - *b);
    a++;b++;
  }
  return result;
}
#ifdef __SSE2__
static int L1Kernel_SSE2(const void *__restrict pVect1, const void *__restrict pVect2, const void *__restrict qty_ptr) 
{
  const uint8_t *a = reinterpret_cast<const uint8_t *>(pVect1);
  const uint8_t *b = reinterpret_cast<const uint8_t *>(pVect2);
  return L1_SSE2(a, b, 128);
}
#endif

#ifdef __AVX2__
static int L1Kernel_AVX2(const void *__restrict pVect1, const void *__restrict pVect2, const void *__restrict qty_ptr) 
{
  const uint8_t *a = reinterpret_cast<const uint8_t *>(pVect1);
  const uint8_t *b = reinterpret_cast<const uint8_t *>(pVect2);
  return L1_AVX2(a, b, 128);
}
#endif

class L1SpaceInteger : public hnswlib::SpaceInterface<int>
{
  hnswlib::DISTFUNC<int> fstdistfunc_;
  size_t data_size_;
  size_t dim_;

public:
  explicit L1SpaceInteger(size_t dim)
  {
    fstdistfunc_ = L1Kernel;
    #ifdef __SSE2__
    //FIXME (RJ): Kernel disabled since there are some troubles on my Linux computer
    /*if(dim == 128) {
      fstdistfunc_ = L1Kernel_SSE2;
    }*/
    #endif
    #ifdef __AVX2__
    if(dim == 128) {
      fstdistfunc_ = L1Kernel_AVX2;
    }
    #endif
    dim_ = dim;
    data_size_ = dim_ * sizeof(uint8_t);
  }

  ~L1SpaceInteger() {}

  size_t get_data_size() override
  {
    return data_size_;
  }

  hnswlib::DISTFUNC<int> get_dist_func() override
  {
    return fstdistfunc_;
  }

  void *get_dist_func_param() override
  {
    return &dim_;
  }
};

#ifdef __AVX2__
static int L2Kernel_AVX2(const void *__restrict pVect1, const void *__restrict pVect2, const void *__restrict qty_ptr) 
{
  const uint8_t *a = reinterpret_cast<const uint8_t *>(pVect1);
  const uint8_t *b = reinterpret_cast<const uint8_t *>(pVect2);
  return L2_AVX2(a, b, 128);
}
#endif

class L2SpaceInteger : public hnswlib::SpaceInterface<int>
{
  hnswlib::DISTFUNC<int> fstdistfunc_;
  size_t data_size_;
  size_t dim_;

public:
  explicit L2SpaceInteger(size_t dim)
  {
    fstdistfunc_ = hnswlib::L2SqrI;
    #ifdef __AVX2__
    if(dim == 128) {
      fstdistfunc_ = L2Kernel_AVX2;
    }
    #endif
    dim_ = dim;
    data_size_ = dim_ * sizeof(uint8_t);
  }

  ~L2SpaceInteger() {}

  size_t get_data_size() override
  {
    return data_size_;
  }

  hnswlib::DISTFUNC<int> get_dist_func() override
  {
    return fstdistfunc_;
  }

  void *get_dist_func_param() override
  {
    return &dim_;
  }
};

} // namespace custom_hnsw
} // namespace matching
} // namespace openMVG
#endif // OPENMVG_MATCHING_METRIC_HNSW_HPP