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
static unsigned int HammingKernel(const void * pVect1, const void * pVect2, const void * qty_ptr) 
{
  constexpr openMVG::matching::Hamming<U> hamming{};
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
  explicit HammingSpace(size_t dim):
    fstdistfunc_(HammingKernel<U>), dim_(dim), data_size_(dim * sizeof(U)) {}

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

static int L1Kernel(const void * pVect1, const void * pVect2, const void * qty_ptr) 
{
  constexpr L1<uint8_t> metricL1{};
  const uint8_t *a = reinterpret_cast<const uint8_t *>(pVect1);
  const uint8_t *b = reinterpret_cast<const uint8_t *>(pVect2);
  return metricL1(a,b,*(reinterpret_cast<const size_t*>(qty_ptr)));
}
#ifdef __SSE2__
static int L1Kernel_SSE2_128(const void * pVect1, const void * pVect2, const void * qty_ptr) 
{
  const uint8_t *a = reinterpret_cast<const uint8_t *>(pVect1);
  const uint8_t *b = reinterpret_cast<const uint8_t *>(pVect2);
  return L1_SSE2(a, b, 128);
}
#endif

class L1SpaceInteger : public hnswlib::SpaceInterface<int>
{
  hnswlib::DISTFUNC<int> fstdistfunc_;
  size_t data_size_;
  size_t dim_;

public:
  explicit L1SpaceInteger(size_t dim): 
    fstdistfunc_(L1Kernel), dim_(dim), data_size_(dim_ * sizeof(uint8_t))
  {
    #ifdef __SSE2__
    if(dim == 128) {
      fstdistfunc_ = L1Kernel_SSE2_128;
    }
    #endif
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
static int L2Kernel_AVX2_128(const void * pVect1, const void * pVect2, const void * qty_ptr) 
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
  explicit L2SpaceInteger(size_t dim):
    fstdistfunc_(hnswlib::L2SqrI), dim_(dim), data_size_(dim_ * sizeof(uint8_t))
  {
    #ifdef __AVX2__
    if(dim == 128) {
      fstdistfunc_ = L2Kernel_AVX2_128;
    }
    #endif
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