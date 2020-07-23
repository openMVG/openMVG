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
  const U *a = static_cast<const U *>(pVect1);
  const U *b = static_cast<const U *>(pVect2);
  return hamming(a, b,*(static_cast<const size_t*>(qty_ptr)));
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
  const uint8_t *a = static_cast<const uint8_t *>(pVect1); // we are sure this was an uint8
  const uint8_t *b = static_cast<const uint8_t *>(pVect2);
  return metricL1(a, b,*(static_cast<const size_t*>(qty_ptr)));
}

class L1SpaceInteger : public hnswlib::SpaceInterface<int>
{
  hnswlib::DISTFUNC<int> fstdistfunc_;
  size_t data_size_;
  size_t dim_;

public:
  explicit L1SpaceInteger(size_t dim): 
    fstdistfunc_(L1Kernel), dim_(dim), data_size_(dim * sizeof(uint8_t))
  {}

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

} // namespace custom_hnsw
} // namespace matching
} // namespace openMVG
#endif // OPENMVG_MATCHING_METRIC_HNSW_HPP