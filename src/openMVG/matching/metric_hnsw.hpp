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
template <typename U>
static unsigned int HammingHnsw(const void *__restrict pVect1, const void *__restrict pVect2, const void *__restrict qty_ptr) 
{
  const U *a = reinterpret_cast<const U *>(pVect1);
  const U *b = reinterpret_cast<const U *>(pVect2);
  const openMVG::matching::Hamming<U> hamming;
  return hamming(a, b,*(reinterpret_cast<const size_t*>(qty_ptr)));
}

template <typename U>
class HammingSpace : public hnswlib::SpaceInterface<unsigned int>
{
  hnswlib::DISTFUNC<unsigned int> fstdistfunc_;
  size_t data_size_;
  size_t dim_;

public:
  HammingSpace(size_t dim)
  {
    fstdistfunc_ = HammingHnsw<U>;
    dim_ = dim;
    data_size_ = dim_ * sizeof(U);
  }

  ~HammingSpace() {}

  size_t get_data_size()
  {
    return data_size_;
  }

  hnswlib::DISTFUNC<unsigned int> get_dist_func()
  {
    return fstdistfunc_;
  }

  void *get_dist_func_param()
  {
    return &dim_;
  }
};

static int L1Hnsw(const void *__restrict pVect1, const void *__restrict pVect2, const void *__restrict qty_ptr) 
{
  size_t size = *(reinterpret_cast<const size_t*>(qty_ptr));
  size = size >> 2;
  uint8_t *a = (uint8_t *)(pVect1); // discard const
  uint8_t *b = (uint8_t *)(pVect2); // discard const

    int result = 0;

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
  return result * result;
}

class L1Space : public hnswlib::SpaceInterface<int>
{
  hnswlib::DISTFUNC<int> fstdistfunc_;
  size_t data_size_;
  size_t dim_;

public:
  L1Space(size_t dim)
  {
    fstdistfunc_ = L1Hnsw;
    dim_ = dim;
    data_size_ = dim_ * sizeof(uint8_t);
  }

  ~L1Space() {}

  size_t get_data_size()
  {
    return data_size_;
  }

  hnswlib::DISTFUNC<int> get_dist_func()
  {
    return fstdistfunc_;
  }

  void *get_dist_func_param()
  {
    return &dim_;
  }
};
} // namespace matching
} // namespace openMVG
#endif // OPENMVG_MATCHING_METRIC_HNSW_HPP