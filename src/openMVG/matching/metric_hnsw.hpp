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
unsigned int HammingHnsw(const void *__restrict pVect1, const void *__restrict pVect2, const void *__restrict qty_ptr) 
{
  size_t qty = *(reinterpret_cast<const size_t*>(qty_ptr));
  const uint8_t *a = reinterpret_cast<const uint8_t *>(pVect1);
  const uint8_t *b = reinterpret_cast<const uint8_t *>(pVect2);
  openMVG::matching::Hamming<uint8_t> hamming;
  return hamming(a, b, qty);
}

class HammingSpace : public hnswlib::SpaceInterface<unsigned int>
{
  hnswlib::DISTFUNC<unsigned int> fstdistfunc_;
  size_t data_size_;
  size_t dim_;

public:
  HammingSpace(size_t dim)
  {
    fstdistfunc_ = HammingHnsw;
    dim_ = dim;
    data_size_ = dim_ * sizeof(uint8_t);
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

#endif // OPENMVG_MATCHING_METRIC_HNSW_HPP