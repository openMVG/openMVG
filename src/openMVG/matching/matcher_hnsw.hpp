// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Romain Janvier and Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_MATCHER_HNSW_HPP
#define OPENMVG_MATCHING_MATCHER_HNSW_HPP

#include <memory>
#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif
#include <vector>

#include "openMVG/matching/matching_interface.hpp"
#include "openMVG/matching/metric.hpp"

#include "third_party/hnswlib/hnswlib.h"

using namespace hnswlib;

namespace openMVG {
namespace matching {

// By default compute square(L2 distance).
template <typename Scalar = float, typename Metric = L2<Scalar>>
class HNSWMatcher: public ArrayMatcher<Scalar, Metric>
{
  public:
  using DistanceType = typename Metric::ResultType;

  HNSWMatcher() = default;
  virtual ~HNSWMatcher()= default;

  /**
   * Build the matching structure
   *
   * \param[in] dataset   Input data.
   * \param[in] nbRows    The number of component.
   * \param[in] dimension Length of the data contained in the dataset.
   *
   * \return True if success.
   */
  bool Build
  (
    const Scalar * dataset,
    int nbRows,
    int dimension
  ) override
  {
    if (nbRows < 1)
    {
      HNSWmetric.reset(nullptr);
      HNSWmatcher.reset(nullptr);  
      return false;
    }

    dimension_ = dimension;
    
    // Here this is tricky since there is no specialization
    if(typeid(DistanceType)== typeid(int)) {
      HNSWmetric.reset(dynamic_cast<SpaceInterface<DistanceType> *>(new L2SpaceI(dimension)));
    } else
    if (typeid(DistanceType) == typeid(float))  {
      HNSWmetric.reset(dynamic_cast<SpaceInterface<DistanceType> *>(new L2Space(dimension)));
    } else {
      std::cerr << "HNSW matcher: this type of distance is not handled Yet" << std::endl;
    }
    
    HNSWmatcher.reset(new HierarchicalNSW<DistanceType>(HNSWmetric.get(), nbRows, 16, 100) );
    HNSWmatcher->setEf(16);
    
    // add first point..
    HNSWmatcher->addPoint((void *)(dataset), (size_t) 0);
    //...and the other in //
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i = 1; i < nbRows; i++) {
        HNSWmatcher->addPoint((void *) (dataset + dimension * i), (size_t) i);
    }

    return true;
  };

  /**
   * Search the nearest Neighbor of the scalar array query.
   *
   * \param[in]   query     The query array.
   * \param[out]  indice    The indice of array in the dataset that.
   *  have been computed as the nearest array.
   * \param[out]  distance  The distance between the two arrays.
   *
   * \return True if success.
   */
  bool SearchNeighbour
  (
    const Scalar * query,
    int * indice,
    DistanceType * distance
  ) override
  {
    if (HNSWmatcher.get() == nullptr)
      return false;
    auto result = HNSWmatcher->searchKnn(query, 1).top();
    *indice = result.second;
    *distance =  result.first;
    return true;
  }

  /**
   * Search the N nearest Neighbor of the scalar array query.
   *
   * \param[in]   query     The query array.
   * \param[in]   nbQuery   The number of query rows.
   * \param[out]  indices   The corresponding (query, neighbor) indices.
   * \param[out]  distances The distances between the matched arrays.
   * \param[in]  NN        The number of maximal neighbor that will be searched.
   *
   * \return True if success.
   */
  bool SearchNeighbours
  (
    const Scalar * query, int nbQuery,
    IndMatches * pvec_indices,
    std::vector<DistanceType> * pvec_distances,
    size_t NN
  ) override
  {
    if (HNSWmatcher.get() == nullptr)
    {
      return false;
    }
    pvec_indices->reserve(nbQuery * NN);
    pvec_distances->reserve(nbQuery * NN);
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i = 0; i < nbQuery; i++) {
      auto result = HNSWmatcher->searchKnn((const void *) (query + dimension_ * i), NN,
        [](const std::pair<DistanceType, size_t> &a, const std::pair<DistanceType, size_t> &b) -> bool {
          return a.first < b.first;
      });
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
      #endif
      {
      for (const auto & res : result)
      {
        pvec_indices->emplace_back(i, res.second);
        pvec_distances->emplace_back(res.first);
      }
      }
    }   
    return true;
  };

private:
  int dimension_;
  std::unique_ptr<SpaceInterface<DistanceType>> HNSWmetric;
  std::unique_ptr<HierarchicalNSW<DistanceType>> HNSWmatcher;
};

}  // namespace matching
}  // namespace openMVG

#endif  // OPENMVG_MATCHING_MATCHER_HNSW_HPP