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
#include "openMVG/matching/metric_hnsw.hpp"

#include "third_party/hnswlib/hnswlib.h"

using namespace hnswlib;

namespace openMVG {
namespace matching {

enum HNSWMETRIC {
  L2_HNSW,
  L1_HNSW,
  HAMMING_HNSW
};

// By default compute square(L2 distance).
template <typename Scalar = float, typename Metric = L2<Scalar>, HNSWMETRIC MetricType = HNSWMETRIC::L2_HNSW>
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
      HNSW_metric_.reset(nullptr);
      HNSW_matcher_.reset(nullptr);  
      return false;
    }

    dimension_ = dimension;

    // Here this is tricky since there is no specialization
    switch (MetricType)
    {
    case  HNSWMETRIC::L1_HNSW:
      if (typeid(DistanceType) == typeid(int)) {
        HNSW_metric_.reset(dynamic_cast<SpaceInterface<DistanceType> *>(new custom_hnsw::L1SpaceInteger(dimension)));
      } else {
        std::cerr << "HNSW matcher: this type of distance is not handled yet for L1 metric" << std::endl;
        return false;
      }
      break;
    case  HNSWMETRIC::L2_HNSW:
      if (typeid(DistanceType) == typeid(int)) {
        HNSW_metric_.reset(dynamic_cast<SpaceInterface<DistanceType> *>(new L2SpaceI(dimension)));
      } else
      if (typeid(DistanceType) == typeid(float)) {
        HNSW_metric_.reset(dynamic_cast<SpaceInterface<DistanceType> *>(new L2Space(dimension)));
      } else {
        std::cerr << "HNSW matcher: this type of distance is not handled yet for L2 metric" << std::endl;
        return false;
      }
      break;
    case  HNSWMETRIC::HAMMING_HNSW:
      if (typeid(DistanceType) == typeid(unsigned int)) {
        HNSW_metric_.reset(dynamic_cast<SpaceInterface<DistanceType> *>(new custom_hnsw::HammingSpace<uint8_t>(dimension)));
      } else {
        std::cerr << "HNSW matcher: this type of distance is not handled yet for Hamming distance" << std::endl;
        return false;
      }
      break;
    default:
        std::cerr << "HNSW matcher: this type of distance is not handled yet" << std::endl;
        return false;  
      break;
    }

    HNSW_matcher_.reset(new HierarchicalNSW<DistanceType>(HNSW_metric_.get(), nbRows, 16, 100) );
    HNSW_matcher_->setEf(16);
    
    // add a first point...
    HNSW_matcher_->addPoint((const void *)(dataset), (size_t)0);
    //...and the others in parallel
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int vector_id = 1; vector_id < nbRows; ++vector_id) {
        HNSW_matcher_->addPoint((const void *)(dataset + dimension * vector_id), (size_t)vector_id);
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
    if (!HNSW_matcher_)
      return false;
    const auto result = HNSW_matcher_->searchKnn(query, 1).top();
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
    if (!HNSW_matcher_)
      return false;
    pvec_indices->reserve(nbQuery * NN);
    pvec_distances->reserve(nbQuery * NN);
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int query_id = 0; query_id < nbQuery; query_id++) {
      const auto result = HNSW_matcher_->searchKnn((const void *)(query + dimension_ * query_id), NN,
        [](const std::pair<DistanceType, size_t> &a, const std::pair<DistanceType, size_t> &b) -> bool {
          return a.first < b.first;
      });
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
      #endif
      {
      for (const auto & res : result)
      {
        pvec_indices->emplace_back(query_id, res.second);
        pvec_distances->emplace_back(res.first);
      }
      }
    }   
    return true;
  };

private:
  int dimension_;
  std::unique_ptr<SpaceInterface<DistanceType>> HNSW_metric_;
  std::unique_ptr<HierarchicalNSW<DistanceType>> HNSW_matcher_;
};

}  // namespace matching
}  // namespace openMVG

#endif  // OPENMVG_MATCHING_MATCHER_HNSW_HPP