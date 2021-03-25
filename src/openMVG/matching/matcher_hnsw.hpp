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
#include <typeindex>
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
    HNSW_metric_.reset();
    HNSW_matcher_.reset();

    if (nbRows < 1)
    {
      return false;
    }

    dimension_ = dimension;

    const std::type_index distance_type(typeid(DistanceType));
    // Here this is tricky since there is no specialization
    switch (MetricType)
    {
    case  HNSWMETRIC::L1_HNSW:
      if (distance_type == typeid(int)) {
        HNSW_metric_.reset(dynamic_cast<SpaceInterface<DistanceType> *>(new custom_hnsw::L1SpaceInteger(dimension)));
      }
      break;
    case  HNSWMETRIC::L2_HNSW:
      if (distance_type == typeid(int)) {
        HNSW_metric_.reset(dynamic_cast<SpaceInterface<DistanceType> *>(new L2SpaceI(dimension)));
      } else
      if (distance_type == typeid(float)) {
        HNSW_metric_.reset(dynamic_cast<SpaceInterface<DistanceType> *>(new L2Space(dimension)));
      }
      break;
    case  HNSWMETRIC::HAMMING_HNSW:
      if (distance_type == typeid(unsigned int)) {
        HNSW_metric_.reset(dynamic_cast<SpaceInterface<DistanceType> *>(new custom_hnsw::HammingSpace<uint8_t>(dimension)));
      }
      break;
    }

    if (!HNSW_metric_) {
      OPENMVG_LOG_ERROR << "HNSW matcher: this type of distance is not handled yet";
      return false;
    }

    HNSW_matcher_.reset(new HierarchicalNSW<DistanceType>(HNSW_metric_.get(), nbRows, 16, 100));

    // add a first point...
    HNSW_matcher_->addPoint(static_cast<const void *>(dataset), static_cast<size_t>(0));
    //...and the others in parallel
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int vector_id = 1; vector_id < nbRows; ++vector_id) {
        HNSW_matcher_->addPoint(static_cast<const void *>(dataset + dimension * vector_id), static_cast<size_t>(vector_id));
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
    HNSW_matcher_->setEf(16); //here we stay conservative but it could probably be lowered in this case (first NN)
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
    // EfSearch parameter could not be < NN.
    // -
    // For vectors with dimensionality of approx. 64-128 and for 2 NNs,
    // EfSearch = 16 produces good results in conjunction with other parameters fixed in this file (EfConstruct = 16, M = 100).
    // But nothing has been evaluated on our side for lower / higher dimensionality and for a higher number of NNs.
    // So for now and for NN > 2, EfSearch is fixed to 2 * NNs without a good a priori knowledge.
    // A good value for EfSearch could really depends on the two other parameters (EfConstruct / M).
    if (NN <= 2) {
      HNSW_matcher_->setEf(16);
    } else {
      HNSW_matcher_->setEf(std::max(NN*2, static_cast<size_t>(nbQuery)));
    }
    pvec_indices->resize(nbQuery * NN);
    pvec_distances->resize(nbQuery * NN);
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int query_id = 0; query_id < nbQuery; query_id++) {
      auto result = HNSW_matcher_->searchKnn(static_cast<const void *>(query + dimension_ * query_id), NN);
      size_t result_id = NN - 1;
      while(!result.empty())
      {
        const auto & res = result.top();
        const size_t match_id = query_id * NN + result_id;
        pvec_indices->operator[](match_id) = {static_cast<IndexT>(query_id), static_cast<IndexT>(res.second)};
        pvec_distances->operator[](match_id) = res.first;
        result.pop();
        result_id--;
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
