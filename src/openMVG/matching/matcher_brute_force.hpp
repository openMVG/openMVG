// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_MATCHER_BRUTE_FORCE_HPP
#define OPENMVG_MATCHING_MATCHER_BRUTE_FORCE_HPP

#include <algorithm>
#include <memory>
#include <future>
#include <thread>
#include <vector>

#include "openMVG/numeric/numeric.h"
#include "openMVG/matching/matching_interface.hpp"
#include "openMVG/matching/metric.hpp"
#include "openMVG/stl/indexed_sort.hpp"

namespace openMVG {
namespace matching {

// By default compute square(L2 distance).
template < typename Scalar = float, typename Metric = L2<Scalar> >
class ArrayMatcherBruteForce : public ArrayMatcher<Scalar, Metric>
{
  public:
  using DistanceType = typename Metric::ResultType;

  ArrayMatcherBruteForce() = default;
  virtual ~ArrayMatcherBruteForce()= default;

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
      memMapping.reset(nullptr);
      return false;
    }
    memMapping.reset(new Eigen::Map<BaseMat>( (Scalar*)dataset, nbRows, dimension));
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
    if (!memMapping || memMapping->rows() < 1)
      return false;

    IndMatches vec_index(1);
    std::vector<DistanceType> dist(1);
    SearchNeighbours_func(query, 0, 1, &vec_index, &dist, 1);
    indice[0] = vec_index[0].j_;
    distance[0] = dist[0];
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
    if (!memMapping ||
        NN > memMapping->rows() ||
        nbQuery < 1)
    {
      return false;
    }

    pvec_distances->resize(nbQuery * NN);
    pvec_indices->resize(nbQuery * NN);

    const int nb_thread = static_cast<int>(std::thread::hardware_concurrency());
    // Compute ranges
    std::vector<int> range;
    SplitRange((int)0 , (int)nbQuery , nb_thread , range);

    std::vector<std::future<void>> fut;
    for (size_t i = 1; i < range.size(); ++i)
    {
      fut.push_back(
        std::async(
          std::launch::async,
          &ArrayMatcherBruteForce<Scalar, Metric>::SearchNeighbours_func,
          this,
          query,
          range[i-1],
          range[i],
          pvec_indices,
          pvec_distances,
          NN));
    }

    for (const auto & fut_it : fut)
    {
      fut_it.wait();
    }
    return true;
  };

private:
  using BaseMat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  /// Use a memory mapping in order to avoid memory re-allocation
  std::unique_ptr< Eigen::Map<BaseMat> > memMapping;

  /**
     * Search the N nearest Neighbor for a section of index of the scalar array query.
     *
     * \param[in]   query     The query array [query_start_index, query_stop_index[.
     * \param[in]   query_start_index  Start of range of index to handle.
     * \param[in]   query_stop_index  End of range to index to handle.
     * \param[out]  indices   The corresponding (query, neighbor) indices (updated for the range).
     * \param[out]  distances The distances between the matched arrays (update for the range).
     * \param[in]  NN        The number of maximal neighbor that will be searched.
     *
     * \return True if success.
     */
  void SearchNeighbours_func
  (
    const Scalar * query,
    size_t query_start_index,
    size_t query_stop_index,
    IndMatches * pvec_indices,
    std::vector<DistanceType> * pvec_distances,
    size_t NN
  )
  {
    // Compute the corresponding nearest neighbor(s) for the
    //  [query_start_index,query_stop_index[ range.
    Metric metric;
    std::vector<DistanceType> vec_distance(memMapping->rows(), 0.0);
    for (size_t queryIndex = query_start_index; queryIndex < query_stop_index; ++queryIndex)
    {
      std::fill(vec_distance.begin(), vec_distance.end(), DistanceType(0));
      const Scalar * queryPtr = query + queryIndex * memMapping->cols();
      for (typename BaseMat::Index i = 0; i < memMapping->rows(); ++i)
      {
        vec_distance[i] = metric(
          queryPtr,
          (*memMapping).data() + i * memMapping->cols(),
          memMapping->cols());
      }

      // Find the N minimum distances
      const int maxMinFound = static_cast<int>(std::min(size_t(NN), vec_distance.size()));
      std::vector<stl::indexed_sort::sort_index_packet_ascend<DistanceType, int>> packet_vec(vec_distance.size());
      stl::indexed_sort::sort_index_helper(packet_vec, &vec_distance[0], maxMinFound);

      for (int i = 0; i < maxMinFound; ++i)
      {
        (*pvec_distances)[queryIndex * NN + i] = packet_vec[i].val;
        (*pvec_indices)[queryIndex * NN + i] = IndMatch(queryIndex, packet_vec[i].index);
      }
    }
  }
};

}  // namespace matching
}  // namespace openMVG

#endif  // OPENMVG_MATCHING_MATCHER_BRUTE_FORCE_HPP