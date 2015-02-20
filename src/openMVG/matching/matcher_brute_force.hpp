// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_ARRAYMATCHER_BRUTE_FORCE_H
#define OPENMVG_MATCHING_ARRAYMATCHER_BRUTE_FORCE_H

#include "openMVG/numeric/numeric.h"
#include "openMVG/matching/matching_interface.hpp"
#include "openMVG/matching/metric.hpp"
#include "openMVG/stl/indexed_sort.hpp"
#include <memory>
#include <iostream>

namespace openMVG {
namespace matching {

// By default compute square(L2 distance).
template < typename Scalar = float, typename Metric = L2_Simple<Scalar> >
class ArrayMatcherBruteForce  : public ArrayMatcher<Scalar, Metric>
{
  public:
  typedef typename Metric::ResultType DistanceType;

  ArrayMatcherBruteForce()   {}
  virtual ~ArrayMatcherBruteForce() {
    memMapping.reset();
  }

  /**
   * Build the matching structure
   *
   * \param[in] dataset   Input data.
   * \param[in] nbRows    The number of component.
   * \param[in] dimension Length of the data contained in the dataset.
   *
   * \return True if success.
   */
  bool Build(const Scalar * dataset, int nbRows, int dimension) {
    if (nbRows < 1) {
      memMapping = auto_ptr< Eigen::Map<BaseMat> >(NULL);
      return false;
    }
    memMapping = auto_ptr< Eigen::Map<BaseMat> >
      (new Eigen::Map<BaseMat>( (Scalar*)dataset, nbRows, dimension) );
    return true;
  };

  /**
   * Search the nearest Neighbor of the scalar array query.
   *
   * \param[in]   query     The query array
   * \param[out]  indice    The indice of array in the dataset that
   *  have been computed as the nearest array.
   * \param[out]  distance  The distance between the two arrays.
   *
   * \return True if success.
   */
  bool SearchNeighbour( const Scalar * query,
                        int * indice, DistanceType * distance)
  {
    if (memMapping.get() != NULL)  {
      //matrix representation of the input data;
      Eigen::Map<BaseMat> mat_query((Scalar*)query, 1, (*memMapping).cols() );
      Metric metric;
      vector<DistanceType> vec_dist((*memMapping).rows(), 0.0);
      for (int i = 0; i < (*memMapping).rows(); ++i)  {
        // Compute Distance Metric
        vec_dist[i] = metric( (Scalar*)query, (*memMapping).row(i).data(), (*memMapping).cols() );
      }
      if (!vec_dist.empty())
      {
        // Find the minimum distance :
        typename vector<DistanceType>::const_iterator min_iter =
          min_element( vec_dist.begin(), vec_dist.end());
        *indice =std::distance(
          typename vector<DistanceType>::const_iterator(vec_dist.begin()),
          min_iter);
        *distance = static_cast<DistanceType>(*min_iter);
      }
      return true;
    }
    else  {
      return false;
    }
  }


/**
   * Search the N nearest Neighbor of the scalar array query.
   *
   * \param[in]   query     The query array
   * \param[in]   nbQuery   The number of query rows
   * \param[out]  indice    The indices of arrays in the dataset that
   *  have been computed as the nearest arrays.
   * \param[out]  distance  The distances between the matched arrays.
   * \param[out]  NN        The number of maximal neighbor that will be searched.
   *
   * \return True if success.
   */
  bool SearchNeighbours( const Scalar * query, int nbQuery,
                          vector<int> * pvec_indice,
                          vector<DistanceType> * pvec_distance,
                          size_t NN)
  {
    if (memMapping.get() == NULL)  {
      return false;
    }
    
    if (NN > (*memMapping).rows() || nbQuery < 1) {
      std::cerr << "Too much asked nearest neighbors" << std::endl;
      return false;
    }

    //matrix representation of the input data;
    Eigen::Map<BaseMat> mat_query((Scalar*)query, nbQuery, (*memMapping).cols());
    Metric metric;

    std::vector<DistanceType> vec_distance((*memMapping).rows(), 0.0);
    for (int queryIndex=0; queryIndex < nbQuery; ++queryIndex) {

      for (int i = 0; i < (*memMapping).rows(); ++i)  {
        vec_distance[i] = metric( mat_query.row(queryIndex).data(),
          (*memMapping).row(i).data(), (*memMapping).cols() );
      }

      // Find the N minimum distances:
      const int maxMinFound = (int) min( size_t(NN), vec_distance.size());
      using namespace indexed_sort;
      vector< sort_index_packet_ascend< DistanceType, int> > packet_vec(vec_distance.size());
      sort_index_helper(packet_vec, &vec_distance[0], maxMinFound);

      for (int i = 0; i < maxMinFound; ++i) {
        pvec_distance->push_back(packet_vec[i].val);
        pvec_indice->push_back(packet_vec[i].index);
      }
    }
    return true;
  };

private:
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> BaseMat;
  /// Use a memory mapping in order to avoid memory re-allocation
  auto_ptr< Eigen::Map<BaseMat> > memMapping;
};

}  // namespace matching
}  // namespace openMVG

#endif  // OPENMVG_MATCHING_ARRAYMATCHER_BRUTE_FORCE_H
