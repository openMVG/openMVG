// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_MATCHING_INTERFACE_HPP
#define OPENMVG_MATCHING_MATCHING_INTERFACE_HPP

#include <vector>

#include "openMVG/matching/indMatch.hpp"

namespace openMVG {
namespace matching {

template < typename Scalar, typename Metric >
class ArrayMatcher
{
  public:
  using ScalarT = Scalar;
  using DistanceType = typename Metric::ResultType;
  using MetricT = Metric;

  ArrayMatcher() = default;
  virtual ~ArrayMatcher() = default;

  /**
   * Build the matching structure
   *
   * \param[in] dataset   Input data.
   * \param[in] nbRows    The number of component.
   * \param[in] dimension Length of the data contained in the dataset.
   *
   * \return True if success.
   */
  virtual bool Build( const Scalar * dataset, int nbRows, int dimension)=0;

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
  virtual bool SearchNeighbour( const Scalar * query,
                                int * indice, DistanceType * distance)=0;


/**
   * Search the N nearest Neighbor of the scalar array query.
   *
   * \param[in]   query     The query array
   * \param[in]   nbQuery   The number of query rows
   * \param[out]  indices   The corresponding (query, neighbor) indices
   * \param[out]  distances The distances between the matched arrays.
   * \param[out]  NN        The number of maximal neighbor that could
   *  will be searched.
   *
   * \return True if success.
   */
  virtual bool SearchNeighbours( const Scalar * query, int nbQuery,
                                  IndMatches * indices,
                                  std::vector<DistanceType> * distances,
                                  size_t NN)=0;
};

}  // namespace matching
}  // namespace openMVG

#endif  // OPENMVG_MATCHING_MATCHING_INTERFACE_HPP
