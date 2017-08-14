// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_MATCHER_KDTREE_FLANN_HPP
#define OPENMVG_MATCHING_MATCHER_KDTREE_FLANN_HPP

#include <memory>
#include <vector>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

#include "openMVG/matching/matching_interface.hpp"

#include <flann/flann.hpp>

namespace openMVG {
namespace matching  {

/// Implement ArrayMatcher as a FLANN KDtree matcher.
// http://www.cs.ubc.ca/~mariusm/index.php/FLANN/FLANN
// David G. Lowe and Marius Muja
//
// By default use squared L2 metric (flann::L2<Scalar>)
// sqrt is monotonic so for performance reason we do not compute it.

template < typename Scalar = float, typename  Metric = flann::L2<Scalar> >
class ArrayMatcher_Kdtree_Flann : public ArrayMatcher<Scalar, Metric>
{
  public:
  using DistanceType = typename Metric::ResultType;

  ArrayMatcher_Kdtree_Flann() = default;

  virtual ~ArrayMatcher_Kdtree_Flann() = default;

  /**
   * Build the matching structure
   *
   * \param[in] dataset   Input data.
   * \param[in] nbRows    The number of component.
   * \param[in] dimension Length of the data contained in the each
   *  row of the dataset.
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
    if (nbRows > 0)
    {
      dimension_ = dimension;
      //-- Build Flann Matrix container (map to already allocated memory)
      datasetM_.reset(
          new flann::Matrix<Scalar>((Scalar*)dataset, nbRows, dimension));

      //-- Build FLANN index
      index_.reset(
          new flann::Index<Metric> (*datasetM_, flann::KDTreeIndexParams(4)));
      index_->buildIndex();

      return true;
    }
    return false;
  }

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
  bool SearchNeighbour
  (
    const Scalar * query,
    int * indice,
    DistanceType * distance
  ) override
  {
    if (index_.get() != nullptr)
    {
      int * indicePTR = indice;
      DistanceType * distancePTR = distance;
      flann::Matrix<Scalar> queries((Scalar*)query, 1, dimension_);

      flann::Matrix<int> indices(indicePTR, 1, 1);
      flann::Matrix<DistanceType> dists(distancePTR, 1, 1);
      // do a knn search, using 128 checks
      return (index_->knnSearch(queries, indices, dists, 1, flann::SearchParams(128)) > 0);
    }
    else
    {
      return false;
    }
  }


/**
   * Search the N nearest Neighbor of the scalar array query.
   *
   * \param[in]   query           The query array
   * \param[in]   nbQuery         The number of query rows
   * \param[out]  indices   The corresponding (query, neighbor) indices
   * \param[out]  pvec_distances  The distances between the matched arrays.
   * \param[in]  NN              The number of maximal neighbor that will be searched.
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
    if (index_.get() != nullptr && NN <= datasetM_->rows)
    {
      std::vector<DistanceType> vec_distances(nbQuery * NN);
      DistanceType * distancePTR = &(vec_distances[0]);
      flann::Matrix<DistanceType> dists(distancePTR, nbQuery, NN);

      std::vector<int> vec_indices(nbQuery * NN, -1);
      flann::Matrix<int> indices(&(vec_indices[0]), nbQuery, NN);

      flann::Matrix<Scalar> queries((Scalar*)query, nbQuery, dimension_);
      // do a knn search, using 128 checks
      flann::SearchParams params(128);
#ifdef OPENMVG_USE_OPENMP
      params.cores = omp_get_max_threads();
#endif
      if (index_->knnSearch(queries, indices, dists, NN, params)>0)
      {
        // Save the resulting found indices
        pvec_indices->reserve(nbQuery * NN);
        pvec_distances->reserve(nbQuery * NN);
        for (size_t i = 0; i < nbQuery; ++i)
        {
          for (size_t j = 0; j < NN; ++j)
          {
            pvec_indices->emplace_back(i, vec_indices[i*NN+j]);
            pvec_distances->emplace_back(vec_distances[i*NN+j]);
          }
        }
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  private:

  std::unique_ptr< flann::Matrix<Scalar> > datasetM_;
  std::unique_ptr< flann::Index<Metric> > index_;
  std::size_t dimension_;
};

} // namespace matching
} // namespace openMVG

#endif // OPENMVG_MATCHING_MATCHER_KDTREE_FLANN_HPP
