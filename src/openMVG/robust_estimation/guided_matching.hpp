// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_GUIDED_MATCHING_H_
#define OPENMVG_ROBUST_ESTIMATION_GUIDED_MATCHING_H_

#include "openMVG/numeric/numeric.h"
#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::matching;

#include "openMVG/features/features.hpp"

#include <vector>

namespace openMVG{

/// Guided Matching (features only):
///  Use a model to find valid correspondences:
///   Keep the best corresponding points for the given model under the
///   user specified distance.
template<
  typename ModelArg, // The used model type
  typename ErrorArg> // The metric to compute distance to the model
void GuidedMatching(
  const ModelArg & mod, // The model
  const Mat & xLeft,    // The left data points
  const Mat & xRight,   // The right data points
  double errorTh,       // Maximal authorized error threshold
  std::vector<IndMatch> & vec_corresponding_index) // Ouput corresponding index
{
  assert(xLeft.rows() == xRight.rows());

  // Looking for the corresponding points that have
  //  the smallest distance (smaller than the provided Threshold)

  for (size_t i = 0; i < xLeft.cols(); ++i) {

    double min = std::numeric_limits<double>::max();
    IndMatch match;
    for (size_t j = 0; j < xRight.cols(); ++j) {
      // Compute the geometric error: error to the model
      const double err = ErrorArg::Error(
        mod,  // The model
        xLeft.col(i), xRight.col(j)); // The corresponding points
      // if smaller error update corresponding index
      if (err < errorTh && err < min) {
        min = err;
        match = IndMatch(i,j);
      }
    }
    if (min < errorTh)  {
      // save the best corresponding index
      vec_corresponding_index.push_back(match);
    }
  }

  // Remove duplicates (when multiple points at same position exist)
  IndMatch::getDeduplicated(vec_corresponding_index);
}

// Struct to help filtering of correspondence according update of 
//  two smallest distance.
// -> useful for descriptor distance ratio filtering
template <typename DistT>
struct distanceRatio
{
  DistT bd, sbd; // best and second best distance
  size_t idx; // best corresponding index

  distanceRatio():
    bd(std::numeric_limits<DistT>::max()),
    sbd(std::numeric_limits<DistT>::max()),
    idx(0)
  { }

  // Update match according the provided distance
  inline void update(size_t index, DistT dist)
  {
    if (dist < bd) // best than any previous
    {
      idx = index;
      // update and swap
      sbd = dist;
      std::swap(bd, sbd);
    }
    else if(dist < sbd)
    {
      sbd = dist;
    }
  }

  // Return if the ratio of distance is ok or not
  inline bool isValid(const double distRatio) const{
    // check:
    // - that two best distance have been found
    // - the distance ratio
    return
      (sbd != std::numeric_limits<DistT>::max()
      && bd < distRatio * sbd);
  }
};

/// Guided Matching (features + descriptors with distance ratio):
///  Use a model to find valid correspondences:
///   Keep the best corresponding points for the given model under the
///   user specified distance ratio.
template<
  typename ModelArg, // The used model type
  typename ErrorArg, // The metric to compute distance to the model
  typename DescriptorT, // The descriptor type
  typename MetricT > // The metric to compare two descriptors
void GuidedMatching(
  const ModelArg & mod, // The model
  const Mat & xLeft,    // The left data points
  const std::vector<DescriptorT > & lDescriptors,
  const Mat & xRight,   // The right data points
  const std::vector<DescriptorT > & rDescriptors,
  double errorTh,       // Maximal authorized error threshold
  double distRatio, // Maximal authorized distance ratio
  std::vector<IndMatch> & vec_corresponding_index) // Ouput corresponding index
{
  assert(xLeft.rows() == xRight.rows());
  assert(xLeft.cols() == lDescriptors.size());
  assert(xRight.cols() == rDescriptors.size());

  MetricT metric;
  typename MetricT::ResultType MetricDist_T;

  // Looking for the corresponding points that have
  //  the satisfy:
  //   1. an geometric distance below the provided Threshold
  //   2. a distance ratio between descriptors of valid geometric correspondencess

  for (size_t i = 0; i < xLeft.cols(); ++i) {

    distanceRatio<typename MetricT::ResultType > dR;
    for (size_t j = 0; j < xRight.cols(); ++j) {
      // Compute the geometric error: error to the model
      const double geomErr = ErrorArg::Error(
        mod,  // The model
        xLeft.col(i), xRight.col(j)); // The corresponding points
      const typename MetricT::ResultType descDist =
        metric( lDescriptors[i].getData(), rDescriptors[j].getData(), DescriptorT::static_size );
      // if smaller error update corresponding index
      if (geomErr < errorTh) {
        dR.update(j, descDist);
      }
    }
    // Add correspondence only if if the distance ratio is valid
    if (dR.isValid(distRatio))  {
      // save the best corresponding index
      vec_corresponding_index.push_back(IndMatch(i,dR.idx));
    }
  }

  // Remove duplicates (when multiple points at same position exist)
  IndMatch::getDeduplicated(vec_corresponding_index);
}

} // namespace openMVG
#endif // OPENMVG_ROBUST_ESTIMATION_GUIDED_MATCHING_H_
