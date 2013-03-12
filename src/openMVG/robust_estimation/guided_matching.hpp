// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_GUIDED_MATCHING_H_
#define OPENMVG_ROBUST_ESTIMATION_GUIDED_MATCHING_H_

#include "openMVG/numeric/numeric.h"
#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::matching;

#include <vector>

namespace openMVG{

/// Guided Matching:
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
      // Compute error to the model
      double err = ErrorArg::Error(
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

} // namespace openMVG
#endif // OPENMVG_ROBUST_ESTIMATION_GUIDED_MATCHING_H_
