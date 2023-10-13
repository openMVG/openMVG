// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON, Ricardo FABBRI and Gabriel ANDRADE.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_UTIL_HPP_BASE
#define OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_UTIL_HPP_BASE

// -----------------------------------------------------------------------------
// Internal helper functions for the Sequential SFM pipeline
// These functions work on sfm_data alone
// They are put on a separate file for organization and also for breaking the
// pipeline into small translation units that are fast to compile in a
// development-build cycle
// -----------------------------------------------------------------------------

#include "openMVG/sfm/base/sfm_matches_provider.hpp"

namespace openMVG {
namespace sfm {

struct SfM_Data;

// TODO(trifocal future) GetTripletWithMostMatches
// Get the PairWiseMatches that have the most support point
std::vector<openMVG::matching::PairWiseMatches::const_iterator>
GetPairWithMostMatches(
    const SfM_Data& sfm_data, 
    const openMVG::matching::PairWiseMatches& matches, 
    int clamp_count = 10);

inline void
invert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double px_tgt_coords[2], 
    double normalized_tgt_coords[2])
{
  const double *tp = px_tgt_coords;
  double *t = normalized_tgt_coords;
  t[1] = tp[1]/K[1][1];
  t[0] = (tp[0] - K[0][1]*t[1])/K[0][0];
  // normalize -- works as a cache for angle computations / dot products
  // TODO: check inside minus if we are normalizing
  double n = hypot(t[0],t[1]);
  t[0] /= n; t[1] /= n;
}


} // namespace sfm
} // namespace openMVG


#endif // OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_UTIL_HPP
