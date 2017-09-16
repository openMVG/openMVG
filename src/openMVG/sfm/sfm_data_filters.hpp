// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_FILTERS_HPP
#define OPENMVG_SFM_SFM_DATA_FILTERS_HPP

#include <set>

#include "openMVG/types.hpp"

namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {
namespace sfm {

/// List the view indexes that have valid camera intrinsic and pose.
std::set<IndexT> Get_Valid_Views
(
  const SfM_Data & sfm_data
);

/// Filter a list of pair: Keep only the pairs that are defined by the index list
template <typename IterablePairs, typename IterableIndex>
Pair_Set Pair_filter
(
  const IterablePairs & pairs,
  const IterableIndex & index
)
{
  Pair_Set kept_pairs;
  for (auto& it : pairs)
  {
    if (index.count(it.first) > 0 &&
      index.count(it.second) > 0)
    kept_pairs.insert(it);
  }
  return kept_pairs;
}

// Remove tracks that have a small angle (tracks with tiny angle leads to instable 3D points)
// Return the number of removed tracks
IndexT RemoveOutliers_PixelResidualError
(
  SfM_Data & sfm_data,
  const double dThresholdPixel,
  const unsigned int minTrackLength = 2
);

// Remove tracks that have a small angle (tracks with tiny angle leads to instable 3D points)
// Return the number of removed tracks
IndexT RemoveOutliers_AngleError
(
  SfM_Data & sfm_data,
  const double dMinAcceptedAngle
);

/// Erase pose with insufficient track observations
bool eraseMissingPoses
(
  SfM_Data & sfm_data,
  const IndexT min_points_per_pose = 6
);

/// Erase observations with no defined pose
bool eraseObservationsWithMissingPoses
(
  SfM_Data & sfm_data,
  const IndexT min_points_per_landmark = 2
);

/// Remove unstable content from analysis of the sfm_data structure
bool eraseUnstablePosesAndObservations
(
  SfM_Data & sfm_data,
  const IndexT min_points_per_pose = 6,
  const IndexT min_points_per_landmark = 2
);

/// Tell if the sfm_data structure is one CC or not
bool IsTracksOneCC
(
  const SfM_Data & sfm_data
);

/// Keep the largest connected component of tracks from the sfm_data structure
void KeepLargestViewCCTracks
(
  SfM_Data & sfm_data
);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_FILTERS_HPP
