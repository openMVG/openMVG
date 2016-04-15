
// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/geometry/Similarity3.hpp"

namespace openMVG {
namespace sfm {

/// Apply a similarity to the SfM_Data scene (transform landmarks & camera poses)
void ApplySimilarity
(
  const geometry::Similarity3 & sim,
  SfM_Data & sfm_data
)
{
  // Transform the landmark positions
  for (auto & iterLandMark : sfm_data.structure)
  {
    iterLandMark.second.X = sim(iterLandMark.second.X);
  }

  // Transform the camera positions
  for (auto & iterPose : sfm_data.poses)
  {
    iterPose.second = sim(iterPose.second);
  }
}

} // namespace sfm
} // namespace openMVG

