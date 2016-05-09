
// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_TRANSFORM_HPP
#define OPENMVG_SFM_DATA_TRANSFORM_HPP

namespace openMVG {

namespace geometry
{
  struct Similarity3;
}

namespace sfm {

struct SfM_Data;

/// Apply a similarity to the SfM_Data scene (transform landmarks & camera poses)
void ApplySimilarity
(
  const geometry::Similarity3 & sim,
  SfM_Data & sfm_data
);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_TRANSFORM_HPP
