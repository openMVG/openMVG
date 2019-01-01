// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_SFM_SFM_DATA_COLORIZATION_HPP
#define OPENMVG_SFM_SFM_DATA_COLORIZATION_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {
namespace sfm {

struct SfM_Data;

bool ColorizeTracks(
  const SfM_Data & sfm_data,
  std::vector<Vec3> & vec_3dPoints,
  std::vector<Vec3> & vec_tracksColor);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_COLORIZATION_HPP
