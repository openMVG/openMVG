// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_SFM_SFM_DATA_IO_HPP
#define OPENMVG_SFM_SFM_DATA_IO_HPP

#include <string>

#include "openMVG/types.hpp"

namespace openMVG {
namespace sfm {

struct SfM_Data;

enum ESfM_Data
{
  // Note: Use power of two values in order to use bitwise operators.
  VIEWS           =  1,
  EXTRINSICS      =  2,
  INTRINSICS      =  4,
  STRUCTURE       =  8,
  CONTROL_POINTS  = 16,
  ALL = VIEWS | EXTRINSICS | INTRINSICS | STRUCTURE | CONTROL_POINTS
};

///Check that each pose have a valid intrinsic and pose id in the existing View ids
bool ValidIds(const SfM_Data & sfm_data, ESfM_Data flags_part);

/// Load SfM_Data SfM scene from a file
bool Load(SfM_Data & sfm_data, const std::string & filename, ESfM_Data flags_part);

/// Save SfM_Data SfM scene to a file
bool Save(const SfM_Data & sfm_data, const std::string & filename, ESfM_Data flags_part);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_IO_HPP
