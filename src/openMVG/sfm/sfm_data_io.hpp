
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_SFM_DATA_IO_HPP
#define OPENMVG_SFM_DATA_IO_HPP

#include "openMVG/sfm/sfm_data.hpp"

namespace openMVG {
namespace sfm {

enum ESfM_Data
{
  VIEWS           = 0x01,
  EXTRINSICS      = 0x02,
  INTRINSICS      = 0x04,
  STRUCTURE       = 0x08,
  CONTROL_POINTS  = 0x16,
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

#endif // OPENMVG_SFM_DATA_IO_HPP
