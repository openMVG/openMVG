// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_IO_BAF_HPP
#define OPENMVG_SFM_SFM_DATA_IO_BAF_HPP

#include <string>

#include "openMVG/sfm/sfm_data_io.hpp"

namespace openMVG {
namespace sfm {

/// Save SfM_Data in an ASCII BAF (Bundle Adjustment File).
// --Header
// #Intrinsics
// #Poses
// #Landmarks
// --Data
// Intrinsic parameters [foc ppx ppy, ...]
// Poses [angle axis, camera center]
// Landmarks [X Y Z #observations id_intrinsic id_pose x y ...]
//--
//- Export also a _imgList.txt file with View filename and id_intrinsic & id_pose.
// filename id_intrinsic id_pose
// The ids allow to establish a link between 3D point observations & the corresponding views
//--
// Export missing poses as Identity pose to keep tracking of the original id_pose indexes
bool Save_BAF(
  const SfM_Data & sfm_data,
  const std::string & filename,
  ESfM_Data flags_part);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_IO_BAF_HPP
