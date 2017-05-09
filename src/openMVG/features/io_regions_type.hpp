// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_IO_REGIONS_TYPE_HPP
#define OPENMVG_FEATURES_IO_REGIONS_TYPE_HPP

#include "openMVG/features/regions.hpp"

namespace openMVG {
namespace features {

// Init the regions_type from an image describer file (used for regions loading)
std::unique_ptr<features::Regions> Init_region_type_from_file
(
  const std::string & sImage_describer_file
);

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_IO_REGIONS_TYPE_HPP
