// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_FEATURE_CONTAINER_HPP
#define OPENMVG_FEATURES_FEATURE_CONTAINER_HPP

#include <vector>

namespace openMVG {
namespace features {

class PointFeature;
using PointFeatures = std::vector<PointFeature>;

class SIOPointFeature;
using SIOPointFeatures = std::vector<SIOPointFeature>;

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_FEATURE_CONTAINER_HPP
