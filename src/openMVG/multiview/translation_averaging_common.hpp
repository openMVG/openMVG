// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_HPP
#define OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_HPP

#include <map>
#include <set>
#include <utility>
#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/types.hpp"

namespace openMVG {

/// Relative information [Rij|tij] for a pair
using relativeInfo = std::pair< Pair, std::pair<Mat3,Vec3>>;

using RelativeInfo_Vec = std::vector< relativeInfo >;
using RelativeInfo_Map = std::map< Pair, std::pair<Mat3, Vec3>>;

// List the pairs used by the relative motions
Pair_Set getPairs(const RelativeInfo_Vec & vec_relative);

Pair_Set getPairs(const std::vector<RelativeInfo_Vec> & vec_relative);

// List the index used by the relative motions
std::set<IndexT> getIndexT(const RelativeInfo_Vec & vec_relative);

// List the index used by the relative motions
std::set<IndexT> getIndexT(const std::vector<RelativeInfo_Vec> & vec_relative);

} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_HPP
