// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_H_
#define OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_H_

#include "openMVG/types.hpp"
#include "openMVG/numeric/numeric.h"
#include <utility>

namespace openMVG {

/// Relative information [Rij|tij] for a pair
typedef std::pair< Pair, std::pair<Mat3,Vec3> > relativeInfo;

typedef std::vector< relativeInfo > RelativeInfo_Vec;
typedef std::map< Pair, std::pair<Mat3,Vec3> > RelativeInfo_Map;

} // namespace openMVG

#endif //OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_H_
