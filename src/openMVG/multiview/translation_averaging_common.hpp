// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_H_
#define OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_H_

#include "openMVG/numeric/numeric.h"
#include <utility>

namespace openMVG {

/// Relative information [Rij[tij] for a pair
typedef std::pair<std::pair<std::size_t,std::size_t>, std::pair<Mat3,Vec3> > relativeInfo;

} // namespace openMVG

#endif //OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_H_
