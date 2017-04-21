// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_UTILS_HPP
#define OPENMVG_MATCHING_IND_MATCH_UTILS_HPP

#include <string>

#include "openMVG/matching/indMatch.hpp"

namespace openMVG {
namespace matching {

bool Load
(
  PairWiseMatches & matches,
  const std::string & filename
);

bool Save
(
  const PairWiseMatches & matches,
  const std::string & filename
);

}  // namespace matching
}  // namespace openMVG

#endif // #define OPENMVG_MATCHING_IND_MATCH_UTILS_HPP
