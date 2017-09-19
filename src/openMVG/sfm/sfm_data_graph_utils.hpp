// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015, 2017 WhuAegeanSea, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_GRAPH_UTILS_HPP
#define OPENMVG_SFM_DATA_GRAPH_UTILS_HPP
#include <string>

namespace openMVG {
namespace sfm {

struct SfM_Data;

bool SplitMatchFileIntoMatchFiles(const SfM_Data & sfm_data, const std::string & match_file,
  const std::string & match_file_components, bool bBiEdge, int nMinNode = 3);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_GRAPH_UTILS_HPP
