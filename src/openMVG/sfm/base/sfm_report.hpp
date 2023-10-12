// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_REPORT_HPP
#define OPENMVG_SFM_SFM_REPORT_HPP

#include <string>

namespace openMVG {
namespace sfm {

struct SfM_Data;

bool Generate_SfM_Report
(
  const SfM_Data & sfm_data,
  const std::string & htmlFilename
);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_REPORT_HPP
