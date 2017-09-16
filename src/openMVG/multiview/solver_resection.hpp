// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_RESECTION_HPP
#define OPENMVG_MULTIVIEW_RESECTION_HPP

namespace openMVG {
namespace resection {

enum class SolverType
{
  DLT_6POINTS = 0,
  P3P_KE_CVPR17 = 1,
  P3P_KNEIP_CVPR11 = 2,
};

}  // namespace resection
}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_RESECTION_HPP
