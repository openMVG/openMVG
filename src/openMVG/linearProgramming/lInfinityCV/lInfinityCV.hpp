// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_L_INFINITY_COMPUTER_VISION_H_
#define OPENMVG_L_INFINITY_COMPUTER_VISION_H_

// Structure and motion problem solver
#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_From_xi_Ri.hpp"
#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_From_xi_Ri_noise.hpp"
#include "openMVG/linearProgramming/lInfinityCV/triplet_tijsAndXis_kernel.hpp"

// Pose estimation solver
#include "openMVG/linearProgramming/lInfinityCV/resection.hpp"
#include "openMVG/linearProgramming/lInfinityCV/resection_kernel.hpp"
// N-View Triangulation solver
#include "openMVG/linearProgramming/lInfinityCV/triangulation.hpp"

//-------------
//-- Global SfM
//-------------
// Compute from global translation by using X-views relative translations guess
#include "openMVG/linearProgramming/lInfinityCV/global_translations_fromTij.hpp"

#endif // OPENMVG_L_INFINITY_COMPUTER_VISION_H_
