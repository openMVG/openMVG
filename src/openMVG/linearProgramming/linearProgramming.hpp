// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINEAR_PROGRAMMING_H_
#define OPENMVG_LINEAR_PROGRAMMING_H_

#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#ifdef OPENMVG_HAVE_MOSEK
#include "openMVG/linearProgramming/linearProgrammingMOSEK.hpp"
#endif

#include "openMVG/linearProgramming/bisectionLP.hpp"

// Multiple View Geometry solver that rely on Linear programming formulations
#include "openMVG/linearProgramming/lInfinityCV/lInfinityCV.hpp"

#endif // OPENMVG_LINEAR_PROGRAMMING_H_
