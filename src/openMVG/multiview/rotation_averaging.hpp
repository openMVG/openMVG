// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_ROTATION_AVERAGING_HPP
#define OPENMVG_MULTIVIEW_ROTATION_AVERAGING_HPP

//--
//-- Implementation related to rotation averaging.
// . Compute global rotation from a list of relative estimates.
// - L2 -> See [1]
// - L1 -> See [2]
//
//- [1] "Robust Multiview Reconstruction."
//- Author : Daniel Martinec.
//- Date : July 2, 2008.
//
//- [2] "Efficient and Robust Large-Scale Rotation Averaging"
//- Authors: Avishek Chatterjee and Venu Madhav Govindu
//- Date: December 2013.
//- Conference: ICCV.
//--


#include "openMVG/multiview/rotation_averaging_common.hpp"
#include "openMVG/multiview/rotation_averaging_l1.hpp"
#include "openMVG/multiview/rotation_averaging_l2.hpp"

#endif // OPENMVG_MULTIVIEW_ROTATION_AVERAGING_HPP
