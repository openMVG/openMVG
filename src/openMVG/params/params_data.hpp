// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PARAMS_DATA_HPP
#define OPENMVG_PARAMS_DATA_HPP

#include "openMVG/types.hpp"
#include "openMVG/params/params_Detection.hpp"
#include "openMVG/params/params_Matching.hpp"
#include <openMVG/params/params_SfM.hpp>


namespace openMVG {
namespace params {

// Store parameter values
struct paramsData
{
	// Detection parameters
	paramsDetection detection;
	// Matching parameters
	paramsMatching matching;
	// IncrementalSfM parameters
	paramsIncrementalSfM incrementalSfM;
	// GlobalSfM parameters
	paramsGlobalSfM globalSfM;
};

} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_DATA_HPP
