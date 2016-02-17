// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PARAMS_GLOBALSFM_HPP
#define OPENMVG_PARAMS_GLOBALSFM_HPP

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include <cereal/cereal.hpp> // Serialization

#include "openMVG/types.hpp"
#include "openMVG/cameras/cameras.hpp"

using namespace openMVG::cameras;

namespace openMVG {
namespace params {

struct paramsGlobalSfM
{
	bool valid;
	// RefineIntrinsics
    //		- 0 intrinsic parameters are kept as constant"
	//		- 1 refine intrinsic parameters (default)"
	bool refineIntrinsics;

    // RotationAveraging method:"
	//		- 1: L1 minimization
	// 		- 2: L2 minimization (default)"
	int rotationAveragingMethod;

    // TranslationAveraging method:"
    //		- 1: L1 minimization"
    //		- 2: L2 minimization of sum of squared Chordal distances"
    //		- 3: SoftL1 minimization (default)"
	int translationAveragingMethod;

	// Min observations needed for the track
	size_t min_obs_per_track;
	// Min points per pose
	size_t min_obs_per_pose;

	// Threshold for outlier removal
	double outlier_max_residual_error;
	double outlier_min_angle_triangulation;


	paramsGlobalSfM(
	bool _valid=false,
	bool _refineIntrinsics=1,
	int _rotationAveragingMethod=2,
	int _translationAveragingMethod=3,
	size_t _min_obs_per_track=3,
	size_t _min_obs_per_pose=12,
	double _outlier_max_residual_error=4.0,
	double _outlier_min_angle_triangulation=2.0)
    :valid(_valid),
	 refineIntrinsics(_refineIntrinsics),
	 rotationAveragingMethod(_rotationAveragingMethod),
	 translationAveragingMethod(_translationAveragingMethod),
	 min_obs_per_track(_min_obs_per_track),
	 min_obs_per_pose(_min_obs_per_pose),
	 outlier_max_residual_error(_outlier_max_residual_error),
	 outlier_min_angle_triangulation(_outlier_min_angle_triangulation)
    {}
};

} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_GLOBALSFM_HPP
