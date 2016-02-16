// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PARAMS_INCREMENTALSFM_HPP
#define OPENMVG_PARAMS_INCREMENTALSFM_HPP

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include <cereal/cereal.hpp> // Serialization

#include "openMVG/types.hpp"
#include "openMVG/cameras/cameras.hpp"

using namespace openMVG::cameras;

namespace openMVG {
namespace params {

struct paramsIncrementalSfM
{
	// RefineIntrinsics
    //		- 0 intrinsic parameters are kept as constant"
	//		- 1 refine intrinsic parameters (default)"
	bool refineIntrinsics;

	// Camera model for view with unknown intrinsic
	//	- (PINHOLE_CAMERA) 1: Pinhole
	//	- (PINHOLE_CAMERA_RADIAL1) 2: Pinhole radial 1
	//	- (PINHOLE_CAMERA_RADIAL3) 3: Pinhole radial 3 (default)
	//	- (PINHOLE_CAMERA_BROWN) 4: Pinhole brown 2
	//	- (PINHOLE_CAMERA_FISHEYE) 5: Pinhole with a simple Fish-eye distortion
	int camera_type;

	paramsIncrementalSfM(
	bool _refineIntrinsics=1,
	int _camera_type=PINHOLE_CAMERA_RADIAL3)
    :refineIntrinsics(_refineIntrinsics),
	 camera_type(_camera_type)
    {}

  // Serialization
  template <class Archive>
  void serialize( Archive & ar )
  {
    ar(cereal::make_nvp("refineIntrinsics", refineIntrinsics),
    	cereal::make_nvp("camera_type", camera_type));
  }
};

struct paramsGlobalSfM
{
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

	paramsGlobalSfM(
	bool _refineIntrinsics=1,
	int _rotationAveragingMethod=2,
	int _translationAveragingMethod=3)
    :refineIntrinsics(_refineIntrinsics),
	 rotationAveragingMethod(_rotationAveragingMethod),
	 translationAveragingMethod(_translationAveragingMethod)
    {}

  // Serialization
  template <class Archive>
  void serialize( Archive & ar )
  {
    ar(cereal::make_nvp("refineIntrinsics", refineIntrinsics),
    	cereal::make_nvp("rotationAveragingMethod", rotationAveragingMethod),
	cereal::make_nvp("translationAveraging", translationAveragingMethod));
  }
};

} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_INCREMENTALSFM_HPP
