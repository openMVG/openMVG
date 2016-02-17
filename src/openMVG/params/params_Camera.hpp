// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PARAMS_CAMERA_HPP
#define OPENMVG_PARAMS_CAMERA_HPP

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include <cereal/cereal.hpp> // Serialization

#include "openMVG/types.hpp"
#include "openMVG/cameras/cameras.hpp"

using namespace openMVG::cameras;

namespace openMVG {
namespace params {

struct paramsCamera
{
	bool valid;
	// Focal length (in px)
	double focal_px;	//(default) -1 (none)

	// Intrinsics as "fx;0;px;0;fy;py;0;0;1"
	std::string kMatrix;

	// Camera model
	//	- (PINHOLE_CAMERA) 1: Pinhole
	//	- (PINHOLE_CAMERA_RADIAL1) 2: Pinhole radial 1
	//	- (PINHOLE_CAMERA_RADIAL3) 3: Pinhole radial 3 (default)
	//	- (PINHOLE_CAMERA_BROWN) 4: Pinhole brown 2
	//	- (PINHOLE_CAMERA_FISHEYE) 5: Pinhole with a simple Fish-eye distortion
	int camera_type;

	// Group camera model
	//	- True: view can share some camera intrinsic parameters (true)
	//	- False: each view have it's own camera intrinsic parameters
	bool shared_intrinsics;

	paramsCamera(
	bool _valid=false,
	double _focal_px = -1.0,
	std::string _kMatrix = "",
	int _camera_type=PINHOLE_CAMERA_RADIAL3,
	bool _shared_intrinsics=true)
    :valid(_valid),
	 focal_px(_focal_px),
     kMatrix(_kMatrix),
     camera_type(_camera_type),
     shared_intrinsics(_shared_intrinsics)
    {}

};

} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_CAMERA_HPP
