// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PARAMS_DETECTION_HPP
#define OPENMVG_PARAMS_DETECTION_HPP

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include <cereal/cereal.hpp> // Serialization


namespace openMVG {
namespace params {

struct paramsDetection
{
	// Feature used
	// 		- SIFT (default)
	//		- SIFTGPU
	// 		- AKAZE_FLOAT: AKAZE with floating point descriptors
	//		- AKAZE_MLDB:  AKAZE with binary descriptors
	std::string feature_type;

	// Feature detector configuration (simple way)
	//		- NORMAL (default)
	//		- HIGH
	//		- ULTRA !!Can take long time!!
	std::string feature_preset;

	// Descriptor orientation
	bool upright;


	paramsDetection(
	std::string _feature_type = "SIFT",
	bool _upright=false,
	std::string _feature_preset = "NORMAL")
    :feature_type(_feature_type),
     upright(_upright),
     feature_preset(_feature_preset)
    {}

  // Serialization
  template <class Archive>
  void serialize( Archive & ar )
  {
    ar(cereal::make_nvp("feature_type", feature_type),
    	cereal::make_nvp("upright", upright),
    	cereal::make_nvp("feature_preset", feature_preset));
  }
};

} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_DETECTION_HPP
