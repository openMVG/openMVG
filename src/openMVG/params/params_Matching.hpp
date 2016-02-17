// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PARAMS_MATCHING_HPP
#define OPENMVG_PARAMS_MATCHING_HPP

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include <cereal/cereal.hpp> // Serialization


namespace openMVG {
namespace params {

struct paramsMatching
{
	// Distance ratio (d1/d2) to discard non meaningful matches (matches with smaller ratio are discarded)
	float max_matching_dist_ratio;	//default 0.8

	// Geometric model for pairwise correspondences filtering thanks to robust model estimation
	// 	- f: (default) fundamental matrix"
	// 	- e: essential matrix"
	// 	- h: homography matrix"
	std::string geometric_model;

	// Sequence matching with overlap of X images
	//	- -1: (default) Off
	//	-  X: match 0 with (1->X), ...]
	//	-  2: match 0 with (1,2), 1 with (2,3), ...
	//	-  3: match 0 with (1,2,3), 1 with (2,3,4), ...
	int video_mode_matching;

	// Nearest neighbor search method
	// 	- AUTO: auto choice from regions type
	//	- For Scalar based regions descriptor:
	//		- BRUTEFORCEL2: L2 BruteForce matching
	//		- ANNL2: L2 Approximate Nearest Neighbor matching
	//		- CASCADEHASHINGL2: L2 Cascade Hashing matching
	//		- FASTCASCADEHASHINGL2: (default)
	//			L2 Cascade Hashing with precomputed hashed regions
	//			(faster than CASCADEHASHINGL2 but use more memory)
	//	- For Binary based descriptor:
	//		- BRUTEFORCEHAMMING: BruteForce Hamming matching
	std::string nearest_matching_method;

	// Guided matching to use the found model to improve the pairwise correspondences
	bool guided_matching;		//(default) false


	paramsMatching(
	float _max_matching_dist_ratio = 0.8,
	std::string _geometric_model = "f",
	int _video_mode_matching=-1,
	std::string _nearest_matching_method = "AUTO",
	bool _guided_matching = false,
	float _geometric_model_initial_residual_tolerance = 4.0,
	bool _pair_overlap_check = true,
	int _min_geometric_feat_matches=50,
	float _min_geo_photo_ratio=0.3)
	:max_matching_dist_ratio(_max_matching_dist_ratio),
	 geometric_model(_geometric_model),
	 video_mode_matching(_video_mode_matching),
	 nearest_matching_method(_nearest_matching_method),
	 guided_matching(_guided_matching)
	{}

	// Serialization
	template <class Archive>
	void serialize( Archive & ar )
	{
		ar(cereal::make_nvp("max_matching_dist_ratio", max_matching_dist_ratio),
			cereal::make_nvp("geometric_model", geometric_model),
		    cereal::make_nvp("video_mode_matching", video_mode_matching),
		    cereal::make_nvp("nearest_matching_method", nearest_matching_method),
		    cereal::make_nvp("guided_matching", guided_matching));


	}
};

} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_MATCHING_HPP
