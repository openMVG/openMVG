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

	// Geometric model used in matching for pairwise correspondences filtering thanks to robust model estimation
	// 	- f: (default) fundamental matrix"
	// 	- e: essential matrix"
	// 	- h: homography matrix"
	std::string matching_geometric_model;

	/// ----------------------------
	/// Tracks
	/// ----------------------------
	// Min observations needed for the track to exist
	size_t min_obs_per_track = 2;

	/// ----------------------------
	/// Init Pair Selection
	/// ----------------------------
	// Select best conditioned pair between k pairs with highest number of tracks
	unsigned init_pair_best_of_k = 20;
	// Min tracks between a pair to be considered for initial pair
	unsigned init_pair_min_tracks = 100;
	// Min and Max angle between cameras for the pair to be considered as initial
	float init_pair_min_angle = 3.0f;
	float init_pair_max_angle = 60.0f;
	// Initial residual tolerance for ACRANSAC for initial pose estimation
	// The ACRANSAC upperbound is square of this number
	double init_pair_pose_init_residual_tolerance = 4.0f;

	// Min residual error for points to be added to the structure
	double init_pair_min_bound_precision_add_point = 1.0;
	// Min angle between triangulation rays for points to be added to the structure
	double init_pair_min_angle_add_point = 2.0;

	/// ----------------------------
	/// SfM
	/// ----------------------------
	// Add all views in one iteration that have at least ratio*N% of the number of the matches of the best image
	float sfm_threshold_group_insert_ratio = 0.75;
	// Min error for new track (if ACRANSAC is bigger error that will be used)
	double sfm_min_bound_residual_error_add_track = 4.0;
	double sfm_min_angle_add_track = 2.0;

	/// ----------------------------
	/// BA
	/// ----------------------------
	// Min number of views needed to use sparse schur instead of dense
	size_t ba_min_sparse_schur = 100;


	/// ----------------------------
	/// Outlier rejection
	/// ----------------------------

	// Threshold for outlier removal after each iteration of BA
	double outlier_max_residual_error_iter = 4.0;
	double outlier_min_angle_triangulation_iter = 2.0;
	// Threshold for outlier removal after final iteration of BA
	double outlier_max_residual_error_final = 4.0;
	double outlier_min_angle_triangulation_final = 2.0;
	// Re-run BA if at least this number of tracks are removed
	size_t outlier_min_tracks_removed_re_ba = 50;

	paramsIncrementalSfM(
	bool _refineIntrinsics=1,
	int _camera_type=PINHOLE_CAMERA_RADIAL3,
	std::string _matching_geometric_model="f",
	size_t _min_obs_per_track = 2,
	unsigned _init_pair_best_of_k = 20,
	unsigned _init_pair_min_tracks = 100,
	float _init_pair_min_angle = 3.0f,
	float _init_pair_max_angle = 60.0f,
	double _init_pair_pose_init_residual_tolerance = 4.0f,
	double _init_pair_min_bound_precision_add_point = 1.0,
	double _init_pair_min_angle_add_point = 2.0,
	float _sfm_threshold_group_insert_ratio = 0.75,
	double _sfm_min_bound_residual_error_add_track = 4.0,
	double _sfm_min_angle_add_track = 2.0,
	size_t _ba_min_sparse_schur = 100,
	double _outlier_max_residual_error_iter = 4.0,
	double _outlier_min_angle_triangulation_iter = 2.0,
	double _outlier_max_residual_error_final = 4.0,
	double _outlier_min_angle_triangulation_final = 2.0,
	size_t _outlier_min_tracks_removed_re_ba = 50)
    :refineIntrinsics(_refineIntrinsics),
	 camera_type(_camera_type),
	 matching_geometric_model(_matching_geometric_model),
	 min_obs_per_track(_min_obs_per_track),
	 init_pair_best_of_k(_init_pair_best_of_k),
	 init_pair_min_tracks(_init_pair_min_tracks),
	 init_pair_min_angle(_init_pair_min_angle),
	 init_pair_max_angle(_init_pair_max_angle),
	 init_pair_pose_init_residual_tolerance(_init_pair_pose_init_residual_tolerance),
	 init_pair_min_bound_precision_add_point(_init_pair_min_bound_precision_add_point),
	 init_pair_min_angle_add_point(_init_pair_min_angle_add_point),
	 sfm_threshold_group_insert_ratio(_sfm_threshold_group_insert_ratio),
	 sfm_min_bound_residual_error_add_track(_sfm_min_bound_residual_error_add_track),
	 sfm_min_angle_add_track(_sfm_min_angle_add_track),
	 ba_min_sparse_schur(_ba_min_sparse_schur),
	 outlier_max_residual_error_iter(_outlier_max_residual_error_iter),
	 outlier_min_angle_triangulation_iter(_outlier_min_angle_triangulation_iter),
	 outlier_max_residual_error_final(_outlier_max_residual_error_final),
	 outlier_min_angle_triangulation_final(_outlier_min_angle_triangulation_final),
	 outlier_min_tracks_removed_re_ba(_outlier_min_tracks_removed_re_ba)
    {}

};

} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_INCREMENTALSFM_HPP
