// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_PARAMS_IO_CEREAL_HPP
#define OPENMVG_PARAMS_IO_CEREAL_HPP

#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>

#include <cereal/types/map.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <openMVG/params/params_io.hpp>

#include <iomanip>
#include <fstream>

namespace openMVG {
namespace params {

// --------------------
// Camera data
// --------------------
template <
// JSONInputArchive/ ...
typename archiveType
>
bool Load_Cereal(
  paramsCamera & data,
  const std::string & filename)
{
  const bool bBinary = stlplus::extension_part(filename) == "bin";

  //Create the stream and check it is ok
  std::ifstream stream(filename.c_str(), std::ios::binary | std::ios::in);
  if (!stream.is_open())
    return false;

  // Data serialization
  try
  {
    archiveType archive(stream);

    std::string version;
    archive(cereal::make_nvp("camera_data_version", version));
    archive(cereal::make_nvp("focal_px", data.focal_px));
    archive(cereal::make_nvp("kMatrix", data.kMatrix));
    archive(cereal::make_nvp("camera_type", data.camera_type));
    archive(cereal::make_nvp("shared_intrinsics", data.shared_intrinsics));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }
  stream.close();
  return true;
}

template <
// JSONOutputArchive/ ...
typename archiveType
>
bool Save_Cereal(
  const paramsCamera & data,
  const std::string & filename)
{
  //Create the stream and check it is ok
  std::ofstream stream(filename.c_str(), std::ios::binary | std::ios::out);
  if (!stream.is_open())
    return false;

  // Data serialization
  {
    archiveType archive(stream);
    const std::string version = "0.2";
    archive(cereal::make_nvp("camera_data_version", version));
    archive(cereal::make_nvp("focal_px", data.focal_px));
    archive(cereal::make_nvp("kMatrix", data.kMatrix));
    archive(cereal::make_nvp("camera_type", data.camera_type));
    archive(cereal::make_nvp("shared_intrinsics", data.shared_intrinsics));
  }
  stream.close();
  return true;
}


// --------------------
// Detection data
// --------------------
template <
// JSONInputArchive/ ...
typename archiveType
>
bool Load_Cereal(
  paramsDetection & data,
  const std::string & filename)
{
  const bool bBinary = stlplus::extension_part(filename) == "bin";

  //Create the stream and check it is ok
  std::ifstream stream(filename.c_str(), std::ios::binary | std::ios::in);
  if (!stream.is_open())
    return false;

  // Data serialization
  try
  {
    archiveType archive(stream);

    std::string version;
    archive(cereal::make_nvp("detection_data_version", version));
    archive(cereal::make_nvp("feature_type", data.feature_type));
    archive(cereal::make_nvp("upright", data.upright));
    archive(cereal::make_nvp("feature_preset", data.feature_preset));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }
  stream.close();
  return true;
}

template <
// JSONOutputArchive/ ...
typename archiveType
>
bool Save_Cereal(
  const paramsDetection & data,
  const std::string & filename)
{
  //Create the stream and check it is ok
  std::ofstream stream(filename.c_str(), std::ios::binary | std::ios::out);
  if (!stream.is_open())
    return false;

  // Data serialization
  {
    archiveType archive(stream);
    const std::string version = "0.2";
    archive(cereal::make_nvp("detection_data_version", version));
	archive(cereal::make_nvp("feature_type", data.feature_type));
	archive(cereal::make_nvp("upright", data.upright));
	archive(cereal::make_nvp("feature_preset", data.feature_preset));
  }
  stream.close();
  return true;
}


// --------------------
// Matching data
// --------------------
template <
// JSONInputArchive/ ...
typename archiveType
>
bool Load_Cereal(
  paramsMatching & data,
  const std::string & filename)
{
  const bool bBinary = stlplus::extension_part(filename) == "bin";

  //Create the stream and check it is ok
  std::ifstream stream(filename.c_str(), std::ios::binary | std::ios::in);
  if (!stream.is_open())
    return false;

  // Data serialization
  try
  {
    archiveType archive(stream);

    std::string version;
    archive(cereal::make_nvp("matching_data_version", version));
    archive(cereal::make_nvp("max_matching_dist_ratio", data.max_matching_dist_ratio));
    archive(cereal::make_nvp("geometric_model", data.geometric_model));
    archive(cereal::make_nvp("video_mode_matching", data.video_mode_matching));
    archive(cereal::make_nvp("nearest_matching_method", data.nearest_matching_method));
    archive(cereal::make_nvp("guided_matching", data.guided_matching));
    archive(cereal::make_nvp("geometric_model_initial_residual_tolerance", data.geometric_model_initial_residual_tolerance));
    archive(cereal::make_nvp("min_geometric_feat_matches", data.min_geometric_feat_matches));
    archive(cereal::make_nvp("min_geometric_photo_ratio", data.min_geometric_photo_ratio));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }
  stream.close();
  return true;
}

template <
// JSONOutputArchive/ ...
typename archiveType
>
bool Save_Cereal(
  const paramsMatching & data,
  const std::string & filename)
{
  //Create the stream and check it is ok
  std::ofstream stream(filename.c_str(), std::ios::binary | std::ios::out);
  if (!stream.is_open())
    return false;

  // Data serialization
  {
    archiveType archive(stream);
    const std::string version = "0.2";
    archive(cereal::make_nvp("matching_data_version", version));
    archive(cereal::make_nvp("max_matching_dist_ratio", data.max_matching_dist_ratio));
    archive(cereal::make_nvp("geometric_model", data.geometric_model));
    archive(cereal::make_nvp("video_mode_matching", data.video_mode_matching));
    archive(cereal::make_nvp("nearest_matching_method", data.nearest_matching_method));
    archive(cereal::make_nvp("guided_matching", data.guided_matching));
    archive(cereal::make_nvp("geometric_model_initial_residual_tolerance", data.geometric_model_initial_residual_tolerance));
    archive(cereal::make_nvp("min_geometric_feat_matches", data.min_geometric_feat_matches));
    archive(cereal::make_nvp("min_geometric_photo_ratio", data.min_geometric_photo_ratio));
  }
  stream.close();
  return true;
}


// --------------------
// Incremental SfM data
// --------------------
template <
// JSONInputArchive/ ...
typename archiveType
>
bool Load_Cereal(
  paramsIncrementalSfM & data,
  const std::string & filename)
{
  const bool bBinary = stlplus::extension_part(filename) == "bin";

  //Create the stream and check it is ok
  std::ifstream stream(filename.c_str(), std::ios::binary | std::ios::in);
  if (!stream.is_open())
    return false;

  // Data serialization
  try
  {
    archiveType archive(stream);

    std::string version;
    archive(cereal::make_nvp("incrementalSfM_data_version", version));
    archive(cereal::make_nvp("refineIntrinsics", data.refineIntrinsics));
    archive(cereal::make_nvp("camera_type", data.camera_type));
    archive(cereal::make_nvp("matching_geometric_model", data.matching_geometric_model));
    archive(cereal::make_nvp("min_obs_per_track", data.min_obs_per_track));
    archive(cereal::make_nvp("init_pair_best_of_k", data.init_pair_best_of_k));
    archive(cereal::make_nvp("init_pair_min_tracks", data.init_pair_min_tracks));
    archive(cereal::make_nvp("init_pair_min_angle", data.init_pair_min_angle));
    archive(cereal::make_nvp("init_pair_max_angle", data.init_pair_max_angle));
    archive(cereal::make_nvp("init_pair_pose_init_residual_tolerance", data.init_pair_pose_init_residual_tolerance));
    archive(cereal::make_nvp("init_pair_min_bound_precision_add_point", data.init_pair_min_bound_precision_add_point));
    archive(cereal::make_nvp("init_pair_min_angle_add_point", data.init_pair_min_angle_add_point));
    archive(cereal::make_nvp("sfm_threshold_group_insert_ratio", data.sfm_threshold_group_insert_ratio));
    archive(cereal::make_nvp("sfm_min_bound_residual_error_add_track", data.sfm_min_bound_residual_error_add_track));
    archive(cereal::make_nvp("sfm_min_angle_add_track", data.sfm_min_angle_add_track));
    archive(cereal::make_nvp("ba_min_sparse_schur", data.ba_min_sparse_schur));
    archive(cereal::make_nvp("outlier_max_residual_error_iter", data.outlier_max_residual_error_iter));
    archive(cereal::make_nvp("outlier_min_angle_triangulation_iter", data.outlier_min_angle_triangulation_iter));
    archive(cereal::make_nvp("outlier_max_residual_error_final", data.outlier_max_residual_error_final));
    archive(cereal::make_nvp("outlier_min_angle_triangulation_final", data.outlier_min_angle_triangulation_final));
    archive(cereal::make_nvp("outlier_min_tracks_removed_re_ba", data.outlier_min_tracks_removed_re_ba));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }
  stream.close();
  return true;
}

template <
// JSONOutputArchive/ ...
typename archiveType
>
bool Save_Cereal(
  const paramsIncrementalSfM & data,
  const std::string & filename)
{
  //Create the stream and check it is ok
  std::ofstream stream(filename.c_str(), std::ios::binary | std::ios::out);
  if (!stream.is_open())
    return false;

  // Data serialization
  {
    archiveType archive(stream);
    const std::string version = "0.2";
    archive(cereal::make_nvp("incrementalSfM_data_version", version));
    archive(cereal::make_nvp("refineIntrinsics", data.refineIntrinsics));
    archive(cereal::make_nvp("camera_type", data.camera_type));
    archive(cereal::make_nvp("matching_geometric_model", data.matching_geometric_model));
    archive(cereal::make_nvp("min_obs_per_track", data.min_obs_per_track));
    archive(cereal::make_nvp("init_pair_best_of_k", data.init_pair_best_of_k));
    archive(cereal::make_nvp("init_pair_min_tracks", data.init_pair_min_tracks));
    archive(cereal::make_nvp("init_pair_min_angle", data.init_pair_min_angle));
    archive(cereal::make_nvp("init_pair_max_angle", data.init_pair_max_angle));
    archive(cereal::make_nvp("init_pair_pose_init_residual_tolerance", data.init_pair_pose_init_residual_tolerance));
    archive(cereal::make_nvp("init_pair_min_bound_precision_add_point", data.init_pair_min_bound_precision_add_point));
    archive(cereal::make_nvp("init_pair_min_angle_add_point", data.init_pair_min_angle_add_point));
    archive(cereal::make_nvp("sfm_threshold_group_insert_ratio", data.sfm_threshold_group_insert_ratio));
    archive(cereal::make_nvp("sfm_min_bound_residual_error_add_track", data.sfm_min_bound_residual_error_add_track));
    archive(cereal::make_nvp("sfm_min_angle_add_track", data.sfm_min_angle_add_track));
    archive(cereal::make_nvp("ba_min_sparse_schur", data.ba_min_sparse_schur));
    archive(cereal::make_nvp("outlier_max_residual_error_iter", data.outlier_max_residual_error_iter));
    archive(cereal::make_nvp("outlier_min_angle_triangulation_iter", data.outlier_min_angle_triangulation_iter));
    archive(cereal::make_nvp("outlier_max_residual_error_final", data.outlier_max_residual_error_final));
    archive(cereal::make_nvp("outlier_min_angle_triangulation_final", data.outlier_min_angle_triangulation_final));
    archive(cereal::make_nvp("outlier_min_tracks_removed_re_ba", data.outlier_min_tracks_removed_re_ba));
  }
  stream.close();
  return true;
}


// --------------------
// Global SfM data
// --------------------
template <
// JSONInputArchive/ ...
typename archiveType
>
bool Load_Cereal(
  paramsGlobalSfM & data,
  const std::string & filename)
{
  const bool bBinary = stlplus::extension_part(filename) == "bin";

  //Create the stream and check it is ok
  std::ifstream stream(filename.c_str(), std::ios::binary | std::ios::in);
  if (!stream.is_open())
    return false;

  // Data serialization
  try
  {
    archiveType archive(stream);

    std::string version;
    archive(cereal::make_nvp("globalSfM_data_version", version));
    archive(cereal::make_nvp("refineIntrinsics", data.refineIntrinsics));
    archive(cereal::make_nvp("rotationAveragingMethod", data.rotationAveragingMethod));
    archive(cereal::make_nvp("translationAveragingMethod", data.translationAveragingMethod));
    archive(cereal::make_nvp("min_obs_per_track", data.min_obs_per_track));
    archive(cereal::make_nvp("min_obs_per_pose", data.min_obs_per_pose));
    archive(cereal::make_nvp("outlier_max_residual_error", data.outlier_max_residual_error));
    archive(cereal::make_nvp("outlier_min_angle_triangulation", data.outlier_min_angle_triangulation));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }
  stream.close();
  return true;
}

template <
// JSONOutputArchive/ ...
typename archiveType
>
bool Save_Cereal(
  const paramsGlobalSfM & data,
  const std::string & filename)
{
  //Create the stream and check it is ok
  std::ofstream stream(filename.c_str(), std::ios::binary | std::ios::out);
  if (!stream.is_open())
    return false;

  // Data serialization
  {
    archiveType archive(stream);
    const std::string version = "0.2";
    archive(cereal::make_nvp("globalSfM_data_version", version));
    archive(cereal::make_nvp("refineIntrinsics", data.refineIntrinsics));
    archive(cereal::make_nvp("rotationAveragingMethod", data.rotationAveragingMethod));
    archive(cereal::make_nvp("translationAveragingMethod", data.translationAveragingMethod));
    archive(cereal::make_nvp("min_obs_per_track", data.min_obs_per_track));
    archive(cereal::make_nvp("min_obs_per_pose", data.min_obs_per_pose));
    archive(cereal::make_nvp("outlier_max_residual_error", data.outlier_max_residual_error));
    archive(cereal::make_nvp("outlier_min_angle_triangulation", data.outlier_min_angle_triangulation));
  }
  stream.close();
  return true;
}


} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_IO_CEREAL_HPP
