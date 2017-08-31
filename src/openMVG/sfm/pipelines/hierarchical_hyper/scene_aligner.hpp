// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#pragma once

#include "openMVG/numeric/numeric.h"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"

namespace ceres {
  class Problem;
}

namespace openMVG {
namespace sfm {

struct SfM_Data;

class SceneAligner;

/**
 * @brief merges two scenes together using the common 3d points of those scenes
 * @param destination_sfm_data where the merged scenes are written
 * @param sfm_data_first the first scene to be merged
 * @param sfm_data_second the second scene to be merged
 * @param common_track_ids the ids of the tracks that are commonly reconstructed in both scenes
 * @param smap_aligner a SceneAligner pointer to choose with which method the scenes should be aligned
 * @return
 */
bool MergeScenesUsingCommonTracks(SfM_Data & destination_sfm_data,
    const SfM_Data & sfm_data_first, // first submap scene
    const SfM_Data & sfm_data_second, // second submap scene
    const std::vector<size_t> & common_track_ids, // which tracks id are commonly reconstructed in both submaps
    SceneAligner *smap_aligner);

/**
 * @brief The SceneAligner class, used to find the transformation between two sfm_data scenes
 * using their common reconstructed 3D points.
 * @note This class does not merge scenes together ! It only finds the transformation between
 * two scenes and the optimized position of the separator landmarks
 */
class SceneAligner
{
public:

  bool computeTransformAndCommonLandmarks
  (Landmarks &destination_landmarks,
    const SfM_Data & sfm_data_first, // first submap scene
    const SfM_Data & sfm_data_second, // second submap scene
    std::vector<double> & second_base_node_pose, // second base node pose with respect to the first base node (what we're looking for)
    double & scaling_factor, // scaling factor to complete the transformation
    const std::vector<size_t> & common_track_ids // which tracks id are commonly reconstructed in both submaps
  );

  explicit SceneAligner(Bundle_Adjustment_Ceres::BA_Ceres_options options = Bundle_Adjustment_Ceres::BA_Ceres_options());

protected:
  Bundle_Adjustment_Ceres::BA_Ceres_options ceres_options_;

  virtual void configureProblem(ceres::Problem & problem,
    Landmarks &destination_landmarks,
    const SfM_Data &sfm_data_first,
    const SfM_Data & sfm_data_second,
    std::vector<double> & second_base_node_pose,
    double & scaling_factor,
    const std::vector<size_t> & common_track_ids
    );
};

/**
 * @brief transforms the content of a sfm data scene into another using a
 * rotation, a translation, and a scale factor
 * @param the original sfm data scene
 * @return the destination sfm data scene
 * @param the rotation matrix
 * @param a translation vector
 * @param a scale factor
 * @note both poses and landmarks are transformed, note that they will be overwritten
 * into the destination scene !
 */
void transformSfMDataScene(SfM_Data & destination_sfm_data,
    const SfM_Data & original_sfm_data,
    const Mat3 & rotation,
    const Vec3 & translation,
    const double scaling_factor);

} // namespace sfm
} // namespace openMVG
