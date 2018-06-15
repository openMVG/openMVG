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

#include "openMVG/geometry/Similarity3.hpp"

namespace ceres {
  class Problem;
}

namespace openMVG {
namespace sfm {

class SceneAligner;

/**
 * @brief finds the common landmarks that are RECONSTRUCTED IN BOTH SCENES. Search is
 * constrained to track ids contained in track_ids. Returns a copy of those landmarks
 * without observation data.
 * @param sfm_data_1 : note that we are in this scene's referential !
 * @param sfm_data_2 : the second scene
 * @param track_ids : the track ids in which the search is constrained.
 * @return common landmarks in the referential of the first scene, but WITHOUT OBSERVATIONS
 */
Landmarks getCommonReconstructedLandmarks(
    const SfM_Data& sfm_data_1,
    const SfM_Data& sfm_data_2,
    const std::set<IndexT>& track_ids);

/**
 * @brief merges two scenes together using the common 3d points of those scenes
 * @param destination_sfm_data where the merged scenes are written
 * @param sfm_data_first the first scene to be merged
 * @param sfm_data_second the second scene to be merged
 * @param common_track_ids the ids of the tracks that are commonly reconstructed in both scenes
 * @param smap_aligner a SceneAligner pointer to choose with which method the scenes should be aligned
 *
 * When merging two scenes together, we compute simultanouesly the similarity transform between the two scenes
 * and the position of the landmarks common to both scenes (reconstructed separator tracks).
 * The way it optimizes the similarity transform and the common points positions is by minimizing the 
 * resulting distance (i.e. the error) between points in a submap and their counterparts in the 
 * other submap.
 *
 * The merged 3d scene (SfMData object) is created in this way :
 * - Views : There are no common views between two submaps so we just copy views from both children in the merged scene.
 * - Intrinsics :
 *   - Intrinsics which are exclusive to only one of the children submaps get copied directly to the scene.
 *   - Intrinsics which are common to both children submaps can have different values,
 *     typically when submap reconstruction was done with intrinsics optimization.
 *     In that case, we copy the intrinsic which minimizes the overall reprojection error in the merged scene.
 * - Structure :
 *   - Landmarks which are exclusive to a single child submap get copied to the scene,
 *     we apply a similarity transform to them so that they are in the referential of
 *     the merged scene.
 *   - Landmarks which are common to both submaps (separator landmarks) are NOT copied from the children submaps.
 *     The position of these landmarks is computed at the same time we compute the similarity transform
 *     so we already know their position.
 * - Poses : Poses are exclusive to each submaps so we just copy their value, up to a similarity transform to
 *           put them in the referential of the merged scene.
 */
bool MergeScenesUsingCommonTracks(SfM_Data & destination_sfm_data,
    const SfM_Data & sfm_data_fst, // first submap scene
    const SfM_Data & sfm_data_snd, // which tracks id are commonly reconstructed in both submaps
    const std::set<IndexT> & separator_track_ids,
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
    geometry::Similarity3 &sim);

  explicit SceneAligner(Bundle_Adjustment_Ceres::BA_Ceres_options options = Bundle_Adjustment_Ceres::BA_Ceres_options());

protected:
  Bundle_Adjustment_Ceres::BA_Ceres_options ceres_options_;

  virtual bool checkScenesAreAlignable(
      const SfM_Data & sfm_data_first,
      const SfM_Data & sfm_data_second,
      const Landmarks & common_landmarks);

  virtual void configureProblem(ceres::Problem & problem,
    Landmarks &destination_landmarks,
    const SfM_Data &sfm_data_first,
    const SfM_Data & sfm_data_second,
    std::vector<double> &second_base_node_pose,
    double &scaling_factor);
};

/**
 * @brief takes a SfM_Data scene and copies it into a destination scene, with a geometric transformation.
 * @note mostly used for merging submaps together
 * @warning root path will be overwritten (at least in current version of the function)
 * @warning intrinsics copying is a deep copy but not views copy... for now at least
 * @param destination_sfm_data : the destination scene, can already contain views, landmarks ... etc
 * @param original_sfm_data : the scene to be transformed into the destination scene
 * @param sim : the similarity transformation to be applied to the original scene
 */
void copySfMDataSceneInto(SfM_Data & destination_sfm_data,
    const SfM_Data & original_sfm_data,
    const geometry::Similarity3 &sim);

} // namespace sfm
} // namespace openMVG
