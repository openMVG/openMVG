// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/test_utility_functions.hpp"
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap.hpp"

void printBaseNodeCoords(const openMVG::Vec3 Translation_gt, openMVG::Vec3 base_node_R, openMVG::Vec3 base_node_t, const openMVG::Vec3 Rotation_gt, double scaling_factor, const double scaling_factor_gt)
{
  std::cout << " new base node coords : \n"
    << base_node_t << std::endl << std::endl
    << base_node_R << std::endl
    << " corresponding angle : " << base_node_R.norm() << std::endl
    << " scale : " << scaling_factor << std::endl;
  std::cout << " ground truth : \n"
    << -Translation_gt << std::endl << std::endl
    << -Rotation_gt << std::endl
    << " corresponding angle : " << Rotation_gt.norm() << std::endl
    << " scale : " << scaling_factor_gt << std::endl;
}

struct TwoScenesConfig
{
  openMVG::sfm::SfM_Data parent_sfm_data;
  openMVG::sfm::SfM_Data sfm_data_A;
  openMVG::sfm::SfM_Data sfm_data_B;
  openMVG::Vec3 Rotation;
  openMVG::Vec3 Translation;
  double scaling_factor_;
  std::vector<size_t> common_track_ids;
};

inline TwoScenesConfig createTwoScenesWithTransformationAndNoise(const int n_poses = 0)
{
  TwoScenesConfig two_scenes;

  openMVG::sfm::Landmarks & landmarks_A = two_scenes.sfm_data_A.structure;
  openMVG::sfm::Landmarks & landmarks_B = two_scenes.sfm_data_B.structure;
  openMVG::sfm::Landmarks & landmarks_parent = two_scenes.parent_sfm_data.structure;
  openMVG::sfm::Poses & poses_A = two_scenes.sfm_data_A.poses;
  openMVG::sfm::Poses & poses_B = two_scenes.sfm_data_B.poses;
  openMVG::sfm::Poses & poses_parent = two_scenes.parent_sfm_data.poses;

  if (n_poses != 0)
  {
    const openMVG::sfm::Poses generated_poses = generateRandomPoses(n_poses);
    poses_A = generated_poses;
    poses_B = generated_poses;
    poses_parent = generated_poses;
  }

  // generate some simple (i.e. non-cluttered together) 3d data
  landmarks_A[0].X = {1,0,0};
  landmarks_A[1].X = {1,1,0};
  landmarks_A[2].X = {-2,0,2};
  landmarks_A[3].X = {0,3,-1};

  // initialize 3d measurement with the same values as the 3d data (those are the common tracks)
  landmarks_B = landmarks_A;
  for (const auto & landmark : landmarks_A)
  {
    two_scenes.common_track_ids.push_back(landmark.first);
  }

  // add some landmarks exclusive to each scene (they should not have the same id)
  int last_id = landmarks_A.size();
  landmarks_A[last_id].X = randomVector(5);
  landmarks_A[last_id+1].X = randomVector(5);
  landmarks_B[last_id+2].X = randomVector(5);
  landmarks_B[last_id+3].X = randomVector(5);
  landmarks_parent = landmarks_A;
  landmarks_parent.insert(landmarks_B.begin(), landmarks_B.end());

  // transformation parameters
  const openMVG::Vec3 rotation_ground_truth = randomVector(M_PI/(2.0 * 1.8)); // 1.8 > sqrt(3)
  two_scenes.Rotation = rotation_ground_truth;
  const openMVG::Vec3 translation_ground_truth = randomVector(2.0);
  two_scenes.Translation = translation_ground_truth;
  const double scaling_factor_ground_truth = randomPositiveDouble(10.0);
  two_scenes.scaling_factor_ = scaling_factor_ground_truth;

  // transform landmarks
  for (auto & lmk_B : landmarks_B)
  {
    // translation
    lmk_B.second.X += translation_ground_truth;

    // rotation
    const openMVG::Vec3 point = lmk_B.second.X;
    ceres::AngleAxisRotatePoint((const double*)rotation_ground_truth.data(), (const double*)point.data(), (double*)lmk_B.second.X.data());

    // scaling
    lmk_B.second.X *= scaling_factor_ground_truth;

    //add some noise
    lmk_B.second.X += randomVector(1e-4);
  }

  // transform poses
  for (auto & mapped_pose : poses_B)
  {
    openMVG::geometry::Pose3 & pose = mapped_pose.second;

    // translation
    openMVG::Vec3 new_center;
    new_center = pose.center() + translation_ground_truth;

    // rotation on z axis
    const openMVG::Vec3 point = new_center;
    ceres::AngleAxisRotatePoint((const double*)rotation_ground_truth.data(), (const double*)point.data(), (double*)new_center.data());
    openMVG::Mat3 new_rotation;
    ceres::AngleAxisToRotationMatrix((const double*)rotation_ground_truth.data(), (double*)new_rotation.data());
    new_rotation = pose.rotation() * new_rotation.inverse();

    // scaling
    new_center *= scaling_factor_ground_truth;

    //add some noise
    new_center += randomVector(1e-4);

    pose = openMVG::geometry::Pose3(new_rotation, new_center);
  }

  // put half of the poses on one scene, the other half on the other scene
  int n_total_cams = poses_B.size();
  for (int i(0); i < n_total_cams; i += 2)
  {
    poses_A.erase(i);
    poses_B.erase(i + 1);
  }

  return two_scenes;
}

inline openMVG::sfm::SfM_Data initializeDestinationSfMData(const TwoScenesConfig & two_scenes)
{
  openMVG::sfm::SfM_Data destination_sfm_data;
  destination_sfm_data.intrinsics = two_scenes.sfm_data_A.intrinsics;
  for (const auto & intrinsic : two_scenes.sfm_data_B.GetIntrinsics())
  {
    if (destination_sfm_data.intrinsics.find(intrinsic.first) == destination_sfm_data.intrinsics.end())
      destination_sfm_data.intrinsics[intrinsic.first] = intrinsic.second;
  }
  destination_sfm_data.poses = two_scenes.sfm_data_A.poses;
  destination_sfm_data.structure = two_scenes.sfm_data_A.structure;

  return destination_sfm_data;
}
