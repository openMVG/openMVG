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

void printBaseNodeCoords(const openMVG::geometry::Similarity3& sim_gt, const openMVG::geometry::Similarity3& sim)
{
  std::cout << "new base node coords : \n"
            << sim.pose_.center() << std::endl << std::endl
            << sim.pose_.rotation() << std::endl
            << "scale : " << sim.scale_ << std::endl;
  std::cout << "ground truth : \n"
            << sim_gt.pose_.center() << std::endl << std::endl
            << sim_gt.pose_.rotation() << std::endl
            << "scale : " << sim_gt.scale_ << std::endl;
}

struct TwoScenesConfig
{
  openMVG::sfm::SfM_Data parent_sfm_data;
  openMVG::sfm::SfM_Data sfm_data_A;
  openMVG::sfm::SfM_Data sfm_data_B;
  openMVG::geometry::Similarity3 sim;
  std::set<openMVG::IndexT> common_track_ids;
};

void addLandmarksToTwoScenes(TwoScenesConfig & two_scenes)
{
  auto& landmarks_parent = two_scenes.parent_sfm_data.structure;
  auto& landmarks_A      = two_scenes.sfm_data_A.structure;
  auto& landmarks_B      = two_scenes.sfm_data_B.structure;

  landmarks_A[0].X = {1,0,0};
  landmarks_A[1].X = {1,1,0};
  landmarks_A[2].X = {-2,0,2};
  landmarks_A[3].X = {0,3,-1};

  // initialize 3d measurement with the same values as the 3d data (those are the common tracks)
  landmarks_B = landmarks_A;
  for (const auto & landmark : landmarks_A)
  {
    two_scenes.common_track_ids.insert(landmark.first);
  }

  // add some landmarks exclusive to each scene (they should not have the same id)
  int last_id = landmarks_A.size();
  landmarks_A[last_id].X = randomVector(5);
  landmarks_A[last_id+1].X = randomVector(5);
  landmarks_B[last_id+2].X = randomVector(5);
  landmarks_B[last_id+3].X = randomVector(5);
  landmarks_parent = landmarks_A;
  landmarks_parent.insert(landmarks_B.begin(), landmarks_B.end());
}

void addNoiseToTwoScenes(TwoScenesConfig& two_scenes)
{
  auto& landmarks_A = two_scenes.sfm_data_A.structure;
  auto& poses_A     = two_scenes.sfm_data_A.poses;
  auto& landmarks_B = two_scenes.sfm_data_B.structure;
  auto& poses_B     = two_scenes.sfm_data_B.poses;
  for (auto & lmk : landmarks_A)
  {
    lmk.second.X += randomVector(1e-4);
  }
  for (auto & lmk : landmarks_B)
  {
    lmk.second.X += randomVector(1e-4);
  }
  for (auto & mapped_pose : poses_A)
  {
    mapped_pose.second.center() += randomVector(1e-4);
  }
  for (auto & mapped_pose : poses_B)
  {
    mapped_pose.second.center() += randomVector(1e-4);
  }
}

void dispatchPosesInChildrenScenes(TwoScenesConfig& two_scenes)
{
  // put half of the poses on one scene, the other half on the other scene
  auto& poses_A = two_scenes.sfm_data_A.poses;
  auto& poses_B = two_scenes.sfm_data_B.poses;
  int n_total_cams = poses_B.size();
  for (int i(0); i < n_total_cams; i += 2)
  {
    poses_A.erase(i);
    poses_B.erase(i + 1);
  }
}

inline TwoScenesConfig createTwoScenesWithTransformationAndNoise(const int n_poses = 0)
{
  TwoScenesConfig two_scenes;

  auto & poses_A           = two_scenes.sfm_data_A.poses;
  auto & poses_B           = two_scenes.sfm_data_B.poses;
  auto & poses_parent      = two_scenes.parent_sfm_data.poses;
  auto & views_A           = two_scenes.sfm_data_A.views;
  auto & views_B           = two_scenes.sfm_data_B.views;
  auto & views_parent      = two_scenes.parent_sfm_data.views;
  auto & intrinsics_A      = two_scenes.sfm_data_A.intrinsics;
  auto & intrinsics_B      = two_scenes.sfm_data_B.intrinsics;
  auto & intrinsics_parent = two_scenes.parent_sfm_data.intrinsics;

  if (n_poses != 0)
  {
    const openMVG::sfm::Poses generated_poses = generateRandomPoses(n_poses);
    poses_A = generated_poses;
    poses_B = generated_poses;
    poses_parent = generated_poses;
    for (const auto & pose : generated_poses)
    {
      const openMVG::IndexT id = pose.first;
      openMVG::sfm::View v("v_" + std::to_string(id),
                           id, 0, id);
      std::shared_ptr<openMVG::sfm::View> v_ptr =
          std::make_shared<openMVG::sfm::View>(v);
      views_A[id] = v_ptr;
      views_B[id] = v_ptr;
      views_parent[id] = v_ptr;
    }
  }

  std::shared_ptr<openMVG::cameras::IntrinsicBase> intrinsic_ptr =
      std::make_shared<openMVG::cameras::Pinhole_Intrinsic>
            (1000, 1000, 1000, 500, 500);
  intrinsics_A[0] = intrinsic_ptr;
  intrinsics_B[0] = intrinsic_ptr;
  intrinsics_parent[0] = intrinsic_ptr;

  // generate some simple (i.e. non-cluttered together) 3d data
  addLandmarksToTwoScenes(two_scenes);

  // transform scene B
  const openMVG::geometry::Similarity3 sim_ground_truth = easySimilarity();
  two_scenes.sim = sim_ground_truth;
  ApplySimilarity(sim_ground_truth, two_scenes.sfm_data_B);

  addNoiseToTwoScenes(two_scenes);

  dispatchPosesInChildrenScenes(two_scenes);

  const std::string common_root_path = "/some/path";
  two_scenes.parent_sfm_data.s_root_path = common_root_path;
  two_scenes.sfm_data_A.s_root_path = common_root_path;
  two_scenes.sfm_data_B.s_root_path = common_root_path;

  return two_scenes;
}

inline openMVG::sfm::Landmarks findCommonLandmarks(const openMVG::sfm::SfM_Data & sfm_data_1, const openMVG::sfm::SfM_Data& sfm_data_2)
{
  openMVG::sfm::Landmarks common_landmarks;
  std::copy_if(sfm_data_1.structure.cbegin(), sfm_data_1.structure.cend(),
               std::inserter(common_landmarks, common_landmarks.begin()),
               [&sfm_data_2](const std::pair<openMVG::IndexT, openMVG::sfm::Landmark>& lmk_pair)
                {return sfm_data_2.structure.count(lmk_pair.first) > 0;});
  return common_landmarks;
}
