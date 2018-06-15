// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner.hpp" // what we're testing
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_test_utilities.hpp"
#include "openMVG/sfm/pipelines/pipelines_test.hpp"

#include "openMVG/test_utility_functions.hpp"
#include "third_party/ceres-solver/include/ceres/rotation.h"
#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;
using namespace std;

TEST(MergeScenesUsingCommonTracks, normalScene_destinationSceneContainsEverything)
{
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise(3);
  SfM_Data destination_sfm_data;

  cout << " Bundle Adjustment ... " << endl;
  SceneAligner smap_aligner;
  MergeScenesUsingCommonTracks(destination_sfm_data, two_scenes.sfm_data_A, two_scenes.sfm_data_B, two_scenes.common_track_ids, &smap_aligner);

  EXPECT_EQ(destination_sfm_data.GetLandmarks().size(), two_scenes.parent_sfm_data.GetLandmarks().size());
  EXPECT_EQ(destination_sfm_data.GetPoses().size(), two_scenes.parent_sfm_data.GetPoses().size());
  EXPECT_EQ(destination_sfm_data.GetIntrinsics().size(), two_scenes.parent_sfm_data.GetIntrinsics().size());
  EXPECT_EQ(destination_sfm_data.GetViews().size(), two_scenes.parent_sfm_data.GetViews().size());
}

TEST(MergeScenesUsingCommonTracks, normalScene_firstSceneOnlyLandmarksKeptConstant)
{
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise(3);
  SfM_Data destination_sfm_data;

  cout << " Bundle Adjustment ... " << endl;
  SceneAligner smap_aligner;
  MergeScenesUsingCommonTracks(
      destination_sfm_data, two_scenes.sfm_data_A, two_scenes.sfm_data_B, two_scenes.common_track_ids, &smap_aligner);

  for (const auto & mapped_landmark : two_scenes.sfm_data_A.GetLandmarks())
  {
    const auto landmark_id = mapped_landmark.first;
    const auto landmark = mapped_landmark.second;

    // separator landmarks are not kept constant by definition
    if (two_scenes.common_track_ids.count(landmark_id) > 0)
    {
      EXPECT_FALSE((landmark.X == destination_sfm_data.structure.at(landmark_id).X));
      continue;
    }

    EXPECT_MATRIX_EQ(landmark.X, destination_sfm_data.structure.at(landmark_id).X);
  }
}

TEST(MergeScenesUsingCommonTracks, checkSeparatorLandmarksGetMergedCorrectly)
{
  SfM_Data scene1, scene2, destination_scene;
  const Observations observations_1_0 = {{0, Observation({0,0}, 1000)},{1, Observation({1,1}, 1001)},
                                        {2, Observation({2,2}, 1002)},{3, Observation({3,3}, 1003)}};
  const Observations observations_2_0 = {{4, Observation({4,4}, 1004)},{5, Observation({5,5}, 1005)}};
  const Observations observations_1_1 = {{6, Observation({6,6}, 1006)},{7, Observation({7,7}, 1007)}};
  const Observations observations_2_2 = {{8, Observation({8,8}, 1008)},{9, Observation({9,9}, 1009)}};

  // in this example, track 0 is the separator track
  const IndexT separatorTrackId = 0;
  const std::set<IndexT> common_track_ids = {separatorTrackId};

  auto x0 = randomVector();
  scene1.structure[0] = {x0, observations_1_0};
  scene1.structure[1] = {randomVector(), observations_1_1};
  scene2.structure[0] = {randomVector(), observations_2_0};
  scene2.structure[2] = {randomVector(), observations_2_2};

  SceneAligner smap_aligner;
  MergeScenesUsingCommonTracks(
        destination_scene, scene1, scene2, common_track_ids, &smap_aligner);

  EXPECT_EQ(destination_scene.structure.at(separatorTrackId).obs.size(),
            observations_1_0.size() + observations_2_0.size());
  EXPECT_MATRIX_NEAR(destination_scene.structure[separatorTrackId].X, x0, 1e-2); // check we are in the right referential
  Observations observations_all_0 = observations_1_0;
  observations_all_0.insert(observations_2_0.cbegin(), observations_2_0.cend());
  for (const auto & obs_pair : observations_all_0)
  {
    const auto view_id = obs_pair.first;
    const auto & destination_obs = destination_scene.structure.at(separatorTrackId).obs.at(view_id);
    EXPECT_EQ(obs_pair.second.x, destination_obs.x);
    EXPECT_EQ(obs_pair.second.id_feat, destination_obs.id_feat);
  }
}

TEST(SceneAligner, sceneWithLandmarksOnly_SubmapAlignmentWorks)
{
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise();

  const Similarity3 sim_gt = two_scenes.sim;
  Similarity3 sim;

  Landmarks separator_landmarks = findCommonLandmarks(two_scenes.sfm_data_A, two_scenes.sfm_data_B);

  cout << " Bundle Adjustment ... " << endl;
  std::unique_ptr<SceneAligner> smap_aligner =
      std::unique_ptr<SceneAligner>(new SceneAligner);
  smap_aligner->computeTransformAndCommonLandmarks(
      separator_landmarks, two_scenes.sfm_data_A, two_scenes.sfm_data_B,
      sim);

  printBaseNodeCoords(sim_gt, sim);

  // test if the bundle adjustment did a good job on evaluating the transformation
  EXPECT_NEAR(sim_gt.scale_, sim.scale_, 1e-2);
  EXPECT_MATRIX_NEAR(sim_gt.pose_.center(), sim.pose_.center(), 1e-2);
  EXPECT_MATRIX_NEAR(sim_gt.pose_.rotation(), sim.pose_.rotation(), 1e-2);
}

TEST(SceneAligner, fullyCalibratedSceneWithLandmarksAndPoses_SubmapRegistrationWorks)
{
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise(3);
  Landmarks separator_landmarks = findCommonLandmarks(two_scenes.sfm_data_A, two_scenes.sfm_data_B);
  const Similarity3 sim_gt = two_scenes.sim;
  Similarity3 sim;

  std::unique_ptr<SceneAligner> smap_aligner =
      std::unique_ptr<SceneAligner>(new SceneAligner);

  smap_aligner->computeTransformAndCommonLandmarks(
      separator_landmarks, two_scenes.sfm_data_A, two_scenes.sfm_data_B,
      sim);

  printBaseNodeCoords(sim_gt, sim);

  // test if the bundle adjustment did a good job on evaluating the transformation
  EXPECT_NEAR(sim_gt.scale_, sim.scale_, 1e-2);
  EXPECT_MATRIX_NEAR(sim_gt.pose_.center(), sim.pose_.center(), 1e-2);
  EXPECT_MATRIX_NEAR(sim_gt.pose_.rotation(), sim.pose_.rotation(), 1e-2);
}

TEST(SceneAligner, twoScenesWithNoCommonTracks_returnFalse)
{
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise(10);

  // remove common tracks from each scene
  for (auto track_it = two_scenes.common_track_ids.begin(); track_it != two_scenes.common_track_ids.end();)
  {
    two_scenes.sfm_data_A.structure.erase(*track_it);
    two_scenes.sfm_data_B.structure.erase(*track_it);
    track_it = two_scenes.common_track_ids.erase(track_it);
  }

  std::unique_ptr<SceneAligner> smap_aligner =
      std::unique_ptr<SceneAligner>(new SceneAligner);

  Landmarks separator_landmarks; // no common tracks
  Similarity3 sim;

  EXPECT_FALSE(smap_aligner->computeTransformAndCommonLandmarks(
      separator_landmarks, two_scenes.sfm_data_A, two_scenes.sfm_data_B,
      sim));
}

// checks that distances are conserved (up to a scale) when
// transforming a sfm data scene
TEST(copySfMDataSceneInto, genericScene_transformationLeavesDistancesConserved)
{
  const Similarity3 sim = randomSimilarity();

  const openMVG::IndexT n_poses(20), n_landmarks(1000);
  const openMVG::sfm::SfM_Data original_sfm_data = generate_random_poses_and_landmarks_in_scene(n_poses, n_landmarks);
  openMVG::sfm::SfM_Data destination_sfm_data;

  const openMVG::Mat distances_before = computeDistancesBetweenPosesAndLandmarks(original_sfm_data);

  openMVG::sfm::copySfMDataSceneInto(destination_sfm_data, original_sfm_data, sim);

  openMVG::Mat distances_after = computeDistancesBetweenPosesAndLandmarks(destination_sfm_data);

  std::cout << sim.scale_ << std::endl;
  distances_after *= sim.scale_;
  EXPECT_MATRIX_NEAR(distances_before, distances_after, 1e-7);
}

TEST(copySfMDataSceneInto, genericScene_dataGetsCopiedAsExpected)
{
  int n_poses = 20;
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise(n_poses);

  const SfM_Data & sceneA     = two_scenes.sfm_data_A;
  const SfM_Data & sceneB     = two_scenes.sfm_data_B;
  const SfM_Data & sceneTotal = two_scenes.parent_sfm_data;
  SfM_Data destination_sfm_data;

  openMVG::sfm::copySfMDataSceneInto(destination_sfm_data, sceneA, Similarity3());

  EXPECT_EQ(sceneA.views.size(), destination_sfm_data.views.size());
  EXPECT_EQ(sceneA.intrinsics.size(), destination_sfm_data.intrinsics.size());
  EXPECT_EQ(sceneA.poses.size(), destination_sfm_data.poses.size());
  EXPECT_EQ(sceneA.structure.size(), destination_sfm_data.structure.size());
  EXPECT_EQ(sceneA.control_points.size(), destination_sfm_data.control_points.size());
  EXPECT_EQ(sceneA.s_root_path, destination_sfm_data.s_root_path);

  const Similarity3 sim = two_scenes.sim;
  copySfMDataSceneInto(destination_sfm_data, sceneB, sim);

  EXPECT_EQ(sceneTotal.views.size(), destination_sfm_data.views.size());
  EXPECT_EQ(sceneTotal.intrinsics.size(), destination_sfm_data.intrinsics.size());
  EXPECT_EQ(sceneTotal.poses.size(), destination_sfm_data.poses.size());
  EXPECT_EQ(sceneTotal.structure.size(), destination_sfm_data.structure.size());
  EXPECT_EQ(sceneTotal.control_points.size(), destination_sfm_data.control_points.size());
  EXPECT_EQ(sceneTotal.s_root_path, destination_sfm_data.s_root_path);
}

TEST(getCommonReconstructedLandmarks, findingCommonLandmarksWorksAsExpected)
{
  SfM_Data scene1, scene2;
  const Observations observations_1_0 = {{0, Observation()},{1, Observation()},
                                        {2, Observation()},{3, Observation()}};
  const Observations observations_2_0 = {{4, Observation()},{5, Observation()}};
  const Observations observations_1_1 = {{6, Observation()},{7, Observation()}};
  const Observations observations_2_2 = {{8, Observation()},{9, Observation()}};

  auto x0 = randomVector();

  scene1.structure[0] = {x0, observations_1_0};
  scene1.structure[1] = {randomVector(), observations_1_1};
  scene2.structure[0] = {randomVector(), observations_2_0};
  scene2.structure[2] = {randomVector(), observations_2_2};

  Landmarks commonLandmarks = getCommonReconstructedLandmarks(scene1, scene2, {0});

  EXPECT_EQ(commonLandmarks.size(), 1);

  // check that we stay in first scene's referential
  EXPECT_EQ(commonLandmarks.at(0).X, x0);
  // check that observations are not copied in the commonlandmarks
  EXPECT_EQ(commonLandmarks.at(0).obs.size(), 0);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
