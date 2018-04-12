// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner.hpp" // what we're testing
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_test_utilities.hpp"

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
      continue;

    EXPECT_MATRIX_EQ(landmark.X, destination_sfm_data.structure.at(landmark_id).X);
  }
}

TEST(SceneAligner, sceneWithLandmarksOnly_SubmapAlignmentWorks)
{
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise();

  const Similarity3 sim_gt = two_scenes.sim;
  Similarity3 sim;

  SfM_Data destination_sfm_data = initializeDestinationSfMData(two_scenes);
  Landmarks & destination_landmarks = destination_sfm_data.structure;

  cout << " Bundle Adjustment ... " << endl;
  std::unique_ptr<SceneAligner> smap_aligner =
      std::unique_ptr<SceneAligner>(new SceneAligner);
  smap_aligner->computeTransformAndCommonLandmarks(
      destination_landmarks, two_scenes.sfm_data_A, two_scenes.sfm_data_B,
      sim, two_scenes.common_track_ids);

  printBaseNodeCoords(sim_gt, sim);

  // test if the bundle adjustment did a good job on evaluating the transformation
  EXPECT_NEAR(sim_gt.scale_, sim.scale_, 1e-2);
  EXPECT_MATRIX_NEAR(sim_gt.pose_.center(), sim.pose_.center(), 1e-2);
  EXPECT_MATRIX_NEAR(sim_gt.pose_.rotation(), sim.pose_.rotation(), 1e-2);
}

TEST(SceneAligner, fullyCalibratedSceneWithLandmarksAndPoses_SubmapRegistrationWorks)
{
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise(3);
  SfM_Data destination_sfm_data = initializeDestinationSfMData(two_scenes);
  Landmarks & destination_landmarks = destination_sfm_data.structure;
  const Similarity3 sim_gt = two_scenes.sim;
  Similarity3 sim;

  std::unique_ptr<SceneAligner> smap_aligner =
      std::unique_ptr<SceneAligner>(new SceneAligner);

  smap_aligner->computeTransformAndCommonLandmarks(
      destination_landmarks, two_scenes.sfm_data_A, two_scenes.sfm_data_B,
      sim, two_scenes.common_track_ids);

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
  for (auto & track_id : two_scenes.common_track_ids)
  {
    two_scenes.sfm_data_A.structure.erase(track_id);
    two_scenes.sfm_data_B.structure.erase(track_id);
    two_scenes.common_track_ids.erase(track_id);
  }

  std::unique_ptr<SceneAligner> smap_aligner =
      std::unique_ptr<SceneAligner>(new SceneAligner);

  SfM_Data destination_sfm_data = initializeDestinationSfMData(two_scenes);
  Landmarks & destination_landmarks = destination_sfm_data.structure;
  Similarity3 sim;

  EXPECT_FALSE(smap_aligner->computeTransformAndCommonLandmarks(
      destination_landmarks, two_scenes.sfm_data_A, two_scenes.sfm_data_B,
      sim, two_scenes.common_track_ids));
}

TEST(SceneAligner, twoScenesWithNoCommonRECONSTRUCTEDTracks_returnFalse)
{
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise(10);

  // remove common tracks from each scene but not from common tracks !
  for (const auto & track_id : two_scenes.common_track_ids)
  {
    two_scenes.sfm_data_A.structure.erase(track_id);
    two_scenes.sfm_data_B.structure.erase(track_id);
  }

  std::unique_ptr<SceneAligner> smap_aligner =
      std::unique_ptr<SceneAligner>(new SceneAligner);

  SfM_Data destination_sfm_data = initializeDestinationSfMData(two_scenes);
  Landmarks & destination_landmarks = destination_sfm_data.structure;
  Similarity3 sim;

  EXPECT_FALSE(smap_aligner->computeTransformAndCommonLandmarks(
      destination_landmarks, two_scenes.sfm_data_A, two_scenes.sfm_data_B,
      sim, two_scenes.common_track_ids));
}

// checks that distances are conserved (up to a scale) when
// transforming a sfm data scene
TEST(transformSfMDataSceneInto, genericScene_transformationLeavesDistancesConserved)
{
  const Similarity3 sim = randomSimilarity();

  const openMVG::IndexT n_poses(20), n_landmarks(1000);
  const openMVG::sfm::SfM_Data original_sfm_data = generate_random_poses_and_landmarks_in_scene(n_poses, n_landmarks);
  openMVG::sfm::SfM_Data destination_sfm_data;

  const openMVG::Mat distances_before = computeDistancesBetweenPosesAndLandmarks(original_sfm_data);

  openMVG::sfm::transformSfMDataSceneInto(destination_sfm_data, original_sfm_data, sim);

  openMVG::Mat distances_after = computeDistancesBetweenPosesAndLandmarks(destination_sfm_data);

  std::cout << sim.scale_ << std::endl;
  distances_after *= sim.scale_;
  EXPECT_MATRIX_NEAR(distances_before, distances_after, 1e-7);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
