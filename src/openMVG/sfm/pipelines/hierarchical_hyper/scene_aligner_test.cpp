#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner.hpp" // what we're testing
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_test_utilities.hpp"

#include "openMVG/test_utility_functions.hpp"

#include "third_party/ceres-solver/include/ceres/rotation.h"

#include "testing/testing.h"
#include "openMVG/cameras/Camera_Common.hpp"

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
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
    if (std::find(two_scenes.common_track_ids.begin(), two_scenes.common_track_ids.end(), landmark_id)
          != two_scenes.common_track_ids.end())
      continue;

    EXPECT_MATRIX_EQ(landmark.X, destination_sfm_data.structure.at(landmark_id).X);
  }
}

TEST(SceneAligner, sceneWithLandmarksOnly_SubmapAlignmentWorks)
{
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise();
  const double scaling_factor_gt = two_scenes.scaling_factor_;
  const Vec3 Rotation_gt = two_scenes.Rotation;
  const Vec3 Translation_gt = two_scenes.Translation;

  std::vector<double> base_node_coords = {0,0,0,0,0,0};
  double scaling_factor(1.0);

  SfM_Data destination_sfm_data = initializeDestinationSfMData(two_scenes);

  cout << " Bundle Adjustment ... " << endl;
  std::unique_ptr<SceneAligner> smap_aligner =
      std::unique_ptr<SceneAligner>(new SceneAligner);
  smap_aligner->computeTransformAndDestinationSeparators(
      destination_sfm_data, two_scenes.sfm_data_A, two_scenes.sfm_data_B,
      base_node_coords, scaling_factor, two_scenes.common_track_ids);

  Vec3 base_node_R = Vec3(base_node_coords[0],base_node_coords[1],base_node_coords[2]);
  Vec3 base_node_t = Vec3(base_node_coords[3],base_node_coords[4],base_node_coords[5]);

  printBaseNodeCoords(Translation_gt, base_node_R, base_node_t, Rotation_gt, scaling_factor, scaling_factor_gt);

  // test if the bundle adjustment did a good job on evaluating the transformation
  EXPECT_NEAR(scaling_factor_gt, scaling_factor, 1e-2);
  const double modulo_rotation = fmod(base_node_R.norm(), M_PI);
  EXPECT_NEAR(modulo_rotation, Rotation_gt.norm(), 1e-2);
  Vec3 inverse_Rotation_gt = -Rotation_gt;
  Vec3 inverse_translation_gt = -Translation_gt;
  EXPECT_MATRIX_NEAR(base_node_R, inverse_Rotation_gt, 1e-2);
  EXPECT_MATRIX_NEAR(base_node_t, inverse_translation_gt, 1e-2);
}

TEST(SceneAligner, fullyCalibratedSceneWithLandmarksAndPoses_SubmapRegistrationWorks)
{
  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise(3);
  SfM_Data destination_sfm_data = initializeDestinationSfMData(two_scenes);
  const double scaling_factor_gt = two_scenes.scaling_factor_;
  const Vec3 Rotation_gt = two_scenes.Rotation;
  const Vec3 Translation_gt = two_scenes.Translation;

  std::vector<double> second_base_node_pose = {0,0,0,0,0,0};
  double scaling_factor(1.0);

  std::unique_ptr<SceneAligner> smap_aligner =
      std::unique_ptr<SceneAligner>(new SceneAligner);

  smap_aligner->computeTransformAndDestinationSeparators(
      destination_sfm_data, two_scenes.sfm_data_A, two_scenes.sfm_data_B,
      second_base_node_pose, scaling_factor, two_scenes.common_track_ids);

  Vec3 base_node_R = Vec3(second_base_node_pose[0],second_base_node_pose[1],second_base_node_pose[2]);
  Vec3 base_node_t = Vec3(second_base_node_pose[3],second_base_node_pose[4],second_base_node_pose[5]);

  printBaseNodeCoords(Translation_gt, base_node_R, base_node_t, Rotation_gt, scaling_factor, scaling_factor_gt);

  // test if the bundle adjustment did a good job on evaluating the transformation
  Vec3 inverse_Rotation_gt = -Rotation_gt;
  Vec3 inverse_translation_gt = -Translation_gt;
  const double modulo_rotation = fmod(base_node_R.norm(), M_PI);
  EXPECT_NEAR(modulo_rotation, Rotation_gt.norm(), 1e-2);
  EXPECT_MATRIX_NEAR(base_node_t, inverse_translation_gt, 1e-3);
  EXPECT_MATRIX_NEAR(base_node_R, inverse_Rotation_gt, 1e-3);
  EXPECT_NEAR(scaling_factor_gt, scaling_factor, 1e-3);
}

// checks that distances are conserved (up to a scale) when
// transforming a sfm data scene
TEST(transformSfMDataScene, genericScene_transformationLeavesDistancesConserved)
{
  const openMVG::Mat3 rotation = randomRotationMatrix();
  const openMVG::Vec3 translation = randomVector();
  const double scaling_factor = randomPositiveDouble();

  const double n_poses(20), n_landmarks(1000);
  const openMVG::sfm::SfM_Data original_sfm_data = generate_random_poses_and_landmarks_in_scene(n_poses, n_landmarks);
  openMVG::sfm::SfM_Data destination_sfm_data;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> distances_before = computeDistances(original_sfm_data);

  openMVG::sfm::transformSfMDataScene(destination_sfm_data, original_sfm_data,
      rotation, translation, scaling_factor);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> distances_after = computeDistances(destination_sfm_data);

  std::cout << scaling_factor << std::endl;
  distances_after *= scaling_factor;
  EXPECT_MATRIX_NEAR(distances_before, distances_after, 1e-7);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
