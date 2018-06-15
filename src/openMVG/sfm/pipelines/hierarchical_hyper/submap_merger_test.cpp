// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/pipelines_test.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_merger.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_test_utilities.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/geometry/Similarity3.hpp"

#include "testing/testing.h"

using namespace openMVG::sfm;
using namespace openMVG::cameras;

void multiplyFocalLength(std::shared_ptr<IntrinsicBase> & modified_intrinsic, double corruption_factor)
{
  Pinhole_Intrinsic_Radial_K3 * cast_intrinsic = dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(modified_intrinsic.get());
  modified_intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>(
        cast_intrinsic->w(), cast_intrinsic->h(),
        corruption_factor*cast_intrinsic->focal(), // "corrupt" focal length
        cast_intrinsic->principal_point()[0],
        cast_intrinsic->principal_point()[1],
        cast_intrinsic->getParams()[0],
        cast_intrinsic->getParams()[1],
        cast_intrinsic->getParams()[2]
      );
}

std::set<openMVG::IndexT> extractTrackIds(const openMVG::sfm::SfM_Data & sfm_data)
{
  std::set<openMVG::IndexT> track_ids;
  std::transform(sfm_data.GetLandmarks().cbegin(),
                 sfm_data.GetLandmarks().cend(),
                 std::inserter(track_ids, track_ids.begin()),
                 [](const std::pair<openMVG::IndexT, Landmark> & l){return l.first;});
  return track_ids;
}

// This function creates a parent and two children submaps from a single scene.
// The parent submap and the first submap are kept in the same frame of reference
// as the input scene, while the second submap is transformed to another arbitrary frame of reference
HsfmSubmaps createSubmapsFromInitialScene(const SfM_Data sfm_data_common)
{
  HsfmSubmap parent_submap;
  parent_submap.separator = extractTrackIds(sfm_data_common);
  parent_submap.track_ids = extractTrackIds(sfm_data_common);
  parent_submap.is_parent = true;
  parent_submap.children_submaps = {1,2};
  parent_submap.sfm_data = sfm_data_common;
  // a parent submap does not contain poses nor structure
  parent_submap.sfm_data.poses.clear();
  parent_submap.sfm_data.structure.clear();

  SfM_Data sfm_data_A = sfm_data_common;
  SfM_Data sfm_data_B = sfm_data_common;
  // put half of the views/poses on one scene, the other half on the other scene
  const int n_total_cams = sfm_data_common.poses.size();
  for (int i(0); i < n_total_cams; i += 2)
  {
    sfm_data_A.views.erase(i);
    sfm_data_B.views.erase(i + 1);
    sfm_data_A.poses.erase(i);
    sfm_data_B.poses.erase(i + 1);
  }

  // make some of the landmarks exclusive to A or B
  // we arbitrarily choose the first 80% the landmarks to
  // exclusive to one of the submaps and the last 20% to be separators (common)
  const int n_total_landmarks = sfm_data_common.structure.size();
  for (int i(0); i < 4*n_total_landmarks/5; i +=2)
  {
    sfm_data_A.structure.erase(i);
    sfm_data_B.structure.erase(i + 1);
    parent_submap.separator.erase(i);
    parent_submap.separator.erase(i+1);
  }

  // submap B is not in the same referential as A
  const geometry::Similarity3 arbitrary_transformation = easySimilarity();
  ApplySimilarity(arbitrary_transformation, sfm_data_B);

  HsfmSubmap submap_A;
  submap_A.track_ids = extractTrackIds(sfm_data_A);
  submap_A.sfm_data = sfm_data_A;
  HsfmSubmap submap_B;
  submap_B.track_ids = extractTrackIds(sfm_data_B);
  submap_B.sfm_data = sfm_data_B;

  const HsfmSubmaps submaps =
  {
    {0, parent_submap},
    {1, submap_A},
    {2, submap_B}
  };

  return submaps;
}

TEST(SubmapMerger, MergingTwoSubmaps_DifferentIntrinsics_GetMerged)
{
  // In this test, we check that the merging process merges the intrinsics
  // as well

  const int nviews = 20;
  const int npoints = 500;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data_common = getInputScene(d, config, PINHOLE_CAMERA_RADIAL3);

  HsfmSubmaps submaps = createSubmapsFromInitialScene(sfm_data_common);
  SfM_Data & second_smap_scene = submaps.at(2).sfm_data;
  SfM_Data & parent_smap_scene = submaps.at(0).sfm_data;

  // make the second submap intrinsic a different intrinsic from the first
  // submap.
  for (auto & view : second_smap_scene.GetViews())
  {
    const openMVG::IndexT view_id = view.first;
    view.second->id_intrinsic = 1;
    parent_smap_scene.views.at(view_id)->id_intrinsic = 1;
  }
  // parent submap contains both intrinsics
  parent_smap_scene.intrinsics[1] = parent_smap_scene.intrinsics[0];

  // the second submap contains only the second intrinsic
  second_smap_scene.intrinsics[1] = second_smap_scene.intrinsics[0];
  second_smap_scene.intrinsics.erase(0);

  // modify the actual values of the intrinsic. This simulates the fact
  // that the intrinsics in the submap might get optimized by the bundle adjustment
  // after clustering
  multiplyFocalLength(second_smap_scene.intrinsics.at(1), 1.5);
  assert(second_smap_scene.intrinsics.size() == 1);
  assert(submaps.at(0).sfm_data.intrinsics.size() == 2);
  assert(submaps.at(1).sfm_data.intrinsics.size() == 1);

  // we don't optimize intrinsics here, we just want to test if they merge well.
  const auto intrinsics_refinement_param = openMVG::cameras::Intrinsic_Parameter_Type::NONE;

  SubmapMerger smap_merger(submaps, intrinsics_refinement_param);
  smap_merger.Merge();
  const HsfmSubmaps merged_submaps = smap_merger.getSubmaps();

  // check merged intrinsics size 2, submaps size 1, etc..
  EXPECT_EQ(merged_submaps.at(0).sfm_data.intrinsics.size(), 2);
  EXPECT_EQ(merged_submaps.at(1).sfm_data.intrinsics.size(), 1);
  EXPECT_EQ(merged_submaps.at(2).sfm_data.intrinsics.size(), 1);

  // check that the intrinsics were correctly propagated back to the
  // parent submap (by checking value of the focal length)
  EXPECT_EQ(dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(merged_submaps.at(0).sfm_data.intrinsics.at(0).get())->focal(),
            dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(merged_submaps.at(1).sfm_data.intrinsics.at(0).get())->focal());
  EXPECT_EQ(dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(merged_submaps.at(0).sfm_data.intrinsics.at(1).get())->focal(),
            dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(merged_submaps.at(2).sfm_data.intrinsics.at(1).get())->focal());
  EXPECT_TRUE(dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(merged_submaps.at(0).sfm_data.intrinsics.at(1).get())->focal() !=
              dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(submaps.at(0).sfm_data.intrinsics.at(1).get())->focal());

  // check that the intrinsics pointers are pointing to different objects (we do a deep copy)
  EXPECT_TRUE(merged_submaps.at(0).sfm_data.intrinsics.at(0)
              != merged_submaps.at(1).sfm_data.intrinsics.at(0));
  EXPECT_TRUE(merged_submaps.at(0).sfm_data.intrinsics.at(1)
              != merged_submaps.at(2).sfm_data.intrinsics.at(1));
}

TEST(SubmapMerger, MergingTwoSubmaps_DivergentIntrinsics_ChoosesTheBestOne)
{
  // In this test, we check that the merging process selects the best intrinsics
  // between the cameras common to both sibling submaps (intrinsics can diverge when being
  // optimized separately in each submap)

  const int nviews = 20;
  const int npoints = 500;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data_common = getInputScene(d, config, PINHOLE_CAMERA_RADIAL3);

  HsfmSubmaps submaps = createSubmapsFromInitialScene(sfm_data_common);

  // alter the intrinsics in the first submap
  HsfmSubmap & smap_child_corrupt = submaps.at(1);
  HsfmSubmap & smap_child_clean = submaps.at(2);
  auto & corrupted_intrinsic = smap_child_corrupt.sfm_data.intrinsics.at(0);
  multiplyFocalLength(corrupted_intrinsic, 1.5);

  // we don't optimize intrinsics here, we just want to test if they merge well.
  const auto intrinsics_refinement_param = openMVG::cameras::Intrinsic_Parameter_Type::NONE;

  SubmapMerger smap_merger(submaps, intrinsics_refinement_param);
  smap_merger.Merge();
  const HsfmSubmaps merged_submaps = smap_merger.getSubmaps();

  // the merged intrinsics should be clean
  EXPECT_TRUE(dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(merged_submaps.at(0).sfm_data.intrinsics.at(0).get())->focal() !=
              dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(smap_child_corrupt.sfm_data.intrinsics.at(0).get())->focal());

  EXPECT_TRUE(dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(merged_submaps.at(0).sfm_data.intrinsics.at(0).get())->focal() ==
              dynamic_cast<Pinhole_Intrinsic_Radial_K3*>(smap_child_clean.sfm_data.intrinsics.at(0).get())->focal());

}

TEST(SubmapMerger, MergingTwoSubmaps_GivesComparableResultsToSingleReconstruction)
{
  // In this test, we check that merging two submaps will give a
  // similar result than a full submap containing everything.
  // i.e. -> we check that the merging works

  // We can do this because by default, SubmapMerger will put everything in
  // the coordinate referential of the first submap it is given. Comparing
  // with a full submap is therefore quite easy since there are no additional
  // rigid transformation to apply to the scene to compare it.

  const int nviews = 20;
  const int npoints = 500;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data_common = getInputScene(d, config, PINHOLE_CAMERA);

  const HsfmSubmaps submaps = createSubmapsFromInitialScene(sfm_data_common);

  SubmapMerger smap_merger(submaps, openMVG::cameras::Intrinsic_Parameter_Type::ADJUST_ALL);
  smap_merger.Merge();
  const HsfmSubmaps merged_submaps = smap_merger.getSubmaps();

  // check merged scene is comparable to the initial common scene
  for (const auto & mapped_pose : merged_submaps.at(0).sfm_data.poses)
  {
    const IndexT pose_id = mapped_pose.first;
    const geometry::Pose3 & merged_pose = mapped_pose.second;
    const geometry::Pose3 & gt_pose = sfm_data_common.poses.at(pose_id);

    EXPECT_MATRIX_NEAR(merged_pose.center(), gt_pose.center(), 1e-10);
    EXPECT_MATRIX_NEAR(merged_pose.rotation(), gt_pose.rotation(), 1e-10);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
