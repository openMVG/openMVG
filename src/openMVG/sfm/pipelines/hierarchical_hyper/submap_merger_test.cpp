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

std::set<openMVG::IndexT> extractTrackIds(const openMVG::sfm::SfM_Data & sfm_data)
{
  std::set<openMVG::IndexT> track_ids;
  std::transform(sfm_data.GetLandmarks().cbegin(),
                 sfm_data.GetLandmarks().cend(),
                 std::inserter(track_ids, track_ids.begin()),
                 [](const std::pair<openMVG::IndexT, Landmark> & l){return l.first;});
  return track_ids;
}

TEST(SubmapMerger, MergingTwoSubmaps_GivesComparableResultsToSingleReconstruction)
{
  // here we check that merging two submaps will give a
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
  // exclusive to one of the submaps and the 20% to be separators (common)
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

  SubmapMerger smap_merger(submaps);
  smap_merger.Merge();
  HsfmSubmaps merged_submaps = smap_merger.getSubmaps();

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
