// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/hierarchical_hyper/hypercluster.hpp" // what we're testing

#include "openMVG/test_utility_functions.hpp"

#include "openMVG/sfm/sfm.hpp"
#include "testing/testing.h"

void generate_sfm_data_and_tracks(openMVG::sfm::SfM_Data & sfm_data, openMVG::tracks::STLMAPTracks & map_tracks, const int n_views, const int n_points);

// test constructor
TEST(HYPERCLUSTER, Basis)
{
  openMVG::sfm::SfM_Data sfm_data;
  openMVG::tracks::STLMAPTracks map_tracks;

  const int n_views(4), n_tracks(50);

  // create some input sfm data and tracks
  generate_sfm_data_and_tracks(sfm_data, map_tracks, n_views, n_tracks);

  // test constructor
  std::unique_ptr<openMVG::sfm::SubmapThresholdChecker> threshold_checker(
    new openMVG::sfm::SubmapTracksThresholdChecker(n_tracks));
  openMVG::sfm::HyperCluster hyper_cluster(sfm_data, map_tracks, std::move(threshold_checker));
}

// test partitioning
TEST(HYPERCLUSTER, simple_partitioning_successful)
{
  openMVG::sfm::SfM_Data sfm_data;
  openMVG::tracks::STLMAPTracks map_tracks;

  const int n_views(20), n_tracks(200);

  // create some input sfm data and tracks
  generate_sfm_data_and_tracks(sfm_data, map_tracks, n_views, n_tracks);

  // create clusterer
  std::unique_ptr<openMVG::sfm::SubmapThresholdChecker> threshold_checker(
    new openMVG::sfm::SubmapTracksThresholdChecker(3*n_tracks/4));
  openMVG::sfm::HyperCluster hyper_cluster(sfm_data, map_tracks, std::move(threshold_checker));

  // test recursive partitioning
  EXPECT_TRUE(hyper_cluster.recursivePartitioning());

  // test submaps
  openMVG::sfm::HsfmSubmaps submaps = hyper_cluster.getSubMaps();

  for (const auto & submap : submaps)
  {
    const auto & smap = submap.second;
    if (!smap.is_parent)
      continue;
    const auto fst_child_id = smap.children_submaps.first;
    const auto snd_child_id = smap.children_submaps.second;

    // the sum of the views of the children submaps should be equal to their parent's number of views
    EXPECT_EQ(submaps.at(fst_child_id).sfm_data.GetViews().size()
        + submaps.at(snd_child_id).sfm_data.GetViews().size(),
        smap.sfm_data.GetViews().size());

    // the sum of the tracks of the children submaps should be equal to their parent's number of tracks
    // (we have to subtract once the separator tracks otherwise they are counted twice)
    EXPECT_EQ(submaps.at(fst_child_id).track_ids.size()
        + submaps.at(snd_child_id).track_ids.size() - smap.separator.size(),
        smap.track_ids.size());

    // check that siblings submaps are made of ALL the tracks from the parent submap
    std::set<openMVG::IndexT> union_siblings_track_ids;
    std::set_union(submaps.at(fst_child_id).track_ids.begin(), submaps.at(fst_child_id).track_ids.end(),
        submaps.at(snd_child_id).track_ids.begin(), submaps.at(snd_child_id).track_ids.end(), std::inserter(union_siblings_track_ids, union_siblings_track_ids.begin()));

    const auto & parent_track_ids = smap.track_ids;

    EXPECT_TRUE(union_siblings_track_ids == parent_track_ids);

    // check that new intrinsics are actually created for the clustered submaps
    // this is important because optimizing intrinsics later on in a submap should
    // not influence intriniscs in other submaps.
    const auto parent_intrinsics = smap.sfm_data.GetIntrinsics().at(0);
    const auto fst_intrinsics = submaps.at(fst_child_id).sfm_data.GetIntrinsics().at(0);
    const auto snd_intrinsics = submaps.at(snd_child_id).sfm_data.GetIntrinsics().at(0);
    EXPECT_TRUE(fst_intrinsics != snd_intrinsics);
    EXPECT_TRUE(fst_intrinsics != parent_intrinsics);
    EXPECT_TRUE(snd_intrinsics != parent_intrinsics);
  }
}

void generate_sfm_data_and_tracks(openMVG::sfm::SfM_Data & sfm_data, openMVG::tracks::STLMAPTracks & map_tracks, const int n_views, const int n_points)
{
  // generate sfm_data with views + set of tracks distributed
  // in the views.

  // Views
  std::set<openMVG::IndexT> all_view_ids;
  for (int i(0); i < n_views; i++)
  {
    const openMVG::IndexT id_view = i, id_pose = i, id_intrinsic = 0; //(shared intrinsics)
    sfm_data.views[i] = std::make_shared<openMVG::sfm::View>("", id_view, id_intrinsic, id_pose, 2000, 2000);
    all_view_ids.insert(i);
  }

  sfm_data.intrinsics[0] = std::make_shared<openMVG::cameras::Pinhole_Intrinsic>();

  // tracks
  for (int i(0); i < n_points; i++)
  {
    openMVG::tracks::submapTrack s_track;
    // generate random submaptrack
    const int n_views_in_track = randomPositiveInteger(n_views);
    std::set<openMVG::IndexT> remaining_view_ids = all_view_ids;
    auto it = remaining_view_ids.cbegin();
    for (int j(0); j<n_views_in_track; j++)
    {
      it = remaining_view_ids.cbegin();
      const int next_id = randomPositiveInteger(remaining_view_ids.size() - 1);
      std::advance(it, next_id);
      s_track[*it] = randomPositiveInteger(RAND_MAX);// complete random feature id...it is not used anyway TODO < true ?
      it = remaining_view_ids.erase(it);
    }
    map_tracks[i] = s_track;
  }
}

/* ************************************************************************* */
int main() {
 TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
