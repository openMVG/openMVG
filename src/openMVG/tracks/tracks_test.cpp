// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/tracks/tracks.hpp"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <vector>
#include <utility>

using namespace openMVG::tracks;
using namespace openMVG::matching;


TEST(Tracks, Simple) {

  // Create some tracks for image (A,B,C)
  // {A,B,C} imageId will be {0,1,2}
  // For those image link some features id depicted below
  //A    B    C
  //0 -> 0 -> 0
  //1 -> 1 -> 6
  //2 -> 3

  // Create the input pairwise correspondences
  PairWiseMatches map_pairwisematches;

  const std::vector<IndMatch> ab = {IndMatch(0,0), IndMatch(1,1), IndMatch(2,3)};
  const std::vector<IndMatch> bc = {IndMatch(0,0), IndMatch(1,6)};
  const int A = 0;
  const int B = 1;
  const int C = 2;
  map_pairwisematches[ {A,B} ] = ab;
  map_pairwisematches[ {B,C} ] = bc;

  //-- Build tracks using the interface tracksbuilder
  TracksBuilder trackBuilder;
  trackBuilder.Build( map_pairwisematches );

  STLMAPTracks map_tracks;
  trackBuilder.ExportToSTL(map_tracks);

  //-------------------
  // Unit Test check
  //-------------------

  //0, {(0,0) (1,0) (2,0)}
  //1, {(0,1) (1,1) (2,6)}
  //2, {(0,2) (1,3)}
  const STLMAPTracks GT_Tracks =
  {
    {0, {{0,0}, {1,0}, {2,0}}},
    {1, {{0,1}, {1,1}, {2,6}}},
    {2, {{0,2}, {1,3}}},
  };
  // Check that computed tracks are the desired one
  CHECK(GT_Tracks == map_tracks);
}

TEST(Tracks, filter_3viewAtLeast) {

  //
  //A    B    C
  //0 -> 0 -> 0
  //1 -> 1 -> 6
  //2 -> 3
  //

  // Create the input pairwise correspondences
  PairWiseMatches map_pairwisematches;

  std::vector<IndMatch> ab = {IndMatch(0,0), IndMatch(1,1), IndMatch(2,3)};
  std::vector<IndMatch> bc = {IndMatch(0,0), IndMatch(1,6)};
  const int A = 0;
  const int B = 1;
  const int C = 2;
  map_pairwisematches[ {A,B} ] = ab;
  map_pairwisematches[ {B,C} ] = bc;

  //-- Build tracks using the interface tracksbuilder
  TracksBuilder trackBuilder;
  trackBuilder.Build(map_pairwisematches);
  CHECK_EQUAL(3, trackBuilder.NbTracks());
  trackBuilder.Filter(3);
  CHECK_EQUAL(2, trackBuilder.NbTracks());
}

TEST(Tracks, Conflict) {

  //
  //A    B    C
  //0 -> 0 -> 0
  //1 -> 1 -> 6
  //{2 -> 3 -> 2
  //      3 -> 8 } This track must be deleted, index 3 appears two times
  //

  // Create the input pairwise correspondences
  PairWiseMatches map_pairwisematches;

  const IndMatch testAB[] = {IndMatch(0,0), IndMatch(1,1), IndMatch(2,3)};
  const IndMatch testBC[] = {IndMatch(0,0), IndMatch(1,6), IndMatch(3,2), IndMatch(3,8)};

  std::vector<IndMatch> ab(testAB, testAB+3);
  std::vector<IndMatch> bc(testBC, testBC+4);
  const int A = 0;
  const int B = 1;
  const int C = 2;
  map_pairwisematches[ {A,B} ] = ab;
  map_pairwisematches[ {B,C} ] = bc;

  //-- Build tracks using the interface tracksbuilder
  TracksBuilder trackBuilder;
  trackBuilder.Build( map_pairwisematches );

  CHECK_EQUAL(3, trackBuilder.NbTracks());
  trackBuilder.Filter(); // Key feature tested here to kill the conflicted track
  CHECK_EQUAL(2, trackBuilder.NbTracks());

  STLMAPTracks map_tracks;
  trackBuilder.ExportToSTL(map_tracks);

  //-------------------
  // Unit Test check
  //-------------------

  //0, {(0,0) (1,0) (2,0)}
  //1, {(0,1) (1,1) (2,6)}
  const STLMAPTracks GT_Tracks =
  {
    {0, {{0,0},{1,0},{2,0}}},
    {1, {{0,1},{1,1},{2,6}}}
  };
  // Check that computed tracks are the desired one
  CHECK(GT_Tracks == map_tracks);
}


TEST(Tracks, TracksInImages) {

  //
  // Test "TracksInImages".
  // Function that allows to retrieve the tracks image observations related to one or many view Ids
  //
  const STLMAPTracks tracks_in =
  {
    {0, {{0,0},{1,1}}}, // Track Id 0: with image observations in view 0 and 1
    {1, {{0,0},{1,1}}}, // Track Id 1: with image observations in view 0 and 1
    {2, {{0,0},{2,2}}}, // Track Id 2: with image observations in view 0 and 2
    {3, {{0,0},{1,1}}}  // Track Id 3: with image observations in view 0 and 1
  };

  {
    STLMAPTracks tracks_out_image0;

    EXPECT_TRUE(TracksUtilsMap::GetTracksInImages({0}, tracks_in, tracks_out_image0));
    EXPECT_EQ(4, tracks_out_image0.size());

    EXPECT_TRUE(TracksUtilsMap::GetTracksInImages({1}, tracks_in, tracks_out_image0));
    EXPECT_EQ(3, tracks_out_image0.size());

    EXPECT_TRUE(TracksUtilsMap::GetTracksInImages({2}, tracks_in, tracks_out_image0));
    EXPECT_EQ(1, tracks_out_image0.size());

    EXPECT_TRUE(TracksUtilsMap::GetTracksInImages({0,1}, tracks_in, tracks_out_image0));
    EXPECT_EQ(3, tracks_out_image0.size());

    EXPECT_TRUE(TracksUtilsMap::GetTracksInImages({0,2}, tracks_in, tracks_out_image0));
    EXPECT_EQ(1, tracks_out_image0.size());

    // Border case (ask tracks for an image id that is not listed in the tracks)
    EXPECT_FALSE(TracksUtilsMap::GetTracksInImages({99}, tracks_in, tracks_out_image0));
    EXPECT_EQ(0, tracks_out_image0.size());
    EXPECT_FALSE(TracksUtilsMap::GetTracksInImages({0,99}, tracks_in, tracks_out_image0));
    EXPECT_EQ(0, tracks_out_image0.size());
  }

  // Test the same behavior but with the class that precompute the track id list per view
  {
    openMVG::tracks::SharedTrackVisibilityHelper shared_track_visibility_helper(tracks_in);

    STLMAPTracks tracks_out_image0;

    EXPECT_TRUE(shared_track_visibility_helper.GetTracksInImages({0}, tracks_out_image0));
    EXPECT_EQ(4, tracks_out_image0.size());

    EXPECT_TRUE(shared_track_visibility_helper.GetTracksInImages({1}, tracks_out_image0));
    EXPECT_EQ(3, tracks_out_image0.size());

    EXPECT_TRUE(shared_track_visibility_helper.GetTracksInImages({2}, tracks_out_image0));
    EXPECT_EQ(1, tracks_out_image0.size());

    EXPECT_TRUE(shared_track_visibility_helper.GetTracksInImages({0,1}, tracks_out_image0));
    EXPECT_EQ(3, tracks_out_image0.size());

    EXPECT_TRUE(shared_track_visibility_helper.GetTracksInImages({0,2}, tracks_out_image0));
    EXPECT_EQ(1, tracks_out_image0.size());

    // Border case (ask tracks for an image id that is not listed in the tracks)
    EXPECT_FALSE(shared_track_visibility_helper.GetTracksInImages({99}, tracks_out_image0));
    EXPECT_EQ(0, tracks_out_image0.size());
    EXPECT_FALSE(shared_track_visibility_helper.GetTracksInImages({0,99}, tracks_out_image0));
    EXPECT_EQ(0, tracks_out_image0.size());
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
