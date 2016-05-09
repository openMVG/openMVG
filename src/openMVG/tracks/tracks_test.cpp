
// Copyright (c) 2012, 2013 Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include "openMVG/tracks/tracks.hpp"
#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::tracks;
using namespace openMVG::matching;

#include <vector>
#include <utility>

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

  const IndMatch testAB[] = {IndMatch(0,0), IndMatch(1,1), IndMatch(2,3)};
  const IndMatch testBC[] = {IndMatch(0,0), IndMatch(1,6)};

  const std::vector<IndMatch> ab(testAB, testAB+3);
  const std::vector<IndMatch> bc(testBC, testBC+2);
  const int A = 0;
  const int B = 1;
  const int C = 2;
  map_pairwisematches[ std::make_pair(A,B) ] = ab;
  map_pairwisematches[ std::make_pair(B,C) ] = bc;

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
  const std::pair<uint32_t,uint32_t> GT_Tracks[] =
  {
    std::make_pair(0,0), std::make_pair(1,0), std::make_pair(2,0),
    std::make_pair(0,1), std::make_pair(1,1), std::make_pair(2,6),
    std::make_pair(0,2), std::make_pair(1,3)
  };

  CHECK_EQUAL(3,  map_tracks.size());
  size_t cpt = 0, i = 0;
  for (STLMAPTracks::const_iterator iterT = map_tracks.begin();
    iterT != map_tracks.end();
    ++iterT, ++i)
  {
    CHECK_EQUAL(i, iterT->first);
    for (submapTrack::const_iterator iter = iterT->second.begin();
      iter != iterT->second.end();
      ++iter)
    {
      CHECK( GT_Tracks[cpt] == std::make_pair(iter->first, iter->second));
      ++cpt;
    }
  }
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

  IndMatch testAB[] = {IndMatch(0,0), IndMatch(1,1), IndMatch(2,3)};
  IndMatch testBC[] = {IndMatch(0,0), IndMatch(1,6)};


  std::vector<IndMatch> ab(testAB, testAB+3);
  std::vector<IndMatch> bc(testBC, testBC+2);
  const int A = 0;
  const int B = 1;
  const int C = 2;
  map_pairwisematches[ std::make_pair(A,B) ] = ab;
  map_pairwisematches[ std::make_pair(B,C) ] = bc;

  //-- Build tracks using the interface tracksbuilder
  TracksBuilder trackBuilder;
  trackBuilder.Build( map_pairwisematches );
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
  map_pairwisematches[ std::make_pair(A,B) ] = ab;
  map_pairwisematches[ std::make_pair(B,C) ] = bc;

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
  const std::pair<uint32_t,uint32_t> GT_Tracks[] =
    {std::make_pair(0,0), std::make_pair(1,0), std::make_pair(2,0),
     std::make_pair(0,1), std::make_pair(1,1), std::make_pair(2,6)};

  CHECK_EQUAL(2,  map_tracks.size());
  size_t cpt = 0, i = 0;
  for (STLMAPTracks::const_iterator iterT = map_tracks.begin();
    iterT != map_tracks.end();
    ++iterT, ++i)
  {
    CHECK_EQUAL(i, iterT->first);
    for (submapTrack::const_iterator iter = iterT->second.begin();
      iter != iterT->second.end();
      ++iter)
    {
      CHECK( GT_Tracks[cpt] == std::make_pair(iter->first, iter->second));
      ++cpt;
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
