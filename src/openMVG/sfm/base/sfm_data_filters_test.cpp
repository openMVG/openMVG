// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"

#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;


void init_scene
(
  SfM_Data &sfm_data,
  const IndexT viewsCount
)
{
  for (IndexT i = 0; i < viewsCount; ++i)
  {
    // Add views
    std::ostringstream os;
    os << "dataset/" << i << ".jpg";
    const IndexT id_view = i, id_pose = i;
    const IndexT id_intrinsic = 0; //(shared or not intrinsics)
    sfm_data.views[id_view] = std::make_shared<View>(os.str(),id_view, id_intrinsic, id_pose, 1000, 1000);

    // Add poses
    sfm_data.poses[i] = Pose3();

    // Add intrinsics
    if (sfm_data.intrinsics.count(0) == 0)
      sfm_data.intrinsics[0] = std::make_shared<Pinhole_Intrinsic>();
  }
 }

TEST(SFM_DATA_FILTERS, LargestViewTrackCC)
{
  // Init a scene with 5 Views & poses
  SfM_Data sfm_data;
  init_scene(sfm_data, 5);

  // Fill with some tracks
  //
  //- View 0,1,2 are sharing 3 tracks
  // TrackId 0 -> 0-1-2
  // TrackId 1 -> 0-1-2
  // TrackId 2 -> 0-1-2
  //
  // While view 3-4 are sharing 1 track
  // TrackId 3 -> 3-4

  Observations obs;
  obs[0] = Observation( Vec2(10,20), 0);
  obs[1] = Observation( Vec2(30,10), 1);
  obs[2] = Observation( Vec2(30,10), 1);
  sfm_data.structure[0].obs = obs;
  sfm_data.structure[0].X = Vec3::Random();

  sfm_data.structure[1].obs = obs;
  sfm_data.structure[1].X = Vec3::Random();
  sfm_data.structure[2].obs = obs;
  sfm_data.structure[2].X = Vec3::Random();

  obs.clear();
  obs[3] = Observation( Vec2(10,20), 0);
  obs[4] = Observation( Vec2(10,20), 0);
  sfm_data.structure[3].obs = obs;
  sfm_data.structure[3].X = Vec3::Random();

  // Track of the SfM_Data scene contains two connected component
  // One for view (0,1,2)
  // One for the view (3,4)
  // Test that function to keep the largest one is working
  EXPECT_EQ (4, sfm_data.GetLandmarks().size());
  EXPECT_FALSE (IsTracksOneCC(sfm_data));
  KeepLargestViewCCTracks(sfm_data);
  EXPECT_EQ (3, sfm_data.GetLandmarks().size());
}

TEST(SFM_DATA_FILTERS, eraseMissingPoses)
{
  // Init a scene with 6 Views & poses
  SfM_Data sfm_data;
  init_scene(sfm_data, 6);

  // Define observation for views:{0,1,2,4,5}
  for (unsigned char i = 0; i < 6; ++i)
  {
    Observations obs;
    obs[i] = Observation( Vec2(10,20), 0);
    sfm_data.structure[i].obs = obs;
    sfm_data.structure[i].X = Vec3::Random();
  }
  // Since all poses are seen at least one time, there is no need to remove poses
  EXPECT_FALSE(eraseMissingPoses(sfm_data, 1));

  // Erase one observation
  sfm_data.structure.erase(5);
  // Since there is no observation for the view 5, pose 5 must be deleted
  EXPECT_TRUE(eraseMissingPoses(sfm_data, 1));
}

TEST(SFM_DATA_FILTERS, eraseObservationsWithMissingPoses)
{
  // Init a scene with 6 Views & poses
  SfM_Data sfm_data;
  init_scene(sfm_data, 6);

  // Define observation for views:{0,1,2,4,5}
  for (unsigned char i = 0; i < 6; ++i)
  {
    Observations obs;
    obs[i] = Observation( Vec2(10,20), 0);
    sfm_data.structure[i].obs = obs;
    sfm_data.structure[i].X = Vec3::Random();
  }
  // Since all poses are seen at least one time, there is no need to remove poses
  EXPECT_FALSE(eraseObservationsWithMissingPoses(sfm_data, 1));

  // Erase one pose
  sfm_data.poses.erase(5);
  // Since there is no pose for observation of the view 5,
  //  all observations that belongs to view 5 must be removed
  EXPECT_EQ(6, sfm_data.structure.size());
  EXPECT_TRUE(eraseObservationsWithMissingPoses(sfm_data, 1));
  EXPECT_EQ(5, sfm_data.structure.size());
}

TEST(SFM_DATA_FILTERS, eraseUnstablePosesAndObservations)
{
  // Init a scene with 6 Views & poses
  SfM_Data sfm_data;
  init_scene(sfm_data, 6);

  // Define observation for views:{0,1,2,4,5}
  for (unsigned char i = 0; i < 6; ++i)
  {
    Observations obs;
    obs[i] = Observation( Vec2(10,20), 0);
    sfm_data.structure[i].obs = obs;
    sfm_data.structure[i].X = Vec3::Random();
  }
  // Since all poses are seen at least one time, there is no need to remove poses
  EXPECT_FALSE(eraseUnstablePosesAndObservations(sfm_data, 1, 1));

  // Erase one pose & one observation
  sfm_data.poses.erase(5);
  sfm_data.structure.erase(4);

  // Pose id 4 & obversation id 5 must be removed
  EXPECT_FALSE(eraseUnstablePosesAndObservations(sfm_data, 1, 1));
  EXPECT_EQ(4, sfm_data.poses.size());
  EXPECT_EQ(4, sfm_data.structure.size());
  EXPECT_EQ(0, sfm_data.poses.count(4));
  EXPECT_EQ(0, sfm_data.structure.count(5));
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
