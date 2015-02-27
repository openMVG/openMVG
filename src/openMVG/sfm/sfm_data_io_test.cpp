
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "testing/testing.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <sstream>

using namespace openMVG;

// Create a SfM scene with desired count of views & poses & intrinsic (shared or not)
// Add a 3D point with observation in 2 view (just in order to have non empty data)
SfM_Data create_test_scene(IndexT viewsCount, bool bSharedIntrinsic)
{
  SfM_Data sfm_data;
  sfm_data.s_root_path = "./";

  for(IndexT i = 0; i < viewsCount; ++i)
  {
    // Add views
    View view;
    view.id_view = i;
    view.id_intrinsic = bSharedIntrinsic ? 0 : i;
    view.id_pose = i;
    std::ostringstream os;
    os << "dataset/" << i << ".jpg";
    view.s_Img_path = os.str();
    sfm_data.views[view.id_view] = view;

    // Add poses
    sfm_data.poses[i] = Pose3();

    // Add intrinsics
    if (bSharedIntrinsic)
    {
      if (i == 0)
        sfm_data.intrinsics[0] = std::make_shared<Intrinsic>();
    }
    else
    {
      sfm_data.intrinsics[i] = std::make_shared<Intrinsic>();
    }
  }

  // Fill with not meaningful tracks
  Observations obs;
  obs[0] = Observation( Vec2(10,20), 0);
  obs[1] = Observation( Vec2(30,10), 1);
  sfm_data.structure[0].obs = obs;
  sfm_data.structure[0].X = Vec3(11,22,33);
  return sfm_data;
}

TEST(SfM_Data_IO, SAVE_LOAD_JSON) {

  const std::vector<std::string> ext_Type = {"json", "bin", "xml"};

  for (int i=0; i < ext_Type.size(); ++i)
  {
    std::ostringstream os;
    os << "SAVE_LOAD" << "." << ext_Type[i];
    const std::string filename = os.str();
    std::cout << "Testing:" << filename << std::endl;

  // SAVE
  {
    const SfM_Data sfm_data = create_test_scene(2, true);
    EXPECT_TRUE( Save(sfm_data, filename, ALL) );
  }

  // LOAD
  {
    const SfM_Data sfm_data = create_test_scene(2, true);
    EXPECT_TRUE( Save(sfm_data, filename, ALL) );
    SfM_Data sfm_data_load;
    ESfM_Data flags_part = ALL;
    EXPECT_TRUE( Load(sfm_data_load, filename, flags_part) );
    EXPECT_EQ( sfm_data_load.views.size(), sfm_data.views.size());
    EXPECT_EQ( sfm_data_load.poses.size(), sfm_data.poses.size());
    EXPECT_EQ( sfm_data_load.intrinsics.size(), sfm_data.intrinsics.size());
    EXPECT_EQ( sfm_data_load.structure.size(), sfm_data.structure.size());
  }

  // LOAD (only a subpart: VIEWS)
  {
    const SfM_Data sfm_data = create_test_scene(2, true);
    EXPECT_TRUE( Save(sfm_data, filename, ALL) );
    SfM_Data sfm_data_load;
    ESfM_Data flags_part = VIEWS;
    EXPECT_TRUE( Load(sfm_data_load, filename, flags_part) );
    EXPECT_EQ( sfm_data_load.views.size(), sfm_data.views.size());
    EXPECT_EQ( sfm_data_load.poses.size(), 0);
    EXPECT_EQ( sfm_data_load.intrinsics.size(), 0);
    EXPECT_EQ( sfm_data_load.structure.size(), 0);
  }

  // LOAD (only a subpart: POSES)
  {
    const SfM_Data sfm_data = create_test_scene(2, true);
    EXPECT_TRUE( Save(sfm_data, filename, ALL) );
    SfM_Data sfm_data_load;
    ESfM_Data flags_part = EXTRINSICS;
    EXPECT_TRUE( Load(sfm_data_load, filename, flags_part) );
    EXPECT_EQ( sfm_data_load.views.size(), 0);
    EXPECT_EQ( sfm_data_load.poses.size(), sfm_data.poses.size());
    EXPECT_EQ( sfm_data_load.intrinsics.size(), 0);
    EXPECT_EQ( sfm_data_load.structure.size(), 0);
  }

  // LOAD (only a subpart: INTRINSICS)
  {
    const SfM_Data sfm_data = create_test_scene(2, true);
    EXPECT_TRUE( Save(sfm_data, filename, ALL) );
    SfM_Data sfm_data_load;
    ESfM_Data flags_part = INTRINSICS;
    EXPECT_TRUE( Load(sfm_data_load, filename, flags_part) );
    EXPECT_EQ( sfm_data_load.views.size(), 0);
    EXPECT_EQ( sfm_data_load.poses.size(), 0);
    EXPECT_EQ( sfm_data_load.intrinsics.size(), sfm_data.intrinsics.size());
    EXPECT_EQ( sfm_data_load.structure.size(), 0);
  }

  // LOAD (subparts: COMBINED)
  {
    const SfM_Data sfm_data = create_test_scene(2,false); //2 intrinsics group here
    EXPECT_TRUE( Save(sfm_data, filename, ALL) );
    SfM_Data sfm_data_load;
    ESfM_Data flags_part = ESfM_Data(INTRINSICS | EXTRINSICS);
    EXPECT_TRUE( Load(sfm_data_load, filename, flags_part) );
    EXPECT_EQ( sfm_data_load.views.size(), 0);
    EXPECT_EQ( sfm_data_load.poses.size(), sfm_data.poses.size());
    EXPECT_EQ( sfm_data_load.intrinsics.size(), sfm_data.intrinsics.size());
    EXPECT_EQ( sfm_data_load.structure.size(), 0);
  }

  // LOAD (subparts: COMBINED)
  {
    const SfM_Data sfm_data = create_test_scene(2, true);
    EXPECT_TRUE( Save(sfm_data, filename, ALL) );
    SfM_Data sfm_data_load;
    ESfM_Data flags_part = ESfM_Data(VIEWS | INTRINSICS | EXTRINSICS);
    EXPECT_TRUE( Load(sfm_data_load, filename, flags_part) );
    EXPECT_EQ( sfm_data_load.views.size(), sfm_data.views.size());
    EXPECT_EQ( sfm_data_load.poses.size(), sfm_data.poses.size());
    EXPECT_EQ( sfm_data_load.intrinsics.size(), sfm_data.intrinsics.size());
    EXPECT_EQ( sfm_data_load.structure.size(), 0);
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
