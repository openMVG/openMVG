#include "testing/testing.h"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "submap_utilities.hpp"
#include "submap_test_utilities.hpp"

using namespace openMVG::sfm;
using namespace openMVG::cameras;
using namespace openMVG;

// Create a SfM scene with desired count of views & poses & intrinsic (shared or not)
// Add a 3D point with observation in 2 view (just in order to have non empty data)
// (Copied from sfm_data_io_test)
SfM_Data create_test_scene(IndexT viewsCount, bool bSharedIntrinsic)
{
  SfM_Data sfm_data;
  sfm_data.s_root_path = "./";

  for (IndexT i = 0; i < viewsCount; ++i)
  {
    // Add views
    std::ostringstream os;
    os << "dataset/" << i << ".jpg";
    const IndexT id_view = i, id_pose = i;
    const IndexT id_intrinsic = bSharedIntrinsic ? 0 : i; //(shared or not intrinsics)
    sfm_data.views[id_view] = std::make_shared<View>(os.str(),id_view, id_intrinsic, id_pose, 1000, 1000);

    // Add poses
    sfm_data.poses[i] = Pose3();

    // Add intrinsics
    if (bSharedIntrinsic)
    {
      if (i == 0)
        sfm_data.intrinsics[0] = std::make_shared<Pinhole_Intrinsic>();
    }
    else
    {
      sfm_data.intrinsics[i] = std::make_shared<Pinhole_Intrinsic>();
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

TEST(subam_IO, SaveLoad_AllFormats_Works)
{
  const std::vector<std::string> extTypes = {"json", "bin", "xml"};
  const int nViews = 20;
  for (const auto & extType : extTypes)
  {
    const std::string filename = "SAVE_LOAD." + extType;
    HsfmSubmap submap;
    submap.sfm_data = create_test_scene(nViews, false);
    submap.is_parent = true;
    submap.children_submaps = {1,2};
    submap.parent_id = 0;
    submap.separator = {0,1,2};
    submap.track_ids = {0,1,2,3,4,5};

    EXPECT_TRUE( Save(submap, filename));

    HsfmSubmap submap_load;
    EXPECT_TRUE( Load(submap_load, filename));
    EXPECT_EQ(submap.children_submaps.first, submap_load.children_submaps.first);
    EXPECT_EQ(submap.children_submaps.second, submap_load.children_submaps.second);
    EXPECT_EQ(submap.parent_id, submap_load.parent_id);
    EXPECT_EQ(submap.separator.size(), submap_load.separator.size());
    EXPECT_EQ(submap.track_ids.size(), submap_load.track_ids.size());
    EXPECT_EQ(submap.is_parent, submap_load.is_parent);
    EXPECT_EQ(submap.sfm_data.GetViews().size(), submap_load.sfm_data.GetViews().size());
    EXPECT_EQ(submap.sfm_data.GetIntrinsics().size(), submap_load.sfm_data.GetIntrinsics().size());
    EXPECT_EQ(submap.sfm_data.GetPoses().size(), submap_load.sfm_data.GetPoses().size());
    EXPECT_EQ(submap.sfm_data.GetLandmarks().size(), submap_load.sfm_data.GetLandmarks().size());
    EXPECT_EQ(submap.sfm_data.GetControl_Points().size(), submap_load.sfm_data.GetControl_Points().size());
  }
}

TEST(full_submaps_IO, SaveLoad_AllFormats_Works)
{
  const std::vector<std::string> extTypes = {"json", "bin", "xml"};
  const int nPoses = 20;

  TwoScenesConfig two_scenes = createTwoScenesWithTransformationAndNoise(nPoses);

  for (const auto & extType : extTypes)
  {
    const std::string filename = "SAVE_LOAD." + extType;
    HsfmSubmaps submaps;
    submaps[0].is_parent = true;
    submaps[1].is_parent = false;
    submaps[2].is_parent = false;

    submaps[0].sfm_data = two_scenes.parent_sfm_data;
    submaps[1].sfm_data = two_scenes.sfm_data_A;
    submaps[2].sfm_data = two_scenes.sfm_data_B;
    submaps[0].children_submaps = {1,2};
    submaps[0].parent_id = 0;
    submaps[1].parent_id = 0;
    submaps[2].parent_id = 0;
    submaps[0].separator = {4,5};
    submaps[0].track_ids = {1,2,3,4,5,6,7,8,9,10};
    submaps[1].track_ids = {1,2,3,4,5};
    submaps[2].track_ids = {4,5,6,7,8,9,10};


    EXPECT_TRUE(Save(submaps, filename));

    HsfmSubmaps submaps_load;
    EXPECT_TRUE(Load(submaps_load, filename));
    EXPECT_EQ(submaps.size(), submaps_load.size());
    for (const auto & smap : submaps)
    {
      const openMVG::IndexT submapId = smap.first;
      const HsfmSubmap & submap = smap.second;
      const HsfmSubmap & submap_load = submaps_load.at(submapId);

      EXPECT_EQ(submap.children_submaps.first, submap_load.children_submaps.first);
      EXPECT_EQ(submap.children_submaps.second, submap_load.children_submaps.second);
      EXPECT_EQ(submap.parent_id, submap_load.parent_id);
      EXPECT_EQ(submap.separator.size(), submap_load.separator.size());
      EXPECT_EQ(submap.track_ids.size(), submap_load.track_ids.size());
      EXPECT_EQ(submap.is_parent, submap_load.is_parent);
      EXPECT_EQ(submap.sfm_data.GetViews().size(), submap_load.sfm_data.GetViews().size());
      EXPECT_EQ(submap.sfm_data.GetIntrinsics().size(), submap_load.sfm_data.GetIntrinsics().size());
      EXPECT_EQ(submap.sfm_data.GetPoses().size(), submap_load.sfm_data.GetPoses().size());
      EXPECT_EQ(submap.sfm_data.GetLandmarks().size(), submap_load.sfm_data.GetLandmarks().size());
      EXPECT_EQ(submap.sfm_data.GetControl_Points().size(), submap_load.sfm_data.GetControl_Points().size());
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
