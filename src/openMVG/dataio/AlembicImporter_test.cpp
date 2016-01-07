#include "openMVG/dataio/AlembicImporter.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "testing/testing.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"


using namespace openMVG;
using namespace openMVG::sfm;


// Create a SfM scene with desired count of views & poses & intrinsic (shared or not)
// Add a 3D point with observation in 2 view (just in order to have non empty data)
SfM_Data create_test_scene(IndexT viewsCount, IndexT pointCount, bool bSharedIntrinsic)
{
  SfM_Data sfm_data;
  sfm_data.s_root_path = "./";

  std::srand(time(NULL));

  for(IndexT i = 0; i < viewsCount; ++i)
  {
    // Add views
    std::ostringstream os;
    os << "dataset/" << i << ".jpg";
    const IndexT id_view = i, id_pose = i;
    const IndexT id_intrinsic = bSharedIntrinsic ? 0 : i; //(shared or not intrinsics)
    sfm_data.views[id_view] = std::make_shared<View>(os.str(),id_view, id_intrinsic, id_pose, 1000, 1000);

    // Add poses
    const Mat3 r = Mat3::Random();
    const Vec3 c = Vec3::Random();
    sfm_data.poses[i] = geometry::Pose3(r, c);
    // sfm_data.poses[i] = Pose3();

    // Add intrinsics
    if (bSharedIntrinsic)
    {
      if (i == 0)
        sfm_data.intrinsics[0] = std::make_shared<cameras::Pinhole_Intrinsic>(1000, 1000, 36.0, std::rand()%10000, std::rand()%10000);
    }
    else
    {
      sfm_data.intrinsics[i] = std::make_shared<cameras::Pinhole_Intrinsic>(1000, 1000, 36.0, std::rand()%10000, std::rand()%10000);
    }
  }

  // Fill with not meaningful tracks
  for(IndexT i = 0; i < pointCount; ++i)
  {

    // Add structure
    Observations obs;
    obs[0] = Observation( Vec2(std::rand() % 10000, std::rand() % 10000), 0);
    obs[1] = Observation( Vec2(std::rand() % 10000, std::rand() % 10000), 1);
    sfm_data.structure[i].obs = obs;
    sfm_data.structure[i].X = Vec3(std::rand() % 10000, std::rand() % 10000, std::rand() % 10000);

    // Add control points    
  }

  return sfm_data;
}

//-----------------
// Test summary:
//-----------------
// - Create a random scene (.json)
// - Export to Alembic
// - Import to Alembic
// - Import to .json
//-----------------
TEST(AlembicImporter, importExport) {

    int flags = ALL;

    // Create a random scene
    SfM_Data sfm_data = create_test_scene(5, 50, true);
    

    /*********************************************************************/
    /*****************              JSON -> JSON           ***************/
    /*********************************************************************/

    // Export as JSON
    std::string jsonFile = "jsonToJson.json";
    {
        EXPECT_TRUE(Save(
        sfm_data,
        jsonFile.c_str(),
        ESfM_Data(flags)));
    }

    // Reload
    SfM_Data sfmJsonToJson;
    {
        EXPECT_TRUE(Load(sfmJsonToJson, jsonFile, ESfM_Data(flags)));
        EXPECT_TRUE(sfm_data == sfmJsonToJson);
    }

    // Export as JSON
    std::string jsonFile2 = "jsonToJson2.json";
    {
        EXPECT_TRUE(Save(
        sfmJsonToJson,
        jsonFile2.c_str(),
        ESfM_Data(flags)));
    }
    
    /*********************************************************************/
    /*****************              ABC -> ABC           *****************/
    /*********************************************************************/

    // Export as ABC
    std::string abcFile = "abcToAbc.abc";
    {
        EXPECT_TRUE(Save(
        sfm_data,
        abcFile.c_str(),
        ESfM_Data(flags)));
    }

    // Reload
    SfM_Data sfmAbcToAbc;
    {
        EXPECT_TRUE(Load(sfmAbcToAbc, abcFile, ESfM_Data(flags)));
        std::string abcFile2 = "abcToJson.json";
        EXPECT_TRUE(Save(
        sfmAbcToAbc,
        abcFile2.c_str(),
        ESfM_Data(flags)));
        EXPECT_TRUE(sfm_data == sfmAbcToAbc);
    }

    // Export as ABC
    std::string abcFile2 = "abcToAbc2.abc";
    {
        EXPECT_TRUE(Save(
        sfmAbcToAbc,
        abcFile2.c_str(),
        ESfM_Data(flags)));
    }



    /*********************************************************************/
    /****************      JSON -> ABC -> ABC -> JSON       **************/
    /*********************************************************************/

    // Export as JSON
    jsonFile = "jsonToABC.json";
    {
        EXPECT_TRUE(Save(
        sfm_data,
        jsonFile.c_str(),
        ESfM_Data(flags)));
    }

    // Reload
    SfM_Data sfmJsonToABC;
    {
        EXPECT_TRUE(Load(sfmJsonToABC, jsonFile, ESfM_Data(flags)));
        EXPECT_EQ( sfm_data.views.size(), sfmJsonToABC.views.size());
        EXPECT_EQ( sfm_data.poses.size(), sfmJsonToABC.poses.size());
        EXPECT_EQ( sfm_data.intrinsics.size(), sfmJsonToABC.intrinsics.size());
        EXPECT_EQ( sfm_data.structure.size(), sfmJsonToABC.structure.size());
        EXPECT_EQ( sfm_data.control_points.size(), sfmJsonToABC.control_points.size());
    }

    // Export as ABC
    abcFile = "jsonToABC.abc";
    {
        EXPECT_TRUE(Save(
        sfmJsonToABC,
        abcFile.c_str(),
        ESfM_Data(flags)));
    }

    // Reload
    SfM_Data sfmJsonToABC2;
    {
        EXPECT_TRUE(Load(sfmJsonToABC2, abcFile, ESfM_Data(flags)));
        EXPECT_EQ( sfm_data.views.size(), sfmJsonToABC2.views.size());
        EXPECT_EQ( sfm_data.poses.size(), sfmJsonToABC2.poses.size());
        EXPECT_EQ( sfm_data.intrinsics.size(), sfmJsonToABC2.intrinsics.size());
        EXPECT_EQ( sfm_data.structure.size(), sfmJsonToABC2.structure.size());
        EXPECT_EQ( sfm_data.control_points.size(), sfmJsonToABC2.control_points.size());
    }

    // Export as ABC
    abcFile2 = "jsonToABC2.abc";
    {
        EXPECT_TRUE(Save(
        sfmJsonToABC2,
        abcFile2.c_str(),
        ESfM_Data(flags)));
    }

    // Reload
    SfM_Data sfmJsonToABC3;
    {
        EXPECT_TRUE(Load(sfmJsonToABC3, abcFile2, ESfM_Data(flags)));
        EXPECT_EQ( sfm_data.views.size(), sfmJsonToABC3.views.size());
        EXPECT_EQ( sfm_data.poses.size(), sfmJsonToABC3.poses.size());
        EXPECT_EQ( sfm_data.intrinsics.size(), sfmJsonToABC3.intrinsics.size());
        EXPECT_EQ( sfm_data.structure.size(), sfmJsonToABC3.structure.size());
        EXPECT_EQ( sfm_data.control_points.size(), sfmJsonToABC3.control_points.size());
    }

    // Export as JSON
    jsonFile2 = "jsonToABC2.json";
    {
        EXPECT_TRUE(Save(
        sfmJsonToABC3,
        jsonFile2.c_str(),
        ESfM_Data(flags)));
    }

}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

